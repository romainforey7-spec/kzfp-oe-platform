"""Panels shared by Figures 4 (ZNF498) and 5 (ZNF18).

Both figures are built from the same three primary panels - volcano with the
peak/no-peak bar, three-stage alluvial, and the target expression heatmap - so
they live here and each figure script only supplies geometry.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram
import rnaseq_panels as RP

CMAP_HM = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
SHORT_BIN = {"<1Kb": "<1", "1 to 10Kb": "1-10", "10 to 50Kb": "10-50"}
LEG_BIN = {"<1Kb": "<1 kb", "1 to 10Kb": "1-10 kb", "10 to 50Kb": "10-50 kb"}
C_CALL = "#B00000"


def volcano(F, rect, bar_rect, g, kz, kb, ymax, claim, leg_xy, xlim=None):
    ax = F.add_axes(rect)
    axb = F.add_axes(bar_rect)
    s = RP.volcano_panel(ax, axb, g, kz, kb, axleg=None, ymax=ymax, claim=claim)
    axb.set_xticklabels(["Peak", "No peak"], rotation=90, fontsize=8)
    axb.set_yticks([0, 50, 100])
    axb.set_yticklabels(["0", "50", "100"], fontsize=8)
    axb.set_ylabel("% genes", fontsize=8, labelpad=1.0)
    axb.tick_params(axis="y", pad=1.0)
    ax.xaxis.labelpad = 1.0
    keys = [("DOWN", True), ("DOWN", False), ("UP", True),
            ("UP", False), ("NS", True), ("NS", False)]
    labs = ["Down, peak", "Down, no peak", "Up, peak",
            "Up, no peak", "n.s., peak", "n.s., no peak"]
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.VOLC[k], edgecolor="none") for k in keys]
    for v in (20, 40, 60):
        hs.append(mlines.Line2D([], [], marker="o", linestyle="none", color="#828282",
                                markeredgecolor="black", markeredgewidth=0.4,
                                markersize=np.sqrt(RP.bubble_area(v)) * 0.75))
        labs.append(f"score {v}")
    if xlim is not None:
        lo, hi = xlim
        ax.set_xlim(lo, hi)
        off = g[(g["l2fc"] > hi) | (g["l2fc"] < lo)]
        if len(off):
            # genes outside the requested axis are drawn as triangles ON the boundary
            # with a printed count, so narrowing the axis never hides data
            yv = np.minimum(off["nl2p"].values.astype(float), ymax)
            xv = np.where(off["l2fc"].values > hi, hi, lo)
            cl = [RP.VOLC[(d if d in ("UP", "DOWN") else "NS", bool(h))]
                  for d, h in zip(off["Direction"], off["haspeak"])]
            ax.scatter(xv, yv, marker=">", s=7, c=cl, edgecolors="none",
                       clip_on=False, zorder=6)
            # top-left corner: empty in every volcano here, and clear of the
            # neighbouring bar axes whose rotated y-label reaches back
            # across this panel's right edge
            ax.text(lo + 0.03 * (hi - lo), ymax * 0.99,
                    f"{len(off)} off scale", fontsize=8, ha="left",
                    va="top", color="#444444")
    ax.legend(hs, labs, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=leg_xy, bbox_transform=F.transFigure,
              handlelength=1.0, handleheight=1.0, handletextpad=0.4,
              labelspacing=0.30, fontsize=8)
    return s


def alluvial(F, rect, s, leg_xy, node_w=0.085):
    ax = F.add_axes(rect)
    RP.alluvial_panel(ax, s, node_w=node_w, gap_frac=0.013)
    ax.set_xticklabels(["Status", "TSS-peak\ndistance (kb)", "Assignation"])
    for t in ax.texts:
        if t.get_rotation() == 90 and t.get_text() in SHORT_BIN:
            t.set_text(SHORT_BIN[t.get_text()])
        elif "_" in t.get_text():
            t.set_text(t.get_text().replace("_", " "))
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.ALLUV[(d, b)], edgecolor="none")
          for d in ("DOWN", "UP") for b in RP.BINS]
    ls = [f"{d[0]}{d[1:].lower()}, {LEG_BIN[b]}" for d in ("DOWN", "UP") for b in RP.BINS]
    ax.legend(hs, ls, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=leg_xy, bbox_transform=F.transFigure,
              handlelength=1.0, handleheight=1.0, handletextpad=0.4,
              labelspacing=0.30, fontsize=8)
    return ax


# HPA gene symbols differ from the names used in the manuscript
HPA_ALIAS = {"ZNF498": "ZSCAN25", "WDR78": "DNAI4"}


def target_heatmap(F, sc, targets, kz, geom, call_genes=(), vmax=3, call_x=0.5,
                   index_col="Cell type",
                   n_cols=6, short=None):
    """Bar of KZFP expression + centred target heatmap + per-column box plot.

    The 34-35 target rows cannot carry 8 pt labels, so rows are unlabelled and
    only the genes named in the text are called out with leaders.
    """
    short = short or {}
    hpa = HPA_ALIAS.get(kz, kz)
    r = sc[sc["Gene name"] == hpa][[index_col, "nTPM"]].copy()
    r = r.sort_values(["nTPM", index_col], ascending=[False, True]).reset_index(drop=True)
    hi = r.head(3)
    lo = r.sort_values(["nTPM", index_col], ascending=[True, True]).head(n_cols - 3)
    lo = lo.sort_values(["nTPM", index_col], ascending=[False, True])
    sel = pd.concat([hi, lo]).drop_duplicates(index_col)
    cols = sel[index_col].tolist()
    barvals = dict(zip(sel[index_col], sel.nTPM))
    m = sc[sc["Gene name"].isin(targets) & sc[index_col].isin(cols)]
    W = m.pivot_table(index="Gene name", columns=index_col, values="nTPM", aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    Z = linkage(C.values, method="complete", metric="euclidean")
    order = dendrogram(Z, no_plot=True)["leaves"]
    C, L = C.iloc[order], L.iloc[order]

    axbar = F.add_axes(geom["bar"])
    axhm = F.add_axes(geom["hm"])
    axbox = F.add_axes(geom["box"])
    axbar.bar(range(len(cols)), [barvals[c] for c in cols], width=0.70,
              color="#4D4D4D", edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5)
    axbar.set_xticks([])
    nm = f"${hpa}$" if hpa != kz else f"${kz}$"
    axbar.set_ylabel(nm + "\nnTPM", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbar.spines[["top", "right"]].set_visible(False)
    axbar.text(1.0, 1.10, f"n = {len(C)} targets", transform=axbar.transAxes,
               ha="right", va="bottom", fontsize=8)

    axhm.imshow(C.values, aspect="auto", cmap=CMAP_HM, vmin=-vmax, vmax=vmax,
                interpolation="nearest")
    axhm.set_xticks([])
    axhm.set_yticks([])
    axhm.set_ylabel("Targets", fontsize=8, labelpad=7)
    for sp in axhm.spines.values():
        sp.set_visible(False)
    present = [(i, gn) for i, gn in enumerate(C.index) if gn in call_genes]
    if present:
        n = len(C)
        # keep each label beside its own row: start at the true y and push apart
        # only as far as the minimum legible gap requires, so leaders stay short.
        # the gap is derived from font metrics, not a heuristic on n: one 8 pt line
        # needs 1.25 * 8 pt of vertical room, and the panel gives
        # (axes height in points / n) points per heatmap row.
        _h_pt = axhm.get_position().height * F.get_size_inches()[1] * 72.0
        min_gap = max(1.0, (8.0 * 1.55) / (_h_pt / n))
        ys = [float(i) for i, _ in present]
        for _ in range(60):
            moved = False
            for j in range(len(ys) - 1):
                d = ys[j + 1] - ys[j]
                if d < min_gap:
                    sh = (min_gap - d) / 2.0
                    ys[j] -= sh
                    ys[j + 1] += sh
                    moved = True
            ys = [min(max(v, 0.0), n - 1.0) for v in ys]
            if not moved:
                break
        for (i, gn), yy in zip(present, ys):
            axhm.annotate(gn, xy=(len(cols) - 0.4, i), xytext=(len(cols) + call_x, yy),
                          fontsize=8, style="italic", color=C_CALL, va="center",
                          ha="left", annotation_clip=False,
                          arrowprops=dict(arrowstyle="-", lw=0.6, color=C_CALL,
                                          shrinkA=0, shrinkB=1))
    cax = F.add_axes(geom["cbar"])
    horiz = geom["cbar"][2] > geom["cbar"][3]
    cb = F.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax), cmap=CMAP_HM),
                    cax=cax, ticks=[-vmax, 0, vmax],
                    orientation="horizontal" if horiz else "vertical")
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    if horiz:
        cb.ax.xaxis.set_ticks_position("bottom")
        cb.set_label("Centred log(nTPM+1)", fontsize=8, labelpad=1.5)
    else:
        cb.ax.yaxis.set_ticks_position("left")
        cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2,
                       linespacing=1.1)
        cax.yaxis.set_label_position("left")

    data = [L[c].dropna().values for c in cols]
    axbox.boxplot(data, positions=range(len(cols)), widths=0.62, patch_artist=True,
                  showfliers=False, whis=0, showcaps=False,
                  medianprops=dict(color="black", lw=0.7),
                  boxprops=dict(facecolor="#D9D9D9", edgecolor="#6E6E6E", lw=0.5),
                  whiskerprops=dict(lw=0))
    for i, v in enumerate(data):
        axbox.plot([i], [np.mean(v)], marker="_", color="red", ms=5.0, mew=0.9)
    axbox.set_xlim(-0.5, len(cols) - 0.5)
    axbox.set_xticks(range(len(cols)))
    axbox.set_xticklabels([short.get(c, c) for c in cols], rotation=90, fontsize=8,
                          ha="center", va="top")
    axbox.tick_params(axis="x", pad=1.5, length=1.8)
    axbox.set_ylabel("log(nTPM+1)", fontsize=8, labelpad=1.5)
    axbox.spines[["top", "right"]].set_visible(False)
    return cols, C
