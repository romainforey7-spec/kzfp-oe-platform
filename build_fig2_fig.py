
"""Figure 2 assembly for PLOS Genetics (7.48 x 8.70 in, all text >= 8 pt).

Layout: A-D stacked in a left-hand column, E as four stacked locus maps on the
right, F (author artwork) at the bottom left.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy import stats as st
import rnaseq_panels as RP
from plos_export import PK  # pins output to the figure's declared size

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
RP.FS.update(star=8.0, leg=8.0, legs=8.0, node=8.0)

OUT = "out/fig2"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70

CMAP_HM = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
TRK_ORD = [("ZNF43", "#2B2B2B", "ZNF43-HA"),
           ("Liver H3K4me3", "#2E7D32", "Liver H3K4me3"),
           ("Thymus H3K4me3", "#66A266", "Thymus H3K4me3"),
           ("Liver H3K9me3", "#5E35B1", "Liver H3K9me3"),
           ("Thymus H3K9me3", "#9C84D4", "Thymus H3K9me3")]
RC2 = {"LTR/ERV1": "#1F6FB4", "SINE/Alu": "#9A9A9A", "LINE/L1": "#7E57C2",
       "SINE/MIR": "#BDBDBD", "LTR/ERVL": "#4FA3D1", "LTR/ERVL-MaLR": "#7FC4E8",
       "Simple_repeat": "#CCCCCC", "Low_complexity": "#DDDDDD"}
RC_OTHER = "#E0E0E0"
HPA_ALIAS = {"WDR78": "DNAI4"}

LX, LW = 0.098, 0.300
EX, EW = 0.545, 0.420
EY = [0.794, 0.602, 0.410, 0.218, 0.026]
EH = 0.182


def LT(F2, x, y, s):
    F2.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def panel_A(F2, g):
    ax = F2.add_axes([LX, 0.830, 0.185, 0.118])
    axb = F2.add_axes([LX + 0.212, 0.830, 0.030, 0.118])
    s = RP.volcano_panel(ax, axb, g, "ZNF43", 10, axleg=None, ymax=8.0, claim="DOWN")
    axb.set_xticklabels(["Peak", "No peak"], rotation=90, fontsize=8)
    ax.xaxis.labelpad = 1.0
    off = g[g.nl2p > 8.0]
    if len(off):
        ax.scatter(off.l2fc, np.full(len(off), 8.0 * 0.985), marker="^", s=8.0,
                   facecolor=[RP.VOLC[(c, h)] for c, h in zip(off.cat, off.haspeak)],
                   edgecolor="black", linewidths=0.25, zorder=5, clip_on=False)
        x0 = float(off.l2fc.iloc[0])
        ax.annotate("1 above axis", xy=(x0, 8.0 * 0.985), xytext=(x0 + 0.20, 8.0 * 0.88),
                    fontsize=8, color="#404040", ha="left", va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="#808080",
                                    shrinkA=0, shrinkB=2))
    keys = [("DOWN", True), ("DOWN", False), ("UP", True),
            ("UP", False), ("NS", True), ("NS", False)]
    labs = ["Down, peak", "Down, no peak", "Up, peak",
            "Up, no peak", "n.s., peak", "n.s., no peak"]
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.VOLC[k], edgecolor="none") for k in keys]
    import matplotlib.lines as _ml
    for v in (20, 40, 60):
        hs.append(_ml.Line2D([], [], marker="o", linestyle="none", color="#828282",
                             markeredgecolor="black", markeredgewidth=0.4,
                             markersize=np.sqrt(RP.bubble_area(v)) * 0.75))
        labs.append(f"score {v}")
    ax.legend(hs, labs, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=(1.33, 1.04), handlelength=1.0, handleheight=1.0,
              handletextpad=0.4, labelspacing=0.30, fontsize=8)
    # peak-score key: outside the data area, under the % bar
    LT(F2, 0.012, 0.996, "A")
    return s


def panel_B(F2, s):
    ax = F2.add_axes([LX, 0.645, 0.232, 0.120])
    RP.alluvial_panel(ax, s, node_w=0.085, gap_frac=0.013)
    ax.set_xticklabels(["Status", "Distance\nTSS - peak", "Assignation"])
    for t in ax.texts:
        if t.get_rotation() == 90:
            t.set_text(t.get_text().replace("Kb", "kb").replace(" to ", "-"))
    SHORT = {"<1Kb": "<1 kb", "1 to 10Kb": "1-10 kb", "10 to 50Kb": "10-50 kb"}
    hs = [Rectangle((0, 0), 1, 1, facecolor=RP.ALLUV[(d, b)], edgecolor="none")
          for d in ("DOWN", "UP") for b in RP.BINS]
    ls = [f"{d[0]}{d[1:].lower()}, {SHORT[b]}" for d in ("DOWN", "UP") for b in RP.BINS]
    ax.legend(hs, ls, ncol=1, frameon=False, loc="upper left",
              bbox_to_anchor=(1.03, 1.05), handlelength=1.0, handleheight=1.0,
              handletextpad=0.4, labelspacing=0.30, fontsize=8)
    LT(F2, 0.012, 0.790, "B")


def panel_C(F2, expr, targets, gene="ZNF43", vmax=4):
    targets = [HPA_ALIAS.get(t, t) for t in targets]
    r = expr[expr["Gene name"] == gene][["Tissue", "nTPM"]].copy()
    r = r.sort_values(["nTPM", "Tissue"], ascending=[False, True]).reset_index(drop=True)
    hi = r.head(3)
    lo = r.sort_values(["nTPM", "Tissue"], ascending=[True, True]).head(3)
    lo = lo.sort_values(["nTPM", "Tissue"], ascending=[False, True])
    sel = pd.concat([hi, lo]).drop_duplicates("Tissue")
    cols = sel["Tissue"].tolist()
    barvals = dict(zip(sel["Tissue"], sel.nTPM))
    m = expr[expr["Gene name"].isin(targets) & expr["Tissue"].isin(cols)]
    W = m.pivot_table(index="Gene name", columns="Tissue", values="nTPM", aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    Z = linkage(C.values, method="complete", metric="euclidean")
    order = dendrogram(Z, no_plot=True)["leaves"]
    C, L = C.iloc[order], L.iloc[order]

    x0, w = 0.215, 0.140
    axbar = F2.add_axes([x0, 0.578, w, 0.019])
    axden = F2.add_axes([x0 - 0.019, 0.474, 0.016, 0.098])
    axhm = F2.add_axes([x0, 0.474, w, 0.098])
    axbox = F2.add_axes([x0, 0.426, w, 0.034])

    axbar.bar(range(len(cols)), [barvals[c] for c in cols], width=0.70,
              color="#4D4D4D", edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5)
    axbar.set_xticks([])
    axbar.set_ylabel("$ZNF43$\n(nTPM)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbar.spines[["top", "right"]].set_visible(False)

    with plt.rc_context({"lines.linewidth": 0.5}):
        dendrogram(Z, ax=axden, orientation="left", no_labels=True,
                   link_color_func=lambda *a: "black")
    axden.invert_yaxis()
    axden.set_xticks([]); axden.set_yticks([])
    axden.axis("off")

    axhm.imshow(C.values, aspect="auto", cmap=CMAP_HM, vmin=-vmax, vmax=vmax,
                interpolation="nearest")
    axhm.set_xticks([])
    axhm.set_yticks(range(len(C)))
    axhm.set_yticklabels(C.index, fontsize=8, style="italic")
    axhm.yaxis.tick_right()
    axhm.tick_params(axis="y", length=0, pad=1.5)
    for sp in axhm.spines.values():
        sp.set_visible(False)

    cax = F2.add_axes([x0 - 0.056, 0.482, 0.009, 0.044])
    cb = F2.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax), cmap=CMAP_HM),
                     cax=cax, ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2, linespacing=1.1)
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
    axbox.set_xticklabels(cols, rotation=45, fontsize=8, ha="right", va="top",
                          rotation_mode="anchor")
    axbox.tick_params(axis="x", pad=1.5, length=1.8)
    axbox.set_ylabel("Targets\nlog(nTPM+1)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbox.spines[["top", "right"]].set_visible(False)
    LT(F2, 0.012, 0.620, "C")
    return cols


def panel_D(F2, M, rho, pv):
    ax = F2.add_axes([LX + 0.022, 0.248, 0.268, 0.092])
    ax.scatter(M.x, M.y, s=22, facecolor="#B0B0B0", edgecolor="#707070",
               linewidth=0.5, zorder=2)
    sl, ic, rr, pp, se = st.linregress(M.x, M.y)
    xx = np.linspace(M.x.min(), M.x.max(), 50)
    ax.plot(xx, ic + sl * xx, color="black", lw=1.0, ls="--", zorder=3)
    PL = {"liver": (0.055, 0.0, "#1F4E9C", "left"),
          "bone marrow": (0.075, 0.075, "#C0392B", "left"),
          "thymus": (0.075, -0.075, "#C0392B", "left")}
    for t, (dx, dy, col, ha) in PL.items():
        x0, y0 = M.loc[t, "x"], M.loc[t, "y"]
        ax.scatter([x0], [y0], s=30, facecolor=col, edgecolor="black",
                   linewidth=0.5, zorder=4)
        ax.annotate(t, xy=(x0, y0), xytext=(x0 + dx, y0 + dy), fontsize=8,
                    fontweight="bold", color=col, ha=ha, va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=col, shrinkA=0,
                                    shrinkB=4) if t != "liver" else None)
    ax.text(0.97, 0.97, f"Spearman $\\rho$ = {rho:.2f}\np = {pv:.1e}\nn = {len(M)} tissues",
            transform=ax.transAxes, fontsize=8, ha="right", va="top", color="#404040",
            linespacing=1.15)
    ax.set_xlabel("$ZNF43$ expr.  log$_{10}$(nTPM+1)", labelpad=1.0)
    ax.set_ylabel("Mean target expr.\nlog$_{10}$(nTPM+1)", linespacing=1.2)
    ax.set_xlim(M.x.min() - 0.10, M.x.max() + 0.40)
    ax.set_ylim(M.y.min() - 0.10, M.y.max() + 0.12)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F2, 0.012, 0.360, "D")


def hspread(xs_, lo, hi, gap):
    o = np.argsort(xs_)
    out = np.array(xs_, float)
    cur = lo
    for i in o:
        out[i] = max(out[i], cur)
        cur = out[i] + gap
    if out.max() > hi:
        out = out - (out.max() - hi)
    return out


def locus_block(F2, g, x0, y0, w, h, WIN2, BW, SUM, TRK, RPALL, TB, label, nb=1000):
    W = WIN2[g]
    ch, a, b, pk, key = W["chrom"], W["a"], W["b"], W["pk"], W["own"]
    xs = np.linspace(a, b, nb)
    sig = {k: np.nan_to_num(BW[k].values(ch, a, b, nbins=nb, agg="mean")) / SUM[k]
           for k in TRK}
    ym = {"ZNF43": max(sig["ZNF43"].max() * 1.12, 1),
          "K4": max(max(sig["Liver H3K4me3"].max(), sig["Thymus H3K4me3"].max()) * 1.12, 1),
          "K9": max(max(sig["Liver H3K9me3"].max(), sig["Thymus H3K9me3"].max()) * 1.12, 1)}
    hann, hgap, htitle = 0.042, 0.0024, 0.015
    hsig = (h - hann - 5 * hgap - htitle) / 5
    t0 = int(TB[TB.symbol == key].tss.iloc[0])
    axes = []
    for i, (k, c, short) in enumerate(TRK_ORD):
        yy = y0 + hann + (4 - i) * (hsig + hgap)
        ax = F2.add_axes([x0, yy, w, hsig])
        grp = "ZNF43" if k == "ZNF43" else ("K4" if "K4" in k else "K9")
        ax.fill_between(xs, 0, sig[k], color=c, lw=0, zorder=2, rasterized=True)
        ax.axvspan(pk[0], pk[1], color="#FFD24D", alpha=0.45, lw=0, zorder=1)
        ax.set_xlim(a, b)
        ax.set_ylim(0, ym[grp])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.008, 0.90, short, transform=ax.transAxes, ha="left", va="top",
                fontsize=8, color=c, zorder=6)
        ax.text(0.995, 0.90, f"{ym[grp]:.0f}", transform=ax.transAxes, ha="right",
                va="top", fontsize=8, color="#404040", zorder=6)
        if "K4" in k:
            ax.axvline(t0, color="#2E7D32", lw=0.5, ls=":", zorder=3)
        ax.spines[["top", "right"]].set_visible(False)
        axes.append(ax)
    axes[0].set_title(f"{label}   {ch}:{a/1e6:.3f}\u2013{b/1e6:.3f} Mb (hg19)",
                      fontsize=8, pad=1.8)

    # Annotation row: three bands, each with its own row label at the left edge
    # so nothing can overlap - peak (top), repeat class (middle), TSS (bottom).
    axn = F2.add_axes([x0, y0, w, hann])
    rp = RPALL[(RPALL.rchr == ch) & RPALL.rstart.between(a, b)][
        ["rname", "rstart", "rend", "rclass"]].drop_duplicates()
    span = b - a
    BAND = {"peak": 2.5, "rep": 1.5, "tss": 0.5}
    HB = 0.17
    axn.add_patch(plt.Rectangle((pk[0], BAND["peak"] - HB),
                                max(pk[1] - pk[0], span * 0.006), 2 * HB,
                                facecolor="#B00000", edgecolor="none"))
    for r_ in rp.itertuples():
        axn.add_patch(plt.Rectangle((r_.rstart, BAND["rep"] - HB),
                                    max(r_.rend - r_.rstart, span * 0.004), 2 * HB,
                                    facecolor=RC2.get(r_.rclass, RC_OTHER),
                                    edgecolor="none"))
    inpk = rp[(rp.rstart < pk[1]) & (rp.rend > pk[0])].copy()
    cls, ccol = None, "#555555"
    if len(inpk):
        inpk["ov"] = np.minimum(inpk.rend, pk[1]) - np.maximum(inpk.rstart, pk[0])
        erv = inpk[inpk.rclass == "LTR/ERV1"]
        e = (erv.sort_values("ov", ascending=False).iloc[0] if len(erv)
             else inpk.sort_values("ov", ascending=False).iloc[0])
        cls, ccol = e.rclass, RC2.get(e.rclass, "#555555")
    axn.text(a, BAND["peak"], "peak", ha="left", va="center", fontsize=8,
             color="#B00000")
    if cls:
        axn.text(a, BAND["rep"], cls, ha="left", va="center", fontsize=8, color=ccol)
    axn.text(a, BAND["tss"], "TSS", ha="left", va="center", fontsize=8, color="#B00000")
    # TSS of the target gene(s): a tick at the exact coordinate, no strand arrow
    tb = W["tss"]
    want = ("APOL1", "APOL2") if g == "APOL1" else (key,)
    hits = tb[tb.symbol.isin(want)].drop_duplicates("symbol").sort_values("tss")
    if len(hits):
        xs_ = hspread(list(hits.tss + span * 0.015), a + span * 0.12, b, span * 0.19)
        for (r_, xt) in zip(hits.itertuples(), xs_):
            axn.plot([r_.tss, r_.tss],
                     [BAND["tss"] - HB * 1.6, BAND["tss"] + HB * 1.6],
                     lw=1.1, color="#B00000", solid_capstyle="butt")
            nm = "DNAI4" if r_.symbol == "WDR78" else r_.symbol
            axn.annotate(nm, xy=(r_.tss, BAND["tss"] + HB * 1.6),
                         xytext=(xt, BAND["tss"]), fontsize=8, style="italic",
                         color="#B00000", ha="left", va="center",
                         arrowprops=dict(arrowstyle="-", lw=0.5, color="#B00000",
                                         shrinkA=1, shrinkB=1))
    axn.set_xlim(a, b)
    axn.set_ylim(0, 3.0)
    axn.set_yticks([])
    axn.set_xticks([])
    axn.set_xticks([]); axn.set_yticks([])
    axn.axis("off")
    return axes[0]


def panel_E(F2, WIN2, BW, SUM, TRK, RPALL, TB, loci):
    for i, (lab, g) in enumerate(loci):
        locus_block(F2, g, EX, EY[i], EW, EH, WIN2, BW, SUM, TRK, RPALL, TB, lab)
    F2.text(0.498, 0.505, "ChIP signal (fold over genome mean)", rotation=90,
            va="center", ha="center", fontsize=8)
    LT(F2, 0.462, 0.996, "E")


def panel_F(F2, img=None):
    """Model drawn natively (model_fig2.py); img is accepted but unused."""
    import model_fig2 as MDL
    ax = F2.add_axes([0.020, 0.010, 0.432, 0.198])
    MDL.draw_model(ax)
    LT(F2, 0.012, 0.228, "F")
    return ax


def build(g, s, tc, targets, MD, rho, pv, WIN2, BW, SUM, TRK, RPALL, TB, loci, modelimg):
    plt.close("all")
    F2 = plt.figure(figsize=(SHEET_W, SHEET_H))
    stA = panel_A(F2, g)
    panel_B(F2, s)
    cols = panel_C(F2, tc, targets)
    panel_D(F2, MD, rho, pv)
    panel_E(F2, WIN2, BW, SUM, TRK, RPALL, TB, loci)
    panel_F(F2, modelimg)
    return F2, stA, cols


def export(F2, out=OUT):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F2.savefig(f"{out}/Figure_2_PLOS.png", dpi=600, **PK(F2))
    F2.savefig(f"{out}/Figure_2_300.png", dpi=300, **PK(F2))
    F2.savefig(f"{out}/Figure_2_PLOS.svg", **PK(F2))
    _im = _PILI.open(f"{out}/Figure_2_300.png").convert("RGB")
    _im.save(f"{out}/Figure_2_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
