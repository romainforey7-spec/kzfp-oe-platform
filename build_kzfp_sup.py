"""Generic KZFP supplement sheet - S6 Fig (ZNF498) and S7 Fig (ZNF18).

S7 Fig (ZNF18) carries five panels: A anti-HA blot, B viability, C K562
expression density, D the KZFP across 1,206 cell lines, E the targets across
the HPA single-cell types.

S6 Fig (ZNF498) carries six, because the target qPCR sits between D and the
heatmap, which is therefore lettered F rather than E:
  A blot | B viability | C density | D cell lines
  E  qPCR of 16 ZNF498 target genes, 3 biological replicates
  F  targets across the single-cell types
The re-lettering is deliberate: panel letters must run in reading order, so the
qPCR takes E and pushes the heatmap to F.

Note on what is NOT here. The published Fig S4 had two extra panels, AcTub(K40)
ciliation micrographs (panel F) and an induction qPCR (panel G). The micrographs
now live in main Fig 4E alongside the ciliation quantification and are not
duplicated here; the induction qPCR is dropped and panel E is the 16-gene target
qPCR instead, recovered from the chart object in the published deck.

Panel E/F extends main Fig 4C / 5C, which show only the six highest-expressing
cell types. ZNF498 is named ZSCAN25 in the Human Protein Atlas tables; the
lookup goes through HPA_ALIAS.
"""
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram

import build_figS2_fig as BS
from build_figS4_fig import SHORT_CT, LT, CMAP_HM
from plos_export import PK  # pins output to the figure's declared size
from presto_screen import series_for

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

HPA_ALIAS = {"ZNF498": "ZSCAN25", "WDR78": "DNAI4"}
N_CT = 40
C_HI = "#1F6FB4"

#: per-KZFP sheet configuration
CFG = {
    "ZNF498": dict(
        sheet="S6", out="out/figS6", figno="4",
        viab=None,   # from presto_screen.series_for("ZNF498")
        call=("STMN3", "DYRK1A", "PAFAH1B1", "TAF1"),
        hi=("Excitatory neurons", "Inhibitory neurons",
            "Oligodendrocyte precursor cells", "Oligodendrocytes",
            "Astrocytes", "Horizontal cells", "Bipolar cells"),
        hi_label="neural",
        height=8.70, reserve_cilia=True),
    "ZNF18": dict(
        sheet="S7", out="out/figS7", figno="5",
        viab=None,   # from presto_screen.series_for("ZNF18")
        call=("SPATA33", "PCDH1", "TENT5C", "DYNC2I2", "BANF1"),
        hi=("Late spermatids", "Early spermatids", "Spermatocytes",
            "Spermatogonia"),
        hi_label="germ-cell",
        height=8.70, reserve_cilia=True),
}


def panel_E(F, sc, targets, kz, cfg, vmax=3, n_ct=N_CT, rects=None):
    """Targets across the single-cell types ranked by the KZFP."""
    hpa = HPA_ALIAS.get(kz, kz)
    r = sc[sc["Gene name"] == hpa][["Cell type", "nTPM"]].copy()
    r = r.sort_values("nTPM", ascending=False)
    cols = r["Cell type"].tolist()[:n_ct]
    barv = dict(zip(r["Cell type"], r.nTPM))
    tg = [HPA_ALIAS.get(t, t) for t in targets]
    m = sc[sc["Gene name"].isin(tg) & sc["Cell type"].isin(cols)]
    W = m.pivot_table(index="Gene name", columns="Cell type", values="nTPM",
                      aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    if len(C) > 2:
        Z = linkage(C.values, method="complete", metric="euclidean")
        C = C.iloc[dendrogram(Z, no_plot=True)["leaves"]]
    x0, w = 0.145, 0.730
    rb, rh, rx = rects or ([x0, 0.566, w, 0.044], [x0, 0.272, w, 0.288],
                           [x0, 0.164, w, 0.100])
    axbar, axhm, axbox = F.add_axes(rb), F.add_axes(rh), F.add_axes(rx)
    axbar.bar(range(len(cols)), [barv[c] for c in cols], width=0.74,
              color="#4D4D4D", edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5)
    axbar.set_xticks([])
    _sym = f"${kz}$" if hpa == kz else f"${kz}$\n(${hpa}$)"
    axbar.set_ylabel(f"{_sym}\n(nTPM)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbar.spines[["top", "right"]].set_visible(False)
    axbar.text(1.0, 1.10, f"{len(cols)} highest-${kz}$ cell types of "
               f"{sc['Cell type'].nunique()}; {len(C)} targets",
               transform=axbar.transAxes, ha="right", va="bottom", fontsize=8)
    axhm.imshow(C.values, aspect="auto", cmap=CMAP_HM, vmin=-vmax, vmax=vmax,
                interpolation="nearest")
    axhm.set_xticks([])
    axhm.set_yticks([])
    axhm.set_ylabel(f"${kz}$ targets", fontsize=8, labelpad=7)
    for sp in axhm.spines.values():
        sp.set_visible(False)
    call = [(i, g) for i, g in enumerate(C.index) if g in cfg["call"]]
    n = len(C)
    ys = (np.linspace(n * 0.06, n * 0.94, len(call)) if len(call) > 1
          else [n * 0.5])
    for (i, g), yy in zip(call, ys):
        axhm.annotate(g, xy=(len(cols) - 0.4, i), xytext=(len(cols) + 1.2, yy),
                      fontsize=8, style="italic", color=C_HI, va="center",
                      ha="left", annotation_clip=False,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=C_HI,
                                      shrinkA=0, shrinkB=1))
    cax = F.add_axes([x0 - 0.064, rh[1] + 0.030, 0.009, 0.060])
    cb = F.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax),
                                          cmap=CMAP_HM), cax=cax,
                    ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2,
                   linespacing=1.1)
    cax.yaxis.set_label_position("left")
    data = [L[c].dropna().values for c in cols]
    axbox.boxplot(data, positions=range(len(cols)), widths=0.62,
                  patch_artist=True, showfliers=False, whis=0, showcaps=False,
                  medianprops=dict(color="black", lw=0.6),
                  boxprops=dict(facecolor="#D9D9D9", edgecolor="#6E6E6E", lw=0.4),
                  whiskerprops=dict(lw=0))
    for i, v in enumerate(data):
        axbox.plot([i], [np.mean(v)], marker="_", color="red", ms=3.6, mew=0.8)
    axbox.set_xlim(-0.5, len(cols) - 0.5)
    axbox.set_xticks(range(len(cols)))
    axbox.set_xticklabels([SHORT_CT.get(c, c) for c in cols], rotation=90,
                          fontsize=8)
    for tl, c in zip(axbox.get_xticklabels(), cols):
        if c in cfg["hi"]:
            tl.set_color(C_HI)
            tl.set_fontweight("bold")
    axbox.tick_params(axis="x", pad=1.0, length=1.6)
    axbox.set_ylabel("Targets\nlog(nTPM+1)", fontsize=8, linespacing=1.1,
                     labelpad=1.5)
    axbox.spines[["top", "right"]].set_visible(False)
    return cols, C


def panel_qpcr(F, QP, rect, letter_xy, letter, kz="ZNF498"):
    """RPE1 qPCR of candidate targets, three biological replicates.

    Values are fold change over the GFP control, which is normalised to 1 in
    every replicate, so the test is a one-sample t-test on log2 fold against 0.
    Plotted on a log2 axis because the range spans 0.70 to 106-fold.
    """
    q = QP[QP.gene != "TBP"].copy()
    rc = [c for c in q.columns if re.fullmatch(r"rep\d+", c)]
    R = (q[rc].astype(float).values if rc
         else np.empty((len(q), 0)))                    # genes x replicates
    ax = F.add_axes(rect)
    x = np.arange(len(q))
    lg = np.log2(q.mean_fold.values)
    up = lg > 0
    ax.bar(x, lg, width=0.66, color=np.where(up, "#C0392B", "#2E6DA4"),
           edgecolor="black", linewidth=0.5, zorder=3)
    for i in range(len(q)):
        v = np.log2(R[i][np.isfinite(R[i])])
        if not len(v):
            continue
        ax.plot(np.full(len(v), i) + np.linspace(-0.16, 0.16, len(v)), v, "o",
                ms=1.9, mfc="white", mec="black", mew=0.45, linestyle="none",
                zorder=5)
    for i, r in enumerate(q.itertuples()):
        if r.mark == "ns":
            continue
        v = np.log2(R[i][np.isfinite(R[i])])
        top = max(v.max() if len(v) else -np.inf, np.log2(r.mean_fold))
        ax.text(i, top + 0.30, r.mark, ha="center", va="bottom", fontsize=8,
                zorder=6)
    ax.axhline(0, lw=0.8, color="black", zorder=4)
    ax.set_xticks(x)
    ax.set_xticklabels(q.gene, fontsize=8, rotation=45, ha="right",
                       rotation_mode="anchor", style="italic")
    ax.set_xlim(-0.7, len(q) - 0.3)
    ax.set_ylabel("$\\log_2$ fold\nvs GFP", fontsize=8, labelpad=1.0)
    ax.tick_params(labelsize=8, length=2.4, width=0.8, pad=1.5)
    ax.spines[["top", "right"]].set_visible(False)
    # Some qPCR tables carry only the per-gene mean and SD (no rep1..repN
    # columns), in which case there are no individual points to bound the axis
    # and the SD sets the range instead.
    fin = np.log2(R[np.isfinite(R)]) if R.size else np.array([])
    if fin.size:
        lo = min(fin.min(), lg.min())
        hi = max(fin.max(), lg.max())
    else:
        sd = np.log2(1.0 + (q.sd_fold.values / np.maximum(q.mean_fold.values, 1e-9)))
        lo = float(np.nanmin(lg - sd))
        hi = float(np.nanmax(lg + sd))
    # headroom solved from font metrics so the tallest asterisk clears the axis
    h_pt = rect[3] * F.get_size_inches()[1] * 72.0
    need = 9.6                                    # 8 pt glyph + leading
    span0 = hi - lo + 0.6
    head = (0.30 + need * span0 / h_pt) / max(1e-6, 1.0 - need / h_pt)
    ax.set_ylim(lo - 0.6, hi + max(head, 0.8))
    nmin, nmax = int(q.n_rep.min()), int(q.n_rep.max())
    nstr = f"n = {nmax}" if nmin == nmax else f"n = {nmin}-{nmax}"
    ax.text(0.0, 1.04, f"hTERT-RPE1, {kz} induction; {nstr} biological "
            f"replicates, BH over {len(q)} genes", transform=ax.transAxes,
            fontsize=8, ha="left", va="bottom", color="#6E6E6E")
    LT(F, *letter_xy, letter)


def _reserve(F, rect, title, body, letter_xy, letter):
    ax = F.add_axes(rect)
    ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#F4F4F4",
                           edgecolor="#BBBBBB", lw=0.6))
    ax.text(0.5, 0.92, title, ha="center", va="top", fontsize=8,
            color="#5A5A5A", fontweight="bold")
    ax.text(0.5, 0.40, body, ha="center", va="center", fontsize=8,
            color="#707070", linespacing=1.35)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, *letter_xy, letter)


VIAB = {}


def build(kz, S45, SS, blot, QP=None):
    cfg = CFG[kz]
    plt.close("all")
    F = plt.figure(figsize=(7.48, cfg["height"]))

    ax = F.add_axes([0.085, 0.838, 0.250, 0.118])
    ax.imshow(blot)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.012, 0.986, "A")

    axb = F.add_axes([0.520, 0.846, 0.240, 0.100])
    rep, wells = VIAB.get(kz, (None, None))
    if rep is None:
        days = [0, 4, 7, 9]
        V = cfg["viab"]
        if V is None:                      # read the recomputed series
            r0 = series_for(kz)[0]
            V = {"D0": 1.0, **{f"D{int(d)}": float(m) for d, m in zip(r0.day, r0["mean"])}}
        axb.plot(days, [V[f"D{d}"] for d in days], marker="o", ms=4,
                 lw=1.1, color="#B00000")
        note = "single induction series"
    else:
        axb.errorbar([0] + list(rep.day), [1.0] + list(rep["mean"]),
                     yerr=[0.0] + list(rep.err), marker="o", ms=4, lw=1.1,
                     color="#B00000", ecolor="#B00000", elinewidth=0.8,
                     capsize=2.0, capthick=0.8, zorder=3)
        if wells is not None and len(wells) > len(rep):
            axb.plot(wells.day, wells.value, "o", ms=2.0, mfc="none",
                     mec="#B00000", mew=0.5, alpha=0.7, zorder=2, clip_on=False)
        nmax, kind = int(rep.n.max()), rep.kind.iloc[0]
        if kind == "plates":
            note = f"mean $\\pm$ SD of {nmax} plates"
        else:
            note = ("one plate; the technical SD over the three Dox and three\n"
                    "No-Dox wells is smaller than the marker")
    axb.axhline(0.85, color="#808080", lw=0.7, ls="--")
    axb.text(9.5, 0.86, "toxicity cut-off", ha="right", va="bottom", fontsize=8,
             color="#808080")
    axb.set_xticks([0, 4, 7, 9])
    axb.set_xlabel("Days after Dox induction", labelpad=1.0)
    axb.set_ylabel(f"Relative viability\n({kz} / GFP)", linespacing=1.15)
    axb.set_ylim(0, 1.15)
    axb.set_xlim(-0.6, 9.8)
    axb.text(0.0, 1.03, note, transform=axb.transAxes, ha="left", va="bottom",
             fontsize=8, color="#6E6E6E", linespacing=1.15)
    axb.spines[["top", "right"]].set_visible(False)
    LT(F, 0.440, 0.986, "B")

    BS.panel_C(F, SS["ALLG"], SS["KZF"], SS["SEL"], gene=kz,
               rect=[0.105, 0.726, 0.280, 0.098], letter=(0.012, 0.836))
    BS.panel_D(F, SS["CL"], SS["SEL"], gene=kz, header=False,
               rect=[0.560, 0.726, 0.280, 0.098], letter=(0.440, 0.836))

    if cfg["reserve_cilia"]:
        panel_qpcr(F, QP, [0.135, 0.606, 0.790, 0.066], (0.012, 0.706), "E", kz)
        rects = ([0.145, 0.482, 0.730, 0.040], [0.145, 0.242, 0.730, 0.234],
                 [0.145, 0.136, 0.730, 0.100])
        cols, C = panel_E(F, S45["sc"], S45[kz]["targets"], kz, cfg, rects=rects)
        LT(F, 0.012, 0.545, "F")
    else:
        cols, C = panel_E(F, S45["sc"], S45[kz]["targets"], kz, cfg)
        LT(F, 0.012, 0.636, "E")
    return F, cols, C


def export(F, kz, name=None):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    cfg = CFG[kz]
    out = cfg["out"]
    os.makedirs(out, exist_ok=True)
    name = name or f"Figure_{cfg['sheet']}"
    F.savefig(f"{out}/{name}_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/{name}_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/{name}_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/{name}_300.png").convert("RGB")
    _im.save(f"{out}/{name}_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
