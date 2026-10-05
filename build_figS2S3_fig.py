
"""The published Figure S2 spans two pages (panels A-I) and cannot be laid out
legibly at 8 pt on one PLOS sheet: G, H and I each need the full 7.48 in width.
It is therefore split.

S2 Fig  (7.48 x 5.60 in) - A anti-HA blot · B viability · C K562 expression
    density · D the same KZFPs across 1,206 cell lines · E peak-to-TSS distance
    by DE class · F TE composition and the LTR/ERV1 binomial test.
S3 Fig  (7.48 x 8.70 in) - G nine K562 chromatin tracks at the seven target
    peaks · H all 50 HPA tissues ranked by ZNF43 with the eight targets ·
    I liver versus thymus H3K4me3 and H3K9me3 at the seven peaks.
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
import build_figS2_fig as BS

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

OUT = "out/figS2"
os.makedirs(OUT, exist_ok=True)
S2_W, S2_H = 7.48, 5.60
S3_W, S3_H = 7.48, 8.70
CMAP_SIG = LinearSegmentedColormap.from_list(
    "sig", ["#2B0B4F", "#2E6E8E", "#3FA68C", "#9FD44F", "#F7F75A"])
CMAP_HM = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
C_LIV, C_THY = "#C2743A", "#4F79B8"
TRK9 = [("ZNF43-HA", False), ("KAP1", False), ("SETDB1", False),
        ("H3K9me3", True), ("H3K27me3", True), ("H3K4me3", True),
        ("H3K4me1", True), ("H3K9ac", True), ("H3K27ac", True)]
HPA_ALIAS = {"WDR78": "DNAI4"}


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


# ============================================================ S2 Fig (A-F)
def build_S2(S, blot):
    plt.close("all")
    F = plt.figure(figsize=(S2_W, S2_H))
    ax = F.add_axes([0.085, 0.700, 0.283, 0.245])
    ax.imshow(blot)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.012, 0.978, "A")

    axb = F.add_axes([0.580, 0.700, 0.280, 0.188])
    R = S["viab_rep"].sort_values("day")
    x = [0] + list(R.day)
    y = [1.0] + list(R["mean"])
    e = [0.0] + list(R["sd"].fillna(0.0))
    axb.errorbar(x, y, yerr=e, marker="o", ms=4, lw=1.1, color="#B00000",
                 ecolor="#B00000", elinewidth=0.8, capsize=2.0, capthick=0.8,
                 zorder=3)
    # the individual well-sets behind the mean
    W = S["viab_wells"]
    axb.plot(W.day, W.value, "o", ms=2.0, mfc="none", mec="#B00000",
             mew=0.5, alpha=0.7, zorder=2, clip_on=False)
    axb.axhline(0.85, color="#808080", lw=0.7, ls="--")
    axb.text(9.5, 0.86, "toxicity cut-off", ha="right", va="bottom", fontsize=8,
             color="#808080")
    axb.set_xticks([0, 4, 7, 9])
    axb.set_xlabel("Days after Dox induction", labelpad=1.0)
    axb.set_ylabel("Relative viability\n(ZNF43 / GFP)", linespacing=1.15)
    axb.set_ylim(0, 1.15)
    axb.set_xlim(-0.6, 9.8)
    nmax = int(R["n"].max())
    axb.text(0.0, 1.03, f"mean $\\pm$ SD of {nmax} plates "
             "from 2 independent experiments", transform=axb.transAxes,
             ha="left", va="bottom", fontsize=8, color="#6E6E6E")
    axb.spines[["top", "right"]].set_visible(False)
    LT(F, 0.480, 0.978, "B")

    BS.panel_C(F, S["ALLG"], S["KZF"], S["SEL"],
               rect=[0.105, 0.400, 0.300, 0.180], letter=(0.012, 0.628))
    BS.panel_D(F, S["CL"], S["SEL"],
               rect=[0.580, 0.400, 0.280, 0.180], letter=(0.480, 0.628))
    BS.panel_E(F, S["G43"], rect=[0.105, 0.120, 0.300, 0.168],
               letter=(0.012, 0.336))
    BS.panel_F(F, S["dn_class_counts"], S["dn_rows"], S["erv1_obs"], S["erv1_exp"],
               S["n_peaks_rep"], S["erv1_fold"], S["erv1_p"],
               rects=([0.620, 0.212, 0.250, 0.080], [0.620, 0.118, 0.250, 0.040]),
               letter=(0.480, 0.336))
    return F


# ============================================================ S3 Fig (G-I)
def panel_G3(F, PK7, BWP, SUMP, BWH, half=10000, nb=120):
    """Nine K562 tracks at the seven target peaks, one heatmap per track."""
    n = len(TRK9)
    x0, W = 0.142, 0.836
    w = (W - (n - 1) * 0.006) / n
    rows = list(PK7.genes)
    out = {}
    for i, (k, raw) in enumerate(TRK9):
        bw = BWH[k] if raw else BWP[k]
        M = np.vstack([np.nan_to_num(bw.values(r.chrom, r.summit - half,
                                               r.summit + half, nbins=nb,
                                               agg="mean"))
                       for r in PK7.itertuples()])
        if not raw:
            M = M / SUMP[k]
        out[k] = M
        xx = x0 + i * (w + 0.006)
        axh = F.add_axes([xx, 0.760, w, 0.126])
        vmax = float(np.nanpercentile(M, 99.5)) or 1.0
        axh.imshow(M, aspect="auto", cmap=CMAP_SIG, vmin=0, vmax=vmax,
                   extent=[-half / 1000, half / 1000, len(rows), 0],
                   interpolation="nearest")
        if i == 0:
            axh.set_xticks([-10, 0, 10])
            axh.tick_params(axis="x", labelsize=8, length=1.8, pad=1.0)
        else:
            axh.set_xticks([0])
            axh.set_xticklabels(["0"], fontsize=8)
            axh.tick_params(axis="x", labelsize=8, length=1.8, pad=1.0)
        axh.set_title(k, fontsize=8, pad=2.0)
        if i == 0:
            axh.set_yticks(np.arange(len(rows)) + 0.5)
            axh.set_yticklabels(rows, fontsize=8, style="italic")
            axh.tick_params(axis="y", length=0, pad=1.5)
            axh.set_ylabel("Target locus", fontsize=8, labelpad=7)
        else:
            axh.set_yticks([])
        for sp in axh.spines.values():
            sp.set_visible(False)
        # colour bar sits ABOVE the track name, not on top of it
        cax = F.add_axes([xx, 0.918, w, 0.007])
        cb = F.colorbar(plt.cm.ScalarMappable(cmap=CMAP_SIG), cax=cax,
                        orientation="horizontal", ticks=[])
        cb.outline.set_linewidth(0.4)
        cax.text(0, 1.35, "0", transform=cax.transAxes, ha="left", va="bottom",
                 fontsize=8)
        cax.text(1, 1.35, f"{vmax:.0f}" if vmax >= 1 else f"{vmax:.1f}",
                 transform=cax.transAxes, ha="right", va="bottom", fontsize=8)
    F.text(x0 + W / 2, 0.736, "Distance to ZNF43 peak summit (kb); each panel spans "
           "$\\pm$10 kb", ha="center", va="top", fontsize=8)
    F.text(x0, 0.948, "fold over genome mean (ZNF43-HA, KAP1, SETDB1) / "
           "normalised ratio (histone marks)", ha="left", va="bottom", fontsize=8,
           color="#5A5A5A")
    LT(F, 0.012, 0.992, "A")
    return out


def panel_H3(F, tc, targets, gene="ZNF43", vmax=4):
    """All 50 HPA consensus tissues ranked by ZNF43 expression."""
    tg = [HPA_ALIAS.get(t, t) for t in targets]
    r = tc[tc["Gene name"] == gene][["Tissue", "nTPM"]].copy()
    r = r.sort_values("nTPM", ascending=False)
    cols = r["Tissue"].tolist()
    barv = dict(zip(r["Tissue"], r.nTPM))
    m = tc[tc["Gene name"].isin(tg) & tc["Tissue"].isin(cols)]
    W = m.pivot_table(index="Gene name", columns="Tissue", values="nTPM",
                      aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    Z = linkage(C.values, method="complete", metric="euclidean")
    order = dendrogram(Z, no_plot=True)["leaves"]
    C, L = C.iloc[order], L.iloc[order]
    x0, w = 0.150, 0.790
    axbar = F.add_axes([x0, 0.628, w, 0.046])
    axden = F.add_axes([x0 - 0.026, 0.496, 0.022, 0.126])
    axhm = F.add_axes([x0, 0.496, w, 0.126])
    axbox = F.add_axes([x0, 0.378, w, 0.110])
    axbar.bar(range(len(cols)), [barv[c] for c in cols], width=0.72,
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
    cax = F.add_axes([x0 - 0.072, 0.520, 0.009, 0.060])
    cb = F.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax),
                                          cmap=CMAP_HM), cax=cax,
                    ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2, linespacing=1.1)
    cax.yaxis.set_label_position("left")
    data = [L[c].dropna().values for c in cols]
    axbox.boxplot(data, positions=range(len(cols)), widths=0.64, patch_artist=True,
                  showfliers=False, whis=0, showcaps=False,
                  medianprops=dict(color="black", lw=0.6),
                  boxprops=dict(facecolor="#D9D9D9", edgecolor="#6E6E6E", lw=0.4),
                  whiskerprops=dict(lw=0))
    for i, v in enumerate(data):
        axbox.plot([i], [np.mean(v)], marker="_", color="red", ms=3.6, mew=0.8)
    axbox.set_xlim(-0.5, len(cols) - 0.5)
    axbox.set_xticks(range(len(cols)))
    axbox.set_xticklabels(cols, rotation=90, fontsize=8)
    axbox.tick_params(axis="x", pad=1.0, length=1.6)
    axbox.set_ylabel("Targets\nlog(nTPM+1)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbar.text(1.0, 1.10, f"{len(cols)} HPA consensus tissues, ranked by $ZNF43$",
               transform=axbar.transAxes, ha="right", va="bottom", fontsize=8)
    axbox.spines[["top", "right"]].set_visible(False)
    LT(F, 0.012, 0.700, "B")
    return cols, C


def panel_I3(F, PK7, BWL, SUML, half=10000, nb=120):
    """Liver versus thymus at the seven target peaks, both marks."""
    rows = list(PK7.genes)
    res = {}
    for j, mk in enumerate(("H3K4me3", "H3K9me3")):
        xb = 0.088 + j * 0.455
        axl = F.add_axes([xb + 0.072, 0.095, 0.112, 0.150])
        axt = F.add_axes([xb + 0.190, 0.095, 0.112, 0.150])
        axv = F.add_axes([xb + 0.354, 0.095, 0.056, 0.150])
        mats = {}
        for tag, axq in (("Liver", axl), ("Thymus", axt)):
            k = f"{tag} {mk}"
            M = np.vstack([np.nan_to_num(BWL[k].values(r.chrom, r.summit - half,
                                                       r.summit + half, nbins=nb,
                                                       agg="mean")) / SUML[k]
                           for r in PK7.itertuples()])
            mats[tag] = M
        vmax = float(np.nanpercentile(np.vstack(list(mats.values())), 99.5)) or 1.0
        for tag, axq in (("Liver", axl), ("Thymus", axt)):
            axq.imshow(mats[tag], aspect="auto", cmap=CMAP_SIG, vmin=0, vmax=vmax,
                       extent=[-half / 1000, half / 1000, len(rows), 0],
                       interpolation="nearest")
            axq.set_xticks([-10, 0, 10] if tag == "Liver" else [0])
            axq.tick_params(axis="x", labelsize=8, length=1.8, pad=1.0)
            axq.set_title(tag, fontsize=8, pad=1.5)
            if tag == "Liver" and j == 0:
                axq.set_yticks(np.arange(len(rows)) + 0.5)
                axq.set_yticklabels(rows, fontsize=8, style="italic")
                axq.tick_params(axis="y", length=0, pad=1.5)
            else:
                axq.set_yticks([])
            for sp in axq.spines.values():
                sp.set_visible(False)
        lv = mats["Liver"].mean(1)
        th = mats["Thymus"].mean(1)
        t, p = st.ttest_rel(lv, th)
        res[mk] = dict(liver=float(lv.mean()), thymus=float(th.mean()),
                       t=float(t), p=float(p), n=len(lv))
        pv = axv.violinplot([lv, th], positions=[0, 1], widths=0.72,
                            showextrema=False, showmedians=False)
        for i, b in enumerate(pv["bodies"]):
            b.set_facecolor(C_LIV if i == 0 else C_THY)
            b.set_edgecolor("black")
            b.set_linewidth(0.5)
            b.set_alpha(0.55)
        for a, b in zip(lv, th):
            axv.plot([0, 1], [a, b], lw=0.5, color="#9A9A9A", zorder=3)
        axv.scatter(np.zeros(len(lv)), lv, s=12, facecolor=C_LIV,
                    edgecolor="black", linewidths=0.4, zorder=4)
        axv.scatter(np.ones(len(th)), th, s=12, facecolor=C_THY,
                    edgecolor="black", linewidths=0.4, zorder=4)
        axv.set_xticks([0, 1])
        axv.set_xticklabels(["Liver", "Thymus"], rotation=90, fontsize=8)
        axv.set_ylabel("Fold over mean", fontsize=8, labelpad=1.0)
        axv.set_xlim(-0.6, 1.6)
        axv.set_ylim(0, max(lv.max(), th.max()) * 1.30)
        axv.text(0.02, 0.99, f"p = {p:.3f}", transform=axv.transAxes, fontsize=8,
                 ha="left", va="top", color="#5A5A5A")
        axv.spines[["top", "right"]].set_visible(False)
        F.text(xb + 0.190, 0.058, f"{mk}, $\\pm$10 kb around the summit",
               ha="center", va="top", fontsize=8)
    LT(F, 0.012, 0.300, "C")
    return res


def build_S3(S):
    plt.close("all")
    F = plt.figure(figsize=(S3_W, S3_H))
    MG = panel_G3(F, S["PK7"], S["BW_prot"], S["SUM_prot"], S["BW_hist9"])
    cols, C = panel_H3(F, S["tc"], S["targets"])
    RI = panel_I3(F, S["PK7"], S["BW_lt"], S["SUM_lt"])
    return F, dict(G=MG, H_cols=cols, H_rows=list(C.index), I=RI)


def export(F, name, out=OUT):
    from PIL import Image as _PILI
    from matplotlib.transforms import Bbox
    _PILI.MAX_IMAGE_PIXELS = None
    os.makedirs(out, exist_ok=True)
    # explicit bbox: a 'tight' savefig.bbox silently changes the sheet size
    w_in, h_in = F.get_size_inches()
    kw = dict(bbox_inches=Bbox([[0, 0], [w_in, h_in]]), pad_inches=0.0,
              facecolor="white")
    F.savefig(f"{out}/{name}_PLOS.png", dpi=600, **kw)
    F.savefig(f"{out}/{name}_300.png", dpi=300, **kw)
    F.savefig(f"{out}/{name}_PLOS.svg", **kw)
    _im = _PILI.open(f"{out}/{name}_300.png").convert("RGB")
    _im.save(f"{out}/{name}_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
