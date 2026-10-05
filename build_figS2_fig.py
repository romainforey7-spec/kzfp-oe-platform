"""Panel library for the KZFP supplement sheets.

This module is no longer a sheet builder. It supplies the four panels that every
KZFP supplement shares, each taking an explicit rect and letter position:

  panel_C  K562 expression density, endogenous and transgene marked
  panel_D  the KZFP across 1,206 cell lines, ranked
  panel_E  TSS-to-peak distance by differential-expression class, with a rug
  panel_F  repeat-class composition and the LTR/ERV1 expected-vs-observed bar

The sheets that call them are build_figS2S3_fig.py (S2 Fig and S3 Fig, the split
of the published two-page Figure S2), build_figS4_fig.py (ZNF257) and
build_kzfp_sup.py (ZNF498 and ZNF18).
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats as st
from plos_export import PK  # pins output to the figure's declared size

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
SHEET_W, SHEET_H = 7.48, 8.70
LX, RX = 0.090, 0.555
C_ERV = "#1F6FB4"
C_ALU = "#9A9A9A"
C_TGT = "#B00000"
C_LIV, C_THY = "#C2743A", "#4F79B8"
CMAP_SIG = LinearSegmentedColormap.from_list(
    "sig", ["#FFFFFF", "#D9E6F2", "#8FB8DC", "#3D7FB8", "#13406E"])
TRK_G = [("ZNF43-HA", "#2B2B2B"), ("KAP1", "#1F4E9C"), ("SETDB1", "#2E7D32"),
         ("H3K9me3", "#5E35B1"), ("H3K4me3", "#2E7D32"), ("H3K27ac", "#C77C39")]


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


# ------------------------------------------------------------------ A
def panel_A(F, img):
    ax = F.add_axes([LX, 0.796, 0.283, 0.150])
    if img is None:
        ax.text(0.5, 0.5, "anti-HA western blot\n(carry-over panel)",
                ha="center", va="center", fontsize=8, color="#808080")
    else:
        ax.imshow(img)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.012, 0.968, "A")
    return ax


# ------------------------------------------------------------------ B
def panel_B(F, viab):
    ax = F.add_axes([LX + 0.030, 0.655, 0.230, 0.098])
    days = [0, 4, 7, 9]
    vals = [viab[f"D{d}"] if f"D{d}" in viab else np.nan for d in days]
    ax.plot(days, vals, marker="o", ms=4, lw=1.1, color=C_TGT)
    ax.axhline(0.85, color="#808080", lw=0.7, ls="--")
    ax.text(9, 0.87, "toxicity cut-off", ha="right", va="bottom", fontsize=8,
            color="#808080")
    ax.set_xticks(days)
    ax.set_xlabel("Days after Dox induction", labelpad=1.0)
    ax.set_ylabel("Relative viability\n(ZNF43 / GFP)", linespacing=1.15)
    ax.set_ylim(0, 1.15)
    ax.set_xlim(-0.6, 9.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.5, 1.02, "one induction series, no replicate error bars",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=8,
            color="#8A8A8A")
    LT(F, 0.012, 0.775, "B")


# ------------------------------------------------------------------ C
def panel_C(F, ALLG, KZF, SEL, gene="ZNF43", rect=None, letter=(0.012, 0.625)):
    ax = F.add_axes(rect or [LX + 0.030, 0.505, 0.230, 0.098])
    va = pd.to_numeric(ALLG["K562_GFP_local"], errors="coerce").dropna()
    vk = pd.to_numeric(KZF["K562_GFP_local"], errors="coerce").dropna()
    for v, col, lab in [(va, "#BFBFBF", f"all genes (n = {len(va):,})"),
                        (vk, "#6E8FCB", f"KZFPs (n = {len(vk):,})")]:
        x = np.log10(v.clip(lower=0) + 1)
        kd = st.gaussian_kde(x, bw_method=0.35)
        xs = np.linspace(0, x.max() * 1.05, 400)
        ax.fill_between(xs, 0, kd(xs), color=col, alpha=0.75, lw=0, label=lab)
    s = SEL[SEL.symbol == gene]
    for _, r in s.iterrows():
        xv = np.log10(max(float(r.value), 0) + 1)
        col = C_TGT if "transgene" in str(r.type).lower() else "#1F4E9C"
        ax.axvline(xv, color=col, lw=1.0)
        lab = str(r.type).lower()
        ax.text(xv, 0.86, lab, rotation=0,
                ha="right" if lab.startswith("endo") else "left",
                va="top", fontsize=8, color=col,
                transform=ax.get_xaxis_transform(), bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=6)
    ax.set_xlabel("log$_{10}$(normalised expression + 1), K562", labelpad=1.0)
    ax.set_ylabel("Density")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, 1.16),
              ncol=2, handlelength=1.0, handleheight=1.0, handletextpad=0.4,
              columnspacing=1.0, labelspacing=0.22, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F, letter[0], letter[1], "C")
    return len(va), len(vk)


# ------------------------------------------------------------------ D
def panel_D(F, CL, SEL, gene="ZNF43", rect=None, letter=(0.012, 0.475),
            header=True):
    ax = F.add_axes(rect or [LX + 0.030, 0.355, 0.230, 0.098])
    d = CL[CL.symbol == gene].copy()
    col = "expression_norm_all_genes" if "expression_norm_all_genes" in d else "expression_raw"
    y = np.log10(pd.to_numeric(d[col], errors="coerce").clip(lower=0) + 1).dropna()
    y = y.sort_values(ascending=False).values
    ax.plot(np.arange(len(y)), y, lw=0.9, color="#4D4D4D")
    tg = SEL[(SEL.symbol == gene) & SEL.type.astype(str).str.contains("transgene", case=False)]
    if len(tg):
        yv = np.log10(max(float(tg.value.iloc[0]), 0) + 1)
        ax.axhline(yv, color=C_TGT, lw=1.0, ls="--")
        ax.text(len(y) * 0.98, yv, "transgene level", ha="right", va="top",
                fontsize=8, color=C_TGT, bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=6)
    ax.set_xlabel(f"{len(y):,} DepMap cell lines, ranked by ${gene}$",
                  labelpad=1.0)
    ax.set_ylabel(f"log$_{{10}}$(${gene}$\n+ 1)", labelpad=1.0,
                  linespacing=1.15)
    if header:
        ax.text(0.0, 1.04, f"Endogenous ${gene}$ never reaches the induced level",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=8,
                color="#6E6E6E")
    ax.set_xlim(0, len(y))
    ax.spines[["top", "right"]].set_visible(False)
    LT(F, letter[0], letter[1], "D")
    return len(y)


# ------------------------------------------------------------------ E
DE_COL = {"DOWN": "#377EB8", "UNCHANGED": "#9A9A9A", "UP": "#B00000"}


def panel_E(F, g, rect=None, letter=(0.012, 0.325)):
    """Peak-to-TSS distance by DE class, with a rug, as in the published panel."""
    ax = F.add_axes(rect or [LX + 0.030, 0.205, 0.230, 0.098])
    for k, lab in (("UNCHANGED", "Unch."), ("DOWN", "Down"), ("UP", "Up")):
        v = g.loc[g.cat == ("NS" if k == "UNCHANGED" else k), "dist"].dropna().values
        if len(v) < 5:
            continue
        x = np.log10(np.asarray(v, float) + 1)
        kd = st.gaussian_kde(x, bw_method=0.30)
        xs = np.linspace(0, np.log10(5.2e4), 400)
        ax.plot(xs, kd(xs), lw=1.2, color=DE_COL[k], label=f"{lab} ({len(v):,})")
        if k != "UNCHANGED":
            xr = x[x <= np.log10(5.2e4)]      # the rug must stay inside the axis
            ax.plot(xr, np.full(len(xr), -0.035), "|", ms=3.2, mew=0.5,
                    color=DE_COL[k], clip_on=False)
    for kb, lab in [(1e3, "1 kb"), (1e4, "10 kb")]:
        ax.axvline(np.log10(kb + 1), color="#707070", lw=0.6, ls=":")
        ax.text(np.log10(kb + 1), 0.86, lab, rotation=0, ha="center",
                va="top", fontsize=8, color="#707070",
                transform=ax.get_xaxis_transform(), bbox=dict(facecolor="white", edgecolor="none", pad=0.6), zorder=6)
    ax.set_xlabel("log$_{10}$(TSS-to-peak distance + 1) (bp)", labelpad=1.0)
    ax.set_ylabel("Gene density")
    ax.set_xlim(0, np.log10(5.2e4))
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(-0.02, 1.16),
              ncol=3, handlelength=1.2, handletextpad=0.4, columnspacing=1.0,
              labelspacing=0.22, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F, letter[0], letter[1], "E")


# ------------------------------------------------------------------ F
def panel_F(F, counts, n_rows, obs, exp, npk, fold, p, rects=None,
            letter=(0.012, 0.180)):
    rt, rb = rects or ([LX + 0.040, 0.092, 0.200, 0.064],
                       [LX + 0.040, 0.048, 0.200, 0.017])
    axt = F.add_axes(rt)
    axb = F.add_axes(rb)
    cols = [C_ERV if c == "LTR/ERV1" else C_ALU for c in counts.index]
    axt.bar(range(len(counts)), counts.values, width=0.58, color=cols,
            edgecolor="black", lw=0.4)
    for i, v in enumerate(counts.values):
        axt.text(i, v + 0.25, str(int(v)), ha="center", va="bottom", fontsize=8)
    axt.set_xticks(range(len(counts)))
    axt.set_xticklabels(counts.index, fontsize=8)
    axt.set_ylabel("Peak-gene pairs")
    axt.set_ylim(0, max(counts.values) * 1.95)
    axt.text(0.02, 1.04, "Repeat family under the peak, for down-regulated\n"
             f"genes with a peak $\\leq$ 10 kb from the TSS (n = {n_rows} pairs)",
             transform=axt.transAxes, ha="left", va="bottom", fontsize=8,
             linespacing=1.15)
    axt.spines[["top", "right"]].set_visible(False)
    axb.barh([1], [obs], height=0.52, color=C_ERV, edgecolor="black", lw=0.6)
    axb.barh([0], [exp], height=0.52, color="#BFBFBF", edgecolor="black", lw=0.6)
    axb.text(obs * 1.03, 1, f"{obs}", ha="left", va="center", fontsize=8)
    axb.text(exp + obs * 0.03, 0, f"{exp:.0f}", ha="left", va="center", fontsize=8)
    axb.set_yticks([0, 1])
    axb.set_yticklabels(["expected", "observed"], fontsize=8)
    axb.set_xlim(0, obs * 1.45)
    axb.set_xticks([0, 300, 600])
    axb.set_xlabel(f"ZNF43 peaks on an LTR/ERV1 element, observed vs expected\n"
                   f"from the genomic share of LTR/ERV1 (of {npk} repeat-overlapping\n"
                   f"peaks; {fold:.1f}-fold, binomial P < 1e-300)",
                   labelpad=1.0, linespacing=1.15)
    axb.tick_params(axis="y", length=0, pad=1.5)
    axb.spines[["top", "right", "left"]].set_visible(False)
    LT(F, letter[0], letter[1], "F")


# ------------------------------------------------------------------ G
def panel_G(F, PK7, BWP, SUMP, BWH, half=5000, nb=200, flank=(2000, 5000)):
    """Mean profile per track at the seven target peaks, plus summit vs flank.

    The three protein tracks are scaled by their genome mean; the normDMSO
    histone files are normalised ratio tracks (genome means ~0.002, negative
    minima) and are plotted raw.
    """
    xs = np.linspace(-half / 1000, half / 1000, nb)
    GX0, GW, GH = RX + 0.058, 0.168, 0.052
    GROWS = [0.886, 0.800, 0.714]
    prof, summ = {}, []
    for i, (k, col) in enumerate(TRK_G):
        raw = k in BWH
        bw = BWH[k] if raw else BWP[k]
        M = np.vstack([np.nan_to_num(bw.values(r.chrom, r.summit - half,
                                               r.summit + half, nbins=nb,
                                               agg="mean"))
                       for r in PK7.itertuples()])
        if not raw:
            M = M / SUMP[k]
        prof[k] = M
        c = np.abs(xs) <= 0.5                                  # summit +/- 500 bp
        fl = (np.abs(xs) >= flank[0] / 1000) & (np.abs(xs) <= flank[1] / 1000)
        summ.append(dict(track=k, raw=raw, summit=float(M[:, c].mean()),
                         flank=float(M[:, fl].mean()),
                         summit_per_peak=M[:, c].mean(1),
                         flank_per_peak=M[:, fl].mean(1)))
        ax = F.add_axes([GX0 + (i % 2) * (GW + 0.034), GROWS[i // 2], GW, GH])
        ax.plot(xs, M.mean(0), lw=1.0, color=col)
        ax.fill_between(xs, 0, M.mean(0), color=col, alpha=0.25, lw=0)
        ax.axhline(0, color="#999999", lw=0.5)
        ax.set_xlim(xs[0], xs[-1])
        ax.set_xticks([-5, 0, 5])
        ax.set_title(k + ("" if raw else "  (fold)"), fontsize=8, pad=1.5, color=col)
        ax.tick_params(labelsize=8, pad=1.0, length=1.8)
        if i // 2 == 2:
            ax.set_xlabel("kb from summit", labelpad=1.0)
        else:
            ax.set_xticklabels([])
        ax.spines[["top", "right"]].set_visible(False)
    SM = pd.DataFrame(summ)
    # summit vs flank, each of the seven peaks shown as a dot
    axs = F.add_axes([RX + 0.070, 0.566, 0.350, 0.082])
    k = np.arange(len(SM))
    for j, (lab, key, fc) in enumerate([("summit $\\pm$ 0.5 kb", "summit", "#4D4D4D"),
                                        ("flank 2-5 kb", "flank", "#C9C9C9")]):
        dx = -0.19 + j * 0.38
        axs.bar(k + dx, SM[key].values, width=0.34, color=fc, edgecolor="black",
                lw=0.4, label=lab)
        for i in k:
            v = SM[f"{key}_per_peak"].iloc[i]
            axs.plot(np.full(len(v), i + dx) + np.linspace(-0.08, 0.08, len(v)), v,
                     "o", ms=2.0, mfc="white", mec="black", mew=0.35, zorder=5)
    axs.set_xticks(k)
    axs.set_xticklabels(SM.track, rotation=45, ha="right", fontsize=8,
                        rotation_mode="anchor")
    axs.axhline(0, color="#999999", lw=0.5)
    axs.set_ylabel("Mean signal\n(fold / ratio)", linespacing=1.15)
    axs.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.03, 1.14),
               ncol=2, handlelength=1.0, handleheight=1.0, handletextpad=0.4,
               columnspacing=1.0, fontsize=8)
    axs.spines[["top", "right"]].set_visible(False)
    LT(F, 0.500, 0.968, "G")
    return prof, SM


# ------------------------------------------------------------------ H
def panel_H(F, H, mark=("bone marrow", "thymus", "liver")):
    """50 HPA tissues ranked by ZNF43; only the three discussed are labelled."""
    ax = F.add_axes([RX + 0.070, 0.408, 0.340, 0.082])
    k = np.arange(len(H))
    zz = np.log10(H.znf43 + 1).values
    tt = (H.tgt_mean / np.log(10)).values
    ax.plot(k, zz, lw=1.2, color="#4D4D4D", label="$ZNF43$")
    ax.plot(k, tt, lw=1.2, color=C_TGT, label="targets (mean)")
    for t in mark:
        if t not in H.index:
            continue
        i = int(np.where(H.index == t)[0][0])
        ax.plot([i, i], [0, max(zz[i], tt[i])], lw=0.6, ls=":", color="#808080",
                zorder=1)
        ax.scatter([i, i], [zz[i], tt[i]], s=16, zorder=5,
                   facecolor=["#4D4D4D", C_TGT], edgecolor="black", linewidths=0.4)
        j = mark.index(t)
        ax.annotate(t, xy=(i, max(zz[i], tt[i])),
                    xytext=(i + (2.0 if i < len(H) / 2 else -2.0),
                            max(zz.max(), tt.max()) * (1.06 + 0.13 * j)),
                    fontsize=8, ha="left" if i < len(H) / 2 else "right",
                    va="center", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#808080",
                                    shrinkA=1, shrinkB=2))
    ax.set_xlabel(f"{len(H)} HPA consensus tissues, ranked by $ZNF43$", labelpad=1.0)
    ax.set_ylabel("log$_{10}$(nTPM + 1)")
    ax.set_xlim(-0.8, len(H) - 0.2)
    ax.set_ylim(0, max(zz.max(), tt.max()) * 1.52)
    ax.set_xticks([])
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(-0.02, -0.04),
              ncol=2, handlelength=1.2, handletextpad=0.4, columnspacing=1.0,
              fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    LT(F, 0.500, 0.514, "H")


# ------------------------------------------------------------------ I
def panel_I(F, PK7, BWL, SUML, half=10000, nb=200):
    ax = F.add_axes([RX + 0.070, 0.252, 0.150, 0.078])
    axs2 = F.add_axes([RX + 0.288, 0.252, 0.122, 0.078])
    res = {}
    for mk, axq in [("H3K4me3", ax), ("H3K9me3", axs2)]:
        lv, th = [], []
        for r in PK7.itertuples():
            a, b = r.summit - half, r.summit + half
            for tag, store in [("Liver", lv), ("Thymus", th)]:
                k = f"{tag} {mk}"
                v = np.nan_to_num(BWL[k].values(r.chrom, a, b, nbins=nb,
                                                agg="mean")) / SUML[k]
                store.append(float(np.nanmean(v)))
        lv, th = np.array(lv), np.array(th)
        t, p = st.ttest_rel(lv, th)
        res[mk] = dict(liver=lv.mean(), thymus=th.mean(), t=float(t), p=float(p),
                       n=len(lv))
        for i, (v, c, lab) in enumerate([(lv, C_LIV, "Liver"), (th, C_THY, "Thymus")]):
            axq.bar([i], [v.mean()], width=0.55, color=c, edgecolor="black", lw=0.4)
            axq.plot(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v,
                     "o", ms=2.6, mfc="white", mec="black", mew=0.4, zorder=4)
        axq.set_xticks([0, 1])
        axq.set_xticklabels(["Liver", "Thymus"], fontsize=8)
        axq.set_title(mk, fontsize=8, pad=2.0)
        axq.set_ylim(0, max(lv.max(), th.max()) * 1.62)
        axq.text(0.5, 0.99, f"paired $t$-test\np = {p:.3f}  (n = {len(lv)})",
                 transform=axq.transAxes, ha="center", va="top", fontsize=8,
                 color="#5A5A5A", linespacing=1.15)
        axq.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Mean fold over genome mean,\nsummit $\\pm$ 10 kb",
                  fontsize=8, linespacing=1.15)
    LT(F, 0.500, 0.358, "I")
    return res


# ------------------------------------------------------------------ J
def panel_J(F, S):
    """The DNAI4/WDR78 locus - the eighth target, not shown in Figure 2E."""
    import build_fig2_fig as B2
    B2.locus_block(F, "DNAI4", RX - 0.005, 0.030, 0.430, 0.170,
                   S["WIN2"], S["BW"], S["SUM"], S["TRK"], S["RPALL"], S["TB"],
                   "DNAI4 / WDR78")
    LT(F, 0.500, 0.222, "J")


def build(S, blot=None):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    panel_A(F, blot)
    panel_B(F, S["viab"])
    nC = panel_C(F, S["ALLG"], S["KZF"], S["SEL"])
    nD = panel_D(F, S["CL"], S["SEL"])
    panel_E(F, S["G43"])
    panel_F(F, S["dn_class_counts"], S["dn_rows"], S["erv1_obs"], S["erv1_exp"],
            S["n_peaks_rep"], S["erv1_fold"], S["erv1_p"])
    MG, SMG = panel_G(F, S["PK7"], S["BW_prot"], S["SUM_prot"], S["BW_hist"])
    panel_H(F, S["H_rank"])
    RI = panel_I(F, S["PK7"], S["BW_lt"], S["SUM_lt"])
    panel_J(F, S)
    return F, dict(nC=nC, nD=nD, G=MG, SMG=SMG, I=RI)


def export(F, out=OUT):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F.savefig(f"{out}/Figure_S2_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/Figure_S2_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/Figure_S2_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/Figure_S2_300.png").convert("RGB")
    _im.save(f"{out}/Figure_S2_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
