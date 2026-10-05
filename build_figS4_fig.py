
"""S4 Fig - the ZNF257 supplement (published Fig. S3, panels A-E).

A anti-HA blot (carry-over) · B viability · C K562 expression density ·
D the four KZFPs across 1,206 cell lines · E the 70 ZNF257 targets across the
HPA single-cell types, ranked by ZNF257, with the MAGEA rows called out.

Published panels F-H (MAGE positions on chrX, the ZNF257 metagene over MAGE
gene bodies, and four published PWMs for the GAGGCA motif) go to S5 Fig: at
8 pt they cannot share a sheet with E, whose column axis alone needs 4.4 in.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.cluster.hierarchy import linkage, dendrogram
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

OUT = "out/figS4"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70
KZ = "ZNF257"
N_CT = 40                       # highest-ZNF257 cell types that fit at 8 pt
CMAP_HM = LinearSegmentedColormap.from_list(
    "kzfp", ["#0000CD", "#8A7DD1", "#C9C3E3", "#E8E8E8", "#F0C0B4", "#E06B52", "#CC1F0A"])
C_MAG = "#9B2FAE"
SHORT_CT = {
    "Oligodendrocyte precursor cells": "OPC", "Distal tubular cells": "Distal tubular",
    "Proximal tubular cells": "Proximal tubular", "Collecting duct cells": "Collecting duct",
    "Extravillous trophoblasts": "Extravillous troph.",
    "Cytotrophoblasts": "Cytotrophoblasts", "Syncytiotrophoblasts": "Syncytiotroph.",
    "Rod photoreceptor cells": "Rod photorec.", "Cone photoreceptor cells": "Cone photorec.",
    "Pancreatic endocrine cells": "Pancr. endocrine",
    "Endometrial stromal cells": "Endometrial stromal",
    "Ovarian stromal cells": "Ovarian stromal",
    "Basal respiratory cells": "Basal respiratory",
    "Basal squamous epithelial cells": "Basal squamous ep.",
    "Squamous epithelial cells": "Squamous ep.",
    "Glandular and luminal cells": "Glandular/luminal",
    "Undifferentiated cells": "Undifferentiated",
    "Exocrine glandular cells": "Exocrine glandular",
    "Serous glandular cells": "Serous glandular",
    "Mucus glandular cells": "Mucus glandular",
    "Breast myoepithelial cells": "Breast myoepith.",
    "Breast glandular cells": "Breast glandular",
    "Prostatic glandular cells": "Prostatic glandular",
    "Gastric mucus-secreting cells": "Gastric mucus-secr.",
    "Intestinal goblet cells": "Intestinal goblet",
    "Alveolar cells type 1": "Alveolar type 1",
    "Alveolar cells type 2": "Alveolar type 2",
    "Lymphatic endothelial cells": "Lymphatic endoth.",
    "Endothelial cells": "Endothelial",
    "Smooth muscle cells": "Smooth muscle",
    "Peritubular cells": "Peritubular",
    "Sertoli cells": "Sertoli", "Leydig cells": "Leydig",
    "Granulosa cells": "Granulosa", "Hofbacher cells": "Hofbacher",
    "Horizontal cells": "Horizontal", "Bipolar cells": "Bipolar",
    "Muller glia cells": "Muller glia", "Microglial cells": "Microglial",
    "dendritic cells": "Dendritic", "Kupffer cells": "Kupffer",
}


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def panel_E(F, sc, targets, kz=KZ, vmax=3, n_ct=N_CT):
    r = sc[sc["Gene name"] == kz][["Cell type", "nTPM"]].copy()
    r = r.sort_values("nTPM", ascending=False)
    cols = r["Cell type"].tolist()[:n_ct]
    barv = dict(zip(r["Cell type"], r.nTPM))
    m = sc[sc["Gene name"].isin(targets) & sc["Cell type"].isin(cols)]
    W = m.pivot_table(index="Gene name", columns="Cell type", values="nTPM",
                      aggfunc="mean")
    L = np.log(W.reindex(columns=cols).dropna(how="all") + 1)
    C = L.sub(L.mean(axis=1), axis=0)
    Z = linkage(C.values, method="complete", metric="euclidean")
    C = C.iloc[dendrogram(Z, no_plot=True)["leaves"]]
    x0, w = 0.145, 0.690
    axbar = F.add_axes([x0, 0.520, w, 0.042])
    axhm = F.add_axes([x0, 0.300, w, 0.214])
    axbox = F.add_axes([x0, 0.188, w, 0.104])
    axbar.bar(range(len(cols)), [barv[c] for c in cols], width=0.74,
              color="#4D4D4D", edgecolor="none")
    axbar.set_xlim(-0.5, len(cols) - 0.5)
    axbar.set_xticks([])
    axbar.set_ylabel(f"${kz}$\n(nTPM)", fontsize=8, linespacing=1.1, labelpad=1.5)
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
    mag = sorted([(i, g) for i, g in enumerate(C.index) if g.startswith("MAGE")])
    n = len(C)
    ys = np.linspace(n * 0.06, n * 0.94, len(mag)) if len(mag) > 1 else [n * 0.5]
    for (i, g), yy in zip(mag, ys):
        axhm.annotate(g, xy=(len(cols) - 0.4, i), xytext=(len(cols) + 1.2, yy),
                      fontsize=8, style="italic", color=C_MAG, va="center",
                      ha="left", annotation_clip=False,
                      arrowprops=dict(arrowstyle="-", lw=0.6, color=C_MAG,
                                      shrinkA=0, shrinkB=1))
    cax = F.add_axes([x0 - 0.064, 0.330, 0.009, 0.060])
    cb = F.colorbar(plt.cm.ScalarMappable(norm=TwoSlopeNorm(0, -vmax, vmax),
                                          cmap=CMAP_HM), cax=cax,
                    ticks=[-vmax, 0, vmax])
    cb.ax.tick_params(labelsize=8, length=1.8, pad=1.5)
    cb.outline.set_linewidth(0.5)
    cb.ax.yaxis.set_ticks_position("left")
    cax.set_ylabel("Centred\nlog(nTPM+1)", fontsize=8, labelpad=2, linespacing=1.1)
    cax.yaxis.set_label_position("left")
    data = [L[c].dropna().values for c in cols]
    axbox.boxplot(data, positions=range(len(cols)), widths=0.62, patch_artist=True,
                  showfliers=False, whis=0, showcaps=False,
                  medianprops=dict(color="black", lw=0.6),
                  boxprops=dict(facecolor="#D9D9D9", edgecolor="#6E6E6E", lw=0.4),
                  whiskerprops=dict(lw=0))
    for i, v in enumerate(data):
        axbox.plot([i], [np.mean(v)], marker="_", color="red", ms=3.6, mew=0.8)
    axbox.set_xlim(-0.5, len(cols) - 0.5)
    axbox.set_xticks(range(len(cols)))
    lbl = [SHORT_CT.get(c, c) for c in cols]
    axbox.set_xticklabels(lbl, rotation=90, fontsize=8)
    for tl, c in zip(axbox.get_xticklabels(), cols):
        if c in ("Spermatocytes", "Spermatogonia", "Early spermatids",
                 "Late spermatids"):
            tl.set_color(C_MAG)
            tl.set_fontweight("bold")
    axbox.tick_params(axis="x", pad=1.0, length=1.6)
    axbox.set_ylabel("Targets\nlog(nTPM+1)", fontsize=8, linespacing=1.1, labelpad=1.5)
    axbox.spines[["top", "right"]].set_visible(False)
    LT(F, 0.012, 0.600, "E")
    return cols, C


from presto_screen import series_for
VIAB_257 = None   # set by build(); falls back to series_for("ZNF257")


def build(S, SS, blot):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    ax = F.add_axes([0.085, 0.838, 0.250, 0.118])
    ax.imshow(blot)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    LT(F, 0.012, 0.986, "A")

    axb = F.add_axes([0.520, 0.846, 0.240, 0.100])
    R = VIAB_257 if VIAB_257 is not None else series_for("ZNF257")[0]
    axb.errorbar([0] + list(R.day), [1.0] + list(R["mean"]),
                 yerr=[0.0] + list(R["sd"].fillna(0.0) if "sd" in R
                                     else R.get("err", pd.Series(0.0, index=R.index))),
                 marker="o", ms=4, lw=1.1,
                 color="#B00000", ecolor="#B00000", elinewidth=0.8,
                 capsize=2.0, capthick=0.8, zorder=3)
    axb.axhline(0.85, color="#808080", lw=0.7, ls="--")
    axb.text(9.5, 0.86, "toxicity cut-off", ha="right", va="bottom", fontsize=8,
             color="#808080")
    axb.set_xticks([0, 4, 7, 9])
    axb.set_xlabel("Days after Dox induction", labelpad=1.0)
    axb.set_ylabel(f"Relative viability\n({KZ} / GFP)", linespacing=1.15)
    axb.set_ylim(0, 1.15)
    axb.set_xlim(-0.6, 9.8)
    axb.text(0.0, 1.03, "one plate; the technical SD over the three Dox and three\n"
             "No-Dox wells is smaller than the marker", transform=axb.transAxes,
             ha="left", va="bottom", fontsize=8, color="#6E6E6E",
             linespacing=1.15)
    axb.spines[["top", "right"]].set_visible(False)
    LT(F, 0.440, 0.986, "B")

    import build_figS2_fig as BS
    BS.panel_C(F, SS["ALLG"], SS["KZF"], SS["SEL"], gene=KZ,
               rect=[0.105, 0.672, 0.280, 0.110], letter=(0.012, 0.818))
    BS.panel_D(F, SS["CL"], SS["SEL"], gene=KZ, header=False,
               rect=[0.560, 0.672, 0.280, 0.110], letter=(0.440, 0.818))
    cols, C = panel_E(F, S["sc"], S["t257"])
    return F, cols, C


def export(F, out=OUT, name="Figure_S4"):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F.savefig(f"{out}/{name}_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/{name}_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/{name}_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/{name}_300.png").convert("RGB")
    _im.save(f"{out}/{name}_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
