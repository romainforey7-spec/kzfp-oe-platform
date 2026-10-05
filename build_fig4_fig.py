
"""Figure 4 (ZNF498) assembly for PLOS Genetics (7.48 x 8.70 in, text >= 8 pt).

A volcano · B alluvial · C target expression heatmap · D corrected GO
enrichment of the 34 genuine targets, with the microtubule block shown under
two tests · E cilia phenotype (reserved for the author's micrographs).
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import upfig_panels as UP

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
import rnaseq_panels as RP
RP.FS.update(star=8.0, leg=8.0, legs=8.0, node=8.0)

OUT = "out/fig4"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 8.70
KZ = "ZNF498"
CALL4 = ("STMN3", "DYRK1A", "PAFAH1B1", "TAF1")
SHORT_CT = {"Oligodendrocyte precursor cells": "OPC",
            "Pancreatic endocrine cells": "Pancr. endocrine",
            "Rod photoreceptor cells": "Rod photorec.",
            "Cone photoreceptor cells": "Cone photorec.",
            "Excitatory neurons": "Excit. neurons",
            "Inhibitory neurons": "Inhib. neurons",
            "Langerhans cells": "Langerhans",
            "granulocytes": "Granulocytes"}
C_MT = "#B00000"


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


ABBR = [("Regulation of Transcription by RNA Polymerase II", "Reg. of transcr., RNA Pol II"),
        ("Negative Regulation of Small GTPase Mediated Signal Transduction",
         "Neg. reg. small GTPase sig."),
        ("Peptidyl-Serine Modification", "Peptidyl-Ser modification"),
        ("Vesicle Cytoskeletal Trafficking", "Vesicle cytoskel. traffick."),
        ("Reg. of microtubule cytoskeleton org.", "Reg. of MT cytoskel. org."),
        ("Regulation of microtubule cytoskeleton organization", "Reg. of MT cytoskel. org."),
        ("Cytoskeleton organization", "Cytoskel. organization"),
        ("Microtubule cytoskeleton", "Microtubule cytoskel."),
        ("Positive Regulation of", "Pos. reg. of"),
        ("Negative Regulation of", "Neg. reg. of"),
        ("Regulation of", "Reg. of"),
        ("RNA Polymerase II", "RNA Pol II"),
        ("DNA-templated Transcription", "DNA-templated transcr."),
        ("Small GTPase Mediated Signal Transduction", "small GTPase signalling"),
        ("Peptidyl-Serine", "Peptidyl-serine"),
        ("Vesicle Cytoskeletal Trafficking", "Vesicle cytoskel. trafficking"),
        ("Tubulin Binding", "Tubulin binding"),
        ("Transcription by", "transcription by")]


def _short(s):
    s = s.replace(" (GO:", "|").split("|")[0]
    for a, b in ABBR:
        s = s.replace(a, b)
    return s


def panel_D(F, terms, mt, rect=(0.150, 0.330, 0.330, 0.152),
            mrect=(0.615, 0.330, 0.300, 0.152)):
    """Left: semantically simplified GO result. Right: the microtubule shift only.

    `terms` is the collapsed table from Fig4D_terms_semantic_simplified.csv
    (9 Enrichr terms merged to 5 at Jaccard >= 0.50 on their gene sets).
    `mt` keeps only the expression-shift rows; the overrepresentation bars are
    dropped because no microtubule term survives that test.
    """
    ax = F.add_axes(list(rect))
    t = terms.sort_values("FDR", ascending=False).reset_index(drop=True)
    t["lab"] = t.representative.map(_short)
    y = np.arange(len(t))
    x = -np.log10(t.FDR.values)
    sz = 16 + t.k.values * 8.0
    ax.scatter(x, y, s=sz, facecolor="#2B2B2B", edgecolor="black", linewidths=0.4,
               zorder=3)
    ax.axvline(-np.log10(0.05), color="#B00000", lw=0.8, ls="--", zorder=2)
    for i, r_ in t.iterrows():
        if r_.n_collapsed > 1:
            ax.text(x[i] + 0.055, i, f"+{int(r_.n_collapsed)-1}", va="center",
                    ha="left", fontsize=8, color="#777777")
    ax.set_yticks(y)
    ax.set_yticklabels(t.lab, fontsize=8)
    ax.tick_params(axis="y", length=0, pad=2.0)
    ax.set_ylim(-0.62, len(t) - 0.38)
    ax.set_xlim(0.9, max(x.max() * 1.20, 2.25))
    ax.text(-np.log10(0.05), -0.55, " FDR 0.05", color="#B00000", fontsize=8,
            ha="left", va="bottom")
    ax.set_xlabel("-Log$_{10}$(FDR)", labelpad=1.0)
    ax.spines[["top", "right"]].set_visible(False)
    hs = [plt.Line2D([], [], marker="o", linestyle="none", color="#2B2B2B",
                     markeredgecolor="black", markeredgewidth=0.4,
                     markersize=np.sqrt(16 + k * 8.0) * 0.52) for k in (3, 8, 14)]
    ax.legend(hs, ["3", "8", "14"], title="Genes", frameon=False,
              loc="lower right", bbox_to_anchor=(0.99, 0.015), handlelength=1.0,
              handletextpad=0.5, labelspacing=0.40, fontsize=8, title_fontsize=8,
              alignment="left", borderaxespad=0.0)

    axm = F.add_axes(list(mrect))
    m = mt[mt.kind != "targets"].iloc[::-1].reset_index(drop=True)
    ym = np.arange(len(m))
    axm.barh(ym, -np.log10(m.fdr.values), height=0.58, color=C_MT,
             edgecolor="black", lw=0.4)
    axm.axvline(-np.log10(0.05), color="#B00000", lw=0.8, ls="--")
    axm.set_yticks(ym)
    axm.set_yticklabels([_short(s) for s in m.lab], fontsize=8)
    axm.tick_params(axis="y", length=0, pad=2.0)
    for tl in axm.get_yticklabels():
        tl.set_color(C_MT)
    for i, r_ in m.iterrows():
        axm.text(-np.log10(r_.fdr) + 0.04, i, f"  z = {r_.stat:.2f}", va="center",
                 ha="left", fontsize=8)
    axm.set_xlabel("-Log$_{10}$(FDR)", labelpad=1.0)
    axm.set_xlim(0, 3.05)
    axm.set_title("Expression shift", fontsize=8, pad=2.5)
    axm.spines[["top", "right"]].set_visible(False)

    return t


def panel_E(F, imgs=None, cil=None, cilst=None, marks=None, order=("GFP", "ZNF498"),
            img_rects=((0.085, 0.058, 0.272, 0.234), (0.375, 0.058, 0.272, 0.234)),
            bar_rect=(0.718, 0.082, 0.142, 0.200)):
    """hTERT-RPE1 ciliation: representative fields and the quantification.

    imgs  : (ciliated_rgb, not_ciliated_rgb) arrays, AcTub(K40) cyan over DAPI grey
    cil   : tidy frame Experiment / Condition / Ciliated_pct, 3 experiments
    cilst : paired-t results per contrast, from Fig4E_ciliation_stats.csv
    marks : list of (x, y) in image pixel coordinates to arrow as cilia in the
            first field; detected as elongated AcTub-positive objects, so the
            placement is provisional until the author confirms it
    order : conditions shown in the bar chart. ZNF18 is omitted per the author's
            annotation 11 on Figure_4_PLOS.png, "Remove ZNF18 from this graph"; it was
            measured in the same three experiments and also reduces ciliation
            (17.7%, 2.36-fold, paired t P = 0.042), which the legend and
            Fig4E_ciliation_stats.csv both report.
    """
    for k, (rc, lab) in enumerate(zip(img_rects, ("Ciliated", "Not ciliated"))):
        ax = F.add_axes(list(rc))
        if imgs is None:
            ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#F4F4F4",
                                   edgecolor="#BBBBBB", lw=0.6))
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
        else:
            im = imgs[k]
            ax.imshow(im, interpolation="bilinear")
            if k == 0 and marks:
                h, w = im.shape[:2]
                for (mx, my) in marks:
                    dx, dy = (70, -70) if mx < w * 0.75 else (-70, -70)
                    ax.annotate("", xy=(mx, my), xytext=(mx + dx, my + dy),
                                arrowprops=dict(arrowstyle="-|>", lw=0.9,
                                                color="#FFD200", shrinkA=0,
                                                shrinkB=3, mutation_scale=5.0))
            sb = 0.26 * im.shape[1]
            ax.plot([im.shape[1] * 0.70, im.shape[1] * 0.70 + sb],
                    [im.shape[0] * 0.945] * 2, lw=2.0, color="white",
                    solid_capstyle="butt")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
            sp.set_edgecolor("#333333")
        ax.set_title(lab, fontsize=8, pad=2.5)
    F.text(np.mean([img_rects[0][0], img_rects[1][0] + img_rects[1][2]]),
           img_rects[0][1] - 0.010,
           "AcTub(K40) cyan / DNA grey" + ("; arrows, cilia" if marks else ""),
           fontsize=8, ha="center", va="top", color="#333333")

    axb = F.add_axes(list(bar_rect))
    if cil is None:
        axb.add_patch(Rectangle((0, 0), 1, 1, facecolor="#F4F4F4",
                                edgecolor="#BBBBBB", lw=0.6))
        axb.set_xlim(0, 1)
        axb.set_ylim(0, 1)
        axb.set_xticks([]); axb.set_yticks([])
        axb.axis("off")
    else:
        order = list(order)
        col = {"GFP": "#B3B3B3", "ZNF18": "#1F4E9C", "ZNF498": "#E07B39"}
        mk = dict(zip(cilst.contrast, cilst["mark"]))
        mn = [cil.loc[cil.Condition == c, "Ciliated_pct"].mean() for c in order]
        sd = [cil.loc[cil.Condition == c, "Ciliated_pct"].std(ddof=1) for c in order]
        x = np.arange(len(order))
        axb.bar(x, mn, yerr=sd, width=0.62, color=[col[c] for c in order],
                edgecolor="black", linewidth=0.5,
                error_kw=dict(lw=0.6, capsize=1.6, capthick=0.6))
        for i, c in enumerate(order):
            v = cil.loc[cil.Condition == c, "Ciliated_pct"].values
            axb.plot(np.full(len(v), i) + np.linspace(-0.13, 0.13, len(v)), v, "o",
                     ms=2.0, mfc="white", mec="black", mew=0.5, zorder=5,
                     linestyle="none")
        top = max(np.array(mn) + np.array(sd))
        for i, c in enumerate(order[1:], start=1):
            m = mk.get(f"GFP vs {c}", "")
            yy = top * (1.08 + 0.17 * (i - 1))
            axb.plot([0, i], [yy, yy], lw=0.6, color="black")
            axb.text(i / 2, yy, m, ha="center", va="bottom", fontsize=8)
        axb.set_ylim(0, top * 1.52)
        axb.set_xticks(x)
        axb.set_xticklabels(order, fontsize=8)
        axb.set_ylabel("Ciliated (%)", fontsize=8, labelpad=1.0)
        axb.tick_params(labelsize=8, length=2.4, width=0.8, pad=1.5)
        axb.spines[["top", "right"]].set_visible(False)
        axb.text(0.5, 1.03, "n = 3", transform=axb.transAxes, fontsize=8,
                 ha="center", va="bottom", color="#555555")


def _flatten_alluvial_labels(ax):
    """De-rotate the node labels so they are legible (author comments 1 and 4)."""
    import matplotlib.patheffects as pe
    for t in ax.texts:
        if t.get_rotation() == 90:
            t.set_rotation(0)
            t.set_ha("center")
            t.set_va("center")
            t.set_fontsize(8)
            t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
            t.set_zorder(20)


def build(S, terms, mt, imgs=None, cil=None, cilst=None, ymaxA=20.0, marks=None):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    d = S[KZ]
    stA = UP.volcano(F, [0.085, 0.800, 0.218, 0.145], [0.318, 0.800, 0.052, 0.145],
                     d["G"], KZ, 1, ymaxA, "UP", (0.384, 0.950), xlim=(-8, 8))
    LT(F, 0.012, 0.996, "A")
    axB = UP.alluvial(F, [0.085, 0.588, 0.250, 0.150], d["S"], (0.398, 0.700))
    _flatten_alluvial_labels(axB)
    LT(F, 0.012, 0.764, "B")
    geom = dict(bar=[0.672, 0.938, 0.196, 0.024],
                hm=[0.672, 0.812, 0.196, 0.120],
                box=[0.672, 0.686, 0.196, 0.052],
                cbar=[0.628, 0.828, 0.009, 0.070])
    cols, C = UP.target_heatmap(F, S["sc"], d["targets"], KZ, geom,
                                call_genes=CALL4, short=SHORT_CT, call_x=0.5)
    LT(F, 0.624, 0.996, "C")
    t = panel_D(F, terms, mt, rect=(0.240, 0.358, 0.250, 0.178),
                mrect=(0.680, 0.358, 0.238, 0.178))
    LT(F, 0.012, 0.552, "D")
    panel_E(F, imgs, cil, cilst, marks=marks)
    LT(F, 0.012, 0.316, "E")
    return F, stA, cols, C, t


def export(F, out=OUT):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    # bbox_inches=None is load-bearing: a 'tight' bbox (the matplotlib default here)
    # crops to the artists and silently changes the sheet size, which breaks the
    # PLOS width limit of 7.5 in.
    from matplotlib.transforms import Bbox
    kw = dict(bbox_inches=Bbox([[0, 0], [SHEET_W, SHEET_H]]), pad_inches=0.0,
              facecolor="white")
    F.savefig(f"{out}/Figure_4_PLOS.png", dpi=600, **kw)
    F.savefig(f"{out}/Figure_4_300.png", dpi=300, **kw)
    F.savefig(f"{out}/Figure_4_PLOS.svg", **kw)
    _im = _PILI.open(f"{out}/Figure_4_300.png").convert("RGB")
    _im.save(f"{out}/Figure_4_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
