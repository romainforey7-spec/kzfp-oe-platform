
"""Figure 5 (ZNF18) assembly for PLOS Genetics.

Three panels only - A volcano, B alluvial, C target expression heatmap - so the
sheet is 7.48 x 5.10 in rather than a full page: the published full-page
version left the bottom half empty, and PLOS sets no minimum height.
"""
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import upfig_panels as UP
import rnaseq_panels as RP

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})
RP.FS.update(star=8.0, leg=8.0, legs=8.0, node=8.0)

OUT = "out/fig5"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 5.10
KZ = "ZNF18"
CALL5 = ("SPATA33", "PCDH1", "TENT5C", "DYNC2I2", "BANF1")
SHORT_CT = {"Late spermatids": "Late spermatids",
            "Early spermatids": "Early spermatids",
            "Rod photoreceptor cells": "Rod photorec.",
            "Pancreatic endocrine cells": "Pancr. endocrine",
            "Cholangiocytes": "Cholangiocytes",
            "Astrocytes": "Astrocytes"}


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def build(S, ymaxA=20.0):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    d = S[KZ]
    stA = UP.volcano(F, [0.085, 0.730, 0.200, 0.165], [0.316, 0.730, 0.050, 0.165],
                     d["G"], KZ, 1, ymaxA, "UP", (0.388, 0.918), xlim=(-8, 8))
    LT(F, 0.012, 0.992, "A")
    UP.alluvial(F, [0.085, 0.380, 0.250, 0.175], d["S"], (0.388, 0.556))
    LT(F, 0.012, 0.620, "B")
    geom = dict(bar=[0.640, 0.845, 0.205, 0.040],
                hm=[0.640, 0.420, 0.205, 0.410],
                box=[0.640, 0.282, 0.205, 0.105],
                cbar=[0.598, 0.580, 0.009, 0.075])
    cols, C = UP.target_heatmap(F, S["sc"], d["targets"], KZ, geom,
                                call_genes=CALL5, short=SHORT_CT)
    LT(F, 0.500, 0.992, "C")
    return F, stA, cols, C


def export(F, out=OUT):
    from PIL import Image as _PILI
    from matplotlib.transforms import Bbox
    _PILI.MAX_IMAGE_PIXELS = None
    # explicit bbox: a 'tight' savefig.bbox silently changes the sheet size
    kw = dict(bbox_inches=Bbox([[0, 0], [SHEET_W, SHEET_H]]), pad_inches=0.0,
              facecolor="white")
    F.savefig(f"{out}/Figure_5_PLOS.png", dpi=600, **kw)
    F.savefig(f"{out}/Figure_5_300.png", dpi=300, **kw)
    F.savefig(f"{out}/Figure_5_PLOS.svg", **kw)
    _im = _PILI.open(f"{out}/Figure_5_300.png").convert("RGB")
    _im.save(f"{out}/Figure_5_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
