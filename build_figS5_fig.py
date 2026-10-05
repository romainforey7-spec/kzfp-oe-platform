
"""S5 Fig - the MAGE locus panels of the ZNF257 supplement (published F-H).

F  the 34 MAGE genes along chromosome X, with the MAGEA cluster expanded
G  ZNF257 ChIP signal over the MAGE gene bodies, TSS to TES, one row per gene
H  reserved: four published PWMs for the GAGGCA motif. CIS-BP is not on the
   network allowlist, so the matrices cannot be fetched here.
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.colors import LinearSegmentedColormap
from plos_export import PK  # pins output to the figure's declared size

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 8,
    "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "axes.titlesize": 8,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "svg.fonttype": "none", "pdf.fonttype": 42,
})

OUT = "out/figS5"
os.makedirs(OUT, exist_ok=True)
SHEET_W, SHEET_H = 7.48, 6.40
CMAP_SIG = LinearSegmentedColormap.from_list(
    "sig", ["#2B0B4F", "#2E6E8E", "#3FA68C", "#9FD44F", "#F7F75A"])
C_MAG = "#9B2FAE"
REPRESSED = ["MAGEA6", "MAGEA12", "MAGEA3", "MAGEA1"]
ZOOM = (148.0e6, 153.0e6)


def LT(F, x, y, s):
    F.text(x, y, s, fontsize=11, fontweight="bold", va="top")


def _spread(vals, lo, hi, gap):
    o = np.argsort(vals)
    out = np.array(vals, float)
    cur = lo
    for i in o:
        out[i] = max(out[i], cur)
        cur = out[i] + gap
    if out.max() > hi:
        out = out - (out.max() - hi)
    return out


def panel_F(F, MB, chrom_len=155270560):
    """All MAGE genes on chrX, with the 148-153 Mb window expanded."""
    g = MB[MB[0] == "chrX"].copy().sort_values(1).reset_index(drop=True)
    axa = F.add_axes([0.075, 0.120, 0.055, 0.780])
    axb = F.add_axes([0.300, 0.120, 0.055, 0.780])
    for ax, (lo, hi) in ((axa, (0, chrom_len)), (axb, ZOOM)):
        ax.add_patch(Rectangle((0, lo), 1, hi - lo, facecolor="#F2F2F2",
                               edgecolor="#8C8C8C", lw=0.7))
        ax.set_xlim(0, 1)
        ax.set_ylim(lo, hi)
        ax.set_xticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    sel = g[g[1].between(*ZOOM)]
    axa.add_patch(Rectangle((-0.05, ZOOM[0]), 1.10, ZOOM[1] - ZOOM[0],
                            facecolor="#F6C9A8", edgecolor="#C2743A", lw=0.8,
                            zorder=3))
    for ax, sub, lo, hi in ((axa, g[~g[1].between(*ZOOM)], 0, chrom_len),
                            (axb, sel, ZOOM[0], ZOOM[1])):
        ys = _spread(list(sub[1]), lo + (hi - lo) * 0.015,
                     hi - (hi - lo) * 0.015, (hi - lo) * 0.030)
        for (_, r), yy in zip(sub.iterrows(), ys):
            rep = r[3] in REPRESSED
            ax.scatter([0.5], [r[1]], s=14, facecolor=C_MAG if rep else "#4D7FB8",
                       edgecolor="black", linewidths=0.4, zorder=5)
            ax.annotate(r[3], xy=(1.0, r[1]), xytext=(1.5, yy), fontsize=8,
                        style="italic", va="center", ha="left",
                        color=C_MAG if rep else "black",
                        fontweight="bold" if rep else "normal",
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="#9A9A9A",
                                        shrinkA=0, shrinkB=1))
    for yy in ZOOM:
        F.add_artist(ConnectionPatch(xyA=(1.05, yy), coordsA=axa.transData,
                                     xyB=(0, yy), coordsB=axb.transData,
                                     lw=0.6, color="#C2743A", linestyle="--"))
    axa.set_ylabel("Position on chromosome X (Mb, hg19)", fontsize=8, labelpad=2)
    axa.set_yticks(np.arange(0, chrom_len, 50e6))
    axa.set_yticklabels([f"{v/1e6:.0f}" for v in np.arange(0, chrom_len, 50e6)],
                        fontsize=8)
    axa.tick_params(axis="y", length=1.8, pad=1.5)
    axa.spines["left"].set_visible(True)
    axb.set_yticks(np.arange(148e6, 153.1e6, 1e6))
    axb.set_yticklabels([f"{v/1e6:.0f}" for v in np.arange(148e6, 153.1e6, 1e6)],
                        fontsize=8)
    axb.tick_params(axis="y", length=1.8, pad=1.5)
    axb.spines["left"].set_visible(True)
    axb.set_title("148-153 Mb", fontsize=8, pad=2.0, color="#C2743A")
    LT(F, 0.012, 0.968, "A")
    return g, sel


def panel_G(F, bw, bg, MB, nb=120, flank=1000):
    """ZNF257 over the MAGE gene bodies: 1 kb flank, TSS to TES, 1 kb flank."""
    g = MB[MB[0] == "chrX"].copy().sort_values(1).reset_index(drop=True)
    rows, M = [], []
    for _, r in g.iterrows():
        a, b = int(r[1]), int(r[2])
        if b - a < 200:
            continue
        up = bw.values("chrX", a - flank, a, nbins=nb // 4, agg="mean")
        body = bw.values("chrX", a, b, nbins=nb // 2, agg="mean")
        dn = bw.values("chrX", b, b + flank, nbins=nb // 4, agg="mean")
        v = np.nan_to_num(np.concatenate([up, body, dn])) / bg
        if r[5] == "-":
            v = v[::-1]
        rows.append(r[3])
        M.append(v)
    M = np.vstack(M)
    ordr = np.argsort(-M.max(1))
    M, rows = M[ordr], [rows[i] for i in ordr]
    ax = F.add_axes([0.600, 0.120, 0.230, 0.780])
    vmax = float(np.nanpercentile(M, 99.5)) or 1.0
    ax.imshow(M, aspect="auto", cmap=CMAP_SIG, vmin=0, vmax=vmax,
              extent=[0, M.shape[1], len(rows), 0], interpolation="nearest")
    ax.set_yticks(np.arange(len(rows)) + 0.5)
    ax.set_yticklabels(rows, fontsize=8, style="italic")
    for tl, nm in zip(ax.get_yticklabels(), rows):
        if nm in REPRESSED:
            tl.set_color(C_MAG)
            tl.set_fontweight("bold")
    ax.tick_params(axis="y", length=0, pad=1.5)
    ax.set_xticks([nb // 4, nb // 4 + nb // 2])
    ax.set_xticklabels(["TSS", "TES"], fontsize=8)
    ax.tick_params(axis="x", length=1.8, pad=1.0)
    ax.set_xlabel(f"$\\pm${flank/1000:.0f} kb flanks", labelpad=1.0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cax = F.add_axes([0.600, 0.925, 0.230, 0.010])
    cb = F.colorbar(plt.cm.ScalarMappable(cmap=CMAP_SIG), cax=cax,
                    orientation="horizontal", ticks=[])
    cb.outline.set_linewidth(0.4)
    cax.text(0, 1.6, "0", transform=cax.transAxes, ha="left", va="bottom",
             fontsize=8)
    cax.text(1, 1.6, f"{vmax:.0f}", transform=cax.transAxes, ha="right",
             va="bottom", fontsize=8)
    cax.text(0.5, 3.2, "ZNF257 (fold over genome mean)", transform=cax.transAxes,
             ha="center", va="bottom", fontsize=8)
    LT(F, 0.500, 0.968, "B")
    return pd.DataFrame(dict(gene=rows, peak_fold=M.max(1).round(2)))


def panel_H(F, logos=None):
    """The four published GAGGCA energy logos, with their CIS-BP accessions.

    Images and accessions both come from the source figure deck
    (20260727_Supps.pptx, slide 5); CIS-BP itself is not reachable from this
    machine, so the matrices behind them were not re-downloaded and the logos
    are reproduced as published rather than redrawn.
    """
    x0, w = 0.845, 0.150
    if logos is None:
        ax = F.add_axes([x0, 0.120, w, 0.780])
        ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#F4F4F4",
                               edgecolor="#BBBBBB", lw=0.6))
        ax.text(0.5, 0.50, "RESERVED\n\nmotif logos", ha="center", va="center",
                fontsize=8, color="#707070", linespacing=1.40)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.axis("off")
    else:
        top, bot = 0.890, 0.130
        n = len(logos)
        slot = (top - bot) / n
        for k, (img, cite, acc) in enumerate(logos):
            yy = top - slot * (k + 1)
            h = img.shape[0] / img.shape[1] * w * (7.48 / 8.70)
            ax = F.add_axes([x0, yy + (slot - h) * 0.46, w, h])
            ax.imshow(img, interpolation="bilinear")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.text(0.0, 1.02, cite, transform=ax.transAxes, fontsize=8,
                    ha="left", va="bottom")
            ax.text(0.0, -0.04, acc, transform=ax.transAxes, fontsize=8,
                    ha="left", va="top", color="#555555")
    LT(F, 0.845, 0.968, "C")

def build(S, logos=None):
    plt.close("all")
    F = plt.figure(figsize=(SHEET_W, SHEET_H))
    g, sel = panel_F(F, S["MB"])
    pk = panel_G(F, S["bw257"], S["S257BG"], S["MB"])
    panel_H(F, logos)
    return F, g, sel, pk


def export(F, out=OUT, name="Figure_S5"):
    from PIL import Image as _PILI
    _PILI.MAX_IMAGE_PIXELS = None
    F.savefig(f"{out}/{name}_PLOS.png", dpi=600, **PK(F))
    F.savefig(f"{out}/{name}_300.png", dpi=300, **PK(F))
    F.savefig(f"{out}/{name}_PLOS.svg", **PK(F))
    _im = _PILI.open(f"{out}/{name}_300.png").convert("RGB")
    _im.save(f"{out}/{name}_PLOS.tif", compression="tiff_lzw", dpi=(300, 300))
    return _im.size
