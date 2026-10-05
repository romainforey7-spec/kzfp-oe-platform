
"""Panel K of Figure 3 - how the MAGEA cluster came to share one ZNF257 site.

The model is evolutionary, not per-promoter: the GAGGCA element sits inside the
5' end of the transcription unit, so tandem duplication of an ancestral MAGEA6
carried it into every descendant copy and the cluster is now co-regulated by a
single factor. The ages in panel I carry the argument - MAGEA6 at 91 Myr is the
oldest member, its descendants are 18-30 Myr.

Feature placement is from Ensembl GRCh37 and is drawn as measured: the motif is
88 bp into the 94-bp first exon of MAGEA6, 6 bp from the splice junction, and the
peak summit is 118 bp further into intron 1. Audited in
Fig3H_motif_gene_feature.csv across all 28 isoforms of the five genes.

Drawn natively; nothing below 8 pt.
"""
from matplotlib.patches import Ellipse, FancyBboxPatch, Rectangle, FancyArrowPatch

C_DNA = "#3A3A3A"
C_ZNF = "#9B2FAE"
C_MOT = "#C8102E"
C_GENE = "#B00000"
C_ANC = "#6E6E6E"

X0, XW = 0.030, 0.940
Y_ANC, Y_NOW = 0.740, 0.250
#: present-day cluster, left to right along chrX, with the age from panel I
NOW = [("A6", "91"), ("A2B", "18"), ("A12", "18"), ("A2", "18"), ("A3", "30")]


def _exon(ax, x, y, w, col, h=0.052, z=4):
    ax.add_patch(FancyBboxPatch((x, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.006",
                                facecolor=col, edgecolor="none", zorder=z))


def _motif(ax, x, y, w=0.016, h=0.070, z=6):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, facecolor=C_MOT,
                           edgecolor="none", zorder=z))


def _znf257(ax, x, y, s=1.0):
    w, h = 0.082 * s, 0.085 * s
    ax.add_patch(Ellipse((x, y), w * 0.60, h, facecolor=C_ZNF, edgecolor="black",
                         lw=0.5, zorder=6))
    bx = x + w * 0.34
    ax.add_patch(FancyBboxPatch((bx, y - h / 2), w * 0.60, h,
                                boxstyle="round,pad=0,rounding_size=0.010",
                                facecolor="#EADBEE", edgecolor="black", lw=0.5,
                                zorder=6))
    for i in range(3):
        ax.add_patch(Rectangle((bx + w * 0.10 + i * w * 0.165, y - h * 0.26),
                               w * 0.075, h * 0.52, facecolor=C_ZNF,
                               edgecolor="none", zorder=7))


def draw_model(ax):
    """Two rows: the ancestral gene with its element, then the present cluster."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    ax.text(0.0, 0.998, "Tandem duplication spread one ZNF257 site",
            ha="left", va="top", fontsize=8, fontweight="bold")

    # --- ancestral locus: exon 1 - intron 1 - exon 2, element at the junction --
    e1x, e1w = X0 + XW * 0.055, XW * 0.055
    i1w = XW * 0.150
    e2x = e1x + e1w + i1w
    ax.plot([X0, e1x], [Y_ANC, Y_ANC], lw=1.2, color=C_DNA, zorder=2)
    ax.plot([e1x + e1w, e2x], [Y_ANC, Y_ANC], lw=0.8, color=C_ANC, zorder=2)
    _exon(ax, e1x, Y_ANC, e1w, C_ANC)
    _exon(ax, e2x, Y_ANC, XW * 0.075, C_ANC)
    ax.plot([e2x + XW * 0.075, X0 + XW * 0.34], [Y_ANC, Y_ANC], lw=1.2,
            color=C_DNA, zorder=2)
    _motif(ax, e1x + e1w * 0.92, Y_ANC)
    ax.annotate("GAGGCA", xy=(e1x + e1w * 0.92, Y_ANC + 0.040),
                xytext=(X0 + XW * 0.330, Y_ANC + 0.125), fontsize=8, color=C_MOT,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=C_MOT,
                                shrinkA=1, shrinkB=1))
    ax.text(e1x + e1w * 0.5, Y_ANC - 0.052, "ex 1", ha="center", va="top",
            fontsize=8, color=C_ANC)
    ax.text(e1x + e1w + i1w * 0.62, Y_ANC - 0.052, "intron 1", ha="center",
            va="top", fontsize=8, color=C_ANC)
    ax.text(X0 + XW * 0.355, Y_ANC - 0.005, "ancestral $MAGEA6$\n91 Myr",
            ha="left", va="center", fontsize=8, color=C_ANC, linespacing=1.25,
            style="italic")
    ax.text(X0 + XW, 0.545, "each copy inherits the element", ha="right",
            va="center", fontsize=8, color="#6E6E6E")

    # --- present-day cluster ----------------------------------------------
    ax.plot([X0, X0 + XW], [Y_NOW, Y_NOW], lw=1.2, color=C_DNA, zorder=2)
    step = XW / len(NOW)
    for i, (nm, age) in enumerate(NOW):
        cx = X0 + step * (i + 0.5)
        gx, gw = cx - step * 0.30, step * 0.60
        _exon(ax, gx, Y_NOW, gw * 0.16, C_GENE)
        ax.plot([gx + gw * 0.16, gx + gw * 0.56], [Y_NOW, Y_NOW], lw=0.8,
                color=C_GENE, zorder=3)
        _exon(ax, gx + gw * 0.56, Y_NOW, gw * 0.44, C_GENE)
        _motif(ax, gx + gw * 0.145, Y_NOW)
        ax.text(cx + step * 0.02, Y_NOW - 0.044, nm, ha="center", va="top",
                fontsize=8, style="italic", color=C_GENE)
        ax.text(cx + step * 0.30, Y_NOW + 0.052, age, ha="center", va="bottom",
                fontsize=8, color="#6E6E6E")
        _znf257(ax, gx + gw * 0.145 - 0.030, Y_NOW + 0.168, s=0.92)
        ax.plot([gx + gw * 0.145, gx + gw * 0.145], [Y_NOW + 0.124, Y_NOW + 0.042],
                lw=0.6, ls=(0, (1.2, 1.0)), color=C_ZNF, zorder=3)
    ax.text(X0, Y_NOW - 0.125, "$MAGEA$ copies, Myr above \u2192 co-regulation",
            ha="left", va="top", fontsize=8, color="#6E6E6E")
    ax.text(X0, Y_NOW + 0.268, "ZNF257", ha="left", va="bottom", fontsize=8,
            style="italic", fontweight="bold", color=C_ZNF)
