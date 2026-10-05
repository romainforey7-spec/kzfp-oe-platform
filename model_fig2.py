"""Figure 2F: ZNF43 model, drawn natively so every label is 8 pt.

Replaces the lifted published artwork, whose lettering measured ~4.5 pt.
Two states of the same locus:
  left  - tissues where ZNF43 is expressed (thymus, bone marrow): ZNF43 is
          bound to the LTR/ERV1 element, H3K9me3 spreads over the promoter,
          the target is repressed.
  right - liver, where ZNF43 is low: the promoter carries H3K4me3/H3K27ac and
          the targets are transcribed.
"""
import numpy as np
from matplotlib.patches import (Rectangle, FancyBboxPatch, Ellipse, Polygon,
                                FancyArrowPatch)

C_DNA = "#3A3A3A"
C_ERV = "#1F6FB4"
C_GENE = "#B3B3B3"
C_ZNF = "#1F4E9C"
C_OFF = "#B00000"
C_ON = "#2E7D32"
C_K9 = "#5E35B1"
C_K4 = "#2E7D32"
FS = 8.0


def _znf43(ax, x, yD):
    """KZFP cartoon sitting on the DNA at yD: KRAB lobe, ZF box, four fingers.

    Both labels are 8 pt and each shape is sized from the label it carries, so
    the two can never touch.
    """
    wE, hE = 0.135, 0.072           # KRAB lobe - fits "KRAB" at 8 pt
    wB, hB = 0.125, 0.072           # ZF box    - fits "ZNF43" at 8 pt
    yB = yD + 0.014                 # box bottom
    yE = yB + hB + 0.006            # lobe bottom
    ax.add_patch(Ellipse((x, yE + hE / 2), wE, hE, facecolor=C_ZNF,
                         edgecolor="black", lw=0.5, zorder=6))
    ax.add_patch(FancyBboxPatch((x - wB / 2, yB), wB, hB,
                                boxstyle="round,pad=0,rounding_size=0.010",
                                facecolor="#6E8FCB", edgecolor="black", lw=0.5,
                                zorder=6))
    for i in range(4):
        fx = x - wB * 0.33 + i * (wB * 0.22)
        ax.add_patch(Polygon([[fx - wB * 0.055, yB],
                              [fx + wB * 0.055, yB],
                              [fx, yD + 0.001]], closed=True,
                             facecolor="#6E8FCB", edgecolor="black", lw=0.4,
                             zorder=5))
    ax.text(x, yE + hE / 2, "KRAB", ha="center", va="center", fontsize=FS,
            color="white", zorder=7)
    ax.text(x, yB + hB / 2, "ZNF43", ha="center", va="center", fontsize=FS,
            color="white", style="italic", zorder=7)
    return yE + hE


def _marks(ax, xs, y, color, label, lab_x=None):
    for xx in xs:
        ax.add_patch(Ellipse((xx, y), 0.016, 0.030, facecolor=color,
                             edgecolor="black", lw=0.3, zorder=5))
    ax.text(lab_x if lab_x is not None else float(np.mean(xs)), y + 0.052, label,
            ha="center", va="bottom", fontsize=FS, color=color)


def _locus(ax, x0, repressed):
    """One copy of the locus. x0 = left edge of a 0.44-wide half.

    Row budget (axes fraction, 1.0 = 1.72 in; an 8 pt line is 0.065):
      0.975 header  0.905 subtitle  0.825 mark label  0.775 marks
      0.565-0.710 protein   0.555 DNA   0.500 LTR label   0.430 TSS arrow
      0.340-0.020 outcome block
    """
    W = 0.44
    yD = 0.555
    ax.plot([x0 + 0.010, x0 + W - 0.010], [yD, yD], color=C_DNA, lw=1.6,
            solid_capstyle="butt", zorder=3)
    xe0, xe1 = x0 + 0.035, x0 + 0.165            # LTR/ERV1 element
    ax.add_patch(Rectangle((xe0, yD - 0.016), xe1 - xe0, 0.032,
                           facecolor=C_ERV, edgecolor="black", lw=0.4, zorder=4))
    xt = x0 + 0.230                              # TSS
    xg1 = x0 + W - 0.025                         # gene body
    ax.add_patch(Rectangle((xt, yD - 0.014), xg1 - xt, 0.028,
                           facecolor=C_GENE, edgecolor="black", lw=0.4, zorder=4))
    ax.text((xe0 + xe1) / 2, 0.500, "LTR/ERV1", ha="center", va="top",
            fontsize=FS, color=C_ERV, style="italic")
    acol = C_OFF if repressed else C_ON
    ax.plot([xt, xt], [yD - 0.016, 0.430], color=acol, lw=1.0, zorder=4)
    mx = np.linspace(xt - 0.012, xt + 0.078, 4)
    if repressed:
        _znf43(ax, (xe0 + xe1) / 2, yD)
        _marks(ax, mx, 0.700, C_K9, "H3K9me3", lab_x=xt + 0.033)
        # blocked transcription: line ending in a perpendicular bar
        ax.plot([xt, xt + 0.055], [0.430, 0.430], color=C_OFF, lw=1.0, zorder=4)
        ax.plot([xt + 0.055, xt + 0.055], [0.406, 0.454], color=C_OFF, lw=1.4,
                solid_capstyle="butt", zorder=6)
        ax.text(xt + 0.072, 0.430, "OFF", ha="left", va="center", fontsize=FS,
                color=C_OFF, fontweight="bold")
    else:
        _marks(ax, mx, 0.700, C_K4, "H3K4me3 / H3K27ac", lab_x=xt + 0.033)
        ax.add_patch(FancyArrowPatch((xt, 0.430), (xt + 0.080, 0.430),
                                     arrowstyle="-|>", mutation_scale=5,
                                     lw=1.0, color=C_ON, zorder=4))
        xw = np.linspace(xt + 0.085, xt + 0.155, 60)
        ax.plot(xw, 0.430 + 0.011 * np.sin((xw - xt) * 170), color=C_ON, lw=0.9,
                zorder=5)
        ax.add_patch(FancyArrowPatch((xt + 0.155, 0.430), (xt + 0.172, 0.430),
                                     arrowstyle="-|>", mutation_scale=5, lw=0.9,
                                     color=C_ON, zorder=5))
        ax.text(xt + 0.090, 0.478, "ON", ha="left", va="center", fontsize=FS,
                color=C_ON, fontweight="bold")


def draw_model(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.axis("off")
    ax.plot([0.50, 0.50], [0.02, 0.985], color="#BBBBBB", lw=0.6, zorder=1)
    for x0, head, sub, rep in [
            (0.03, "Thymus, bone marrow", "$ZNF43$ high", True),
            (0.53, "Liver", "$ZNF43$ low", False)]:
        ax.text(x0 + 0.22, 0.975, head, ha="center", va="top", fontsize=FS,
                fontweight="bold")
        ax.text(x0 + 0.22, 0.905, sub, ha="center", va="top", fontsize=FS,
                color="#555555")
        _locus(ax, x0, rep)
    ax.text(0.25, 0.350, "Target genes\nrepressed", ha="center", va="top",
            fontsize=FS, color=C_OFF, linespacing=1.20)
    ax.text(0.75, 0.350,
            "$APOL1$, $APOL2$, $HYAL1$, $CLTB$,\n"
            "$ECI1$, $GSTO1$, $NFU1$, $DNAI4$\n"
            "detoxification,\nfatty-acid metabolism",
            ha="center", va="top", fontsize=FS, color=C_ON, linespacing=1.20)
    return ax
