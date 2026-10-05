"""Peak-caller-independent scan of the ZNF498 coverage track on chr18.

Answers the question 'are there really no ZNF498 peaks on chr18?' without using the
peak caller at all: every 200-bp window (the median called-peak width is 207 bp) of
ZNF498_n_noid.bw is scored as mean coverage / genome mean, and every window reaching
the weakest called peak anywhere in the genome (5.68x) is reported.

Panel A  chr18, the 18 regions that reach that threshold, coloured by the repeat class
         they overlap (17 of 18 overlap one), with the 37 upregulated chr18 gene TSSs
         as a rug below.
Panel B  the same scan on three size-matched chromosomes that DO carry called peaks,
         as a positive control: the scan recovers 25 of their 28 called peaks.

build(RG, SUM, ctrl, up_tss, thr, peri) -> matplotlib Figure
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

RC = {"font.family": "Arial", "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
      "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
      "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
      "xtick.major.size": 2.4, "ytick.major.size": 2.4,
      "svg.fonttype": "none", "pdf.fonttype": 42}

C_GRP = {"Satellite": "#C0392B", "Simple repeat": "#E08A3C",
         "Transposon": "#7D6FB0", "No repeat": "#1F4E9C"}
ORDER = ["Satellite", "Simple repeat", "Transposon", "No repeat"]
SHEET_W, SHEET_H = 7.48, 3.60


def _letter(F, x, y, s):
    F.text(x, y, s, fontsize=10, fontweight="bold", va="top", ha="left")


def panel_A(F, RG, up_tss, thr, peri, rect, rug_rect):
    ax = F.add_axes(rect)
    ax.axvspan(peri[0] / 1e6, peri[1] / 1e6, color="#EAEAEA", lw=0, zorder=0)
    ax.text(np.mean(peri) / 1e6, thr * 0.95, "\u03b1-satellite", ha="center", va="bottom",
            color="#777777", fontsize=8, zorder=2)
    ax.axhline(thr, color="#444444", lw=0.8, ls="--", zorder=1)
    ax.text(30.0, thr * 1.25,
            f"weakest called peak genome-wide ({thr:.2f}\u00d7)",
            ha="left", va="bottom", fontsize=8, color="#444444")
    for g in ORDER:
        d = RG[RG.grp == g]
        if not len(d):
            continue
        x = d.mid.values / 1e6
        ax.vlines(x, thr, d.maxfold.values, color=C_GRP[g], lw=1.0, zorder=3)
        ax.plot(x, d.maxfold.values, "o", ms=3.2, mfc=C_GRP[g], mec="white",
                mew=0.4, zorder=4, clip_on=False)
    for r, dx, ha in [(RG.nlargest(1, "maxfold").iloc[0], 1.6, "left"),
                      (RG.nlargest(2, "maxfold").iloc[1], 1.6, "left")]:
        ax.annotate(r.rep_name.split("|")[0], (r.mid / 1e6 + dx, r.maxfold),
                    ha=ha, va="center", fontsize=8, color=C_GRP[r.grp])
    near = RG.loc[RG.dist_UPgene.idxmin()]
    ax.annotate(f"(CGGGG)n, {int(near.dist_UPgene)} bp from {near.nearest_gene}",
                xy=(near.mid / 1e6, near.maxfold), xytext=(12.0, 90),
                fontsize=8, color="#333333", ha="left", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.6, color="#888888",
                                shrinkA=0, shrinkB=2, relpos=(0.0, 0.0)))
    ax.set_yscale("log")
    ax.set_ylim(thr * 0.80, 5000)
    ax.set_xlim(-0.6, 78.6)
    ax.set_ylabel("Coverage /\ngenome mean", labelpad=1.5, linespacing=1.2)
    ax.set_title("ZNF498 coverage track scanned on chr18 in 200-bp windows, "
                 "independently of the peak caller", fontsize=8, pad=3.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelbottom=False, length=0)

    h = [Line2D([], [], color=C_GRP[g], lw=1.0, marker="o", ms=3.2, mec="white",
                mew=0.4, label=g) for g in ORDER if (RG.grp == g).any()]
    ax.legend(handles=h, loc="upper left", bbox_to_anchor=(0.285, 1.015), frameon=False,
              ncol=2, handlelength=1.1, handletextpad=0.5, labelspacing=0.3,
              columnspacing=1.3, borderaxespad=0.0,
              title="repeat class overlapping the window", title_fontsize=8,
              alignment="left")

    axr = F.add_axes(rug_rect, sharex=ax)
    axr.vlines(np.asarray(up_tss) / 1e6, 0, 1, color="#1A1A1A", lw=0.7)
    axr.set_ylim(0, 1)
    axr.set_yticks([])
    axr.set_xlabel(f"chr18 position (Mb, hg19); ticks above mark the TSS of the "
                   f"{len(up_tss)} upregulated chr18 genes", labelpad=2.0)
    for s in ("top", "left", "right"):
        axr.spines[s].set_visible(False)
    return ax


def panel_B(F, SUM, rect, THR_LAB):
    ax = F.add_axes(rect)
    chs = list(SUM.chr)
    xs = np.arange(len(chs))
    bot = np.zeros(len(chs))
    keys = [("satellite_simple", "Satellite / simple repeat", "#C0392B"),
            ("transposon", "Transposon", "#7D6FB0"),
            ("no_repeat", "No repeat", "#1F4E9C")]
    for k, lab, col in keys:
        v = (SUM[k] / SUM.Mb).values
        ax.bar(xs, v, 0.62, bottom=bot, color=col, lw=0.4, edgecolor="white", label=lab)
        bot += v
    for i, r in enumerate(SUM.itertuples()):
        ax.text(i, bot[i] + 0.035,
                "no called\npeaks" if r.called_peaks == 0
                else f"{r.recovered}/{r.called_peaks} called\npeaks recovered",
                ha="center", va="bottom", fontsize=8, linespacing=1.15,
                color="#C0392B" if r.called_peaks == 0 else "#333333")
    ax.set_xticks(xs)
    ax.set_xticklabels([c.replace("chr", "chr ") for c in chs])
    ax.set_ylim(0, max(bot) * 1.58)
    ax.set_ylabel(f"Regions reaching {THR_LAB}\u00d7\nper Mb", labelpad=1.5, linespacing=1.2)
    ax.set_title("Same scan on size-matched chromosomes", fontsize=8, pad=3.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.195), frameon=False, ncol=3,
              handlelength=0.9, handleheight=0.9, handletextpad=0.45, columnspacing=1.2,
              borderaxespad=0.0)
    return ax


def build(RG, SUM, ctrl, up_tss, thr, peri):
    with mpl.rc_context(RC):
        F = plt.figure(figsize=(SHEET_W, SHEET_H))
        panel_A(F, RG, up_tss, thr, peri,
                rect=[0.098, 0.590, 0.880, 0.320],
                rug_rect=[0.098, 0.534, 0.880, 0.024])
        panel_B(F, SUM, rect=[0.098, 0.150, 0.300, 0.250], THR_LAB=f"{thr:.2f}")
        _letter(F, 0.008, 0.990, "A")
        _letter(F, 0.008, 0.455, "B")
        F.text(0.465, 0.430,
               "17 of the 18 chr18 regions reaching called-peak signal\n"
               "strength overlap a repeat: pericentromeric and subtelomeric\n"
               "satellite, the telomeric (CCCTAA)n array at the chromosome\n"
               "start, Alu, LTR and simple repeats. An input control removes\n"
               "exactly these, which is why the ZNF498-vs-input peak call\n"
               "returns none. The 18th has no repeat annotation but sits\n"
               "714 bp from one. Only one region lies within 10 kb of an\n"
               "upregulated gene, a (CGGGG)n repeat 737 bp from PPP4R1;\n"
               "chr18 carries 4-fold fewer such regions per Mb than controls.",
               fontsize=8, va="top", ha="left", linespacing=1.42, color="#1A1A1A")
    return F
