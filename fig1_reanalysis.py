#!/usr/bin/env python3
"""
Figure 1 reanalysis - Forey et al., PLOS Genetics PGENETICS-D-26-00836
Regenerates panels 1C, 1D, 1E, 1F and 1H from primary data.

Toxicity definition: day-9 normalised proliferation score <= 0.85 (Methods).
Every panel declares its own denominator, which is set by data availability.

Inputs (paths relative to REVISION_ROOT):
  Supplementary Table S1.xlsx                    per-KZFP D0/D4/D7/D9 scores, 366 KZFPs
  20250110_TE_Fam/OE_ANALYSIS.xlsx               merged_OE_Tables: evolutionary age, SCAN flag
  20250110_TE_Fam/GenomicLocation.txt            per-KZFP enrichment p at promoters and TEs, 300 KZFPs
  20250110_TE_Fam/enrich_kzfp_perFam/Fam.tsv     per-KZFP per-family enrichment, 3 tests

Outputs: fig1/Fig1{C,D,E,F,H}.svg and .tif (600 dpi, LZW)

NOTE ON PANEL F: in the published panel the x-axis labels for LTR/ERVK, LTR/ERVL
and LTR/ERVL-MaLR are rotated relative to the bars. This script labels every bar
with the family it was actually computed from.

PANEL G is not generated: the per-KZFP count of TSS bound is absent from the
source data (NbPeaks is populated for 44 of 223 KZFPs and no TSS table exists).
"""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy import stats

ROOT = os.environ.get("REVISION_ROOT", ".")
TE_DIR = os.path.join(ROOT, "20250110_TE_Fam")
OUT = "fig1"
os.makedirs(OUT, exist_ok=True)

TOXIC_THRESHOLD = 0.85
ENRICH_ALPHA = 0.05
ENRICH_TEST = "padj.hypergeom.reg"   # the test that reproduces the published panel
AGE_BINWIDTH = 10                    # Myr

# colours sampled directly from the published Figure 1 TIFF
NT_C, TX_C = "#7F7FFF", "#DF7F7F"
mpl.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8,
                     "xtick.major.width": 0.8, "ytick.major.width": 0.8,
                     "svg.fonttype": "none"})

# KZFP names that differ between the ChIP tables and Table S1
ALIAS = {"ZFP69": "ZNF69", "ZKSCAN7": "ZNF167", "ZKSCAN8": "ZNF192"}

LEG = [mpl.patches.Patch(fc=NT_C, ec="black", lw=.5, label="Not Toxic"),
       mpl.patches.Patch(fc=TX_C, ec="black", lw=.5, label="Toxic")]


def load_scores():
    s1 = pd.read_excel(os.path.join(ROOT, "Supplementary Table S1.xlsx"))
    s1.columns = ["KZFP", "D0", "D4", "D7", "D9"]
    return s1


def make_toxic_lookup(s1):
    d9 = dict(zip(s1.KZFP, s1.D9))

    def toxic(name):
        """True / False / None - None when the KZFP was not screened."""
        v = d9.get(name, d9.get(ALIAS.get(name, name)))
        return None if v is None or pd.isna(v) else bool(v <= TOXIC_THRESHOLD)

    return toxic


def build_master(s1):
    oe = os.path.join(TE_DIR, "OE_ANALYSIS.xlsx")
    merged = pd.read_excel(oe, sheet_name="merged_OE_Tables").drop_duplicates("Row.Labels")
    age = merged.set_index("Row.Labels")["age_combined"].to_dict()
    scan = merged.set_index("Row.Labels")["PresenceOfSCAN"].to_dict()
    m = s1.copy()
    m["toxic"] = m.D9 <= TOXIC_THRESHOLD
    m["age"] = m.KZFP.map(age)
    m["is_scan"] = m.KZFP.map(scan).fillna("").str.contains("SCAN")
    return m


def save(fig, name):
    fig.savefig(os.path.join(OUT, name + ".svg"), bbox_inches="tight")
    png = os.path.join(OUT, name + ".png")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    Image.open(png).convert("RGB").save(
        os.path.join(OUT, name + ".tif"), format="TIFF",
        compression="tiff_lzw", dpi=(600, 600))


def stacked(ax, ntox, tox, xlabels, fontsize=7.5):
    """100% stacked bars, toxic at the bottom, absolute counts printed inside."""
    idx = np.arange(len(ntox))
    ntox, tox = np.asarray(ntox), np.asarray(tox)
    total = ntox + tox
    pt, pn = 100 * tox / total, 100 * ntox / total
    ax.bar(idx, pt, 0.72, color=TX_C, edgecolor="black", linewidth=0.5)
    ax.bar(idx, pn, 0.72, bottom=pt, color=NT_C, edgecolor="black", linewidth=0.5)
    for i in idx:
        ax.text(i, pt[i] / 2, str(tox[i]), ha="center", va="center", fontsize=fontsize)
        ax.text(i, pt[i] + pn[i] / 2, str(ntox[i]), ha="center", va="center", fontsize=fontsize)
    ax.set_xticks(idx)
    ax.set_xticklabels(xlabels, rotation=45, ha="right", fontsize=fontsize)
    ax.set_ylabel("KZFPs (%)")
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(handles=LEG, frameon=False, fontsize=fontsize, loc="lower center",
              bbox_to_anchor=(0.5, 1.01), ncol=2, handlelength=0.9,
              handleheight=0.9, columnspacing=0.8, borderpad=0.1)


def panel_c(m):
    nt, tx = int((~m.toxic).sum()), int(m.toxic.sum())
    total = nt + tx
    fig, ax = plt.subplots(figsize=(2.5, 2.0))
    ax.pie([nt, tx], colors=[NT_C, TX_C], startangle=90, counterclock=False,
           wedgeprops=dict(linewidth=0.6, edgecolor="black"))
    ax.text(1.12, 0.62, "Toxic : {} ({:.1f}%)".format(tx, 100 * tx / total),
            fontsize=8, ha="left", va="center")
    ax.text(0.55, -0.95, "Not Toxic : {} ({:.1f}%)".format(nt, 100 * nt / total),
            fontsize=8, ha="left", va="center")
    ax.set_xlim(-1.15, 3.0)
    ax.set_ylim(-1.3, 1.25)
    ax.axis("off")
    save(fig, "Fig1C")
    print("C  {} screened: {} toxic, {} non-toxic".format(total, tx, nt))


def panel_d(m):
    ma = m.dropna(subset=["age"])
    fig, ax = plt.subplots(figsize=(3.1, 2.2))
    bins = np.arange(0, ma.age.max() + AGE_BINWIDTH, AGE_BINWIDTH)
    for sel, colour, label in [(~ma.toxic, NT_C, "Not Toxic"), (ma.toxic, TX_C, "Toxic")]:
        ax.hist(ma.age[sel], bins=bins, orientation="horizontal", density=True,
                color=colour, alpha=0.6, edgecolor="black", linewidth=0.4, label=label)
    ax.invert_yaxis()
    ax.set_ylabel("Age Combined (myo)")
    ax.set_xlabel("Fraction of KZFPs (1 = 100%)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", handlelength=0.9,
              handleheight=0.9, borderpad=0.1, labelspacing=0.3)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "Fig1D")
    print("D  {} with an age: {} toxic, {} non-toxic".format(
        len(ma), int(ma.toxic.sum()), int((~ma.toxic).sum())))


def panel_e(toxic):
    gl = pd.read_csv(os.path.join(TE_DIR, "GenomicLocation.txt"), sep="\t")
    gl["tox"] = gl.KZFP.map(toxic)
    gl = gl.dropna(subset=["tox"])
    gl["te"] = gl["TEs.pval"] <= ENRICH_ALPHA
    gl["pr"] = gl["promoter.pval"] <= ENRICH_ALPHA
    gl["cat"] = pd.Categorical(
        np.select([gl.te & ~gl.pr, ~gl.te & gl.pr, gl.te & gl.pr],
                  ["TEs", "Prom.", "Both"], default="None"),
        ["TEs", "Prom.", "Both", "None"], ordered=True)
    ct = pd.crosstab(gl.cat, gl.tox.map({True: "Toxic", False: "Not Toxic"}))
    fig, ax = plt.subplots(figsize=(2.3, 2.3))
    stacked(ax, ct["Not Toxic"].tolist(), ct["Toxic"].tolist(), list(ct.index))
    save(fig, "Fig1E")
    print("E  {} with location data and a score".format(len(gl)))
    print(ct)
    return ct


def clean_name(filename):
    name = re.sub(r"_(pubM|n)?_?vs_293T.*$", "", filename)
    for _ in range(2):
        name = re.sub(r"_(rep\d?|lowconf|n_lowconf|n|pubM|rep_n|rep_pubM)$", "", name)
    return name


def panel_f(toxic):
    fam = pd.read_csv(os.path.join(TE_DIR, "enrich_kzfp_perFam", "Fam.tsv"), sep="\t")
    fam["kzfp"] = fam.filename.map(clean_name)
    fam["tox"] = fam.kzfp.map(toxic)
    sig = fam[(fam[ENRICH_TEST] <= ENRICH_ALPHA) & fam.tox.notna()]
    families = ["LINE/L1", "LTR/ERV1", "LTR/ERVK", "LTR/ERVL",
                "LTR/ERVL-MaLR", "nonTE", "Retroposon/SVA", "SINE/Alu"]
    rows = []
    for family in families:
        g = sig[sig.subfam_name == family]
        rows.append([family,
                     g[~g.tox.astype(bool)].kzfp.nunique(),
                     g[g.tox.astype(bool)].kzfp.nunique()])
    ff = pd.DataFrame(rows, columns=["Family", "NotToxic", "Toxic"])
    fig, ax = plt.subplots(figsize=(3.3, 2.35))
    stacked(ax, ff.NotToxic.tolist(), ff.Toxic.tolist(),
            [f.replace("/", ".") for f in ff.Family])
    save(fig, "Fig1F")
    print("F")
    print(ff.to_string(index=False))
    return ff


def panel_h(m):
    a = m.D9[m.is_scan].dropna()
    b = m.D9[~m.is_scan].dropna()
    _, pval = stats.mannwhitneyu(a, b, alternative="two-sided")
    stars = "***" if pval <= .001 else "**" if pval <= .01 else "*" if pval <= .05 else "ns"
    fig, ax = plt.subplots(figsize=(1.55, 2.3))
    violins = ax.violinplot([a, b], positions=[0, 1], widths=0.8,
                            showextrema=False, showmedians=False)
    for body in violins["bodies"]:
        body.set_facecolor("#9E9E9E")
        body.set_edgecolor("black")
        body.set_linewidth(0.5)
        body.set_alpha(1)
    ax.boxplot([a, b], positions=[0, 1], widths=0.14, patch_artist=True, showfliers=False,
               medianprops=dict(color="black", lw=0.8),
               boxprops=dict(facecolor="white", lw=0.5),
               whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    ymax = max(a.max(), b.max())
    ax.plot([0, 0, 1, 1], [ymax * 1.04, ymax * 1.09, ymax * 1.09, ymax * 1.04],
            lw=0.7, color="black")
    ax.text(0.5, ymax * 1.10, stars, ha="center", va="bottom", fontsize=9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["SCAN", "No\nSCAN"], fontsize=8)
    ax.set_ylabel("Norm. A570/A600\n(Day 9)", fontsize=8)
    ax.set_ylim(0, ymax * 1.22)
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, "Fig1H")
    print("H  SCAN n={} median={:.3f} | No SCAN n={} median={:.3f} | p={:.4g} ({})".format(
        len(a), a.median(), len(b), b.median(), pval, stars))


def main():
    s1 = load_scores()
    toxic = make_toxic_lookup(s1)
    master = build_master(s1)
    panel_c(master)
    panel_d(master)
    panel_e(toxic)
    panel_f(toxic)
    panel_h(master)
    print("\nPanel G not generated: per-KZFP TSS count absent from the source data.")


if __name__ == "__main__":
    main()
