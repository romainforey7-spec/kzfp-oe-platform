"""Rebuild the objects Figures 4 (ZNF498) and 5 (ZNF18) need.

Both KZFPs act as activators in the overexpression screen, so the claimed
class is UP and the peak window is 1 kb.

The chr18 sentinel trap (E-07): bedtools closest -d emits dist = -1 for every
gene on a chromosome that carries no peak. gene_table and alluvial_table both
drop dist < 0, so the 355 ZNF498 sentinel rows (288 chr18, 54 chrY, 13 chrM) no
longer masquerade as the most proximal targets.

Those rows carry no PEAK-based distance, but the question of whether ZNF498 is
bound there has now been settled directly from the coverage track
(ZNF498_n_noid.bw, hg19). At TSS +/- 1 kb, normalised to the genome mean of
6.93: the 34 measured targets sit at a median 2.75x (85% above 2x, Mann-Whitney
vs unchanged genes p = 2.3e-17), whereas the 38 chr18 UP genes sit at a median
1.07x against an unchanged-gene background of 1.18x (p = 0.88, i.e. no
enrichment). chr18 is covered normally - 0.91x the genome mean with 95.8%
non-zero bins, like chr17 1.02x and chr19 1.01x - so this is a measurement, not
a gap. The three genes behind the published microtubule term read ROCK1 1.62x,
MAPRE2 1.90x and SKA1 0.63x, inside the background distribution, while STMN3
(chr20, the one with a called peak) reads 4.07x. The peak caller did not miss
them: they are unbound. Audit in ZNF498_bigwig_promoter_occupancy.csv and
ZNF498_bigwig_check_summary.csv.
"""
import os
import numpy as np
import pandas as pd
import rnaseq_panels as RP
import fig2_data as D2

R, AN, GX = D2.R, D2.AN, D2.GX
COLS = D2.COLS
FIGOF = {"ZNF498": "4", "ZNF18": "5"}
CLAIM = "UP"
WINDOW = 1000


def load_one(kz):
    p = os.path.join(AN, f"analysis_results_{kz}", "DE_TSS_peak_repeat_merged.tsv")
    Z = D2.load_merged(p)
    g0 = RP.gene_table(Z, WINDOW)
    G = g0[g0.tens.notna() & (g0.tens.astype(str) != ".")]
    Sl = RP.group_assign(RP.alluvial_table(Z, maxdist=50000), min_frac=0.045)
    tg = sorted(G[(G.cat == CLAIM) & G.haspeak].symbol_de.dropna().unique())
    # how many rows the chr18-style sentinel would have contributed
    sent = int((pd.to_numeric(Z.dist, errors="coerce") < 0).sum())
    chroms_no_peak = sorted(set(Z.tchr.dropna().unique())
                            - set(Z.loc[Z.pname.notna() & (Z.pname != "."),
                                        "pchr"].dropna().unique()))
    return dict(Z=Z, G=G, S=Sl, targets=tg, n_sentinel=sent,
                chroms_without_peak=[c for c in chroms_no_peak if c != "."])


def load():
    S = {}
    for kz in FIGOF:
        S[kz] = load_one(kz)
    # shared single-cell / tissue tables for the target heatmaps
    S["sc"] = pd.read_csv(os.path.join(GX, "rna_single_cell_type.tsv"), sep="\t")
    S["tc"] = pd.read_csv(os.path.join(GX, "rna_tissue_consensus.tsv"), sep="\t")
    return S
