"""Rebuild every object Figure 3 needs, from the primary tables.

Usage:
    import fig3_data as D3D
    S = D3D.load()          # dict of all panel inputs
"""
import os
import json
import numpy as np
import pandas as pd
import rnaseq_panels as RP
from bigwig import BigWig

R = r"C:\Users\forey\Desktop\20260924_For_Revision"
AN = os.path.join(R, "20250811_AnalyseRNAseq")
GX = os.path.join(R, "20250828_GTEX_ANALYSIS")
OUT = "out/fig3"

# 32 fields of DE_TSS_peak_repeat_merged.tsv (positions 24-31 are the repeat hit)
COLS = ["ensembl", "entrez", "symbol_de", "genename", "mean_cond", "mean_GFP",
        "foldChange", "padj", "diffExpr", "DEflag", "Direction",
        "tchr", "tstart", "tend", "tens", "tscore", "tstrand", "symbol",
        "pchr", "pstart", "pend", "pname", "pscore", "pstrand",
        "rchr", "rstart", "rend", "rname", "rscore", "rstrand", "rclass", "dist"]

STAGES = ["Spermatogonia", "Spermatocytes", "Early spermatids", "Late spermatids"]
GERM = ["Spermatogonia", "Spermatocytes"]
MAGEA_DOWN = ["MAGEA6", "MAGEA12", "MAGEA3", "MAGEA1"]
SHOW3D = ["ZNF257", "MAGEA6", "MAGEA3", "MAGEA12"]
LAB3E = ["MAGEA6", "MAGEA12", "MAGEA3"]
ORD3I = ["MAGEA2", "MAGEA12", "MAGEA6", "MAGEA2B", "MAGEA3"]
AGE3G = [("MAGEA6", 91.0), ("MAGEA3", 29.6), ("MAGEA1", 29.6),
         ("MAGEA12", 18.0), ("MAGEA2B", 18.0), ("MAGEA2", 18.0)]
K_ROWS = [("MAGEs", ["MAGEA6", "MAGEA3"]),
          ("ZNF257", ["ZNF257"]),
          ("Spermatogonia", ["KDM1B", "MAGEB2", "PIWIL4"]),
          ("Spermatocytes", ["ANKRD31", "RBM44", "TOP2A"]),
          ("Early spermatids", ["CEP55", "KPNA5", "PBK"]),
          ("Late spermatids", ["PRM1", "PRM2", "TNP1"])]
MOTIF = "GAGGCA"
# the five contiguous MAGEA genes whose promoters carry the ZNF257 motif
CLUSTER5 = ["MAGEA6", "MAGEA2B", "MAGEA12", "MAGEA2", "MAGEA3"]
S257BG = 10.62          # genome mean of ZNF257_pubM.bw (totalSummary)


def _closest_peak(tssdf, peaks):
    """Genuine TSS-to-peak distance per chromosome (no bedtools -1 sentinel)."""
    up = peaks[["pchr", "pstart", "pend", "pname", "pscore"]].drop_duplicates().copy()
    res = []
    for c, grp in tssdf.groupby("chr", sort=False):
        pk = up[up.pchr == c]
        if len(pk) == 0:
            continue
        pk = pk.sort_values("pstart").reset_index(drop=True)
        ps, pe = pk.pstart.values, pk.pend.values
        t = grp.start.values
        idx = np.searchsorted(ps, t)
        best = np.zeros(len(t), dtype=int)
        bd = np.full(len(t), np.inf)
        for off in (-1, 0, 1):
            j = np.clip(idx + off, 0, len(pk) - 1)
            d = np.where(t < ps[j], ps[j] - t,
                         np.where(t >= pe[j], t - pe[j] + 1, 0)).astype(float)
            m = d < bd
            bd[m], best[m] = d[m], j[m]
        o = grp.copy()
        o["pname"] = pk.pname.values[best]
        o["dist"] = bd.astype(int)
        res.append(o)
    return pd.concat(res, ignore_index=True)


def build_znf257_merged():
    """ZNF257 has no DE_TSS_peak_repeat_merged.tsv; rebuild it from the parts."""
    tss = pd.read_csv(os.path.join(AN, "1807_hg19_ens_coding_genes_tss_symbol.bed"),
                      sep="\t", header=None,
                      names=["chr", "start", "stop", "ensembl", "score", "strand", "symbol"])
    pk = pd.read_csv(os.path.join(AN, "analysis_results_ZNF257",
                                  "peaks_with_repeat_info.tsv"), sep="\t", header=None)
    pk.columns = ["pchr", "pstart", "pend", "pname", "pscore", "pstrand",
                  "rchr", "rstart", "rend", "rname", "rscore", "rstrand", "rclass"]
    de = pd.read_csv(os.path.join(AN, "analysis_results_ZNF257",
                                  "step2_DE_direction.tsv"), sep="\t")
    de = de.rename(columns={"mean.257": "mean_cond", "mean.GFP": "mean_GFP"})[
        ["ensembl", "entrez", "symbol", "genename", "mean_cond", "mean_GFP",
         "foldChange", "padj", "diffExpr", "DE", "Direction"]]
    de.columns = ["ensembl", "entrez", "symbol_de", "genename", "mean_cond",
                  "mean_GFP", "foldChange", "padj", "diffExpr", "DEflag", "Direction"]
    cp = _closest_peak(tss, pk).rename(
        columns={"chr": "tchr", "start": "tstart", "stop": "tend",
                 "ensembl": "tens", "score": "tscore", "strand": "tstrand"})
    cp = cp.merge(pk[["pname", "pchr", "pstart", "pend", "pscore", "pstrand",
                      "rchr", "rstart", "rend", "rname", "rscore", "rstrand",
                      "rclass"]].drop_duplicates("pname"), on="pname", how="left")
    g = de.merge(cp, left_on="ensembl", right_on="tens", how="right")
    for c in COLS:
        if c not in g.columns:
            g[c] = np.nan
    g = g[COLS].copy()
    for c in ("tstart", "tend", "pstart", "pend", "rstart", "rend", "dist",
              "foldChange", "padj", "pscore", "tscore"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    return g


def germ_table(sc):
    """Per-gene mean nTPM in spermatogonia+spermatocytes vs all other cell types."""
    g = sc[sc["Cell type"].isin(GERM)].groupby("Gene name").nTPM.mean()
    o = sc[~sc["Cell type"].isin(GERM)].groupby("Gene name").nTPM.mean()
    D = pd.concat([g.rename("germ"), o.rename("other")], axis=1).dropna()
    D["lg_germ"] = np.log10(D.germ + 1)
    D["lg_other"] = np.log10(D.other + 1)
    D["delta"] = D.lg_germ - D.lg_other
    D["pct_delta"] = D.delta.rank(pct=True) * 100
    return D


def align_promoters(cache=os.path.join(OUT, "magea_promoters.json")):
    """FAMSA alignment of the five 801 bp MAGEA promoter windows."""
    from pyfamsa import Aligner, Sequence
    raw = json.load(open(cache))
    seqs = {k: v[0] for k, v in raw.items()}       # [sequence, coords, strand]
    coords = {k: (v[1], v[2]) for k, v in raw.items()}
    al = Aligner(guide_tree="upgma").align(
        [Sequence(k.encode(), v.encode()) for k, v in seqs.items()])
    AL = {s.id.decode(): s.sequence.decode() for s in al}
    names = list(AL)
    L = len(AL[names[0]])
    cols = [[AL[n][i] for n in names] for i in range(L)]
    ident = np.array([(sum(1 for x in c if x == c[0]) / len(c)) if c[0] != "-" else 0.0
                      for c in cols])
    consensus = "".join(max(set(c), key=c.count) for c in cols)
    # all alignment columns where every promoter reads the motif identically
    cand = [i for i in range(L - len(MOTIF))
            if all(AL[n][i:i + len(MOTIF)] == MOTIF for n in names)]
    if not cand:
        raise RuntimeError("motif not conserved in all five promoters")
    # the copy sitting in the most conserved neighbourhood
    pos = max(cand, key=lambda i: ident[max(0, i - 20):i + 26].mean())
    return seqs, coords, AL, ident, consensus, int(pos), cand


def load():
    S = {}
    S["R"] = R
    # ---- peaks / DE / TSS ------------------------------------------------
    Z = build_znf257_merged()
    S["Z257"] = Z
    g0 = RP.gene_table(Z, 1000)
    S["G257"] = g0[g0.tens.notna() & (g0.tens.astype(str) != ".")]
    S["S257"] = RP.group_assign(RP.alluvial_table(Z, maxdist=50000), min_frac=0.045)
    # panel E is the whole transcriptome, not only the TSS-annotated subset
    de = pd.read_csv(os.path.join(AN, "analysis_results_ZNF257",
                                  "step2_DE_direction.tsv"), sep="\t")
    de = de.rename(columns={"symbol": "symbol_de"})[
        ["ensembl", "symbol_de", "foldChange", "padj", "Direction"]
    ].drop_duplicates("ensembl")
    S["de257"] = de
    g = S["G257"]
    S["t257"] = sorted(g[(g.cat == "DOWN") & g.haspeak].symbol_de.dropna().unique())
    # ---- HPA single cell -------------------------------------------------
    sc = pd.read_csv(os.path.join(GX, "rna_single_cell_type.tsv"), sep="\t")
    S["sc"] = sc
    S["D3"] = germ_table(sc)
    # ---- bigwig / MAGEA locus -------------------------------------------
    bw = BigWig(os.path.join(R, "20251013_ZNF257PhyloP", "ZNF257_pubM.bw"))
    S["bw257"] = bw
    S["S257BG"] = S257BG
    MB = pd.read_csv(os.path.join(R, "20260115_ZNF257_MAGEA", "MAGE.bed"),
                     sep="\t", header=None)
    S["MB"] = MB
    clu = MB[(MB[0] == "chrX") & MB[3].isin(CLUSTER5)]
    FA, FB = int(clu[1].min()) - 8000, int(clu[2].max()) + 8000
    S["FA"], S["FB"] = FA, FB
    pk = Z[(Z.pchr == "chrX") & Z.pstart.between(FA, FB)].drop_duplicates("pname").copy()
    pk["summit"] = ((pk.pstart + pk.pend) / 2).astype(int)
    S["pk257"] = pk.sort_values("summit")
    # peak summit nearest the MAGEA6 TSS, for panel H
    gi = MB[(MB[0] == "chrX") & (MB[3] == "MAGEA6")].iloc[0]
    tss6 = int(gi[1]) if gi[5] == "+" else int(gi[2])
    S["tss6"] = tss6
    S["summit"] = int(pk.loc[(pk.summit - tss6).abs().idxmin(), "summit"])
    # ---- promoter alignment ---------------------------------------------
    (S["SEQS"], S["COORDS"], S["AL"], S["ident"], S["consensus"],
     S["MSTART"], S["MOTIF_HITS"]) = align_promoters()
    S["ORD3I"] = ORD3I
    S["AGE3G"] = AGE3G
    S["K_ROWS"] = K_ROWS
    S["STAGES"] = STAGES
    S["SHOW3D"] = SHOW3D
    S["LAB3E"] = LAB3E
    return S
