"""Rebuild every object Figure 2 needs, from the primary tables.

Usage:
    import fig2_data as D2D
    S = D2D.load()
"""
import os
import struct
import numpy as np
import pandas as pd
from scipy import stats as st
import rnaseq_panels as RP
from bigwig import BigWig

R = r"C:\Users\forey\Desktop\20260924_For_Revision"
AN = os.path.join(R, "20250811_AnalyseRNAseq")
GX = os.path.join(R, "20250828_GTEX_ANALYSIS")
LTDIR = os.path.join(R, "20250205_RomainAnalysis", "20250205_BASE",
                     "20251212_ZNF43_Liver_Thymus")

COLS = ["ensembl", "entrez", "symbol_de", "genename", "mean_cond", "mean_GFP",
        "foldChange", "padj", "diffExpr", "DEflag", "Direction",
        "tchr", "tstart", "tend", "tens", "tscore", "tstrand", "symbol",
        "pchr", "pstart", "pend", "pname", "pscore", "pstrand",
        "rchr", "rstart", "rend", "rname", "rscore", "rstrand", "rclass", "dist"]

# 8 ZNF43 targets (DOWN & peak <= 10 kb). WDR78 is DNAI4 in the HPA tables.
TGT8 = ["HYAL1", "APOL1", "APOL2", "ECI1", "CLTB", "GSTO1", "NFU1", "DNAI4"]
LOCI5 = [("APOL1 / APOL2", "APOL1"), ("HYAL1", "HYAL1"), ("CLTB", "CLTB"),
         ("ECI1", "ECI1"), ("NFU1", "NFU1")]
TRKFILES = {"ZNF43": os.path.join(LTDIR, "ZNF43_n.bw"),
            "Liver H3K4me3": os.path.join(LTDIR, "liver", "Liver_H3K4me3.bigWig"),
            "Thymus H3K4me3": os.path.join(LTDIR, "Thymus", "Thymus_H3K4me3.bigWig"),
            "Liver H3K9me3": os.path.join(LTDIR, "liver", "Liver_H3K9me3.bigWig"),
            "Thymus H3K9me3": os.path.join(LTDIR, "Thymus", "Thymus_H3K9me3.bigWig")}


def load_merged(path):
    df = pd.read_csv(path, sep="\t", header=None, dtype=str, low_memory=False)
    df = df[~df[0].astype(str).str.startswith("ensembl")]
    df.columns = COLS
    for c in ("tstart", "tend", "pstart", "pend", "rstart", "rend", "dist",
              "foldChange", "padj", "pscore", "tscore"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def tsum(bw):
    """Genome mean of a bigWig, from its totalSummary block."""
    bw.f.seek(bw.totalSummaryOffset)
    vc, mn, mx, sd, ssq = struct.unpack(bw.e + "Qdddd", bw.f.read(40))
    return sd / vc


def load():
    S = {}
    Z = load_merged(os.path.join(AN, "analysis_results_ZNF43",
                                 "DE_TSS_peak_repeat_merged.tsv"))
    S["Z43"] = Z
    g0 = RP.gene_table(Z, 10000)
    S["G43"] = g0[g0.tens.notna() & (g0.tens.astype(str) != ".")]
    S["S43"] = RP.group_assign(RP.alluvial_table(Z, maxdist=50000), min_frac=0.045)

    # ---- tracks ----------------------------------------------------------
    BW = {k: BigWig(v) for k, v in TRKFILES.items()}
    S["BW"], S["TRK"] = BW, list(TRKFILES)
    S["SUM"] = {k: tsum(b) for k, b in BW.items()}

    # ---- TSS bed ---------------------------------------------------------
    TB = pd.read_csv(os.path.join(AN, "1807_hg19_ens_coding_genes_tss_symbol.bed"),
                     sep="\t", header=None,
                     names=["chr", "start", "end", "ensembl", "score", "strand", "symbol"])
    TB["tss"] = TB.start
    S["TB"] = TB
    TBI = {c: gg.sort_values("tss").reset_index(drop=True) for c, gg in TB.groupby("chr")}

    def tss_in(ch, a, b):
        d = TBI.get(ch)
        cols = ["symbol", "tss", "strand"]
        return d[d.tss.between(a, b)][cols] if d is not None else pd.DataFrame(columns=cols)

    # ---- repeats / peaks -------------------------------------------------
    RPALL = Z[Z.rname.notna() & (Z.rname != ".")].copy()
    PKALL = Z[Z.pname.notna() & (Z.pname != ".")].copy()
    S["RPALL"], S["PKALL"] = RPALL, PKALL

    sub = Z[Z.symbol.isin(TGT8) & Z.dist.between(0, 10000)]
    PKT = sub[["symbol", "pchr", "pstart", "pend", "pname", "pscore"]
              ].drop_duplicates("symbol").set_index("symbol")
    # DNAI4 == WDR78: nearest peak to the genuine chr1:67.38 Mb TSS
    t0w = int(TB[TB.symbol == "WDR78"].tss.iloc[0])
    cw = PKALL[PKALL.pchr == "chr1"].drop_duplicates("pname").copy()
    cw["d"] = np.where(t0w < cw.pstart, cw.pstart - t0w,
                       np.where(t0w >= cw.pend, t0w - cw.pend + 1, 0))
    wr = cw.sort_values("d").iloc[0]
    PKT.loc["DNAI4"] = dict(pchr=wr.pchr, pstart=float(wr.pstart),
                            pend=float(wr.pend), pname=wr.pname, pscore=wr.pscore)
    S["PKT"], S["WDR78_dist"] = PKT, int(wr.d)
    T8B = TB[TB.symbol.isin(TGT8 + ["WDR78"])].copy()

    WIN2 = {}
    for gn in TGT8:
        key = "WDR78" if gn == "DNAI4" else gn
        pr = PKT.loc[gn]
        c = pr.pchr
        summit = int((pr.pstart + pr.pend) // 2)
        a, b = summit - 15000, summit + 15000
        own = T8B[T8B.symbol == key]
        if len(own):
            t0 = int(own.tss.iloc[0])
            a, b = min(a, t0 - 3000), max(b, t0 + 3000)
        WIN2[gn] = dict(chrom=c, a=a, b=b, summit=summit,
                        pk=(int(pr.pstart), int(pr.pend)),
                        tss=tss_in(c, a, b).reset_index(drop=True), own=key)
    S["WIN2"], S["LOCI5"] = WIN2, LOCI5

    # ---- panel C: HPA consensus tissues ---------------------------------
    S["tc"] = pd.read_csv(os.path.join(GX, "rna_tissue_consensus.tsv"), sep="\t")
    S["targets"] = ["HYAL1", "APOL1", "APOL2", "ECI1", "CLTB", "GSTO1", "NFU1", "WDR78"]

    # ---- panel D: target expression vs ZNF43 across tissues -------------
    CT = os.path.join(GX, "ConsensusTissusGeneExpression.txt")
    d0 = pd.read_csv(CT, sep="\t", nrows=3)
    use = [c for c in d0.columns if c in ("Gene name", "Tissue", "nTPM")]
    dat = pd.read_csv(CT, sep="\t", usecols=use, low_memory=False)
    dat["nTPM"] = pd.to_numeric(dat["nTPM"], errors="coerce")
    for c in ("Gene name", "Tissue"):
        dat[c] = dat[c].astype(str).str.strip()
    SIG = ["ECI1", "DNAI4", "GSTO1", "APOL1", "APOL2", "HYAL1", "NFU1", "CLTB"]
    a_ = dat[dat["Gene name"].isin(SIG)].groupby("Tissue").nTPM.mean().rename("tgt")
    b_ = dat[dat["Gene name"] == "ZNF43"].groupby("Tissue").nTPM.mean().rename("rep")
    M = pd.concat([a_, b_], axis=1).dropna()
    M["x"] = np.log10(M.rep + 1)
    M["y"] = np.log10(M.tgt + 1)
    S["MD"] = M
    S["rho"], S["pv"] = st.spearmanr(M.x, M.y)
    return S
