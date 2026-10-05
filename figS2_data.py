"""Rebuild every object Figure S2 needs, from the primary tables.

Panels: A anti-HA blot (carry-over image) · B PrestoBlue viability ·
C K562 expression density (endogenous vs transgene) · D the same KZFPs across
1,206 cell lines · E peak-to-TSS distance density · F TE-family composition and
an LTR/ERV1 permutation test · G K562 ChIP at the seven target peaks ·
H HPA tissues ranked by ZNF43 · I liver vs thymus histone marks at the seven
peaks · J the DNAI4/WDR78 locus.

Usage:
    import figS2_data as DS2
    S = DS2.load()
"""
import os
import struct
import numpy as np
import pandas as pd
from scipy import stats as st
from bigwig import BigWig
import fig2_data as D2

R = D2.R
AN = D2.AN
GX = D2.GX
LTDIR = D2.LTDIR
HB = os.path.join(R, "20250827__Heatmap_AP_OE_Fig5_Specific_Integrants")
TGK = os.path.join(R, "20260630_OERevision",
                   "20260630_ExpressionFromTransgene_RNAseqApproach",
                   "K562_KZFP_expression_outputs_v4")
REPBED = os.path.join(AN, "2005_hg19_REPEATS_ucsc_merged_s_20170127.bed")

# K562 tracks. The normDMSO_* files are normalised RATIO tracks (genome means
# ~0.002, minima below zero) so they are plotted RAW; only the three protein
# tracks have positive minima and may be scaled by the genome mean.
K562_PROT = {"ZNF43-HA": os.path.join(HB, "ZNF43", "ZNF43_n.bw"),
             "KAP1": os.path.join(HB, "KAP1_K562_rep1.bw"),
             "SETDB1": os.path.join(HB, "SetDB1_K562_rep1.bw")}
K562_HIST = {f"H3K{m}": os.path.join(HB, f"normDMSO_H3K{m}.bw")
             for m in ("9me3", "4me3", "27ac")}
LT_HIST = {"Liver H3K4me3": os.path.join(LTDIR, "liver", "Liver_H3K4me3.bigWig"),
           "Thymus H3K4me3": os.path.join(LTDIR, "Thymus", "Thymus_H3K4me3.bigWig"),
           "Liver H3K9me3": os.path.join(LTDIR, "liver", "Liver_H3K9me3.bigWig"),
           "Thymus H3K9me3": os.path.join(LTDIR, "Thymus", "Thymus_H3K9me3.bigWig")}

TGT8 = ["HYAL1", "APOL1", "APOL2", "ECI1", "CLTB", "GSTO1", "NFU1", "DNAI4"]
SCREEN_VIAB = {"D4": 0.64, "D7": 0.49, "D9": 0.33}     # ZNF43, OE screen table

# ---- B: viability replicates -----------------------------------------
# The published screen table gives one value per day with no dispersion.
# presto_screen.py recovers the individual plates from the raw reader exports;
# this reads its output so the figure and the supporting table can never
# disagree. If the CSV is absent the panel falls back to the published point
# estimates with no error bars.
PRESTO_REPLICATES = "PrestoBlue_viability_replicates.csv"


def viability_replicates(kzfp="ZNF43", path=PRESTO_REPLICATES):
    """(viab_rep, viab_wells) for one bait, from presto_screen.py output."""
    if not os.path.exists(path):
        rep = pd.DataFrame([dict(day=int(d.lstrip("D")), mean=v, sd=np.nan, n_plates=1)
                            for d, v in SCREEN_VIAB.items()])
        return rep, pd.DataFrame(columns=["day", "value"])
    R = pd.read_csv(path)
    R = R[R["KZFP"].astype(str).str.upper() == kzfp.upper()].copy()
    # the builder plots days on a numeric axis, so "D4" -> 4
    R["day"] = R["day"].astype(str).str.extract(r"(\d+)").astype(int)
    rep = (R.groupby("day")["viab"]
             .agg(mean="mean", sd="std", n="size").reset_index())
    wells = R.rename(columns={"viab": "value"})[["day", "value"]]
    return rep, wells
NPERM = 10000


def _tsum(bw):
    bw.f.seek(bw.totalSummaryOffset)
    vc, mn, mx, sd, ssq = struct.unpack(bw.e + "Qdddd", bw.f.read(40))
    return sd / vc


def target_peaks(Z, PKT):
    """The seven unique peaks carrying the eight DOWN targets (<=10 kb)."""
    P = PKT.reset_index()[["symbol", "pchr", "pstart", "pend", "pname", "pscore"]].copy()
    P["summit"] = ((P.pstart + P.pend) / 2).astype(int)
    grp = P.groupby("pname").agg(chrom=("pchr", "first"), summit=("summit", "first"),
                                 start=("pstart", "first"), end=("pend", "first"),
                                 score=("pscore", "first"),
                                 genes=("symbol", lambda s: "/".join(sorted(s))))
    return grp.reset_index().sort_values("pname").reset_index(drop=True)


def repeat_class_fractions():
    """Genome-wide bp fraction per repeat class, from the merged RepeatMasker bed."""
    tot = {}
    for ch in pd.read_csv(REPBED, sep="\t", header=None, usecols=[1, 2, 6],
                          names=["s", "e", "cls"], chunksize=1_000_000):
        ch["bp"] = ch.e - ch.s
        for k, v in ch.groupby("cls").bp.sum().items():
            tot[k] = tot.get(k, 0) + int(v)
    S = pd.Series(tot).sort_values(ascending=False)
    return S / S.sum()


def perm_ltr_erv1(Z, n_target_peaks, frac_erv1, nperm=NPERM, seed=0):
    """How many of n random repeat-overlapping peaks would fall on LTR/ERV1."""
    rng = np.random.default_rng(seed)
    null = rng.binomial(n_target_peaks, frac_erv1, size=nperm)
    return null


def load():
    S = dict(D2.load())                       # reuse the ZNF43 objects
    Z, PKT = S["Z43"], S["PKT"]
    S["PK7"] = target_peaks(Z, PKT)

    # ---- B: viability ---------------------------------------------------
    S["viab"] = pd.Series({"D0": 1.0, **SCREEN_VIAB})
    S["viab_rep"], S["viab_wells"] = viability_replicates("ZNF43")

    # ---- C / D: transgene vs endogenous expression ----------------------
    S["SEL"] = pd.read_csv(os.path.join(TGK, "selected_KZFP_endogenous_and_transgene_values.csv"))
    S["ALLG"] = pd.read_csv(os.path.join(TGK, "K562_local_all_genes.csv"))
    S["KZF"] = pd.read_csv(os.path.join(TGK, "K562_local_KZFPs.csv"))
    S["CL"] = pd.read_csv(os.path.join(TGK, "selected_KZFPs_across_cell_lines_normalized.csv"))

    # ---- E: peak-to-TSS distances ---------------------------------------
    g = S["G43"]
    S["dist_all"] = g.dist.dropna().values
    S["dist_down"] = g[g.cat == "DOWN"].dist.dropna().values

    # ---- F: repeat composition + binomial enrichment ---------------------
    # composition of the DOWN & 0-10 kb gene/peak/repeat rows (small n: report
    # counts, not percentages, and state n on the panel)
    dn = Z[(Z.Direction == "DOWN") & Z.dist.between(0, 10000)
           & Z.rclass.notna() & (Z.rclass != ".")]
    S["dn_rows"] = int(len(dn))
    S["dn_class_counts"] = dn.rclass.value_counts()
    # peaks on LTR/ERV1: unique peak names carrying at least one LTR/ERV1 repeat
    allpk = Z[Z.pname.notna() & (Z.pname != ".")]
    rp = allpk[allpk.rname.notna() & (allpk.rname != ".")]
    S["n_peaks"] = int(allpk.pname.nunique())
    S["n_peaks_rep"] = int(rp.pname.nunique())
    S["erv1_obs"] = int(rp.loc[rp.rclass == "LTR/ERV1", "pname"].nunique())
    S["genome_frac"] = repeat_class_fractions()
    S["erv1_exp_frac"] = float(S["genome_frac"]["LTR/ERV1"])
    S["erv1_exp"] = S["n_peaks_rep"] * S["erv1_exp_frac"]
    S["erv1_fold"] = S["erv1_obs"] / S["erv1_exp"]
    S["erv1_p"] = float(st.binomtest(S["erv1_obs"], S["n_peaks_rep"],
                                     S["erv1_exp_frac"], alternative="greater").pvalue)

    # ---- G: K562 tracks at the seven peaks -------------------------------
    BWP = {k: BigWig(v) for k, v in K562_PROT.items()}
    S["BW_prot"], S["SUM_prot"] = BWP, {k: _tsum(b) for k, b in BWP.items()}
    S["BW_hist"] = {k: BigWig(v) for k, v in K562_HIST.items()}

    # ---- H: HPA tissues ranked by ZNF43 ----------------------------------
    tc = S["tc"]
    z = tc[tc["Gene name"] == "ZNF43"].set_index("Tissue").nTPM
    HPA_TGT = ["HYAL1", "APOL1", "APOL2", "ECI1", "CLTB", "GSTO1", "NFU1", "WDR78"]
    tgt = tc[tc["Gene name"].isin(HPA_TGT)]
    H = pd.DataFrame({"znf43": z})
    H["tgt_mean"] = np.log(tgt.groupby("Tissue").nTPM.mean() + 1)
    H = H.dropna().sort_values("znf43", ascending=False)
    S["H_rank"] = H
    S["H_tgt"] = tgt

    # ---- I: liver vs thymus at the seven peaks ---------------------------
    BWL = {k: BigWig(v) for k, v in LT_HIST.items()}
    S["BW_lt"], S["SUM_lt"] = BWL, {k: _tsum(b) for k, b in BWL.items()}
    return S
