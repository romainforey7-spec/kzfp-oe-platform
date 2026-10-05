# Known gaps in this repository

Everything in the paper is reproduced by the code here with two exceptions, both
recorded honestly rather than papered over.

## 1. The data-assembly step for Figure 1 and S1 Fig

`build_fig1.py` and `build_figS1_fig.py` are the drawing halves only. They expect
their data objects to already exist in the namespace they are executed in - 23
names for Figure 1, about 20 for S1 Fig - and they raise `NameError` if imported
on their own. Every other figure has a matching data module (`fig2_data.py`,
`fig3_data.py`, `fig45_data.py`, `figS2_data.py`) that rebuilds its objects from
the primary tables; Figures 1 and S1 do not yet.

What the two builders need is documented in their module docstrings. The
quantities themselves are all in the supporting tables: the per-KZFP screen
scores and the toxicity threshold in S1 Table, the per-well plate readings in
S2 Table, the TE-family enrichment calls in the source tables described below.

## 2. The enrichment test behind the TE-family p-values

The p-values rendered in Fig 1F and S1D come from `fam_data.tsv`. That file
carries a single `padj.final` column and no record of the test that produced it.
It is not equal to any of the three corrected p-values in the subfamily-level
`Fam.tsv` (`padj.hypergeom.reg`, `padj.hypergeom.alafisher`, `padj.binomial`),
and it is not reproducible from its own count columns by an upper-tail binomial
or a Poisson tail under either Benjamini-Hochberg or Bonferroni correction.
S1 Fig panel F states the test as a hypergeometric with BH adjustment, which is
the most defensible reading, but the provenance is not closed.

## Superseded files

`build_figS1.py`, `recover_fig1.py`, `recover_figS1.py`, `dnds_panels.py`,
`target_heatmap.py` and `make_figs.py` are earlier drafts and
one-off recovery scripts. They set type well below the 8 pt floor (down to 3.2 pt)
and some of them write to the same filenames as the current builders. The last
three are imported by no builder in this repository, so nothing delivered depends
on them; this was checked rather than assumed. Each carries a SUPERSEDED
banner in its docstring. They are kept for the record and should not be run.
