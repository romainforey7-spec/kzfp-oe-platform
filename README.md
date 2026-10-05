# Analysis code — Forey et al., *PLOS Genetics*

An overexpression platform reveals the functional diversity of human KRAB-Zinc
Finger Proteins in maintaining cellular homeostasis.

Manuscript PGENETICS-D-26-00836. Corresponding authors: Romain Forey and
Didier Trono, School of Life Sciences, EPFL, Lausanne, Switzerland.

## What this repository does

It regenerates every figure, supporting figure and statistical table in the
paper from the primary data. There is no hidden intermediate state: each figure
has a data module that rebuilds its objects from the source tables and a builder
module that draws the sheet.

## Install

```
conda env create -f environment.yml && conda activate kzfp-oe-platform
# or
pip install -r requirements.txt
```

Python 3.11.16. The bigWig coverage tracks are read by
`bigwig.py`, a pure-standard-library reader, so the repository has **no compiled
dependencies**.

## Layout

| file | what it does |
| --- | --- |
| `plos_export.py` | the only sanctioned way to write a figure. Pins the output to the figure's declared size and emits 600 dpi PNG, 300 dpi PNG, live-text SVG and LZW TIFF. `check_sheet()` measures any file against the PLOS bounds. |
| `bigwig.py` | pure-stdlib bigWig region reader, used for every coverage track |
| `presto_screen.py` | PrestoBlue proliferation screen: raw plate-reader exports to per-construct viability, with the plate layout and all exclusions documented |
| `fig2_data.py`, `fig3_data.py`, `fig45_data.py`, `figS2_data.py` | rebuild each figure's objects from the primary tables |
| `build_fig<N>_fig.py`, `build_figS<N>_fig.py`, `build_kzfp_sup.py` | draw the sheets |
| `rnaseq_panels.py`, `upfig_panels.py`, `target_heatmap.py`, `dnds_panels.py`, `model_fig2.py`, `model_fig3.py` | shared panel builders |
| `moderated_t.py` | empirical-Bayes moderated t-statistics (Smyth 2004) |

`build_figS2S3_fig.py` builds both S2 Fig and S3 Fig; `build_kzfp_sup.py` builds
S6 Fig (ZNF498) and S7 Fig (ZNF18) from one generic sheet.

## Figure typography

Every sheet is 7.48 in wide (PLOS full column) and at most 8.70 in tall. All
body text, axis labels, ticks and legends are 8 pt. Mathematical sub- and
superscripts render at 5.6 pt, which is matplotlib's standard 70 % of the base
size.

Do not call `savefig` directly. If `savefig.bbox` is set to `"tight"` anywhere in
the session, matplotlib ignores the declared figure size and crops to content,
changing the sheet width — and because PLOS scales artwork by width, that scales
the type below the legibility floor. `plos_export.py` clears the setting at
import and passes an explicit bounding box on every write.

## Reproducing the PrestoBlue replicates

```
python presto_screen.py /path/to/screen/archive
```

Writes `PrestoBlue_viability_replicates.csv` and
`PrestoBlue_viability_summary.csv`. Across the screen the recomputed values
agree with the published table at Spearman rho = 0.846 over 815 construct-day
pairs (300 constructs). Note that the archive contains mirrored copies of some
day folders; the script deduplicates on experiment/day/plate/construct, without
which replicate counts double and standard deviations collapse to zero.

## Known gaps

Two are documented in `REPO_GAPS.md`: the data-assembly step for Figure 1 and
S1 Fig is not yet in the repository (their builders are the drawing halves only),
and the enrichment test behind the TE-family p-values in Fig 1F and S1D is not
recorded in the file that supplies them. Six files carry a SUPERSEDED banner in
their own docstring and should not be run: `build_figS1.py`, `recover_fig1.py`,
`recover_figS1.py`, `dnds_panels.py`, `target_heatmap.py` and `make_figs.py`.
`build_figS1.py` is the one to be most careful with - it writes to the same
filenames as the current S1 builder and sets type as small as 3.2 pt.

## Data availability

RNA-seq and ChIP-exo data: GEO, accession to be added on acceptance.
AP-MS data: PRIDE, accession to be added on acceptance.
Protocols: protocols.io, DOI to be added on acceptance.

## License

MIT. See `CITATION.cff` for how to cite.
