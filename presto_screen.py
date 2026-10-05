"""PrestoBlue proliferation screen: reprocessing from the raw plate-reader exports.

Why this script exists
----------------------
The published screen table carries one viability value per KZFP per time point with
no dispersion. Reviewer 1 asked for replicates. The per-well readings needed to
recover them were still on disk, in two forms, and this script reads both:

* ``sumup_presto_*_D<n>.xlsx`` -- one workbook per experiment per day. Each
  incubation sheet (1h / 3h / 4h) holds the plate's optical densities AND its own
  plate map, written as a three-column block (row letter, column number,
  construct name) starting at spreadsheet column P. The block is repeated at
  +20, +40 and +60 columns for plates 2, 3 and 4 of the same batch.
* ``*_<1h|3h|4h>_plate<n>.xlsx`` -- the plate reader's own export, used for the
  plates whose sumup sheet was left blank. The readings are on the sheet
  ``Results by well`` under seven header rows, in the columns ``Well``,
  ``Raw OD(570)`` and ``Raw OD(600)``.

Plate layout (decoded from the formula chain in the sumup workbooks, then checked
against the stored blank value in cell P29 of one sheet, which it reproduces
exactly)::

    columns  2-11   one construct each; column 7 is the GFP normaliser
    rows     B,C,D  no-doxycycline triplicate
    rows     E,F,G  doxycycline triplicate
    rows     A,H    medium-only blanks

Viability is computed on the 3 h read as::

    od    = OD570 - OD600
    blank = mean(od) over rows A and H of the same plate
    viab  = mean(od[E,F,G] - blank) / mean(od[B,C,D] - blank)

i.e. each construct is its own no-Dox reference and GFP is NOT divided into the
ratio -- the sumup formula chain applies the GFP normalisation outside this
quantity, not inside it.

Validation. Recomputed values agree with the published screen table at
Spearman rho = 0.846 over 815 construct-day pairs covering 300 constructs. The
four baits of the paper come out at ZNF43 0.644/0.485/0.348, ZNF498
0.766/0.298/0.114, ZNF18 0.211/0.100/0.045 and ZNF257 0.826/0.581/0.257,
against published 0.64/0.49/0.33, 0.72/0.31/0.13, 0.20/0.10/0.04 and
0.88/0.74/0.40. ZNF257 is the one series that does not reproduce; no incubation
time and no normalisation recovers the published 0.74 and 0.40.

Known limitations, carried into the paper
-----------------------------------------
* Experiment CR019 is excluded: it correlates with the published table at
  rho = 0.166 over its 20 constructs, i.e. its plate map cannot be trusted.
* The repeat run ``CR015_2`` is excluded: its layout could not be verified, and
  under the layout above it returns ZNF18 viabilities of 0.94/0.86/0.69 against
  0.21/0.10/0.04 from CR015 itself.
* Per-experiment agreement with the published table ranges from rho = 0.944
  (CR021) down to 0.763 (CR031); see PrestoBlue_per_experiment_agreement.csv.
  CR036 (rho = 0.853) and CR037 (rho = 0.809) rest on only 6 and 8 constructs
  respectively, so their coefficients are weakly determined; CR036 nonetheless
  contributes three of the four plates behind the ZNF43 numbers.
* The archive holds mirrored copies of several day folders, so the same plate
  export is reachable by more than one path. process() deduplicates on
  experiment/day/plate/construct; without that, replicate counts double and
  every standard deviation collapses to zero.
"""
import os
import re
import glob
import numpy as np
import pandas as pd
import openpyxl

ROWS_NODOX = ["B", "C", "D"]
ROWS_DOX = ["E", "F", "G"]
ROWS_BLANK = ["A", "H"]
CONSTRUCT_COLS = list(range(2, 12))
GFP_COL = 7
PLATE_OFFSETS = [0, 20, 40, 60]          # map block for plates 1-4
MAP_COL0 = openpyxl.utils.column_index_from_string("R")  # first name block
MAP_ROWS = range(5, 15)                  # ten names, one per assay-plate column
EXCLUDE_EXPERIMENTS = {"CR019"}          # plate map not trustworthy (rho = 0.166)
EXCLUDE_RUNS = {"CR015_2"}               # unverifiable layout


def read_plate_map(ws, offset):
    """Construct names for one assay plate, as {assay_plate_column: name}.

    Each incubation sheet carries four name blocks, one per assay plate, in
    spreadsheet column R and at +20, +40 and +60 columns from it. A block is ten
    consecutive cells (rows 5-14) listing the constructs in assay-plate column
    order, so the first maps to plate column 2 and the last to plate column 11.
    The LacZ control is always first and the GFP normaliser always sixth, which
    puts GFP in plate column 7 -- the position the sumup formulas divide by.
    Wells recorded as "dead" carried no viable culture and are dropped.
    """
    col = MAP_COL0 + offset
    out = {}
    for i, r in enumerate(MAP_ROWS):
        v = ws.cell(row=r, column=col).value
        if v is None:
            continue
        name = str(v).strip()
        if not name or name.lower() == "dead":
            continue
        out[CONSTRUCT_COLS[i]] = name
    return out


def read_raw_plate(path):
    """One plate-reader export -> {(row_letter, column): OD570 - OD600}."""
    df = pd.read_excel(path, sheet_name="Results by well", header=7)
    cols = {c.strip(): c for c in df.columns if isinstance(c, str)}
    w, c570, c600 = cols["Well"], cols["Raw OD(570)"], cols["Raw OD(600)"]
    grid = {}
    for _, r in df.iterrows():
        m = re.match(r"^([A-H])\s*0?(\d{1,2})$", str(r[w]).strip())
        if not m:
            continue
        od570, od600 = pd.to_numeric(r[c570], errors="coerce"), pd.to_numeric(r[c600], errors="coerce")
        if pd.isna(od570) or pd.isna(od600):
            continue
        grid[(m.group(1), int(m.group(2)))] = float(od570) - float(od600)
    return grid


def viability(grid, names_by_col):
    """Blank-corrected Dox / no-Dox ratio for every construct on one plate."""
    blanks = [grid[(r, c)] for r in ROWS_BLANK for c in CONSTRUCT_COLS if (r, c) in grid]
    if not blanks:
        return {}
    blank = float(np.mean(blanks))
    out = {}
    for col, name in names_by_col.items():
        if col not in CONSTRUCT_COLS:
            continue
        nodox = [grid[(r, col)] - blank for r in ROWS_NODOX if (r, col) in grid]
        dox = [grid[(r, col)] - blank for r in ROWS_DOX if (r, col) in grid]
        if len(nodox) < 2 or len(dox) < 2:
            continue
        denom = float(np.mean(nodox))
        if denom <= 0:
            continue
        out[name] = dict(viab=float(np.mean(dox)) / denom,
                         n_nodox=len(nodox), n_dox=len(dox),
                         is_gfp=(col == GFP_COL), plate_column=col)
    return out


def parse_run(fname):
    """('CR036', 'D4', '3h', 2) from a raw-export filename, or Nones."""
    m = re.search(r"(CR\d+)_?(D\d+)?_(1h|3h|4h)_plate(\d)", fname, re.I)
    if not m:
        return (None, None, None, None)
    return (m.group(1).upper(), (m.group(2) or "").upper(), m.group(3).lower(), int(m.group(4)))


def _sumup_for(path):
    """The sumup workbook that describes a raw export: same folder, else parent."""
    d = os.path.dirname(path)
    for cand in (d, os.path.dirname(d)):
        g = sorted(glob.glob(os.path.join(cand, "sumup_presto*.xls*")))
        if g:
            return g[0]
    return None


def maps_from_sumup(path, incubation="3h"):
    """{plate_index: {assay_column: construct}} for one sumup workbook."""
    wb = openpyxl.load_workbook(path, data_only=True)
    try:
        sheets = [w for w in wb.worksheets if re.fullmatch(rf"\s*{incubation}\s*", w.title, re.I)]
        if not sheets:
            sheets = [w for w in wb.worksheets if re.search(incubation, w.title, re.I)]
        if not sheets:
            return {}
        ws = sheets[0]
        return {i: read_plate_map(ws, off)
                for i, off in enumerate(PLATE_OFFSETS, start=1)}
    finally:
        wb.close()


def process(root, incubation="3h"):
    """Every raw plate export under *root* -> one row per construct per plate.

    Each export is paired with the sumup workbook in its own folder, so the
    construct names can never be taken from a different experiment or day.
    """
    rows = []
    cache = {}
    for path in sorted(glob.glob(os.path.join(root, "**", "*plate?.xls*"), recursive=True)):
        base = os.path.basename(path)
        if base.startswith("~$") or base.lower().startswith("sumup"):
            continue
        exp, day, inc, plate = parse_run(base)
        if exp is None or inc != incubation or exp in EXCLUDE_EXPERIMENTS:
            continue
        if any(bad in base for bad in EXCLUDE_RUNS):
            continue
        sp = _sumup_for(path)
        if sp is None:
            continue
        if sp not in cache:
            try:
                cache[sp] = maps_from_sumup(sp, incubation)
            except Exception:
                cache[sp] = {}
        names = cache[sp].get(plate)
        if not names:
            continue
        try:
            grid = read_raw_plate(path)
        except Exception:
            continue
        for name, v in viability(grid, names).items():
            rows.append(dict(experiment=exp, day=day or "", plate=plate, incubation=inc,
                             construct=name, **v, source=base, sumup=os.path.basename(sp)))
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # The archive holds mirrored copies of several day folders, so the same plate
    # export is reachable by more than one path. Without this the replicate count
    # and the standard deviation are both wrong (n doubles, SD collapses to 0).
    return df.drop_duplicates(subset=["experiment", "day", "plate", "construct"],
                              keep="first").reset_index(drop=True)


def summarise(rep):
    """Mean, SD and plate count per construct per day."""
    g = rep.groupby(["construct", "day"])["viab"]
    out = g.agg(mean="mean", sd="std", n_plates="size").reset_index()
    return out.pivot(index="construct", columns="day",
                     values=["mean", "sd", "n_plates"]).reset_index()


PUBLISHED_SERIES = {                       # the values in the published screen table
    "ZNF43":   {4: 0.64, 7: 0.49, 9: 0.33},
    "ZNF257":  {4: 0.88, 7: 0.74, 9: 0.40},
    "ZNF498":  {4: 0.72, 7: 0.31, 9: 0.13},
    "ZNF18":   {4: 0.20, 7: 0.10, 9: 0.04},
    "ZSCAN25": {4: 0.36, 7: 0.10, 9: 0.05},
}


def series_for(kzfp, path="PrestoBlue_viability_replicates.csv"):
    """Viability series for one bait: day, mean, sd, n, and the per-plate values.

    Reads the output of this script so every figure panel and the supporting
    table come from one source. If the CSV is missing it returns the published
    point estimates with no dispersion, so a figure still builds -- but the
    ``recomputed`` flag on the returned frame is then False.
    """
    if os.path.exists(path):
        R = pd.read_csv(path)
        R = R[R["KZFP"].astype(str).str.upper() == kzfp.upper()].copy()
        if len(R):
            R["day"] = R["day"].astype(str).str.extract(r"(\d+)").astype(int)
            rep = (R.groupby("day")["viab"]
                     .agg(mean="mean", sd="std", n="size").reset_index())
            rep.attrs["recomputed"] = True
            return rep, R.rename(columns={"viab": "value"})[["day", "value"]]
    rep = pd.DataFrame([dict(day=d, mean=v, sd=np.nan, n=1)
                        for d, v in sorted(PUBLISHED_SERIES.get(kzfp.upper(), {}).items())])
    rep.attrs["recomputed"] = False
    return rep, pd.DataFrame(columns=["day", "value"])


if __name__ == "__main__":
    import sys
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    rep = process(root)
    rep.to_csv("PrestoBlue_viability_replicates.csv", index=False)
    summarise(rep).to_csv("PrestoBlue_viability_summary.csv", index=False)
    print(f"{len(rep)} construct-plate measurements from "
          f"{rep['source'].nunique()} plate exports across {rep['experiment'].nunique()} experiments")
