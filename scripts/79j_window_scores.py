#!/usr/bin/env python
"""Carry the per-window (early / late) reversal scores into the ranking tables (D170).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-20: *"we should return both scores separately for a drug that has both."*

WHY. `79c` builds one profile per drug x cell line as the mean of that drug's active signatures
**over doses AND times**, so every drug collapses to one score per disease. For some drugs the two
windows say opposite things and the mean means neither. D152 measured this for anti-estrogens: in
MCF7 the ESR1 programme reads -0.056 at 6 h (correct) and +0.094 at 24 h (inverted), and after the
D156 floor change fulvestrant in BRCA-LumA scores **96.3rd percentile early against 60.2nd late**,
clomifene **97.9 against 12.9**. The window-separated numbers already exist
(`results/systems_validation/integrated/window_reversal.tsv.gz`, D155) and stop in a side table.
This script is the plumbing, not a method change: no score is recomputed, nothing is replaced.

WHAT IT ADDS to the split tables, in place and idempotently:
  rev_early, rev_late          the per-window reversal, as computed by D155 (79e's formula)
  rev_early_pct, rev_late_pct  within-cohort percentile of each, over the drugs that HAVE that
                               window, so the two are each read against their own population
  window_coverage              early_only | late_only | both | (blank = no window data)
  window_agreement             concordant | discordant | not_evaluable
  best_window                  for `both` rows only: which window scores higher, and by how much
  window_gap_pct               rev_early_pct - rev_late_pct (blank unless `both`)

HONEST LIMITS, stated here because the columns invite over-reading:
  * The window table covers **25 of 41 cohorts**. 50,187 of 64,094 split rows get no window data at
    all; those keep blank columns and must not be read as "no window effect".
  * Of the rows that have both windows, **37.4% are discordant** (D155/D156). That is a large
    minority, so a single pooled score is genuinely ambiguous for them - which is the point - but it
    also means a window flag is not a rare special case to be hand-waved.
  * `rev_early_pct` and `rev_late_pct` are percentiles within *different* populations (the drugs
    having that window). They are comparable to each other as rankings, not as magnitudes.
  * Median conditions per window is 1 (D155 run.json), so most window scores rest on a single
    condition. Treat a large gap as a flag for inspection, not as a measured effect size.

Env: SPLIT_SFX (default _v2_tf) selects the tables. Reruns are safe: the columns are recomputed.
"""
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
RES = ROOT / "results/lincs"
WIN = ROOT / "results/systems_validation/integrated/window_reversal.tsv.gz"
SFX = os.environ.get("SPLIT_SFX", "_v2_tf")
COLS = ["rev_early", "rev_late", "rev_early_pct", "rev_late_pct", "window_coverage",
        "window_agreement", "best_window", "window_gap_pct"]


def annotate(path, unit):
    if not path.exists():
        print(f"  skip (absent): {path.name}")
        return None
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d.drop(columns=[c for c in COLS if c in d.columns])          # idempotent
    w = pd.read_csv(WIN, sep="\t", low_memory=False).rename(columns={"cohort": unit})
    w = w[[unit, "pert_id", "rev_early", "rev_late", "window_coverage", "window_agreement"]]
    w = w.drop_duplicates([unit, "pert_id"])
    d = d.merge(w, on=[unit, "pert_id"], how="left")

    # percentile within cohort, over the drugs that HAVE that window
    for side in ("early", "late"):
        col, pct = f"rev_{side}", f"rev_{side}_pct"
        d[pct] = (d.groupby(unit)[col]
                  .transform(lambda v: v.rank(method="max") / v.notna().sum() * 100
                             if v.notna().any() else v))
        d[pct] = d[pct].round(1)

    both = d.window_coverage.eq("both")
    d["best_window"] = pd.NA
    d.loc[both, "best_window"] = (d.loc[both, "rev_early_pct"] >= d.loc[both, "rev_late_pct"]
                                  ).map({True: "early", False: "late"})
    d["window_gap_pct"] = pd.NA
    d.loc[both, "window_gap_pct"] = (d.loc[both, "rev_early_pct"] - d.loc[both, "rev_late_pct"]).round(1)
    d.to_csv(path, sep="\t", index=False)

    n = len(d)
    cov = d.window_coverage.value_counts(dropna=False)
    print(f"  {path.name}: {n:,} rows")
    for k in ("both", "early_only", "late_only"):
        print(f"      {k:<11} {int(cov.get(k, 0)):>7,}")
    print(f"      no window  {int(d.window_coverage.isna().sum()):>7,}  "
          f"({100 * d.window_coverage.isna().mean():.1f}% - 16 of 41 cohorts are not in the window table)")
    if both.any():
        g = d.loc[both, "window_gap_pct"].astype(float)
        dis = d.loc[both, "window_agreement"].eq("discordant")
        print(f"      of 'both': {dis.sum():,} discordant ({100 * dis.mean():.1f}%); "
              f"best_window early {int((d.loc[both,'best_window'] == 'early').sum()):,} / "
              f"late {int((d.loc[both,'best_window'] == 'late').sum()):,}; "
              f"|gap| > 25 pct pts: {(g.abs() > 25).sum():,}")
    return d


def main():
    print(f"79j: carrying per-window scores into the{SFX} split tables (D170)")
    if not WIN.exists():
        print(f"  MISSING {WIN} - run the D155 window step first", file=sys.stderr)
        return 1
    annotate(RES / f"drug_reversal_split{SFX}.tsv", "group_id")
    annotate(RES / f"stratum_reversal_split{SFX}.tsv", "stratum")
    print("  nothing recomputed; no existing column changed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
