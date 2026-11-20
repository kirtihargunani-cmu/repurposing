#!/usr/bin/env python
"""Positive-hit view of the per-window reversal scores (user, 2026-11-20).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

Two decisions from the user, both narrowing scope deliberately:
  * compare EARLY as one group and LATE as another. No dose-response progression is modelled -
    two buckets, nothing fitted.
  * lead with POSITIVE HITS and flag everything else as a consideration rather than a filter.
    A drug that reverses a disease at 24 h but not 6 h is a hit with a timing note, not a reject.

So a drug is ranked by its BEST window, and the timing is carried beside it:
    best_window / best_pct      the hit, and which window produced it
    other_pct                   what the other window said (NaN if absent)
    timing_note                 the consideration to flag

Input: window_reversal.tsv.gz.  Output: window_hits.tsv.gz + per-cohort top table.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/integrated"


def main():
    d = pd.read_csv(OUT / "window_reversal.tsv.gz", sep="\t")
    for w in ("early", "late"):
        d[f"pct_{w}"] = d.groupby("cohort")[f"rev_{w}"].rank(pct=True) * 100
    e, l = d.pct_early, d.pct_late
    d["best_window"] = np.where(e.fillna(-1) >= l.fillna(-1), "early", "late")
    d["best_pct"] = np.where(d.best_window == "early", e, l)
    d["other_pct"] = np.where(d.best_window == "early", l, e)
    d["best_rev"] = np.where(d.best_window == "early", d.rev_early, d.rev_late)
    d["n_conditions_best"] = np.where(d.best_window == "early",
                                      d.n_conditions_early, d.n_conditions_late)
    both = d.window_coverage == "both"
    d["timing_note"] = np.select(
        [d.window_coverage == "early_only",
         d.window_coverage == "late_only",
         both & (d.window_agreement == "discordant"),
         both & (d.best_pct - d.other_pct >= 25),
         both],
        ["early evidence only (no 24 h data)",
         "late evidence only (no 6 h data)",
         "REVERSES between 6 h and 24 h",
         "much stronger in the %s window" % "best",
         "consistent across both windows"], default="")
    d.loc[both & (d.window_agreement == "concordant") & (d.best_pct - d.other_pct >= 25),
          "timing_note"] = ("stronger " + d.best_window + ", same direction in both")
    d = d.sort_values(["cohort", "best_pct"], ascending=[True, False])
    d.to_csv(OUT / "window_hits.tsv.gz", sep="\t", index=False)

    hits = d[d.best_pct >= 90]
    summ = dict(step="positive-hit view of per-window reversal", n_rows=int(len(d)),
                n_cohorts=int(d.cohort.nunique()),
                n_hits_top10pct=int(len(hits)),
                hits_by_best_window=dict(hits.best_window.value_counts().items()),
                hits_by_coverage=dict(hits.window_coverage.value_counts().items()),
                hits_reversing=int((hits.timing_note == "REVERSES between 6 h and 24 h").sum()),
                median_conditions_behind_a_hit=float(hits.n_conditions_best.median()))
    (OUT / "window_hits_run.json").write_text(json.dumps(summ, indent=1, default=str))
    print(f"{len(d):,} (cohort, drug) rows; {len(hits):,} hits at >= 90th percentile\n")
    print(f"  hits by best window : {summ['hits_by_best_window']}")
    print(f"  hits by coverage    : {summ['hits_by_coverage']}")
    print(f"  hits that REVERSE between windows: {summ['hits_reversing']} "
          f"({100*summ['hits_reversing']/max(len(hits),1):.1f}% of hits)")
    print(f"  median conditions behind a hit: {summ['median_conditions_behind_a_hit']:.0f}\n")
    print("Example — BRCA-LumA top 8 by best window:")
    ex = d[d.cohort == "BRCA-LumA"].head(8)
    print(ex[["cmap_name", "best_window", "best_pct", "other_pct", "n_conditions_best",
              "timing_note"]].to_string(index=False, float_format=lambda x: f"{x:6.1f}"))


if __name__ == "__main__":
    sys.exit(main())
