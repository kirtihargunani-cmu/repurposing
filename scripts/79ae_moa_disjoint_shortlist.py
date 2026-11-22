#!/usr/bin/env python
"""MOA-disjoint lead shortlist from the D184 concordance output (D185, 2026-11-22).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

No new statistic. This formalises the filter that was first run by hand, so the shortlist is
rebuilt by the chain instead of being a one-off query against `79ad`'s output.

A (cohort, window, drug) row qualifies when ALL of:
  * set == "T"               top decile by the window's own reversal percentile, and
                             label_phase == 0 for THAT cohort (79ad amendment 2/13)
  * moa_disjoint             annotated, and sharing no `evidence_targets` entry and no
                             `hub_moa` class with any member of the approved set A
  * conc  >  that cell's own median conc over N    (a per-cell bar, never a global one:
                                                    null medians run 0.02 - 0.38)
  * conc_pc1removed > 0      survives removing the cohort's leading shared component, so the
                             row is not carried by the generic proliferation-like axis (D184 s7)

`approved_other_cohort` flags rows whose drug carries label_phase >= 3 in SOME OTHER cohort of
the table. These are NOT discoveries: per-cohort labelling correctly leaves tamoxifen unlabelled
in BRCA-Basal (it is not indicated in basal/TNBC) and lapatinib unlabelled outside Her2, so both
enter T. The flag exists so they are never read as novel.

THIS IS A LEAD LIST, NOT A RESULT. D184 licenses the claim that T as a SET converges on the
approved profile. It licenses nothing about any individual row here, each of which is one cosine
against a centroid built from a median of 4-10 approved drugs.

Output: results/lincs/approved_concordance/moa_disjoint_shortlist.tsv
"""
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
DIR = RES / "approved_concordance"
SRC = DIR / "per_drug.tsv"
OUT = DIR / "moa_disjoint_shortlist.tsv"

COLS = ["cohort", "window", "cmap_name", "pert_id", "conc", "conc_pc1removed",
        "n_med", "conc_over_null", "rev", "approved_other_cohort"]


def main() -> int:
    if not SRC.exists():
        print(f"ABORT: {SRC} missing - run scripts/79ad_approved_concordance.py first.", file=sys.stderr)
        return 1
    d = pd.read_csv(SRC, sep="\t")

    # the highest clinical phase this drug reaches in ANY cohort of the table
    max_phase = d.groupby("cmap_name").label_phase.max().rename("max_phase_any")
    # the per-cell bar: the median concordance of the matched null set N in that same cell
    n_med = d[d.set == "N"].groupby(["cohort", "window"]).conc.median().rename("n_med")

    t = (d[(d.set == "T") & d.moa_disjoint]
         .join(n_med, on=["cohort", "window"])
         .join(max_phase, on="cmap_name"))

    r = t[(t.conc > t.n_med) & (t.conc_pc1removed > 0)].copy()
    r["approved_other_cohort"] = r.max_phase_any >= 3
    r["conc_over_null"] = r.conc - r.n_med
    r = r.sort_values(["cohort", "window", "conc"], ascending=[True, True, False])

    DIR.mkdir(parents=True, exist_ok=True)
    r[COLS].to_csv(OUT, sep="\t", index=False, float_format="%.4f")

    n_cells = r.groupby(["cohort", "window"]).ngroups
    print(f"{len(r)} of {len(t)} MOA-disjoint T rows survive, across {n_cells} (cohort, window) cells")
    print(f"  approved for some other cohort (flagged, not novel): {int(r.approved_other_cohort.sum())}")
    print(f"  wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
