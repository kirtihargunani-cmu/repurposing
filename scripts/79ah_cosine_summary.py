#!/usr/bin/env python
"""Reduce the granular contribution cosine to what a page can carry (D194).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

`79e3 COSINE=1` writes one row per (cohort, drug, window, REGULATOR) - 18.5M rows, 501 MB. That
table is the granular view and stays on disk; this reduces it to two small ones so the window
ranking page can embed them.

    cos_contrib   within one regulon, cosine between the gene-level DISEASE contributions (50b)
                  and the gene-level DRUG contributions (79e3). Both already carry each edge's own
                  tfmode x likelihood, so the per-gene weights are the pipeline's own, not rebuilt
                  from z-scores. NEGATIVE = the drug moves the genes that account for that
                  regulator's dysregulation, in the opposite direction, gene by gene.

Per drug the regulators are combined |NES|-weighted, matching how the reversal score weights them,
so a regulator that barely moves the disease cannot dominate the summary.

Outputs:
  results/lincs/contribution_cosine_per_drug.tsv.gz   cohort, window, pert_id, cmap_name,
                                                      n_regulators, cos_mean, cos_weighted,
                                                      frac_opposing
  results/lincs/contribution_cosine_summary.tsv       per (cohort, window) medians
"""
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SRC = RES / "regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"
OUT_D = RES / "contribution_cosine_per_drug.tsv.gz"
OUT_S = RES / "contribution_cosine_summary.tsv"
COLS = ["cohort", "pert_id", "cmap_name", "window", "regulator", "disease_nes", "cos_contrib"]
CHUNK = 2_000_000


def main() -> int:
    if not SRC.exists():
        print(f"ABORT: {SRC} missing - run 79e3 with COSINE=1 ALL=1 first.", file=sys.stderr)
        return 1
    parts, n_in = [], 0
    for ch in pd.read_csv(SRC, sep="\t", usecols=lambda c: c in COLS, chunksize=CHUNK,
                          low_memory=False):
        if "cos_contrib" not in ch.columns:
            print("ABORT: cos_contrib absent - the source was built without COSINE=1.", file=sys.stderr)
            return 1
        n_in += len(ch)
        ch = ch.dropna(subset=["cos_contrib"])
        ch["w"] = ch.disease_nes.abs().fillna(0.0)
        ch["cw"] = ch.cos_contrib * ch.w
        ch["neg"] = (ch.cos_contrib < 0).astype(float)
        parts.append(ch.groupby(["cohort", "pert_id", "cmap_name", "window"], sort=False)
                     .agg(n=("cos_contrib", "size"), s=("cos_contrib", "sum"),
                          cw=("cw", "sum"), w=("w", "sum"), neg=("neg", "sum")).reset_index())
        print(f"  {n_in:,} rows read", flush=True)

    # chunk boundaries can split a drug, so re-aggregate the partial sums
    A = (pd.concat(parts, ignore_index=True)
         .groupby(["cohort", "pert_id", "cmap_name", "window"], as_index=False)
         [["n", "s", "cw", "w", "neg"]].sum())
    A["cos_mean"] = A.s / A.n
    A["cos_weighted"] = A.cw / A.w.clip(lower=1e-9)
    A["frac_opposing"] = A.neg / A.n
    A = A.rename(columns={"n": "n_regulators"}).drop(columns=["s", "cw", "w", "neg"])
    A.to_csv(OUT_D, sep="\t", index=False, float_format="%.6g")

    S = (A.groupby(["cohort", "window"])
         .agg(n_drugs=("pert_id", "size"), n_regulators=("n_regulators", "median"),
              cos_weighted_median=("cos_weighted", "median"),
              cos_mean_median=("cos_mean", "median"),
              frac_opposing_median=("frac_opposing", "median")).round(5).reset_index())
    S.to_csv(OUT_S, sep="\t", index=False)

    print(f"\n{n_in:,} granular rows -> {len(A):,} (cohort, drug, window) -> {len(S)} (cohort, window)")
    print(f"wrote {OUT_D}\n      {OUT_S}\n")
    print(S.sort_values("cos_weighted_median").head(6).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
