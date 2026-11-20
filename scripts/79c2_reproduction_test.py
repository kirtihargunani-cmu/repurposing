#!/usr/bin/env python
"""Acceptance test 1 for D153: do the condition-level profiles reproduce the production artifact?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_condition_level_profiles.md section 6, test 1 - the GATE. If this fails, the new
artifact is wrong and nothing downstream of it is trustworthy.

The test is on the RAW arm. 79c computes one profile per (pert_id, ModelID) as the mean over its
active SIGNATURES. 79c2 computes a profile per condition as the mean over that condition's
signatures. Averaging conditions back is therefore a mean of means, which equals the original only
when WEIGHTED BY n_sig - so the test uses the n_sig-weighted mean, and a plain unweighted mean is
reported alongside to show how much the weighting matters.

The SPECIFIC arm is deliberately NOT expected to match: 79c removes one generic response per line,
79c2 removes one per (line, timepoint). That difference is the point of the change, and its size is
measured here rather than assumed.

Output: results/systems_validation/observability/d153_reproduction_test.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/systems_validation/observability"
TOL = 1e-4


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    gi = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    nG = len(gi)
    pidx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    cidx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    Praw = np.memmap(PROC / "drug_line_profiles_bing.f32", dtype="float32", mode="r",
                     shape=(len(pidx), nG))
    Craw = np.memmap(PROC / "drug_condition_profiles_raw.f32", dtype="float32", mode="r",
                     shape=(len(cidx), nG))
    Pspec = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype="float32", mode="r",
                      shape=(len(pidx), nG))
    Cspec = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype="float32", mode="r",
                      shape=(len(cidx), nG))
    print(f"production {len(pidx):,} profiles | conditions {len(cidx):,} over "
          f"{cidx.groupby(['pert_id','ModelID']).ngroups:,} pairs | {nG:,} genes")

    pkey = pidx.pert_id + "|" + pidx.ModelID
    ckey = cidx.pert_id + "|" + cidx.ModelID
    missing = sorted(set(pkey) - set(ckey))
    extra = sorted(set(ckey) - set(pkey))
    print(f"pairs in production but not in conditions: {len(missing)}  {missing[:3]}")
    print(f"pairs in conditions but not in production: {len(extra)}")

    pos = {k: i for i, k in enumerate(pkey)}
    wmax = umax = 0.0
    wbad = []
    for key, g in cidx.groupby(ckey.values):
        i = pos.get(key)
        if i is None:
            continue
        rows = g.index.values
        w = g.n_sig.values.astype(np.float64)
        blk = np.asarray(Craw[rows], dtype=np.float64)
        wmean = (blk * w[:, None]).sum(0) / w.sum()
        umean = blk.mean(0)
        tgt = np.asarray(Praw[i], dtype=np.float64)
        dw = float(np.abs(wmean - tgt).max())
        du = float(np.abs(umean - tgt).max())
        wmax = max(wmax, dw); umax = max(umax, du)
        if dw > TOL:
            wbad.append((key, round(dw, 6), len(rows)))

    # how much does the per-timepoint generic response change the specific arm?
    sdiff = []
    for key, g in cidx.groupby(ckey.values):
        i = pos.get(key)
        if i is None:
            continue
        rows = g.index.values
        w = g.n_sig.values.astype(np.float64)
        wmean = (np.asarray(Cspec[rows], dtype=np.float64) * w[:, None]).sum(0) / w.sum()
        tgt = np.asarray(Pspec[i], dtype=np.float64)
        d = wmean - tgt
        n = np.linalg.norm(tgt)
        sdiff.append(float(np.linalg.norm(d) / n) if n > 1e-6 else np.nan)
    sdiff = np.array([x for x in sdiff if x == x])

    res = dict(test="D153 acceptance test 1 - reproduction of the production artifact",
               tolerance=TOL,
               n_production_profiles=int(len(pidx)), n_conditions=int(len(cidx)),
               n_pairs_in_conditions=int(cidx.groupby(["pert_id", "ModelID"]).ngroups),
               pairs_missing_from_conditions=missing, pairs_extra=extra,
               raw_weighted_max_abs_diff=round(wmax, 8),
               raw_unweighted_max_abs_diff=round(umax, 8),
               raw_weighted_pairs_over_tolerance=len(wbad), worst_weighted=wbad[:5],
               specific_relative_l2_diff=dict(
                   median=round(float(np.median(sdiff)), 5),
                   p90=round(float(np.percentile(sdiff, 90)), 5),
                   max=round(float(sdiff.max()), 5)),
               passed=bool(wmax <= TOL and not missing))
    (OUT / "d153_reproduction_test.json").write_text(json.dumps(res, indent=1, default=str))
    print(f"\nRAW arm, n_sig-weighted mean of conditions vs production:")
    print(f"   max |difference| = {wmax:.2e}   (tolerance {TOL})  -> "
          f"{'MATCH' if wmax <= TOL else 'MISMATCH'}")
    print(f"   unweighted mean would give {umax:.2e}  (shows why the weighting is needed)")
    print(f"\nSPECIFIC arm, relative L2 difference (expected non-zero: per-timepoint generic response):")
    print(f"   median {res['specific_relative_l2_diff']['median']:.4f}  "
          f"p90 {res['specific_relative_l2_diff']['p90']:.4f}  max {res['specific_relative_l2_diff']['max']:.4f}")
    print(f"\nGATE: {'PASS' if res['passed'] else 'FAIL'}")
    return 0 if res["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
