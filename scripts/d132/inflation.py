#!/usr/bin/env python
"""D132 step 3: the per-regulator inflation factor, and the corrected activity sandbox.

inflation_r = the SD of regulator r's aREA NES across the null profiles. The analytic normalisation 79d relies on
(`viper(..., method = "none")` with no `nullmodel`, dividing by sqrt(sum likelihood^2)) is correct exactly when
that SD is 1, so the SD it actually has IS the factor by which the assumption is wrong - the same logic D131 used
for `rev`, but here it is a property of each REGULON and therefore varies by regulator.

Writes:
  results/lincs/d132/inflation_<null>.tsv          per (set, net, regulator): mean, sd, and the inflation factor
  data/processed/lincs/drug_viper_specific_v2_tf_d132corr/      corrected activity = production / inflation
  data/processed/lincs/drug_viper_specific_collectri_d132corr/  (sandbox; production folders untouched)
"""
import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs/d132"
SETS = {"tf": "drug_viper_specific_v2_tf", "collectri": "drug_viper_specific_collectri"}
T0 = time.time()
LOG = []


def log(m):
    line = f"[{time.time() - T0:6.1f}s] {m}"
    LOG.append(line)
    print(line, flush=True)


def read_nes(folder, net, n_prof):
    """R wrote `regulators x profiles` column-major, so the file is profile-major: (n_profiles, n_regulators)."""
    regs = (folder / f"{net}.regulators.txt").read_text().split()
    a = np.fromfile(folder / f"{net}.f32", dtype=np.float32)
    if a.size != n_prof * len(regs):
        raise SystemExit(f"{folder/net}: {a.size} != {n_prof} x {len(regs)}")
    return regs, a.reshape(n_prof, len(regs))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nulls", default="replicate,rotation")
    ap.add_argument("--correct-with", default="replicate", help="which null defines the correction applied")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    RES.mkdir(parents=True, exist_ok=True)
    n_prof = len(pd.read_csv(PROC / "d132/replicate_index.tsv", sep="\t"))
    tables = {}

    for null in a.nulls.split(","):
        tag = "_d132" + null[:3]
        rows = []
        for setname, base in SETS.items():
            folder = PROC / (base + tag)
            if not folder.exists():
                log(f"{null}/{setname}: {folder.name} missing, skipped")
                continue
            nets = sorted(p.stem.replace(".regulators", "") for p in folder.glob("*.regulators.txt"))
            for net in nets:
                regs, X = read_nes(folder, net, n_prof)
                sd = X.std(0, ddof=1)
                rows.append(pd.DataFrame(dict(set=setname, net=net, regulator=regs,
                                              null_mean=X.mean(0), inflation=sd)))
            log(f"{null}/{setname}: {len(nets)} networks, "
                f"{sum(len(r) for r in rows if r.set.iloc[0] == setname)} regulator rows")
        D = pd.concat(rows, ignore_index=True)
        D.to_csv(RES / f"inflation_{null}.tsv", sep="\t", index=False)
        tables[null] = D
        q = D.inflation.describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
        iqr = q["75%"] - q["25%"]
        log(f"{null}: {len(D)} regulators; inflation median {q['50%']:.3f} "
            f"(IQR {q['25%']:.3f}-{q['75%']:.3f}, 5-95% {q['5%']:.3f}-{q['95%']:.3f}, "
            f"min {q['min']:.3f}, max {q['max']:.3f})")
        log(f"{null}: IQR / median = {iqr / q['50%']:.4f}  -> "
            f"{'NON-UNIFORM' if iqr / q['50%'] > 0.10 else 'UNIFORM'} by the D132(3) rule (>0.10 = non-uniform); "
            f"CV = {D.inflation.std(ddof=1) / D.inflation.mean():.4f}")

    if len(tables) > 1:
        k = list(tables)
        j = tables[k[0]].merge(tables[k[1]], on=["set", "net", "regulator"], suffixes=("_a", "_b"))
        log(f"agreement between the two nulls on {len(j)} shared regulators: "
            f"Spearman {np.corrcoef(j.inflation_a.rank(), j.inflation_b.rank())[0, 1]:.4f}, "
            f"median ratio {np.median(j.inflation_b / j.inflation_a):.3f} ({k[1]} / {k[0]})")

    # --- corrected activity sandbox --------------------------------------------------------------------------------
    D = tables[a.correct_with].set_index(["set", "net", "regulator"]).inflation
    written = 0
    pidx = len(pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"))
    for setname, base in SETS.items():
        src, dst = PROC / base, PROC / (base + "_d132corr")
        dst.mkdir(parents=True, exist_ok=True)
        for rf in sorted(src.glob("*.regulators.txt")):
            net = rf.stem.replace(".regulators", "")
            regs, X = read_nes(src, net, pidx)
            f = np.array([D.get((setname, net, r), np.nan) for r in regs], dtype=np.float64)
            miss = ~np.isfinite(f) | (f <= 0)
            if miss.any():          # a regulator with no null estimate keeps its production value, and is counted
                f[miss] = 1.0
            Y = (X / f.astype(np.float32)).astype(np.float32)
            Y.tofile(dst / f"{net}.f32")
            shutil.copy(rf, dst / rf.name)
            written += 1
            if miss.any():
                log(f"{setname}/{net}: {int(miss.sum())} of {len(regs)} regulators had no null factor, left uncorrected")
        log(f"{setname}: corrected sandbox written to {dst.name}")

    (RES / "d132_inflation_run.json").write_text(json.dumps(dict(
        stage="D132-inflation", nulls=a.nulls, correction_null=a.correct_with, n_null_profiles=n_prof,
        n_net_files_written=written,
        git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                               text=True).stdout.strip(),
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
