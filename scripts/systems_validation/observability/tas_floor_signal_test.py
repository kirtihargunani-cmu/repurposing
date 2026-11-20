#!/usr/bin/env python
"""D156 acceptance test 1 (REPLACEMENT): do newly-admitted 6 h signatures carry drug-specific signal?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

The spec originally named the 79n knockdown benchmark as the hard-fail guard. That was WRONG: 79n
scores genetic perturbations (trt_sh.cgs, trt_xpr) and never reads a compound signature or the TAS
floor, so it cannot respond to this change at all - it would have passed unchanged however much
noise was admitted.

This is the replacement, and it does respond. A signature carries drug-specific signal if it
resembles the SAME drug measured elsewhere more than it resembles a random other drug. Concretely,
for each 6 h signature of drug d in line l, compare:

    r_self   correlation with d's 24 h profile in the same line l
    r_other  correlation with a random other drug's 24 h profile in line l

If the newly admitted signatures are noise, r_self - r_other goes to zero for them. If they carry
the drug's own response, it stays positive. The test reports this separately for signatures ABOVE
the current 0.4 floor and for the ones a lower anchor would ADMIT, so the question is not "is there
signal" but "does the admitted set still have it".

Output: results/systems_validation/observability/tas_floor_signal_test.{tsv,json}
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/systems_validation/observability"
SEED, N_OTHER = 20261120, 20
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    X = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    cp = sig[sig.pert_type == "trt_cp"].copy()
    cp["row"] = cp.index.values
    early = cp[cp.pert_time <= 6]
    late = cp[(cp.pert_time >= 24) & (cp.tas >= 0.4)]        # reference: always the accepted late set
    lateref = (late.groupby(["pert_id", "ModelID"]).row.apply(list))
    print(f"{len(early):,} 6 h signatures | {len(late):,} accepted 24 h reference signatures")

    rng = np.random.default_rng(SEED)
    by_line = late.groupby("ModelID").row.apply(list).to_dict()
    rows = []
    cand = early[early.ModelID.isin(by_line)]
    cand = cand.sample(min(6000, len(cand)), random_state=SEED)
    for t in cand.itertuples():
        key = (t.pert_id, t.ModelID)
        if key not in lateref.index:
            continue
        v = np.asarray(X[t.row], dtype=np.float64)
        self_rows = lateref.loc[key]
        s = np.asarray(X[self_rows], dtype=np.float64).mean(0)
        pool = [r for r in by_line[t.ModelID] if r not in set(self_rows)]
        if len(pool) < 5 or np.std(v) == 0 or np.std(s) == 0:
            continue
        oth = rng.choice(pool, min(N_OTHER, len(pool)), replace=False)
        r_self = float(np.corrcoef(v, s)[0, 1])
        r_other = float(np.median([np.corrcoef(v, np.asarray(X[o], dtype=np.float64))[0, 1] for o in oth]))
        rows.append(dict(sig_id=t.sig_id, pert_id=t.pert_id, ModelID=t.ModelID, tas=t.tas,
                         pert_time=t.pert_time, r_self=r_self, r_other=r_other,
                         delta=r_self - r_other))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tas_floor_signal_test.tsv", sep="\t", index=False)

    bands = [("above current floor (TAS>=0.400)", df.tas >= 0.400),
             ("ADMITTED at anchor 0.645 (0.331-0.400)", (df.tas >= 0.331) & (df.tas < 0.400)),
             ("ADMITTED at anchor 0.55 (0.284-0.331)", (df.tas >= 0.284) & (df.tas < 0.331)),
             ("still excluded (TAS<0.284)", df.tas < 0.284)]
    print(f"\n{len(df):,} 6 h signatures with a same-drug 24 h reference in the same line\n")
    print(f"{'band':42s} {'n':>6s} {'r_self':>8s} {'r_other':>8s} {'delta':>8s} {'% delta>0':>10s}")
    res = []
    for lab, m in bands:
        g = df[m]
        if not len(g):
            continue
        row = dict(band=lab, n=int(len(g)), r_self=round(float(g.r_self.median()), 4),
                   r_other=round(float(g.r_other.median()), 4),
                   delta=round(float(g.delta.median()), 4),
                   pct_positive=round(100 * float((g.delta > 0).mean()), 1))
        res.append(row)
        print(f"{lab:42s} {row['n']:>6,} {row['r_self']:>8.4f} {row['r_other']:>8.4f} "
              f"{row['delta']:>8.4f} {row['pct_positive']:>9.1f}%")
    (OUT / "tas_floor_signal_test.json").write_text(json.dumps(dict(
        test="D156 acceptance test 1 (replacement): drug-specific signal in admitted 6 h signatures",
        note="79n cannot serve as this guard - it scores genetic perturbations and never reads the compound TAS floor",
        seed=SEED, n_signatures=int(len(df)), bands=res), indent=2))
    print("\nwrote", OUT / "tas_floor_signal_test.json")


if __name__ == "__main__":
    sys.exit(main())
