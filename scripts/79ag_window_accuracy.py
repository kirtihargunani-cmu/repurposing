#!/usr/bin/env python
"""Per-cohort, per-window approval-recovery AUC, matched against the pooled score (D192).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "can we calculate accuracy with the timepoint resolved profiles?"

WHAT THIS IS, AND WHAT IT IS NOT. This is **approval recovery**, not accuracy. `label_phase == 0`
is NOT a negative - the data are positive-unlabelled, and a drug unlabelled for a cohort may work
and simply never have been trialled there. Under [[discovery-is-the-goal]] a LOW value is
consistent with the pipeline doing its job. Approval is weak validation (D127 biomarker pairs is
the primary check), so this number is reported, never optimised against.

THE ONE METHODOLOGICAL POINT. A window score exists only for the drugs measured in that window,
and those drugs are not a random subset - they concentrate in well-assayed cohorts. Comparing a
window AUC over that subset against the pooled AUC over ALL drugs therefore compares two different
populations and flatters whichever has the easier one. This is the same selection effect that
invalidated D184's first run (`layered_pct` populated for 13.8% of rows).

So for every (cohort, window) cell this computes BOTH on exactly the same rows:

    auc_window              AUC of rev_<window>_pct
    auc_pooled_same_rows    AUC of rev_pct, restricted to the rows that HAVE that window
    delta                   auc_window - auc_pooled_same_rows

`delta` is the only number that says whether resolving the timepoint helped.

Standard errors are Hanley-McNeil, which assumes one score per case; with n_pos as low as 6 in
some cells they are indicative, not inferential. Cells with < 10 rows or < 3 positives are dropped.

Output: results/lincs/window_accuracy.tsv  (consumed by 79x2_window_page.py)
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SFX = os.environ.get("SPLIT_SFX", "_v2_tf")
MIN_ROWS, MIN_POS = 10, 3
OUT = RES / "window_accuracy.tsv"


def auc_se(y, s):
    """AUC with a Hanley-McNeil standard error. Returns (auc, se, n_pos, n_neg)."""
    m = ~np.isnan(s)
    y, s = np.asarray(y)[m], np.asarray(s)[m]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < MIN_POS or n0 < MIN_POS:
        return None, None, n1, n0
    r = pd.Series(s).rank().to_numpy()
    a = (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    se = np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n0 - 1) * (q2 - a * a)) / (n1 * n0))
    return float(a), float(se), n1, n0


def main() -> int:
    lab = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    lab = lab[lab.unit_kind == "group"][["unit_id", "pert_id", "label_phase"]]
    sp = pd.read_csv(RES / f"drug_reversal_split{SFX}.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "rev_pct", "rev_early_pct", "rev_late_pct"])
    d = sp.merge(lab, left_on=["group_id", "pert_id"], right_on=["unit_id", "pert_id"])
    d["y"] = (d.label_phase >= 3).astype(int)

    rows = []
    for g, dd in d.groupby("group_id"):
        for win, col in (("early", "rev_early_pct"), ("late", "rev_late_pct")):
            sub = dd[dd[col].notna()]
            if len(sub) < MIN_ROWS:
                continue
            aw, sew, n1, n0 = auc_se(sub.y.values, sub[col].values.astype(float))
            ap, _, _, _ = auc_se(sub.y.values, sub.rev_pct.values.astype(float))
            if aw is None or ap is None:
                continue
            rows.append(dict(cohort=g, window=win, n=len(sub), n_pos=n1,
                             auc_window=round(aw, 4), se=round(sew, 4),
                             auc_pooled_same_rows=round(ap, 4), delta=round(aw - ap, 4),
                             z=round((aw - 0.5) / sew, 2)))
    if not rows:
        print("no computable cells", file=sys.stderr)
        return 1
    R = pd.DataFrame(rows).sort_values(["window", "auc_window"], ascending=[True, False])
    R.to_csv(OUT, sep="\t", index=False)

    print(f"{len(R)} (cohort, window) cells over {R.cohort.nunique()} cohorts -> {OUT}")
    for w, g in R.groupby("window"):
        print(f"  {w:5s} n={len(g):2d}  median AUC window {g.auc_window.median():.3f} vs "
              f"pooled-same-rows {g.auc_pooled_same_rows.median():.3f}  "
              f"median delta {g.delta.median():+.3f}  window better in {int((g.delta > 0).sum())}/{len(g)}")
    better = int((R.delta > 0).sum())
    print(f"\nVERDICT: resolving the timepoint improves approval recovery in {better}/{len(R)} cells "
          f"(median delta {R.delta.median():+.3f}).")
    (RES / "window_accuracy_run.json").write_text(json.dumps(dict(
        script="79ag_window_accuracy.py", split_sfx=SFX, min_rows=MIN_ROWS, min_pos=MIN_POS,
        n_cells=len(R), n_cohorts=int(R.cohort.nunique()),
        window_better_cells=better, median_delta=float(R.delta.median()),
        caveat="approval recovery on positive-unlabelled data; reported, never optimised against"),
        indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
