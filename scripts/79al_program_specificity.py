#!/usr/bin/env python
"""Specificity-decomposed program score, and whether it ranks better (D201).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_program_specificity.md  (pre-registered; read it before changing anything here)

Programs are partitioned by the OBJECTIVE D199 annotations - artefact (plate / RNA quality),
proliferation, clean - and the opposition is summed within each part. Nothing is deleted: all four
sums are carried so proliferation stays available as its own ranking, per the user's instruction
not to remove it yet.

The original design - residualising on ||drug activity|| - was withdrawn before registration
because the premise measured false (median |rho| 0.054). See SPEC section 1.
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
OUT = RES / "programs"
MIN_ROWS, MIN_POS = 10, 3


def auc_se(y, s):
    m = ~np.isnan(s)
    y, s = np.asarray(y)[m], np.asarray(s)[m]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < MIN_POS or n0 < MIN_POS:
        return None, n1
    r = pd.Series(s).rank().to_numpy()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)), n1


def main() -> int:
    D = pd.read_csv(OUT / "program_drug.tsv.gz", sep="\t")
    F = pd.read_csv(OUT / "program_filter.tsv", sep="\t")
    key = ["cohort", "program"]
    cls = F[key].copy()
    cls["part"] = np.where(F.artefact, "art",
                           np.where(F.fail_proliferation, "prolif", "clean"))
    D = D.merge(cls, on=key, how="left")
    D["part"] = D.part.fillna("clean")

    S = (D.pivot_table(index=["cohort", "pert_id", "cmap_name", "window"], columns="part",
                       values="opposition", aggfunc="sum").fillna(0.0).reset_index())
    for c in ("clean", "prolif", "art"):
        if c not in S:
            S[c] = 0.0
    S = S.rename(columns={"clean": "S_clean", "prolif": "S_prolif", "art": "S_art"})
    S["S_total"] = S.S_clean + S.S_prolif + S.S_art
    den = S.S_clean.abs() + S.S_prolif.abs()
    S["specificity"] = np.where(den > 0, S.S_clean / den.replace(0, np.nan), np.nan)
    S.to_csv(OUT / "specificity_scores.tsv.gz", sep="\t", index=False, float_format="%.6g")

    sp = pd.read_csv(RES / "drug_reversal_split_v2_tf.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "rev_pct", "rev_early_pct", "rev_late_pct"])
    lab = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    lab = lab[lab.unit_kind == "group"][["unit_id", "pert_id", "label_phase"]]
    M = (S.merge(sp, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
         .merge(lab, left_on=["cohort", "pert_id"], right_on=["unit_id", "pert_id"]))
    M["y"] = (M.label_phase >= 3).astype(int)
    M["rev_win"] = np.where(M.window == "early", M.rev_early_pct, M.rev_late_pct)

    SCORES = [("B0_rev_pct", "rev_pct"), ("B1_rev_window", "rev_win"), ("B2_S_total", "S_total"),
              ("T1_S_clean", "S_clean"), ("T2_S_prolif", "S_prolif"), ("T3_specificity", "specificity")]
    rows, dropped = [], []
    for (c, w), g in M.groupby(["cohort", "window"]):
        if len(g) < MIN_ROWS or g.y.sum() < MIN_POS:
            dropped.append(dict(cohort=c, window=w, n=len(g), n_pos=int(g.y.sum())))
            continue
        r = dict(cohort=c, window=w, n=len(g), n_pos=int(g.y.sum()))
        for name, col in SCORES:
            a, _ = auc_se(g.y.values, g[col].values.astype(float))
            r[name] = None if a is None else round(a, 4)
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "specificity_auc.tsv", sep="\t", index=False)

    print(f"{len(R)} (cohort, window) cells, {len(dropped)} dropped\n")
    print(f"{'score':18s} {'median AUC':>10s} {'vs B0':>8s} {'better':>9s}")
    base = R.B0_rev_pct
    verdict = {}
    for name, _ in SCORES:
        v = R[name].dropna()
        d = (R[name] - base).dropna()
        print(f"{name:18s} {v.median():10.3f} {d.median():+8.3f} {int((d>0).sum()):5d}/{len(d)}")
        verdict[name] = dict(median_auc=float(v.median()), median_delta_vs_B0=float(d.median()),
                             better=int((d > 0).sum()), n=int(len(d)))
    dc = (R.T1_S_clean - base).dropna()
    c1 = bool(dc.median() >= 0.03 and (dc > 0).mean() >= 0.60)
    c2 = bool(R.T2_S_prolif.median() <= R.T1_S_clean.median())
    print(f"\nACCEPTANCE (SPEC s4)")
    print(f"  1. S_clean beats rev_pct by >=0.03 and in >=60% of cells: "
          f"delta {dc.median():+.3f}, {100*(dc>0).mean():.0f}% -> {'PASS' if c1 else 'FAIL'}")
    print(f"  2. gain not explained by proliferation (S_prolif <= S_clean): "
          f"{R.T2_S_prolif.median():.3f} vs {R.T1_S_clean.median():.3f} -> {'PASS' if c2 else 'FAIL'}")
    for w, g in R.groupby("window"):
        d = (g.T1_S_clean - g.B0_rev_pct).dropna()
        print(f"  by window {w:5s}: S_clean {g.T1_S_clean.median():.3f} vs rev {g.B0_rev_pct.median():.3f} "
              f"(delta {d.median():+.3f}, better {int((d>0).sum())}/{len(d)})")
    (OUT / "specificity_run.json").write_text(json.dumps(dict(
        script="79al_program_specificity.py", spec="docs/SPEC_program_specificity.md",
        n_cells=len(R), dropped=dropped, scores=verdict,
        acceptance=dict(c1_rank_gain=c1, c2_not_proliferation=c2, passed=bool(c1 and c2))), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
