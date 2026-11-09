#!/usr/bin/env python
"""D125 A0e (with A0c / A0d): ceiling and attainment for the reported known-use metrics (docs/SPEC_graph_metric_scoring.md
A0c-A0e; reporting rule R4 in docs/LLM_COORDINATION_NOTES.md).

Per cohort, the reliability used (A0c): A0b's r where A0b covers the cohort (its in-key LINCS lines hold an eligible
pair), else A0a's, flagged r_optimistic; within an estimator, r_pooled cohorts (< R_POOLED_MIN_DRUGS drugs) use that
estimator's median (A0d). |r_A0b - r_A0a| >= R_DIFF_NOTE is stated in every row it touches.
    D = 2 AUC - 1;  AUC_max = (1 + sqrt(r)) / 2;  attainment = D / sqrt(r)
Rows: the pooled metrics of results/lincs/known_use/primary.tsv (79l; r = median of the per-cohort r over the 41 groups,
said so in r_basis) and the per-cohort AUCs of a1_per_cohort.tsv (cohort r). The ceiling applies to scores built from
the LINCS reversal path (spec_mech_pct, spec_mech_pct_pregate: r is measured on the 79i spec_mech_pct, the pregate column;
the post-gate column is a deterministic demotion of it). Other scores keep their raw AUC with ceiling_applies False and the
reason (layered_pct uses DepMap / TCGA / ChEMBL, not LINCS; onc is the class baseline; rev / spec not measured).
Raw AUC and attainment are always in the same row; attainment is never reported alone. r <= 0 -> no ceiling (NaN).
Output: results/lincs/d125/a0/ceiling.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# --- config block (values in a0_config.py only) ---
from a0_config import R_DIFF_NOTE  # noqa: E402
import a0_common as C  # noqa: E402

LINCS_SCORES = {"spec_mech_pct", "spec_mech_pct_pregate"}
REASON = {"layered_pct": "not LINCS-derived (DepMap, TCGA, ChEMBL): the LINCS ceiling does not apply",
          "onc": "class baseline (anticancer 0/1): no ceiling"}


def reliability():
    a = pd.read_csv(C.OUT / "reliability_a0a.tsv", sep="\t")
    fb = C.OUT / "reliability_a0b.tsv"
    b = pd.read_csv(fb, sep="\t") if fb.exists() else None
    t = a[["kind", "cohort", "r_c", "r_used", "r_pooled"]].rename(columns={"r_c": "r_a0a", "r_used": "r_used_a0a",
                                                                          "r_pooled": "r_pooled_a0a"})
    if b is not None:
        t = t.merge(b[["kind", "cohort", "r_c", "r_used", "r_pooled", "r_optimistic", "indicative_only"]].rename(
            columns={"r_c": "r_a0b", "r_used": "r_used_a0b", "r_pooled": "r_pooled_a0b", "r_optimistic": "a0b_uncovered",
                     "indicative_only": "a0b_indicative"}), on=["kind", "cohort"], how="left")
    else:
        t = t.assign(r_a0b=np.nan, r_used_a0b=np.nan, r_pooled_a0b=np.nan, a0b_uncovered=True, a0b_indicative=np.nan)
    use_b = t.a0b_uncovered.eq(False) & t.r_used_a0b.notna()
    t["r"] = np.where(use_b, t.r_used_a0b, t.r_used_a0a)
    t["r_estimator"] = np.where(use_b, "A0b", "A0a")
    t["r_optimistic"] = ~use_b
    t["r_pooled"] = np.where(use_b, t.r_pooled_a0b, t.r_pooled_a0a).astype(bool)
    t["a0b_minus_a0a"] = t.r_a0b - t.r_a0a
    t["r_differ_note"] = (t.a0b_minus_a0a.abs() >= R_DIFF_NOTE).fillna(False)
    return t, b is not None


def ceiling(auc, r):
    d = 2 * auc - 1
    if not np.isfinite(r) or r <= 0:
        return d, np.nan, np.nan
    return d, (1 + np.sqrt(r)) / 2, d / np.sqrt(r)


def main():
    t, has_b = reliability()
    g = t[t.kind == "group"].set_index("cohort")
    rows = []
    prim = pd.read_csv(C.RES / "known_use/primary.tsv", sep="\t")
    for _, p in prim.iterrows():
        app = p.score in LINCS_SCORES
        r = float(g.r.median()) if app else np.nan
        d, amax, att = ceiling(p.auc, r)
        rows.append(dict(level="pooled", source=p.source, subset=p.subset, test=p.test, score=p.score, cohort="",
                         n_cohorts=p.n_cohorts, auc=p.auc, D_obs=round(d, 4), ceiling_applies=app,
                         r=r, r_basis="median of the per-cohort r over the 41 groups" if app else "",
                         r_estimator=("A0b+A0a" if g.r_estimator.eq("A0b").any() and g.r_estimator.eq("A0a").any()
                                      else g.r_estimator.iloc[0]) if app else "",
                         r_optimistic=bool(g.r_optimistic.any()) if app else np.nan,
                         r_optimistic_cohorts=int(g.r_optimistic.sum()) if app else np.nan,
                         r_pooled=bool(g.r_pooled.any()) if app else np.nan,
                         r_a0a=float(g.r_a0a.median()) if app else np.nan,
                         r_a0b=float(g.r_a0b.median()) if app and has_b else np.nan,
                         r_differ_note=bool(g.r_differ_note.any()) if app else np.nan,
                         auc_max=amax, attainment=att, note="" if app else REASON.get(p.score, "")))
    per = pd.read_csv(C.RES / "known_use/a1_per_cohort.tsv", sep="\t")
    for _, p in per[per.subset == "all"].iterrows():
        if p.cohort not in g.index or not np.isfinite(p.auc_spec_mech_pct):
            continue
        x = g.loc[p.cohort]
        d, amax, att = ceiling(p.auc_spec_mech_pct, x.r)
        rows.append(dict(level="cohort", source=p.source, subset=p.subset, test="A1 per cohort (79l)",
                         score="spec_mech_pct", cohort=p.cohort, n_cohorts=1, auc=p.auc_spec_mech_pct, D_obs=round(d, 4),
                         ceiling_applies=True, r=x.r, r_basis="this cohort's r" + (" (estimator median: r_pooled)" if x.r_pooled else ""),
                         r_estimator=x.r_estimator, r_optimistic=bool(x.r_optimistic), r_optimistic_cohorts=np.nan,
                         r_pooled=bool(x.r_pooled), r_a0a=x.r_a0a, r_a0b=x.r_a0b, r_differ_note=bool(x.r_differ_note),
                         auc_max=amax, attainment=att,
                         note=f"n labelled {p.n_labelled}, unlabelled {p.n_unlabelled}"))
    out = pd.DataFrame(rows)
    out["note"] = [n + ("; " if n else "") + f"A0b and A0a differ by >= {R_DIFF_NOTE} in r" if dn == True else n
                   for n, dn in zip(out.note, out.r_differ_note)]
    out["note"] = [n + ("; " if n else "") + "r <= 0: no attainable signal measured" if a and not np.isfinite(m) else n
                   for n, a, m in zip(out.note, out.ceiling_applies, out.auc_max)]
    out.to_csv(C.OUT / "ceiling.tsv", sep="\t", index=False, float_format="%.4f")
    pd.set_option("display.width", 250)
    print(out[out.level == "pooled"][["source", "subset", "test", "score", "auc", "D_obs", "r", "auc_max", "attainment",
                                      "r_optimistic"]].to_string(index=False))


if __name__ == "__main__":
    main()
