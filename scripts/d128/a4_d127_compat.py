#!/usr/bin/env python
"""D128 A4 action 4: replay the exact D127 development pairs, DESCRIPTIVE ONLY.

D127's 49 primary instances (plus its 18 previously-seen ones) are **already seen** development benchmarks. The plan
forbids optimising on them, so this script computes and reports; it chooses nothing and tunes nothing. Every pair,
compound and cohort definition is taken from `results/lincs/d127/instances.tsv` unchanged - the same `drug`
(`cmap_name`), the same high and low cohort - and the only thing added is D128's M and R beside the production lead
percentile and the target-role / reversal baselines that D127 itself reported.

`pert_id` is the compound key (P12): a compound with several `pert_id`s is reported per `pert_id` and as their mean;
the instance-level verdict uses the mean, and the per-`pert_id` spread is kept so a merged verdict can be checked.

Outputs results/lincs/d128/a4/: d127_seen_pairs.tsv, d127_seen_summary.tsv.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402

A4 = C.OUT / "a4"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3", default="a3")
    a = ap.parse_args()
    I = pd.read_csv(C.ROOT / "results/lincs/d127/instances.tsv", sep="\t")
    S = pd.read_csv(C.OUT / a.a3 / "scores_by_cohort.tsv", sep="\t", low_memory=False)

    keep = ["cohort", "pert_id", "cmap_name", "M", "M_pct", "R", "R_pct", "M_n_genes", "M_n_landmark", "R_n_genes",
            "n_units_M", "unit_scored_fraction", "layered_pct", "split_spec_mech_pct_pregate", "target_role_pct",
            "split_rev_pct", "M_abstain_reason", "source", "in_comparison_population"]
    keep = [c for c in keep if c in S.columns]
    rows = []
    for r in I.itertuples(index=False):
        h = S[(S.cohort == r.high) & (S.cmap_name == r.drug)][keep]
        l = S[(S.cohort == r.low) & (S.cmap_name == r.drug)][keep]
        d = dict(pair_set=r.pair_set, drug_class=r.drug_class, high=r.high, low=r.low, drug=r.drug,
                 d127_hi_layered_pct=r.hi_pct, d127_lo_layered_pct=r.lo_pct, d127_gap=r.gap,
                 d127_correct=bool(r.correct),
                 n_pert_ids_high=len(h), n_pert_ids_low=len(l),
                 pert_ids="|".join(sorted(set(h.pert_id) | set(l.pert_id))))
        for side, t in (("hi", h), ("lo", l)):
            for col, name in (("M", "M"), ("M_pct", "M_pct"), ("R", "R"), ("R_pct", "R_pct"),
                              ("layered_pct", "layered_pct"), ("split_spec_mech_pct_pregate", "spec_mech_pregate"),
                              ("target_role_pct", "target_role_pct"), ("M_n_genes", "M_n_genes"),
                              ("M_n_landmark", "M_n_landmark"), ("unit_scored_fraction", "unit_scored_fraction")):
                if col in t.columns:
                    d[f"{side}_{name}"] = float(t[col].mean()) if len(t) and t[col].notna().any() else np.nan
            d[f"{side}_M_spread_over_pert_ids"] = (float(t.M.max() - t.M.min())
                                                   if len(t) > 1 and t.M.notna().any() else np.nan)
            d[f"{side}_abstain"] = ";".join(sorted({x for x in t.M_abstain_reason.dropna() if x})) if len(t) else "no_row"
        d["M_gap"] = d.get("hi_M", np.nan) - d.get("lo_M", np.nan)
        d["M_pct_gap"] = d.get("hi_M_pct", np.nan) - d.get("lo_M_pct", np.nan)
        d["R_gap"] = d.get("hi_R", np.nan) - d.get("lo_R", np.nan)
        d["layered_pct_gap_d128_join"] = d.get("hi_layered_pct", np.nan) - d.get("lo_layered_pct", np.nan)
        d["M_correct"] = bool(d["M_gap"] > 0) if np.isfinite(d["M_gap"]) else None
        d["M_pct_correct"] = bool(d["M_pct_gap"] > 0) if np.isfinite(d["M_pct_gap"]) else None
        d["R_correct"] = bool(d["R_gap"] > 0) if np.isfinite(d["R_gap"]) else None
        d["M_computable"] = bool(np.isfinite(d["M_gap"]))
        rows.append(d)
    D = pd.DataFrame(rows)
    D.to_csv(A4 / "d127_seen_pairs.tsv", sep="\t", index=False)

    def agg(g):
        out = dict(n_instances=len(g), n_M_computable=int(g.M_computable.sum()),
                   d127_share_correct=round(float(g.d127_correct.mean()), 4),
                   d127_median_gap=round(float(g.d127_gap.median()), 2))
        m = g[g.M_computable]
        out.update(M_share_correct=round(float(m.M_correct.mean()), 4) if len(m) else np.nan,
                   M_median_gap=round(float(m.M_gap.median()), 5) if len(m) else np.nan,
                   M_pct_share_correct=round(float(m.M_pct_correct.dropna().mean()), 4) if len(m) else np.nan,
                   M_pct_median_gap=round(float(m.M_pct_gap.median()), 2) if len(m) else np.nan,
                   R_share_correct=round(float(m.R_correct.dropna().mean()), 4) if m.R_correct.notna().any() else np.nan,
                   R_median_gap=round(float(m.R_gap.median()), 5) if len(m) else np.nan,
                   median_hi_M_n_genes=float(m.hi_M_n_genes.median()) if len(m) else np.nan,
                   median_hi_M_n_landmark=float(m.hi_M_n_landmark.median()) if len(m) else np.nan)
        return pd.Series(out)

    SU = pd.concat([D.groupby("pair_set").apply(agg, include_groups=False).reset_index(),
                    D.groupby(["pair_set", "drug_class"]).apply(agg, include_groups=False).reset_index()],
                   ignore_index=True)
    SU["note"] = "DESCRIPTIVE ONLY: these are previously seen development instances; nothing is tuned on them."
    SU.to_csv(A4 / "d127_seen_summary.tsv", sep="\t", index=False)
    print(SU.to_string())
    print(D[["pair_set", "drug_class", "high", "low", "drug", "hi_M", "lo_M", "M_gap", "M_correct",
             "d127_gap", "d127_correct"]].to_string())

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A4 / "d127_replay_run.json").write_text(json.dumps(dict(
        stage="D128-A4 action 4 (D127 replay)", git_sha=sha, descriptive_only=True,
        n_instances=len(D), n_primary=int((D.pair_set == "primary").sum()),
        summary=SU.to_dict("records"), seconds=round(time.time() - time.time(), 1)), indent=1, default=str))


if __name__ == "__main__":
    main()
