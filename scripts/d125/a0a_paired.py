#!/usr/bin/env python
"""D125 A0a: paired-signature reliability, every cohort (docs/SPEC_graph_metric_scoring.md A0a; D125-A0 (2)-(4)).

Stratum = (dose, time) cell of a drug x line pair; a cell's profile = mean of its active level-5 signatures. Pairs with
>= A0A_MIN_STRATA cells: two cells drawn at random (seed A0A_SEED) give pass A and pass B; every other pair keeps its
production profile in both passes. Each pass runs the production path (a0_common: 79c rules with the production g_l,
79d on the changed columns spliced into the production VIPER, 79e groups + subgroups, 79i STRATA=1) in the sandbox.
r_c = per cohort (41 groups, 14 subgroups) Spearman of spec_mech_pct between the passes over gate-passing drugs that draw
on >= 1 changed pair; CI Fisher z (Bonett-Wright); < R_POOLED_MIN_DRUGS drugs -> r_pooled (corrected with the median).
Checks before any pass is read (stop on failure): (1) the rebuilt production profiles of the changed pairs equal 79c's;
(2) 79d on them equals the production VIPER columns; (3) the sandbox scoring of unchanged production inputs equals the
production spec_mech_pct_pregate (groups and subgroups).
Outputs: results/lincs/d125/a0/reliability_a0a.tsv, a0a_strata.tsv (the two cells per pair), a0a_run.json;
         pass tables data/processed/lincs/d125_a0/passes/ (reruns reuse them).
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# --- config block (values in a0_config.py only) ---
from a0_config import (A0A_STRATUM, A0A_MIN_STRATA, A0A_SEED, SCORE_IDENTITY_TOL, THREADS, R_POOLED_MIN_DRUGS)  # noqa: E402
import a0_common as C  # noqa: E402

TAG = "_d125a0a"


def baseline_check():
    """(3): sandbox scoring of production inputs vs production spec_mech_pct_pregate."""
    tag = "_d125a0_base"
    b = C.load_pass(tag)
    if b is None:
        C.cleanup(tag)
        pre = C.stage_baseline(tag)
        b, secs = C.score_pass(tag, pre, THREADS)
        C.save_pass(tag, b)
        C.record_time("score_pass", 0, secs)
    g = pd.read_csv(C.RES / f"{C.GROUP_TABLE}.tsv", sep="\t", low_memory=False,
                    usecols=["group_id", "pert_id", "spec_mech_pct_pregate"]).rename(columns={"group_id": "cohort"})
    s = pd.read_csv(C.RES / f"{C.STRATUM_TABLE}.tsv", sep="\t", low_memory=False,
                    usecols=["stratum", "pert_id", "spec_mech_pct_pregate"]).rename(columns={"stratum": "cohort"})
    prod = pd.concat([g.assign(kind="group"), s.assign(kind="subgroup")])
    m = b.merge(prod, on=["kind", "cohort", "pert_id"], how="outer", indicator=True)
    both = m[m._merge == "both"]
    d = (both.spec_mech_pct - both.spec_mech_pct_pregate).abs()
    nan_mismatch = int((both.spec_mech_pct.isna() != both.spec_mech_pct_pregate.isna()).sum())
    res = dict(rows_sandbox=int((m._merge != "right_only").sum()), rows_production=int((m._merge != "left_only").sum()),
               rows_only_one_side=int((m._merge != "both").sum()), max_abs_diff=float(d.max()),
               rows_diff_gt_tol=int((d > SCORE_IDENTITY_TOL).sum()), nan_mismatch=nan_mismatch)
    ok = res["rows_only_one_side"] == 0 and res["rows_diff_gt_tol"] == 0 and nan_mismatch == 0
    C.log(f"score identity (sandbox vs production pregate): {res} -> {'ok' if ok else 'FAILED'}")
    if not ok:
        raise SystemExit("score identity failed: the sandbox does not reproduce production spec_mech_pct_pregate; "
                         "A0 would not describe the production score. Stopping (report, do not adjust).")
    return res, b


def main():
    C.init_work()
    t0 = time.time()
    act = C.active_signatures()
    P = C.profiles()
    cells = act.groupby(["col", *A0A_STRATUM], dropna=False).row.apply(list).reset_index()
    cells = cells.sort_values(["col", "pert_time", "pert_dose"]).reset_index(drop=True)
    ncell = cells.groupby("col").size()
    cols = np.array(sorted(ncell.index[ncell >= A0A_MIN_STRATA]))
    rng = np.random.default_rng(A0A_SEED)
    pick, rows_a, rows_b = [], {}, {}
    for c, d in cells[cells.col.isin(cols)].groupby("col", sort=True):
        i, j = rng.choice(len(d), 2, replace=False)
        a, b = d.iloc[i], d.iloc[j]
        rows_a[c], rows_b[c] = a.row, b.row
        pick.append(dict(col=c, n_cells=len(d), dose_a=a.pert_dose, time_a=a.pert_time, n_sig_a=len(a.row),
                         dose_b=b.pert_dose, time_b=b.pert_time, n_sig_b=len(b.row)))
    pick = pd.DataFrame(pick).merge(P.idx[["col", "pert_id", "ModelID"]], on="col")
    pick.to_csv(C.OUT / "a0a_strata.tsv", sep="\t", index=False)
    C.log(f"A0a: {len(cols)} pairs with >= {A0A_MIN_STRATA} dose x time cells; "
          f"same time {int((pick.time_a == pick.time_b).sum())}, same dose {int((pick.dose_a == pick.dose_b).sum())}")

    run = dict(stratum=list(A0A_STRATUM), seed=A0A_SEED, n_changed_pairs=len(cols))
    if C.load_pass(TAG + "A") is None or C.load_pass(TAG + "B") is None:
        C.build_sandbox()
        run["profile_identity"] = P.identity_check(act, cols)
        t = time.time()
        sub = C.viper_subset("_d125a0_id", cols, P.spec[cols], THREADS)
        C.record_time("viper_subset", len(cols), time.time() - t)
        run["viper_identity"] = C.check_viper_identity(sub, cols)
        del sub
    else:
        C.log("A0a passes on disk; identity checks were run with them (a0a_run.json of that run)")
    run["score_identity"], base = baseline_check()

    passes = {}
    for k, rows in [("A", rows_a), ("B", rows_b)]:
        tag = TAG + k
        x = C.load_pass(tag)
        if x is None:
            C.cleanup(tag)
            built = [P.build(c, rows[c]) for c in cols]
            raw = np.array([r for r, _ in built]); spec = np.array([s for _, s in built])
            t = time.time()
            sub = C.viper_subset(tag, cols, spec, THREADS)
            C.record_time("viper_subset", len(cols), time.time() - t)
            pre = C.stage_pass(tag, cols, raw, spec, sub)
            del sub
            x, secs = C.score_pass(tag, pre, THREADS)
            C.record_time("score_pass", len(cols), secs)
            C.save_pass(tag, x)
        passes[k] = x
        C.log(f"pass {k}: {len(x)} cohort x drug rows")

    used, _ = C.pairs_used(passes["A"][["kind", "cohort", "parent_group", "pert_id", "source"]])
    r = C.r_by_cohort(passes["A"], passes["B"], cols, used)
    r, med = C.finish_reliability(r, "A0a")
    r["n_splits"], r["se"] = np.nan, np.nan     # A0a: no resampling loop (spec A0a)
    r["r_optimistic"] = True                     # A0a alone; A0e sets the final flag against A0b's coverage
    cols_out = ["kind", "cohort", "estimator", "r_c", "ci_lo", "ci_hi", "n_drugs", "n_gate", "n_splits", "se",
                "r_pooled", "r_used", "r_optimistic"]
    r[cols_out].to_csv(C.OUT / "reliability_a0a.tsv", sep="\t", index=False, float_format="%.4f")
    run.update(median_r_c=med, median_r_c_groups=float(r[r.kind == "group"].r_c.median()),
               median_r_c_subgroups=float(r[r.kind == "subgroup"].r_c.median()),
               n_r_pooled=int(r.r_pooled.sum()), r_pooled_min_drugs=R_POOLED_MIN_DRUGS,
               minutes=round((time.time() - t0) / 60, 1))
    C.write_json(C.OUT / "a0a_run.json", run)
    C.log(f"A0a: median r_c {med:.3f} (groups {run['median_r_c_groups']:.3f}, subgroups {run['median_r_c_subgroups']:.3f}); "
          f"r_pooled {run['n_r_pooled']}")
    print(r[cols_out].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
