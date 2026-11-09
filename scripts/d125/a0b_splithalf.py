#!/usr/bin/env python
"""D125 A0b: split-half reliability on the deep pairs (docs/SPEC_graph_metric_scoring.md A0b; D125-A0).

Scope comes from A0.0 (results/lincs/d125/a0/eligible_pairs_gate.json; run a00_eligible_pairs.py first): the pairs
with >= A00_MIN_REPLICATES active level-5 signatures. A0.0 says cohort-level (>= 500 pairs in >= 3 lines), so this is
the cohort-level path; the line-level path is not implemented (it stops if A0.0 ever says line-level).
Split s (seed A0B_SEED + s): each eligible pair's signatures shuffled and halved (floor(n/2) / rest); half 1 and half 2
are two passes of the production path (a0_common, as A0a); per cohort r_c^s = Spearman of spec_mech_pct between the
halves over gate-passing drugs whose every scoring line is an
eligible pair (D125-A0.1; >= 1 eligible pair reported as r_c_any). Batches of A0B_BATCH splits; after >= A0B_FLOOR splits stop
when SE(median r_c) <= A0B_SE_STABLE (0.01, D125-A0.1), else at A0B_CAP (then the ceiling is indicative only). SE(median r_c) = sd over
splits of the median-across-cohorts r_c^s / sqrt(splits). Per cohort: r_c = mean over splits, SE = sd / sqrt(splits),
CI = split percentiles; covered (r_optimistic False) when the cohort's in-key LINCS lines hold >= 1 eligible pair.

DRY_RUN=1 (default): reports planned passes, estimated wall time and disk (from the A0a timings in
data/processed/lincs/d125_a0/timings.tsv and the sizes on disk) and writes a0b_dryrun.json; nothing is executed.
DRY_RUN=0 executes; refuses while A0B_SE_STABLE is None (no numeric stopping rule is pre-registered).
Outputs: results/lincs/d125/a0/reliability_a0b.tsv, a0b_run.json, a0b_dryrun.json; per-split r in
         data/processed/lincs/d125_a0/a0b_splits.tsv (reruns resume from it).
"""
import json
import os
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# --- config block (values in a0_config.py only) ---
from a0_config import (A0B_BATCH, A0B_FLOOR, A0B_CAP, A0B_SE_STABLE, A0B_SEED, CI_LEVEL, THREADS,  # noqa: E402
                       PASS_PARALLEL)
import a0_common as C  # noqa: E402

TAG = "_d125a0b"
SPLITS = C.WORK / "a0b_splits.tsv"


def gate():
    f = C.OUT / "eligible_pairs_gate.json"
    if not f.exists():
        raise SystemExit("A0.0 has not run (eligible_pairs_gate.json missing): run a00_eligible_pairs.py first")
    g = json.loads(f.read_text())
    if not g["cohort_level"]:
        raise SystemExit(f"A0.0 says line-level bound ({g['n_pairs']} pairs, {g['n_lines']} lines): the line-level "
                         "path of A0b is not implemented; stop and report")
    e = pd.read_csv(C.OUT / "eligible_pairs.tsv", sep="\t")
    return g, np.array(sorted(e.col[e.a0b_eligible]))


def dir_bytes(p):
    return sum(f.stat().st_size for f in Path(p).rglob("*") if f.is_file())


def dry_run(g, cols):
    P = C.production_index()
    n_prof, n_genes = len(P), len(pd.read_csv(C.PROC / "drug_line_profiles_genes.tsv", sep="\t"))
    tf = C.WORK / "timings.tsv"
    t = pd.read_csv(tf, sep="\t") if tf.exists() else pd.DataFrame(columns=["stage", "n_cols", "seconds"])
    v = t[(t.stage == "viper_subset") & (t.n_cols > 0)]
    s = t[t.stage == "score_pass"]
    per_col = float((v.seconds / v.n_cols).median()) if len(v) else np.nan      # at THREADS cores (A0a runs)
    viper_s = per_col * len(cols)
    score_s = float(s.seconds.median()) if len(s) else np.nan
    tpp = max(THREADS // PASS_PARALLEL, 1)
    # VIPER parallelises over networks (mclapply): time ~ cores^-1; scoring is mostly one core (BLAS threads aside)
    pass_s = viper_s * THREADS / tpp + score_s
    prof_b = 2 * n_prof * n_genes * 4
    viper_b = sum(dir_bytes(C.PROC / d) for d in C.VIPER_SETS)
    tables_b = sum((C.RES / f).stat().st_size for f in ["drug_reversal_v2_consensus.tsv", "stratum_reversal_v2_consensus.tsv",
                                                        f"{C.GROUP_TABLE}.tsv", f"{C.STRATUM_TABLE}.tsv"])
    sub_b = len(cols) * n_genes * 4
    kept = [f.stat().st_size for f in C.PASSES.glob("*.tsv.gz")]
    kept_b = float(np.median(kept)) if kept else np.nan
    free = shutil.disk_usage(C.PROC).free
    plan = dict(
        a00=dict(n_pairs=g["n_pairs"], n_lines=g["n_lines"], scope="cohort-level"),
        splits=dict(batch=A0B_BATCH, floor=A0B_FLOOR, cap=A0B_CAP, se_stable=A0B_SE_STABLE,
                    stopping_rule_set=A0B_SE_STABLE is not None),
        passes=dict(per_split=2, min=2 * A0B_FLOOR, max=2 * A0B_CAP, viper_columns_per_pass=len(cols),
                    viper_columns_of_production=n_prof, parallel_passes=PASS_PARALLEL, threads_per_pass=tpp),
        timing_basis=dict(viper_seconds_per_column_at_threads=round(per_col, 4) if np.isfinite(per_col) else None,
                          viper_runs_measured=len(v), score_pass_seconds=round(score_s, 1) if np.isfinite(score_s) else None,
                          score_runs_measured=len(s), source=str(tf)),
        est_minutes=dict(per_pass=round(pass_s / 60, 1) if np.isfinite(pass_s) else None,
                         floor_wall=round(2 * A0B_FLOOR * pass_s / PASS_PARALLEL / 60, 0) if np.isfinite(pass_s) else None,
                         cap_wall=round(2 * A0B_CAP * pass_s / PASS_PARALLEL / 60, 0) if np.isfinite(pass_s) else None),
        disk_gb=dict(transient_per_pass=round((prof_b + viper_b + tables_b + sub_b) / 1e9, 2),
                     transient_peak=round(PASS_PARALLEL * (prof_b + viper_b + tables_b + sub_b) / 1e9, 2),
                     kept_per_pass=round(kept_b / 1e9, 4) if np.isfinite(kept_b) else None,
                     kept_at_cap=round(2 * A0B_CAP * kept_b / 1e9, 2) if np.isfinite(kept_b) else None,
                     written_at_cap=round(2 * A0B_CAP * (prof_b + viper_b + tables_b + sub_b) / 1e12, 2),
                     written_at_cap_unit="TB (streamed: each pass's files are removed once its scores are saved)",
                     free_now=round(free / 1e9, 1)),
        caveats=["VIPER time scaled linearly from the A0a subset runs (same networks, more or fewer columns)",
                 "per-pass scoring time from the A0a / baseline sandbox passes (79e groups + subgroups, 79i)",
                 "estimates exclude profile building (seconds) and splicing I/O (~1.5 GB written per pass)"])
    C.write_json(C.OUT / "a0b_dryrun.json", plan)
    print(json.dumps(plan, indent=1))
    if not plan["splits"]["stopping_rule_set"]:
        print("\nA0b cannot execute: A0B_SE_STABLE is None (the spec gives no numeric stopping criterion). "
              "Record one as a D125-A0 amendment before DRY_RUN=0.")


def halves(act_rows, cols, s):
    rng = np.random.default_rng(A0B_SEED + s)
    h1, h2 = {}, {}
    for c in cols:
        r = np.array(act_rows[c])
        r = r[rng.permutation(len(r))]
        k = len(r) // 2
        h1[c], h2[c] = r[:k], r[k:]
    return h1, h2


def run_half(tag, cols, rows, P, tpp):
    x = C.load_pass(tag)
    if x is not None:
        return x
    C.cleanup(tag)
    built = [P.build(c, rows[c]) for c in cols]
    raw = np.array([r for r, _ in built]); spec = np.array([sp for _, sp in built])
    t = time.time()
    sub = C.viper_subset(tag, cols, spec, tpp)
    C.record_time("viper_subset_a0b", len(cols), time.time() - t)
    pre = C.stage_pass(tag, cols, raw, spec, sub)
    del sub
    x, secs = C.score_pass(tag, pre, tpp)
    C.record_time("score_pass_a0b", len(cols), secs)
    C.save_pass(tag, x)
    return x


def coverage(cols):
    """cohort -> covered: its in-key LINCS lines (79e direct lines) hold >= 1 eligible pair."""
    elig_lines = set(C.production_index().ModelID.values[cols])
    return {k: bool(elig_lines & v) for k, v in C.cohort_inkey_lines().items()}


def execute(g, cols):
    if A0B_SE_STABLE is None:
        raise SystemExit("A0B_SE_STABLE is None: no numeric stopping rule is pre-registered (spec A0b). "
                         "Record one as a D125-A0 amendment first; not executing.")
    act = C.active_signatures()
    P = C.profiles()
    act_rows = act.groupby("col").row.apply(list)
    C.build_sandbox()
    run = dict(n_pairs=len(cols), n_lines=g["n_lines"], se_stable=A0B_SE_STABLE, floor=A0B_FLOOR, cap=A0B_CAP)
    run["profile_identity"] = P.identity_check(act, cols)
    sub = C.viper_subset("_d125a0b_id", cols, P.spec[cols], THREADS)
    run["viper_identity"] = C.check_viper_identity(sub, cols)
    del sub
    tpp = max(THREADS // PASS_PARALLEL, 1)
    done = pd.read_csv(SPLITS, sep="\t") if SPLITS.exists() else pd.DataFrame()
    used, start = None, 0
    while True:
        for s in range(start, min(start + A0B_BATCH, A0B_CAP)):
            if len(done) and (done.split == s).any():
                continue
            h1, h2 = halves(act_rows, cols, s)
            with ThreadPoolExecutor(PASS_PARALLEL) as ex:
                f1 = ex.submit(run_half, f"{TAG}_s{s:03d}h1", cols, h1, P, tpp)
                f2 = ex.submit(run_half, f"{TAG}_s{s:03d}h2", cols, h2, P, tpp)
                a, b = f1.result(), f2.result()
            if used is None:
                used, _ = C.pairs_used(a[["kind", "cohort", "parent_group", "pert_id", "source"]])
            r = C.r_by_cohort(a, b, cols, used).assign(split=s)
            r.to_csv(SPLITS, sep="\t", index=False, mode="a", header=not SPLITS.exists())
            done = pd.concat([done, r], ignore_index=True)
            C.log(f"A0b split {s}: median r_c {r.r_c.median():.3f}")
        start += A0B_BATCH
        n = done.split.nunique()
        med = done.groupby("split").r_c.median()
        se = float(med.std(ddof=1) / np.sqrt(n)) if n > 1 else np.inf
        stable = n >= A0B_FLOOR and se <= A0B_SE_STABLE
        C.log(f"A0b after {n} splits: median r_c {med.mean():.3f}, SE {se:.4f} -> {'stable' if stable else 'continue'}")
        if stable or start >= A0B_CAP:
            break
    n = done.split.nunique()
    q = (1 - CI_LEVEL) / 2
    agg = done.groupby(["kind", "cohort"]).agg(
        r_c=("r_c", "mean"), sd=("r_c", "std"), ci_lo=("r_c", lambda v: v.quantile(q)), ci_hi=("r_c", lambda v: v.quantile(1 - q)),
        n_drugs=("n_drugs", "median"), n_drugs_direct=("n_drugs_direct", "median"), n_gate=("n_gate", "median"),
        n_splits=("split", "nunique"), r_c_any=("r_c_any", "mean"), n_drugs_any=("n_drugs_any", "median")).reset_index()
    agg["se"] = agg.sd / np.sqrt(agg.n_splits)
    agg, medr = C.finish_reliability(agg, "A0b")
    cov = coverage(cols)
    agg["r_optimistic"] = [not cov.get((k, c), False) for k, c in zip(agg.kind, agg.cohort)]
    med = done.groupby("split").r_c.median()
    se_med = float(med.std(ddof=1) / np.sqrt(n))
    agg["stabilised"] = se_med <= A0B_SE_STABLE and n >= A0B_FLOOR
    agg["indicative_only"] = ~agg.stabilised
    cols_out = ["kind", "cohort", "estimator", "r_c", "ci_lo", "ci_hi", "n_drugs", "n_drugs_direct", "n_gate", "n_splits",
                "se", "r_pooled", "r_used", "r_optimistic", "stabilised", "indicative_only", "r_c_any", "n_drugs_any"]
    agg[cols_out].to_csv(C.OUT / "reliability_a0b.tsv", sep="\t", index=False, float_format="%.4f")
    run.update(n_splits=n, se_median_r_c=se_med, median_r_c=medr, stabilised=bool(agg.stabilised.iloc[0]))
    C.write_json(C.OUT / "a0b_run.json", run)


def main():
    C.init_work()
    g, cols = gate()
    C.log(f"A0b scope (A0.0): {len(cols)} pairs in {g['n_lines']} lines, cohort-level")
    if os.environ.get("DRY_RUN", "1") != "0":
        dry_run(g, cols)
    else:
        execute(g, cols)


if __name__ == "__main__":
    main()
