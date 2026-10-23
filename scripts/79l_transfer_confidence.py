#!/usr/bin/env python
"""Confidence that a borrowed drug's effect carries over to a disease it was not tested in (first version, user
2026-10-20). Diagnostic: nothing replaces rev / spec / tiers yet.

Why: borrowed drugs (scripts/79b: not active in the disease's tissue key, active in >= 3 cancer lines elsewhere) are
73% of scored disease x drug rows and get the same weight as drugs measured in the tissue. Some drugs do the same thing
in every line, others do not, and that is measurable.
Calibration data: the leave-one-line-out test (scripts/77, results/lincs/extrapolation_lolo_pairs.tsv): for each deep
line L and drug d, r_mean = Pearson r (978 landmark genes, specific profiles) between d's profile in L and the plain
mean of d's profiles in its other lines, i.e. exactly the borrowed estimate of 79e.
Features (per held-out pair in training; per disease x drug when applied):
  cons     drug consistency: mean r_mean of the drug's OTHER held-out lines (training: leave L out; applied: all
           of the drug's held-out lines)
  log_k    log2 of the number of pairs cons is based on (few pairs = noisy cons)
  log_n    log2 of the number of lines averaged (training: n_other_lines; applied: n_lines_used from 79e)
Model: OLS of r_mean on cons, log_k, cons * log_k, log_n; prediction clipped to [0, 1]. CEILING = r of two
single signatures of the same drug in the same line (extrapolation_generic_check.tsv); averaged profiles can exceed it.
Checked by 5-fold CV grouped by drug and leave-one-held-out-line-out, against an n-lines-only model; calibration also
reported on the pairs whose nearest tested line is in another lineage (the borrowed situation).
Distance to the tested lines (Spearman -0.24 with r_mean) is not used: the committed Celligner tables keep only 10
neighbors per line, so disease-to-tested-line distance cannot be computed for borrowed drugs. CV reports what it would
add (dist_nearest, held-out line to nearest tested line).
conf = min(predicted r / CEILING, 1): 1 = transfers as well as a replicate signature in the same line; direct drugs get
conf = 1 (assumed, not modelled).
rev_conf = rev * conf: a borrowed score shrunk toward 0 by its expected transfer accuracy. rev_conf_pct and
tier_rev_conf as rev_pct / tier_rev in 79e (approved + passes gate + rev_conf > 0 + top TOP_PCT percent).
Outputs: results/lincs/transfer_confidence_model.tsv (coefficients and CV), transfer_confidence_calibration.tsv,
         drug_transfer_confidence.tsv (disease x drug), transfer_confidence_summary.tsv (per disease)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
TOP_PCT, N_FOLDS, SEED, TOP_N = 10, 5, 0, 10


def design(cons, log_k, log_n, dist=None):
    cols = [np.ones_like(cons), cons, log_k, cons * log_k, log_n]
    if dist is not None:
        cols.append(dist)
    return np.column_stack(cols)


def fit(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


def cv(df, X, groups, rng):
    """Out-of-fold predictions, folds grouped by `groups`."""
    g = pd.Series(groups).astype("category").cat.codes.values
    if g.max() + 1 > N_FOLDS:
        fold = rng.permutation(g.max() + 1)[g] % N_FOLDS
    else:
        fold = g
    pred = np.empty(len(df))
    for f in np.unique(fold):
        tr = fold != f
        pred[~tr] = X[~tr] @ fit(X[tr], df.r_mean.values[tr])
    return pred


def metrics(y, p):
    return {"r2": 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(),
            "spearman": pd.Series(y).corr(pd.Series(p), method="spearman")}


def main():
    rng = np.random.default_rng(SEED)
    ceiling = pd.read_csv(RES / "extrapolation_generic_check.tsv", sep="\t").set_index("profiles")
    CEILING = float(ceiling.filter(like="ceiling", axis=0).iloc[0, 0])
    p = pd.read_csv(RES / "extrapolation_lolo_pairs.tsv", sep="\t")
    nb = pd.read_csv(ROOT / "results/celligner/line_neighbors.tsv", sep="\t")
    lineage = nb.drop_duplicates("ModelID").set_index("ModelID").OncotreeLineage

    # training rows: drug held out in >= 2 lines, so cons can leave L out
    g = p.groupby("pert_id").r_mean
    k = g.transform("count")
    p["cons"] = (g.transform("sum") - p.r_mean) / (k - 1)
    p["log_k"] = np.log2((k - 1).clip(lower=1))  # k = 1 rows are dropped below
    p["log_n"] = np.log2(p.n_other_lines)
    p["other_lineage"] = p.nearest_line.map(lineage).values != p.ModelID.map(lineage).values
    t = p[k > 1].reset_index(drop=True)
    y = t.r_mean.values

    X = design(t.cons.values, t.log_k.values, t.log_n.values)
    X_n = np.column_stack([np.ones(len(t)), t.log_n.values])
    X_d = design(t.cons.values, t.log_k.values, t.log_n.values, t.dist_nearest.values / 100)
    rows = []
    for name, XX in [("full", X), ("n_lines_only", X_n), ("full_plus_distance", X_d)]:
        for split, grp in [("by_drug", t.pert_id.values), ("by_heldout_line", t.ModelID.values)]:
            pr = np.clip(cv(t, XX, grp, rng), 0, 1)
            m = metrics(y, pr)
            mo = metrics(y[t.other_lineage.values], pr[t.other_lineage.values])
            rows.append({"model": name, "cv": split, "n_pairs": len(t), "r2": m["r2"], "spearman": m["spearman"],
                         "r2_other_lineage": mo["r2"], "spearman_other_lineage": mo["spearman"]})
        if name == "full":
            oof = np.clip(cv(t, X, t.pert_id.values, rng), 0, 1)
    beta = fit(X, y)
    model = pd.DataFrame(rows)
    coef = pd.DataFrame({"model": "coefficients", "cv": ["intercept", "cons", "log_k", "cons_x_log_k", "log_n"],
                         "r2": beta})
    pd.concat([model, coef]).round(4).to_csv(RES / "transfer_confidence_model.tsv", sep="\t", index=False)

    # calibration (out-of-fold, grouped by drug)
    t["pred"] = oof
    t["bin"] = pd.qcut(t.pred, 10, duplicates="drop")
    cal = (t.groupby(["other_lineage", "bin"], observed=True)
             .agg(n=("r_mean", "size"), pred_mean=("pred", "mean"), obs_median=("r_mean", "median"),
                  obs_mean=("r_mean", "mean"), frac_gt_0_3=("r_mean", lambda v: (v > 0.3).mean()))
             .reset_index())
    cal["bin"] = cal.bin.astype(str)
    cal.round(3).to_csv(RES / "transfer_confidence_calibration.tsv", sep="\t", index=False)

    # apply to every disease x drug
    dc = p.groupby("pert_id").r_mean.agg(cons="mean", k="count")
    r = pd.read_csv(RES / "drug_reversal.tsv", sep="\t", low_memory=False)
    r = r.join(dc, on="pert_id")
    b = r.source == "borrowed"
    assert r.loc[b, "cons"].notna().all(), "borrowed drug without leave-one-line-out data"
    rb = r[b]
    pred = design(rb.cons.values, np.log2(rb.k.values), np.log2(rb.n_lines_used.values)) @ beta
    r["pred_r"] = np.nan
    r.loc[b, "pred_r"] = np.clip(pred, 0, 1)
    r["conf"] = np.where(b, np.minimum(r.pred_r / CEILING, 1), 1.0)
    r["rev_conf"] = r.rev * r.conf
    # rev_pct in 79e = percent of drugs with a lower rev; same here
    r["rev_conf_pct"] = r.groupby("group_id").rev_conf.transform(
        lambda v: (v.rank(method="min") - 1) / len(v) * 100).round(1)
    top = r.rev_conf_pct >= 100 - TOP_PCT
    r["tier_rev_conf"] = np.where(r.approved & r.passes_gate & (r.rev_conf > 0) & top, 1,
                                  np.where(~r.approved & r.passes_gate & (r.rev_conf > 0) & top, 2, 3))
    r["rank_rev_conf"] = r.groupby("group_id").rev_conf.rank(ascending=False, method="first").astype(int)
    cols = ["group_id", "pert_id", "cmap_name", "source", "approved", "passes_gate", "n_lines_used", "cons", "k",
            "pred_r", "conf", "rev", "rev_pct", "rank_rev", "tier_rev", "rev_conf", "rev_conf_pct", "rank_rev_conf",
            "tier_rev_conf"]
    r[cols].round(4).to_csv(RES / "drug_transfer_confidence.tsv", sep="\t", index=False)

    summ = []
    for gid, d in r.groupby("group_id"):
        top_old = set(d[d.rank_rev <= TOP_N].pert_id)
        top_new = set(d[d.rank_rev_conf <= TOP_N].pert_id)
        t1_old, t1_new = d[d.tier_rev == 1], d[d.tier_rev_conf == 1]
        db = d[d.source == "borrowed"]
        summ.append({"group_id": gid, "n_drugs": len(d), "pct_borrowed": round(100 * len(db) / len(d), 1),
                     "conf_median_borrowed": db.conf.median(),
                     f"top{TOP_N}_kept": len(top_old & top_new),
                     f"top{TOP_N}_borrowed_before": int((d[d.rank_rev <= TOP_N].source == "borrowed").sum()),
                     f"top{TOP_N}_borrowed_after": int((d[d.rank_rev_conf <= TOP_N].source == "borrowed").sum()),
                     "tier1_before": len(t1_old), "tier1_after": len(t1_new),
                     "tier1_kept": len(set(t1_old.pert_id) & set(t1_new.pert_id)),
                     "tier1_borrowed_conf_lt_0_5_before": int(((t1_old.source == "borrowed") & (t1_old.conf < 0.5)).sum())})
    pd.DataFrame(summ).round(3).to_csv(RES / "transfer_confidence_summary.tsv", sep="\t", index=False)

    print(model.round(3).to_string(index=False))
    print("coefficients:", dict(zip(["intercept", "cons", "log_k", "cons_x_log_k", "log_n"], beta.round(4))))
    print(f"ceiling {CEILING}; borrowed conf:", r.loc[b, "conf"].describe().round(3).to_dict())


if __name__ == "__main__":
    main()
