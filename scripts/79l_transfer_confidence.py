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
VERSION=v2 (user 2026-10-28, D77; hybrid): run BEFORE 79e, which reads the result.
  excess distance (Celligner 70-PC space, data/processed/celligner/pcs_70.npy): how much further the drug's tested lines
    are than the best available line. Training pair (held-out line L, drug d): dist(L, nearest of d's other lines) -
    dist(L, nearest other DepMap cancer line). Disease x drug: per patient, dist(patient, nearest line the 79e estimate
    uses) - dist(patient, nearest DepMap cancer line); median over the disease's patients (weights_hybrid.tsv.gz).
  model: cons, log_k, cons * log_k, log_n, excess / 100 (CV with and without excess reported); conf as above, applied
    to borrowed drugs (lines = the drug's cancer lines outside the tissue key, as 79e); direct drugs conf = 1.
  low_transfer = borrowed and conf < CONF_MIN (fixed from the conf distribution before any ranking was looked at);
    79e keeps flagged drugs out of tiers 1-2 (as unmatched_lines, D76).
  Outputs: results/lincs/transfer_confidence_model_v2.tsv, transfer_confidence_calibration_v2.tsv,
    drug_transfer_confidence_v2.tsv (group_id, pert_id, source, n_lines, excess_dist, pred_r, conf, low_transfer)
Outputs: results/lincs/transfer_confidence_model.tsv (coefficients and CV), transfer_confidence_calibration.tsv,
         drug_transfer_confidence.tsv (disease x drug), transfer_confidence_summary.tsv (per disease)
"""
from pathlib import Path

import os

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
PROC = ROOT / "data/processed/lincs"
TOP_PCT, N_FOLDS, SEED, TOP_N = 10, 5, 0, 10
V2 = os.environ.get("VERSION", "v1") == "v2"
CONF_MIN = float(os.environ.get("CONF_MIN", "0.5"))   # v2 flag cutoff (D77): fixed 2026-10-28 from the conf distribution (no break; 11.5% of borrowed pairs < 0.5) before any ranking
CEL = ROOT / "data/processed/celligner"


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


def pcs():
    X = np.load(CEL / "pcs_70.npy"); rows = pd.Index(open(CEL / "pcs_70_rows.txt").read().split())
    return X, rows


def nearest(Xa, Xb, chunk=2000):
    """Euclidean distance matrix rows of Xa vs Xb (float32)."""
    nb = (Xb ** 2).sum(1)
    out = np.empty((len(Xa), len(Xb)), dtype=np.float32)
    for i in range(0, len(Xa), chunk):
        a = Xa[i:i + chunk]
        out[i:i + chunk] = np.sqrt(np.maximum((a ** 2).sum(1)[:, None] + nb[None] - 2 * a @ Xb.T, 0))
    return out


def main_v2():
    rng = np.random.default_rng(SEED)
    CEILING = float(pd.read_csv(RES / "extrapolation_generic_check.tsv", sep="\t").set_index("profiles")
                    .filter(like="ceiling", axis=0).iloc[0, 0])
    X, rows = pcs()
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    cancer = models.index[models.OncotreePrimaryDisease != "Non-Cancerous"]
    lines = rows[rows.isin(cancer)]; Xl = X[rows.get_indexer(lines)]
    lpos = pd.Series(np.arange(len(lines)), index=lines)
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    lines_of = idx.groupby("pert_id").ModelID.apply(lambda v: sorted(set(v) & set(lines)))

    # training pairs
    p = pd.read_csv(RES / "extrapolation_lolo_pairs.tsv", sep="\t")
    p = p[p.ModelID.isin(lines)].reset_index(drop=True)
    Dll = nearest(Xl, Xl); np.fill_diagonal(Dll, np.inf)
    best = Dll.min(1)
    ex = np.full(len(p), np.nan)
    for i, (m, d) in enumerate(zip(p.ModelID, p.pert_id)):
        o = [lpos[x] for x in lines_of.get(d, []) if x != m]
        if o:
            ex[i] = Dll[lpos[m], o].min() - best[lpos[m]]
    p["excess"] = ex
    g = p.groupby("pert_id").r_mean
    k = g.transform("count")
    p["cons"] = (g.transform("sum") - p.r_mean) / (k - 1)
    p["log_k"] = np.log2((k - 1).clip(lower=1)); p["log_n"] = np.log2(p.n_other_lines)
    t = p[(k > 1) & p.excess.notna()].reset_index(drop=True); y = t.r_mean.values
    X0 = design(t.cons.values, t.log_k.values, t.log_n.values)
    Xd = design(t.cons.values, t.log_k.values, t.log_n.values, t.excess.values / 100)
    res = []
    for name, XX in [("full", X0), ("full_plus_excess_distance", Xd)]:
        for split, grp in [("by_drug", t.pert_id.values), ("by_heldout_line", t.ModelID.values)]:
            res.append({"model": name, "cv": split, "n_pairs": len(t), **metrics(y, np.clip(cv(t, XX, grp, rng), 0, 1))})
    beta = fit(Xd, y)
    t["pred"] = np.clip(cv(t, Xd, t.pert_id.values, rng), 0, 1)
    t["bin"] = pd.qcut(t.pred, 10, duplicates="drop")
    cal = t.groupby("bin", observed=True).agg(n=("r_mean", "size"), pred_mean=("pred", "mean"), obs_mean=("r_mean", "mean")).reset_index()
    cal["bin"] = cal.bin.astype(str); cal.round(3).to_csv(RES / "transfer_confidence_calibration_v2.tsv", sep="\t", index=False)
    coef = pd.DataFrame({"model": "coefficients", "cv": ["intercept", "cons", "log_k", "cons_x_log_k", "log_n", "excess_per_100"], "r2": beta})
    mdl = pd.DataFrame(res); pd.concat([mdl, coef]).round(4).to_csv(RES / "transfer_confidence_model_v2.tsv", sep="\t", index=False)

    # disease x drug
    dc = p.groupby("pert_id").r_mean.agg(cons="mean", k="count")
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False)
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    wh = pd.read_csv(CEL / "weights_hybrid.tsv.gz", sep="\t", usecols=["sample", "group_id"]).drop_duplicates()
    out = []
    for gid, gd in dd.groupby("group_id"):
        k_ = keys.loc[gid]; km = models[models.index.isin(cancer) & models.OncotreeLineage.isin(str(k_.lineages).split(";"))]
        if isinstance(k_.diseases, str) and k_.diseases:
            km = km[km.OncotreePrimaryDisease.isin(k_.diseases.split(";"))]
        pts = rows[rows.isin(wh[wh.group_id == gid]["sample"])]
        Dp = nearest(X[rows.get_indexer(pts)], Xl); bp = Dp.min(1)
        for d, src in zip(gd.pert_id, gd.source):
            L = [x for x in lines_of.get(d, []) if (x in km.index) == (src == "direct")]
            e = float(np.median(Dp[:, [lpos[x] for x in L]].min(1) - bp)) if L and len(pts) else np.nan
            out.append((gid, d, src, len(L), e))
    r = pd.DataFrame(out, columns=["group_id", "pert_id", "source", "n_lines", "excess_dist"]).join(dc, on="pert_id")
    b = (r.source == "borrowed") & r.cons.notna() & r.excess_dist.notna()
    r["pred_r"] = np.nan
    r.loc[b, "pred_r"] = np.clip(design(r.cons[b].values, np.log2(r.k[b].values), np.log2(r.n_lines[b].clip(lower=1).values),
                                        r.excess_dist[b].values / 100) @ beta, 0, 1)
    r["conf"] = np.where(r.source == "borrowed", np.minimum(r.pred_r / CEILING, 1), 1.0)
    r["low_transfer"] = (r.source == "borrowed") & ~(r.conf >= CONF_MIN) if not np.isnan(CONF_MIN) else False
    r.round(4).to_csv(RES / "drug_transfer_confidence_v2.tsv", sep="\t", index=False)
    print(mdl.round(4).to_string(index=False))
    print("coefficients:", dict(zip(["intercept", "cons", "log_k", "cons_x_log_k", "log_n", "excess_per_100"], beta.round(4))))
    print("excess distance, training pairs:", t.excess.describe().round(1).to_dict())
    print("excess distance, borrowed disease x drug:", r.excess_dist[r.source == "borrowed"].describe().round(1).to_dict())
    print("borrowed conf:", r.conf[r.source == "borrowed"].describe(percentiles=[.05, .1, .25, .5]).round(3).to_dict())
    print("borrowed without leave-one-line-out data (conf NaN):", int(((r.source == "borrowed") & r.cons.isna()).sum()))


if __name__ == "__main__":
    main_v2() if V2 else main()
