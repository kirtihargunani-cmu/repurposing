#!/usr/bin/env python
"""Can a drug's missing timepoint window be imputed? (D183)

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_window_imputation.md  -- PRE-REGISTERED, committed before this script was written.
Every threshold, baseline, split and metric below is from that document. Nothing here is tuned to a result.

Stages (all in one run, in this order, because later stages depend on earlier ones):
  1. RELIABILITY (spec section 5). Split-half within window on the 1,224 both-rows with >= 2 conditions in
     BOTH windows, N_SPLIT random halves each, score each half with the D155 formula, Spearman between
     halves across rows, Spearman-Brown to full length. Gives rel_early / rel_late and the attainable
     ceiling sqrt(rel_target). Computed FIRST and frozen: no model is fitted before it exists.
  2. BASELINES B0-B3 (section 6), per direction, on the grouped folds.
  3. MODELS M1-M4 (section 8), same folds.
  4. ACCEPTANCE (section 8), the four pre-registered criteria, per direction.

Directions are kept separate throughout and never pooled: late->early and early->late have different
coverage and different condition counts, and one pooled number would hide a one-sided failure.

Splits: grouped 5-fold by pert_id (primary) and leave-one-tissue-key-out (secondary, a GATING criterion -
the both-window pool spans 8 tissue keys while the rows to impute span 29 cohorts).

Env: N_SPLIT (20) | SEED (0) | N_CORES (6) | FOLDS (5)
Outputs: results/lincs/window_imputation/{reliability,baselines,models,acceptance,per_cohort}.tsv,
         predictions.tsv.gz, run.json
"""
import json
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parent.parent
PROC, RES = ROOT / "data/processed/lincs", ROOT / "results/lincs"
SFX = "_v2_tf"
VDIR = PROC / f"drug_viper_specific{SFX}_cond"
WIN = ROOT / "results/systems_validation/integrated/window_reversal.tsv.gz"
MR_DIR = ROOT / f"results/step6/aracne_project{SFX}"
OUT = RES / "window_imputation"
N_SPLIT = int(os.environ.get("N_SPLIT", "20"))
SEED = int(os.environ.get("SEED", "0"))
FOLDS = int(os.environ.get("FOLDS", "5"))
DIRS = [("late", "early"), ("early", "late")]          # (observed, target)

# ---- pre-registered acceptance thresholds (SPEC section 8). Do not edit to recover a number.
ACC_BEAT_B1 = 0.10        # Spearman must beat identity by >= this
ACC_CEIL_FRAC = 0.60      # and reach >= this fraction of the reliability ceiling
ACC_SIGN = 0.75           # and sign agreement >= this
ACC_TISSUE_DROP = 0.10    # and lose < this much Spearman leaving a tissue out


def load_common():
    idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
    idx["col"] = np.arange(len(idx))
    t = idx.pert_time.where(idx.pert_time >= 0)
    idx["window"] = np.where(t.isna(), "unknown", np.where(t <= 6, "early", np.where(t >= 24, "late", "other")))
    mr = pd.read_csv(MR_DIR / "group_mr.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id").key
    return idx, mr, lw, keys


def cohort_scorer(cohort, net, mr, lw, idx):
    """Returns score(cols) -> rev, exactly window_reversal's formula on a set of condition rows."""
    regs = pd.Index((VDIR / f"{net}.regulators.txt").read_text().split())
    A = np.memmap(VDIR / f"{net}.f32", dtype=np.float32, mode="r", shape=(len(idx), len(regs)))
    m = mr[(mr.group_id == cohort) & mr.regulator.isin(regs)]
    ci = regs.get_indexer(m.regulator)
    s = m.nes.values.astype(np.float64)
    norm = np.linalg.norm(s)
    wg = lw[(lw.group_id == cohort) & lw.in_key].set_index("ModelID").weight

    def score(sub):
        ww = wg.reindex(sub.ModelID).fillna(0).values
        ww = ww if ww.sum() > 0 else np.ones(len(sub))
        D = (ww[:, None] * np.asarray(A[sub.col.values][:, ci], dtype=np.float64)).sum(0) / ww.sum()
        return float(-(D @ s) / norm)

    return score, len(m)


def reliability(w, idx, mr, lw):
    """SPEC section 5. Split-half within window; Spearman-Brown to full length."""
    pool = w[(w.window_coverage == "both") & (w.n_conditions_early >= 2) & (w.n_conditions_late >= 2)]
    rng = np.random.default_rng(SEED)
    rows = []
    for cohort, g in pool.groupby("cohort"):
        net = g.network.iloc[0]
        if not (VDIR / f"{net}.f32").exists():
            continue
        score, _ = cohort_scorer(cohort, net, mr, lw, idx)
        wgset = set(lw[(lw.group_id == cohort) & lw.in_key].ModelID)
        sub_all = idx[idx.ModelID.isin(wgset)]
        for r in g.itertuples():
            for win in ("early", "late"):
                c = sub_all[(sub_all.pert_id == r.pert_id) & (sub_all.window == win)]
                if len(c) < 2:
                    continue
                for k in range(N_SPLIT):
                    p = rng.permutation(len(c))
                    h = len(c) // 2
                    a, b = c.iloc[p[:h]], c.iloc[p[h:]]        # use every condition; halves may differ by 1
                    if not len(a) or not len(b):
                        continue
                    rows.append(dict(cohort=cohort, pert_id=r.pert_id, window=win, rep=k,
                                     half_a=score(a), half_b=score(b), n_cond=len(c)))
        print(f"  reliability {cohort:22s} rows {len(rows):,}", flush=True)
    R = pd.DataFrame(rows)
    out = []
    for win, g in R.groupby("window"):
        per = [spearmanr(x.half_a, x.half_b).correlation for _, x in g.groupby("rep") if len(x) > 5]
        rh = float(np.nanmean(per))
        full = 2 * rh / (1 + rh) if rh > -1 else np.nan        # Spearman-Brown
        out.append(dict(window=win, n_rows=int(g.pert_id.nunique()), n_obs=int(len(g) / max(1, g.rep.nunique())),
                        r_half=round(rh, 4), reliability=round(full, 4),
                        ceiling=round(float(np.sqrt(max(full, 0))), 4),
                        r_half_sd=round(float(np.nanstd(per)), 4)))
    return pd.DataFrame(out), R


def features(w, keys):
    """Both-window rows with the section-8 feature set. No held-out value is ever a feature."""
    b = w[w.window_coverage == "both"].copy()
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    b["tissue"] = b.cohort.map(keys)
    b["moa"] = b.pert_id.map(ann.hub_moa).fillna("")
    b["n_targets"] = pd.to_numeric(b.pert_id.map(ann.hub_target).str.count(r"\|").add(1), errors="coerce").fillna(0)
    b["anticancer"] = b.pert_id.map(ann.atc).fillna("").str.contains("L01|L02").astype(int)
    nm = pd.read_csv(RES / f"drug_reversal{SFX}.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_mrs", "line_match", "source"]).drop_duplicates(["group_id", "pert_id"])
    b = b.merge(nm.rename(columns={"group_id": "cohort"}), on=["cohort", "pert_id"], how="left")
    b["moa_code"] = pd.factorize(b.moa)[0]
    b["tissue_code"] = pd.factorize(b.tissue)[0]
    return b


def metrics(pred, actual):
    ok = np.isfinite(pred) & np.isfinite(actual)
    if ok.sum() < 10:
        return dict(n=int(ok.sum()), spearman=np.nan, sign_agree=np.nan, mae_over_sd=np.nan)
    p, a = pred[ok], actual[ok]
    return dict(n=int(ok.sum()), spearman=round(float(spearmanr(p, a).correlation), 4),
                sign_agree=round(float((np.sign(p) == np.sign(a)).mean()), 4),
                mae_over_sd=round(float(np.abs(p - a).mean() / (a.std() + 1e-12)), 4))


def run_direction(b, obs, tgt, splits):
    """Baselines and models for one direction, over the given (name, iterable-of-folds) splits."""
    X_cols = [f"rev_{obs}", f"n_conditions_{obs}", f"n_conditions_{tgt}", f"n_doses_{obs}", "n_mrs", "line_match"]
    X3_cols = X_cols + ["n_targets", "anticancer", "moa_code", "tissue_code"]
    d = b.dropna(subset=[f"rev_{obs}", f"rev_{tgt}"]).reset_index(drop=True)
    y = d[f"rev_{tgt}"].values
    res, preds = [], []
    for split_name, folds in splits:
        store = {k: np.full(len(d), np.nan) for k in ("B0", "B1", "B2", "B3", "M1", "M2", "M3")}
        for tr, te in folds(d):
            # --- baselines
            med = d.iloc[tr].groupby("cohort")[f"rev_{tgt}"].median()
            store["B0"][te] = d.iloc[te].cohort.map(med).fillna(d.iloc[tr][f"rev_{tgt}"].median()).values
            store["B1"][te] = d.iloc[te][f"rev_{obs}"].values                       # identity
            lr = LinearRegression().fit(d.iloc[tr][[f"rev_{obs}"]], y[tr])
            store["B2"][te] = lr.predict(d.iloc[te][[f"rev_{obs}"]])
            dm = d.iloc[tr].groupby("pert_id")[f"rev_{tgt}"].mean()                 # drug mean elsewhere
            store["B3"][te] = d.iloc[te].pert_id.map(dm).fillna(d.iloc[tr][f"rev_{tgt}"].median()).values
            # --- models
            store["M1"][te] = lr.predict(d.iloc[te][[f"rev_{obs}"]])
            Xtr = d.iloc[tr][X_cols].astype(float).fillna(0); Xte = d.iloc[te][X_cols].astype(float).fillna(0)
            store["M2"][te] = LinearRegression().fit(Xtr, y[tr]).predict(Xte)
            g = HistGradientBoostingRegressor(random_state=SEED, max_iter=300, early_stopping=True,
                                              validation_fraction=0.2)
            g.fit(d.iloc[tr][X3_cols].astype(float).fillna(0), y[tr])
            store["M3"][te] = g.predict(d.iloc[te][X3_cols].astype(float).fillna(0))
        for k, v in store.items():
            res.append(dict(direction=f"{obs}->{tgt}", split=split_name, method=k, **metrics(v, y)))
        p = d[["cohort", "pert_id", "cmap_name", f"rev_{obs}", f"rev_{tgt}"]].copy()
        p["direction"], p["split"] = f"{obs}->{tgt}", split_name
        for k, v in store.items():
            p[k] = v
        preds.append(p)
    return pd.DataFrame(res), pd.concat(preds, ignore_index=True)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    idx, mr, lw, keys = load_common()
    w = pd.read_csv(WIN, sep="\t", low_memory=False)

    print("STAGE 1: reliability (frozen before any model is fitted)", flush=True)
    rel, raw = reliability(w, idx, mr, lw)
    rel.to_csv(OUT / "reliability.tsv", sep="\t", index=False)
    raw.to_csv(OUT / "reliability_halves.tsv.gz", sep="\t", index=False)
    print(rel.to_string(index=False), flush=True)
    ceil = rel.set_index("window").ceiling.to_dict()

    print("\nSTAGE 2-3: baselines and models", flush=True)
    b = features(w, keys)
    gk = GroupKFold(n_splits=FOLDS)
    def primary(d):
        return list(gk.split(d, groups=d.pert_id))
    def tissue_out(d):
        out = []
        for t in sorted(d.tissue.dropna().unique()):
            te = np.flatnonzero((d.tissue == t).values)
            tr = np.flatnonzero((d.tissue != t).values)
            if len(te) >= 20 and len(tr) >= 50:
                out.append((tr, te))
        return out
    splits = [("grouped_by_drug", primary), ("leave_one_tissue_out", tissue_out)]

    allres, allpred = [], []
    for obs, tgt in DIRS:
        r, p = run_direction(b, obs, tgt, splits)
        allres.append(r); allpred.append(p)
        print(r.to_string(index=False), flush=True)
    R = pd.concat(allres, ignore_index=True)
    P = pd.concat(allpred, ignore_index=True)
    R.to_csv(OUT / "models.tsv", sep="\t", index=False)
    P.to_csv(OUT / "predictions.tsv.gz", sep="\t", index=False, float_format="%.6g")

    print("\nSTAGE 4: acceptance (SPEC section 8, pre-registered)", flush=True)
    acc = []
    for obs, tgt in DIRS:
        dn = f"{obs}->{tgt}"
        pr = R[(R.direction == dn) & (R.split == "grouped_by_drug")].set_index("method")
        to = R[(R.direction == dn) & (R.split == "leave_one_tissue_out")].set_index("method")
        b1 = pr.loc["B1", "spearman"]
        best = pr.loc[["M1", "M2", "M3"], "spearman"].astype(float)
        bm = best.idxmax()
        c = ceil.get(tgt, np.nan)
        row = dict(direction=dn, best_model=bm, spearman=best.max(), b1_spearman=b1,
                   ceiling=c, sign_agree=pr.loc[bm, "sign_agree"],
                   spearman_tissue_out=to.loc[bm, "spearman"] if bm in to.index else np.nan)
        row["c1_beats_B1"] = bool(row["spearman"] - b1 >= ACC_BEAT_B1)
        row["c2_ceiling_frac"] = bool(c and np.isfinite(c) and row["spearman"] >= ACC_CEIL_FRAC * c)
        row["c3_sign"] = bool(row["sign_agree"] >= ACC_SIGN)
        row["c4_tissue"] = bool(np.isfinite(row["spearman_tissue_out"]) and
                                row["spearman"] - row["spearman_tissue_out"] < ACC_TISSUE_DROP)
        row["VERDICT"] = "PASS" if all(row[k] for k in ("c1_beats_B1", "c2_ceiling_frac", "c3_sign", "c4_tissue")) \
            else ("RESTRICTED" if row["c1_beats_B1"] and row["c3_sign"] and row["c2_ceiling_frac"] else "CLOSED")
        acc.append(row)
    A = pd.DataFrame(acc)
    A.to_csv(OUT / "acceptance.tsv", sep="\t", index=False)
    print(A.to_string(index=False), flush=True)

    json.dump(dict(spec="docs/SPEC_window_imputation.md", seed=SEED, n_split=N_SPLIT, folds=FOLDS,
                   thresholds=dict(beat_b1=ACC_BEAT_B1, ceiling_frac=ACC_CEIL_FRAC, sign=ACC_SIGN,
                                   tissue_drop=ACC_TISSUE_DROP),
                   reliability=rel.to_dict("records"), verdicts=A[["direction", "VERDICT"]].to_dict("records"),
                   n_both_rows=int((w.window_coverage == "both").sum()), runtime_s=round(time.time() - t0, 1)),
              open(OUT / "run.json", "w"), indent=1)
    print(f"\n-> {OUT}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
