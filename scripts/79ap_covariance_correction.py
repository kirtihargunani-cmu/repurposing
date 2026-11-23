#!/usr/bin/env python
"""D205: rev + reliability-weighted gene covariance. Registered test.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_covariance_correction.md (pre-registered; read it before changing anything here)

    cov_i   = rev_gene_contrib_i - rev_contrib_i / n_i      ( = -n_i cov(d_i,C_i)/||s|| )
    rev_cov = rev + sum_i w_i * cov_i

w_i is MEASURED, not tuned: the median split-half cosine of that regulon's gene vector
(D202), clipped to [0,1], falling back to the cell median. A regulon whose gene vector does not
reproduce contributes nothing; one measured from 10+ conditions (reliability 0.79) contributes
nearly fully. There is no lambda and no threshold.

Judged on the D120 target-role components, which use no approval data and share no inputs with
`rev` (which carries no target annotation at all). Approval AUC is reported and NOT gated: D192
measured it bimodal and D196 showed the proliferation axis beats mechanism where approved sets are
cytotoxic.
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
CDIR = RES / "regulon_contributions_v2_tf"
OUT = RES / "interpretation"
MIN_ROWS, MIN_POS = 50, 3


def auc(y, s):
    m = ~np.isnan(s)
    y, s = np.asarray(y)[m], np.asarray(s)[m]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < MIN_POS or n0 < MIN_POS:
        return np.nan
    r = pd.Series(s).rank().to_numpy()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- w_i: the MEASURED reliability of each regulon's gene vector (SPEC s2)
    SH = pd.read_csv(CDIR / "splithalf_reliability.tsv.gz", sep="\t", low_memory=False,
                     usecols=["cohort", "window", "regulator", "cos_half"])
    W = (SH.dropna(subset=["cos_half"])
         .groupby(["cohort", "window", "regulator"], as_index=False).cos_half.median()
         .rename(columns={"cos_half": "w"}))
    W["w"] = W.w.clip(0, 1)
    cellmed = W.groupby(["cohort", "window"], as_index=False).w.median().rename(columns={"w": "w_cell"})
    W.to_csv(OUT / "d205_weights.tsv", sep="\t", index=False, float_format="%.5g")
    print(f"reliability weights: {len(W):,} (cohort, window, regulator); "
          f"median w {W.w.median():.3f}, cell-median fallback for the rest")

    # ---- the scores, streamed
    cols = ["cohort", "pert_id", "cmap_name", "window", "regulator", "n_genes_regulon",
            "rev_contrib", "rev_gene_contrib"]
    acc = []
    for ch in pd.read_csv(CDIR / "regulator_contributions_by_window_all.tsv.gz", sep="\t",
                          usecols=cols, chunksize=3_000_000, low_memory=False):
        ch["covt"] = ch.rev_gene_contrib - ch.rev_contrib / ch.n_genes_regulon.clip(lower=1)
        ch = ch.merge(W, on=["cohort", "window", "regulator"], how="left") \
               .merge(cellmed, on=["cohort", "window"], how="left")
        ch["w"] = ch.w.fillna(ch.w_cell).fillna(0.0)
        ch["wcov"] = ch.w * ch.covt
        k = ["cohort", "pert_id", "cmap_name", "window"]
        acc.append(ch.groupby(k, sort=False).agg(
            rev=("rev_contrib", "sum"), rev_gene=("rev_gene_contrib", "sum"),
            covt=("covt", "sum"), wcov=("wcov", "sum")).reset_index())
    S = (pd.concat(acc, ignore_index=True)
         .groupby(["cohort", "pert_id", "cmap_name", "window"], as_index=False)
         [["rev", "rev_gene", "covt", "wcov"]].sum())
    del acc
    S["rev_cov"] = S.rev + S.wcov                 # the registered score
    S["rev_cov_unw"] = S.rev + S.covt             # diagnostic only (SPEC s2)
    S.to_csv(OUT / "d205_rev_cov.tsv.gz", sep="\t", index=False, float_format="%.6g")

    # ---- anchors
    lab = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    lab = lab[lab.unit_kind == "group"][["unit_id", "pert_id", "label_phase"]]
    TR = pd.read_csv(RES / "target_role/drug_target_role.tsv", sep="\t", low_memory=False)
    TR["k"] = TR.drug.str.lower().str.strip()
    M = S.merge(lab, left_on=["cohort", "pert_id"], right_on=["unit_id", "pert_id"], how="left")
    M["y"] = (M.label_phase >= 3).astype(int)
    M["k"] = M.cmap_name.str.lower().str.strip()
    M = M.merge(TR[["cohort", "k", "T", "D", "A", "O"]], on=["cohort", "k"], how="left")

    # proliferation-program score, for criterion 4
    P = pd.read_csv(RES / "programs/program_drug.tsv.gz", sep="\t")
    F = pd.read_csv(RES / "programs/program_filter.tsv", sep="\t")
    pro = set(map(tuple, F[F.fail_proliferation][["cohort", "program"]].values))
    P = P[[(c, p) in pro for c, p in zip(P.cohort, P.program)]]
    PR = P.groupby(["cohort", "pert_id", "window"], as_index=False).opposition.sum() \
          .rename(columns={"opposition": "prolif"})
    M = M.merge(PR, on=["cohort", "pert_id", "window"], how="left")

    SCORES = ["rev", "rev_gene", "covt", "rev_cov", "rev_cov_unw"]
    rows, dropped = [], []
    for (c, w), g in M.groupby(["cohort", "window"]):
        if len(g) < MIN_ROWS or g.y.sum() < MIN_POS:
            dropped.append(dict(cohort=c, window=w, n=len(g), n_pos=int(g.y.sum())))
            continue
        r = dict(cohort=c, window=w, n=len(g), n_pos=int(g.y.sum()))
        for nm in SCORES:
            r["auc_" + nm] = auc(g.y.values, g[nm].values.astype(float))
            for comp in ("T", "D", "A", "O"):
                gg = g.dropna(subset=[comp])
                r[f"{comp}_{nm}"] = (gg[nm].corr(gg[comp], method="spearman")
                                     if len(gg) >= 20 and gg[nm].nunique() > 4 else np.nan)
            gp = g.dropna(subset=["prolif"])
            r["prolif_" + nm] = (gp[nm].corr(gp.prolif, method="spearman")
                                 if len(gp) >= 20 else np.nan)
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "d205_acceptance.tsv", sep="\t", index=False)

    def rep(col, base="rev"):
        i = R[[f"{col}_{base}", f"{col}_rev_cov"]].dropna()
        d = i[f"{col}_rev_cov"] - i[f"{col}_{base}"]
        p = binomtest(int((d > 0).sum()), len(d), 0.5).pvalue if len(d) else np.nan
        return i[f"{col}_{base}"].median(), i[f"{col}_rev_cov"].median(), d.median(), int((d > 0).sum()), len(d), p

    print(f"\n{len(R)} cells ({len(dropped)} dropped)\n")
    print(f"{'score':14s} " + " ".join(f"{c:>9s}" for c in ("AUC", "T", "D", "A", "O")))
    for nm in SCORES:
        print(f"{nm:14s} " + " ".join(
            f"{R[f'{p}_{nm}'].median():9.3f}" for p in ("auc", "T", "D", "A", "O")))

    print("\nACCEPTANCE (SPEC s4), rev_cov vs rev")
    v = {}
    for crit, col, bar in ((1, "T", 0.05), (2, "D", 0.03)):
        b, a, d, k, n, p = rep(col)
        ok = bool(d >= bar and k / n >= 0.60 and p < 0.05)
        v[f"c{crit}"] = ok
        print(f"  {crit}. {col}: {b:+.3f} -> {a:+.3f}  delta {d:+.3f} (bar {bar:+.2f}), "
              f"{k}/{n} ({100*k/n:.0f}%, bar 60%), p={p:.4f} -> {'PASS' if ok else 'FAIL'}")
    print("  3. both windows:")
    okw = True
    for w, g in R.groupby("window"):
        i = g[["T_rev", "T_rev_cov"]].dropna()
        d = (i.T_rev_cov - i.T_rev)
        pos = bool(len(d) and d.median() > 0)
        okw &= pos
        print(f"     {w:5s} n={len(i):2d}  T {i.T_rev.median():+.3f} -> {i.T_rev_cov.median():+.3f} "
              f"(delta {d.median():+.3f}) -> {'ok' if pos else 'NOT POSITIVE'}")
    v["c3"] = okw
    i = R[["prolif_rev", "prolif_rev_cov"]].dropna()
    dp = (i.prolif_rev_cov - i.prolif_rev).median()
    v["c4"] = bool(dp <= 0.10)
    print(f"  4. proliferation guard: {i.prolif_rev.median():+.3f} -> {i.prolif_rev_cov.median():+.3f} "
          f"(delta {dp:+.3f}, bar <= +0.10) -> {'PASS' if v['c4'] else 'FAIL'}")
    v["passed"] = bool(all(v.values()))
    print(f"\n  D205 {'PASSES' if v['passed'] else 'FAILS'} as registered")
    bu, au, du, ku, nu, _ = rep("T", "rev")
    iu = R[["T_rev", "T_rev_cov_unw"]].dropna()
    print(f"  (diagnostic, not a result) unweighted rev+cov: T {iu.T_rev_cov_unw.median():+.3f} "
          f"vs weighted {au:+.3f} — weighting {'helps' if au > iu.T_rev_cov_unw.median() else 'does not help'}")

    (OUT / "d205_run.json").write_text(json.dumps(dict(
        script="79ap_covariance_correction.py", spec="docs/SPEC_covariance_correction.md",
        n_cells=int(len(R)), dropped=dropped, acceptance=v), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
