#!/usr/bin/env python
"""C1: does one NES value mean several things for a master regulator? (D202)

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_interpretation_validation.md section 4 (pre-registered; read it before changing
anything here)

THE PREMISE THE WHOLE GENE-LEVEL LAYER RESTS ON. D179/D194/D195 keep a per-gene breakdown inside
every regulator because the same activity value was believed reachable by different routes. That
belief has never been tested directly. If it is false - if equal activity implies equal gene
pattern - then the nested arrays are machinery, not information, and this script says so.

DESIGN, as registered.

  unit                (cohort, window, regulator) with >= MIN_DRUGS drugs scored
  d_act               |activity_i - activity_j| for a pair of drugs on that regulator
  s_gene              cosine(contribution vector_i, contribution vector_j) over the regulon's
                      targets - the gene-level pattern similarity
  matched / unmatched d_act in the LOWEST / HIGHEST decile FOR THAT REGULATOR (a relative
                      tolerance, so regulators with different activity spreads are treated alike)

THE CEILING COMES FIRST. A vector measured once cannot be matched better than it matches itself.
`splithalf_reliability.tsv.gz` (79e3 SPLITHALF=1) splits each drug's conditions in half and takes
the per-regulon cosine between the halves; its median is R. Section 4.2: if R is low the test is
INCONCLUSIVE, not confirmatory - the instrument would simply be unable to resolve gene pattern.

THE CONFOUND, NAMED IN THE SPEC BEFORE RUNNING. Both contribution vectors carry the same edge
weights tfmode x likelihood, so their product carries w^2 and EVERY pair shows elevated cosine.
That is why nothing here is judged on an absolute cosine: the criteria are matched-vs-ceiling and
the rank association, both relative.

ACCEPTANCE (SPEC 4.4) - C1 is SUPPORTED when BOTH hold:
  1. median s_gene among activity-matched pairs < 0.60 * R
  2. |Spearman(d_act, s_gene)| < 0.30, pooled as the median over regulators
Failing either means the gene layer is largely redundant with the NES.

Output: results/lincs/interpretation/c1_nes_ambiguity.tsv (+ _pairs sample, run.json)
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
CDIR = RES / "regulon_contributions_v2_tf"
OUT = RES / "interpretation"

MIN_DRUGS = int(os.environ.get("MIN_DRUGS", "8"))
DECILE = 0.10
MAX_PAIRS = int(os.environ.get("MAX_PAIRS", "4000"))   # per unit, sampled with a fixed seed
SEED = 0


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- the ceiling, before anything else (SPEC 4.2)
    shf = CDIR / "splithalf_reliability.tsv.gz"
    if not shf.exists():
        print(f"ABORT: {shf} missing - run 79e3 with SPLITHALF=1 first.", file=sys.stderr)
        return 1
    SH = pd.read_csv(shf, sep="\t", low_memory=False)
    R_all = float(SH.cos_half.dropna().median())
    R_by = SH.dropna(subset=["cos_half"]).groupby(["cohort", "window"]).cos_half.median()
    print(f"CEILING R (median split-half gene-vector cosine) = {R_all:.4f} "
          f"over {len(SH):,} (drug, window, regulator) rows")

    # ---- the comparison
    cols = ["cohort", "window", "regulator", "pert_id", "cmap_name", "target", "activity_contrib"]
    T = pd.read_csv(CDIR / "target_contributions_by_window.tsv.gz", sep="\t",
                    usecols=lambda c: c in cols, low_memory=False)
    # the regulator-level activity each gene vector sums to
    A = (pd.read_csv(CDIR / "regulator_contributions_by_window.tsv.gz", sep="\t",
                     usecols=["cohort", "pert_id", "window", "regulator", "drug_activity"],
                     low_memory=False).drop_duplicates(["cohort", "pert_id", "window", "regulator"]))

    rng = np.random.default_rng(SEED)
    rows, pair_sample = [], []
    for (coh, win, reg), g in T.groupby(["cohort", "window", "regulator"], sort=False):
        drugs = g.pert_id.unique()
        if len(drugs) < MIN_DRUGS:
            continue
        M = g.pivot_table(index="pert_id", columns="target", values="activity_contrib")
        M = M.dropna(axis=1, how="any")
        if M.shape[1] < 15:
            continue
        act = (A[(A.cohort == coh) & (A.window == win) & (A.regulator == reg)]
               .set_index("pert_id").drug_activity.reindex(M.index))
        if act.isna().any():
            continue
        X = M.to_numpy(dtype=np.float64)
        nrm = np.linalg.norm(X, axis=1)
        ok = nrm > 0
        if ok.sum() < MIN_DRUGS:
            continue
        X, nrm, a = X[ok], nrm[ok], act.to_numpy()[ok]
        S = (X @ X.T) / np.outer(nrm, nrm)                       # gene-pattern similarity
        Dv = np.abs(a[:, None] - a[None, :])                     # activity difference
        iu = np.triu_indices(len(a), k=1)
        d_act, s_gene = Dv[iu], S[iu]
        if len(d_act) > MAX_PAIRS:
            sel = rng.choice(len(d_act), MAX_PAIRS, replace=False)
            d_act, s_gene = d_act[sel], s_gene[sel]
        lo, hi = np.quantile(d_act, DECILE), np.quantile(d_act, 1 - DECILE)
        m, u = s_gene[d_act <= lo], s_gene[d_act >= hi]
        if len(m) < 5 or len(u) < 5:
            continue
        rho = pd.Series(d_act).corr(pd.Series(s_gene), method="spearman")
        rows.append(dict(cohort=coh, window=win, regulator=reg, n_drugs=int(ok.sum()),
                         n_genes=int(M.shape[1]), n_pairs=int(len(d_act)),
                         ceiling=float(R_by.get((coh, win), R_all)),
                         s_matched=float(np.median(m)), s_unmatched=float(np.median(u)),
                         rho_dact_sgene=float(rho) if rho == rho else np.nan))
        if len(pair_sample) < 200000:
            k = min(50, len(d_act))
            s2 = rng.choice(len(d_act), k, replace=False)
            pair_sample.append(pd.DataFrame(dict(cohort=coh, window=win, regulator=reg,
                                                 d_act=d_act[s2], s_gene=s_gene[s2])))

    C = pd.DataFrame(rows)
    if C.empty:
        print("no usable units", file=sys.stderr)
        return 1
    C["s_matched_over_ceiling"] = C.s_matched / C.ceiling
    C.to_csv(OUT / "c1_nes_ambiguity.tsv", sep="\t", index=False, float_format="%.5g")
    pd.concat(pair_sample, ignore_index=True).to_csv(
        OUT / "c1_pairs_sample.tsv.gz", sep="\t", index=False, float_format="%.5g")

    med_matched = float(C.s_matched.median())
    med_ratio = float(C.s_matched_over_ceiling.median())
    med_rho = float(C.rho_dact_sgene.median())
    c1 = bool(med_ratio < 0.60)
    c2 = bool(abs(med_rho) < 0.30)

    print(f"\n{len(C):,} units (cohort, window, regulator), "
          f"{C.cohort.nunique()} cohorts, median {int(C.n_drugs.median())} drugs/unit")
    print(f"  median s_gene, activity-MATCHED pairs   {med_matched:.4f}")
    print(f"  median s_gene, activity-UNMATCHED pairs {C.s_unmatched.median():.4f}")
    print(f"  matched / ceiling                        {med_ratio:.4f}")
    print(f"  median Spearman(d_act, s_gene)           {med_rho:+.4f}")
    print(f"\nACCEPTANCE (SPEC 4.4)")
    print(f"  1. matched < 0.60 x ceiling : {med_ratio:.3f} -> {'PASS' if c1 else 'FAIL'}")
    print(f"  2. |rho| < 0.30             : {abs(med_rho):.3f} -> {'PASS' if c2 else 'FAIL'}")
    print(f"  C1 {'SUPPORTED' if (c1 and c2) else 'NOT SUPPORTED'}")
    if R_all < 0.3:
        print("  NOTE: ceiling is low - SPEC 4.2 says the test is INCONCLUSIVE, not confirmatory.")

    (OUT / "c1_run.json").write_text(json.dumps(dict(
        script="79am_c1_nes_ambiguity.py", spec="docs/SPEC_interpretation_validation.md#4",
        min_drugs=MIN_DRUGS, decile=DECILE, max_pairs=MAX_PAIRS, seed=SEED,
        ceiling_R=R_all, n_units=int(len(C)),
        median_s_matched=med_matched, median_s_unmatched=float(C.s_unmatched.median()),
        matched_over_ceiling=med_ratio, median_rho=med_rho,
        acceptance=dict(c1_matched_below_ceiling=c1, c2_weak_rho=c2,
                        supported=bool(c1 and c2)),
        inconclusive_low_ceiling=bool(R_all < 0.3)), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
