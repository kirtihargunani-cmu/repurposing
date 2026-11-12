#!/usr/bin/env python
"""Drug-side dimension diagnostic: how many effective dimensions does the regulator space have across DRUGS?

REPORTED ONLY. No threshold, no pass/fail, no decision number, nothing in production is read or written. This
answers one descriptive question that has never been measured: the patient-side effective dimension is known
(scripts/86_subgroup_variability.R -> results/step10_v2/variability.tsv, median 6.51 over the 41 cohorts), but the
same quantity across drugs is not. It decides whether any condensing proposal is worth a decision at all.

Method mirrors 86 exactly so the two numbers are comparable:
  PCA on the drugs x regulators matrix, regulators centered;
  effective_dim = participation ratio of the eigenvalues, (sum lambda)^2 / sum(lambda^2);
  coherence = median over drugs of spearman(drug, leave-one-out mean of the others)   [86's loo_rho, drugs for patients]

Two arms, because cohorts differ in both drug count (1,144-1,871) and regulator count (1,003-1,223), and both
shift an eigenvalue spectrum on their own:
  own     each cohort's full matrix, in its own space -- descriptive, not comparable across cohorts
  common  the 1,144 pert_ids and 571 regulators present in ALL 55 cohorts -- identical shape everywhere, so
          differences are biology and not sample size. 86 needed a fixed-size subsample for this; here the
          common frame is exact, so no subsampling and no RNG.

Input is the 79e consensus drug activity dumped by D126-B1 (D126_ACT_DIR), i.e. the same quantity the score uses.

Outputs (new paths, nothing production collides):
  results/lincs/diagnostics/drug_space_dimension.tsv      one row per cohort x arm
  results/lincs/diagnostics/drug_space_pc1_loadings.tsv   top regulators by |PC1 loading|, common arm
  results/lincs/diagnostics/drug_space_pc1_alignment.tsv  pairwise |cos| between cohorts' PC1, common arm
  results/lincs/diagnostics/drug_space_run.json
"""
import os

# --- config block: every constant here and nowhere else -------------------------------------------------------
THREADS = 8                      # CPU rule: well under half of the 64 available
VAR_LEVELS = (0.5, 0.8, 0.9, 0.95)   # report the PC count needed to reach each share of variance
TOP_N = 25                       # regulators listed per cohort by |PC1 loading|
CENTER = True                    # center regulators across drugs, as 86's prcomp(center = TRUE)
SCALE = False                    # no scaling: the activities are already on one NES-like scale
WATCH = ("ESR1", "AR", "ESR2")   # receptors whose PC1 rank is worth printing; reported, nothing depends on them

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = str(THREADS)

import glob  # noqa: E402
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import rankdata  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
ACT_DIR = ROOT / "data/processed/lincs/d126_act"
OUT = ROOT / "results/lincs/diagnostics"
PATIENT = ROOT / "results/step10_v2/variability.tsv"


def load(g):
    """The cohort's consensus drug activity: drugs x regulators, with its drug and regulator labels."""
    regs = (ACT_DIR / f"{g}.regulators.txt").read_text().split("\n")
    drugs = pd.read_csv(ACT_DIR / f"{g}.drugs.tsv", sep="\t")
    D = np.fromfile(ACT_DIR / f"{g}.f32", dtype=np.float32).reshape(-1, len(regs)).astype(np.float64)
    return D, np.array(regs), drugs


def spectrum(D):
    """Eigenvalues of the regulator covariance across drugs, via SVD. Returns (lambda, PC1 loading, PC1 score)."""
    X = D - D.mean(0) if CENTER else D.copy()
    if SCALE:
        s = X.std(0)
        X = X / np.where(s > 0, s, 1.0)
    U, s, Vt = np.linalg.svd(X, full_matrices=False)
    lam = s ** 2 / max(len(X) - 1, 1)
    return lam, Vt[0], U[:, 0] * s[0]


def srho(a, b):
    a, b = rankdata(a) - (len(a) + 1) / 2, rankdata(b) - (len(b) + 1) / 2
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def coherence(D):
    """86's loo_rho with drugs in place of patients: each drug vs the mean of all the others."""
    tot = D.sum(0)
    n = len(D)
    out = np.empty(n)
    for j in range(n):
        ref = (tot - D[j]) / (n - 1)
        a, b = rankdata(D[j]), rankdata(ref)
        a -= a.mean()
        b -= b.mean()
        out[j] = float(a @ b / np.sqrt((a @ a) * (b @ b)))
    return out


def stats(D, lam, score):
    tot = lam.sum()
    ve = lam / tot
    cum = np.cumsum(ve)
    row = dict(n_drugs=D.shape[0], n_regulators=D.shape[1],
               effective_dim=round(float(tot ** 2 / (lam ** 2).sum()), 2),
               pc1_var=round(float(ve[0]), 3), pc2_var=round(float(ve[1]), 3), pc3_var=round(float(ve[2]), 3))
    for q in VAR_LEVELS:
        row[f"pcs_to_{int(q * 100)}"] = int(np.searchsorted(cum, q) + 1)
    coh = coherence(D)
    row["coherence"] = round(float(np.median(coh)), 3)
    row["coherence_q10"] = round(float(np.quantile(coh, 0.1)), 3)
    # is PC1 a direction, or just "how hard the drug hits"? |PC1 score| vs the drug's own activity magnitude
    row["pc1_amplitude_rho"] = round(srho(np.abs(score), np.linalg.norm(D, axis=1)), 3)
    return row


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    groups = sorted(os.path.basename(f)[:-len(".regulators.txt")]
                    for f in glob.glob(str(ACT_DIR / "*.regulators.txt")))
    print(f"{len(groups)} cohorts from {ACT_DIR.relative_to(ROOT)}")

    # the common frame: pert_ids and regulators present in every cohort
    shared_r, shared_d = None, None
    for g in groups:
        r = set((ACT_DIR / f"{g}.regulators.txt").read_text().split("\n"))
        d = set(pd.read_csv(ACT_DIR / f"{g}.drugs.tsv", sep="\t").pert_id)
        shared_r = r if shared_r is None else shared_r & r
        shared_d = d if shared_d is None else shared_d & d
    shared_r = sorted(shared_r)
    shared_d = sorted(shared_d)
    print(f"common frame: {len(shared_d)} drugs x {len(shared_r)} regulators in all {len(groups)} cohorts")

    rows, loads, pc1 = [], [], {}
    for g in groups:
        D, regs, drugs = load(g)
        lam, v1_own, sc_own = spectrum(D)
        own_row = dict(cohort=g, arm="own", **stats(D, lam, sc_own))
        for w in WATCH:          # where the receptor sits on PC1 in the cohort's OWN space, if it is scored there
            j = np.where(regs == w)[0]
            own_row[f"{w.lower()}_pc1_rank"] = (int(np.where(np.argsort(-np.abs(v1_own)) == j[0])[0][0]) + 1
                                                if len(j) else None)
        rows.append(own_row)

        ri = pd.Index(regs).get_indexer(shared_r)
        di = pd.Index(drugs.pert_id).get_indexer(shared_d)
        assert (ri >= 0).all() and (di >= 0).all(), g
        C = D[np.ix_(di, ri)]
        lam_c, v1, sc_c = spectrum(C)
        rows.append(dict(cohort=g, arm="common", **stats(C, lam_c, sc_c)))

        if v1[np.argmax(np.abs(v1))] < 0:        # sign is arbitrary; orient each PC1 by its largest loading
            v1 = -v1
        pc1[g] = v1
        o = np.argsort(-np.abs(v1))[:TOP_N]
        loads += [dict(cohort=g, rank=k + 1, regulator=shared_r[i], loading=round(float(v1[i]), 4))
                  for k, i in enumerate(o)]
        rk = {w: (int(np.where(np.argsort(-np.abs(v1)) == shared_r.index(w))[0][0]) + 1
                  if w in shared_r else None) for w in WATCH}
        print(f"  {g:30s} own eff_dim {rows[-2]['effective_dim']:6.2f}  common eff_dim {rows[-1]['effective_dim']:6.2f}"
              f"  pc1 {rows[-1]['pc1_var']:.3f}  " + " ".join(f"{w} #{rk[w]}" for w in WATCH if rk[w]), flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "drug_space_dimension.tsv", sep="\t", index=False)
    pd.DataFrame(loads).to_csv(OUT / "drug_space_pc1_loadings.tsv", sep="\t", index=False)

    # how much do the cohorts share a leading axis
    M = np.array([pc1[g] for g in groups])
    M = M / np.linalg.norm(M, axis=1, keepdims=True)
    cos = np.abs(M @ M.T)
    iu = np.triu_indices(len(groups), 1)
    pd.DataFrame(dict(cohort_a=[groups[i] for i in iu[0]], cohort_b=[groups[j] for j in iu[1]],
                      abs_cos=np.round(cos[iu], 4))).to_csv(OUT / "drug_space_pc1_alignment.tsv",
                                                            sep="\t", index=False)

    com = out[out.arm == "common"]
    own = out[out.arm == "own"]
    summary = dict(
        n_cohorts=len(groups), common_drugs=len(shared_d), common_regulators=len(shared_r),
        center=CENTER, scale=SCALE,
        own=dict(effective_dim_median=float(own.effective_dim.median()),
                 effective_dim_min=float(own.effective_dim.min()), effective_dim_max=float(own.effective_dim.max()),
                 pc1_var_median=float(own.pc1_var.median())),
        common=dict(effective_dim_median=float(com.effective_dim.median()),
                    effective_dim_min=float(com.effective_dim.min()), effective_dim_max=float(com.effective_dim.max()),
                    pc1_var_median=float(com.pc1_var.median()),
                    pcs_to_80_median=float(com.pcs_to_80.median()),
                    pcs_to_90_median=float(com.pcs_to_90.median()),
                    coherence_median=float(com.coherence.median()),
                    pc1_amplitude_rho_median=float(com.pc1_amplitude_rho.median())),
        own_pc1_amplitude_rho_median=float(own.pc1_amplitude_rho.median()),
        esr1_pc1_rank=dict(n_cohorts_scoring_esr1=int(own.esr1_pc1_rank.notna().sum()),
                           median_rank=(float(own.esr1_pc1_rank.median())
                                        if own.esr1_pc1_rank.notna().any() else None),
                           best_rank=(float(own.esr1_pc1_rank.min())
                                      if own.esr1_pc1_rank.notna().any() else None)),
        pc1_alignment=dict(median_abs_cos=round(float(np.median(cos[iu])), 4),
                           min_abs_cos=round(float(cos[iu].min()), 4),
                           frac_pairs_above_0_8=round(float((cos[iu] > 0.8).mean()), 4)))
    if PATIENT.exists():                        # the comparison this diagnostic exists to make
        v = pd.read_csv(PATIENT, sep="\t")
        c = v[v.level == "cohort"]
        summary["patient_side_for_comparison"] = dict(
            source=str(PATIENT.relative_to(ROOT)), n=int(len(c)),
            effective_dim_median=float(c.effective_dim.median()),
            eff_dim_n30_median=float(c.eff_dim_n30.median()), pc1_var_median=float(c.pc1_var.median()),
            note="across patients within a cohort, 41 cohorts; the drug-side numbers are across drugs")
    (OUT / "drug_space_run.json").write_text(json.dumps(summary, indent=1))

    pd.set_option("display.width", 200)
    print("\n=== common frame (same 1,144 drugs x 571 regulators in every cohort) ===")
    print(com[["effective_dim", "pc1_var", "pcs_to_50", "pcs_to_80", "pcs_to_90", "pcs_to_95",
               "coherence", "pc1_amplitude_rho"]].describe().round(3).to_string())
    print("\nhighest and lowest effective dimension:")
    print(pd.concat([com.nsmallest(5, "effective_dim"), com.nlargest(5, "effective_dim")])
          [["cohort", "effective_dim", "pc1_var", "pcs_to_80", "coherence"]].to_string(index=False))
    print("\nPC1 agreement across cohorts: median |cos| "
          f"{summary['pc1_alignment']['median_abs_cos']}, min {summary['pc1_alignment']['min_abs_cos']}, "
          f"{summary['pc1_alignment']['frac_pairs_above_0_8']:.0%} of pairs above 0.8")
    top = pd.DataFrame(loads).groupby("regulator").size().sort_values(ascending=False)
    print(f"\nregulators in the top {TOP_N} of PC1 most often (of {len(groups)} cohorts):")
    print(", ".join(f"{r} {n}" for r, n in top.head(15).items()))
    e = own[own.esr1_pc1_rank.notna()]
    print(f"\nESR1 is a scored regulator in {len(e)} of {len(groups)} cohorts; its rank on PC1 by |loading| "
          f"(own space): median {e.esr1_pc1_rank.median():.0f}, best {e.esr1_pc1_rank.min():.0f}, "
          f"worst {e.esr1_pc1_rank.max():.0f}")
    br = e[e.cohort.str.startswith("BRCA")]
    if len(br):
        print("  breast cohorts: " + ", ".join(f"{r.cohort} #{int(r.esr1_pc1_rank)}/{r.n_regulators}"
                                               for r in br.itertuples()))
    print(f"\n|PC1 score| vs the drug's own activity magnitude: median Spearman "
          f"{own.pc1_amplitude_rho.median():.3f} (own), {com.pc1_amplitude_rho.median():.3f} (common)")
    if "patient_side_for_comparison" in summary:
        p = summary["patient_side_for_comparison"]
        print(f"\npatient side for comparison: effective_dim median {p['effective_dim_median']}, "
              f"pc1_var median {p['pc1_var_median']} (across patients, {p['n']} cohorts)")
    print(f"\nwritten to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
