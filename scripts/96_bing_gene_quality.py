#!/usr/bin/env python
"""Per-gene quality table for the LINCS BING gene set: reproducibility, redundancy with the landmarks, effective
independent information and variance compression (runs in `celligner`).

DIAGNOSTIC ONLY. Nothing is filtered, dropped, re-weighted or curated here; no scoring code is read or written; no
gene set is chosen. The filter design is a separate, later decision by the user.

Why (user, 2026-11-16): D128's landmark sensitivity found that 25% of scored units reverse the sign of their
mechanism score depending only on whether inferred genes are included. Two explanations have to be told apart:
(a) some inferred genes are simply noisy, (b) inferred genes are near-deterministic functions of the 978 landmarks
they were predicted from, so a 478-gene footprint carries far less independent information than 478, which would
both dilute scores toward zero and make permutation nulls falsely narrow. `geneinfo_beta.txt` carries only the
coarse three-tier `feature_space` and no per-gene inference-accuracy column, so every statistic here is derived
empirically from our own signatures.

EVERY STATISTIC IS LABEL-FREE: computed from the expression matrix and the experiment metadata (compound id, cell
line, dose, time) alone. No drug target annotation, no disease contrast, no M, no R, no D128 output is read.
All definitions and constants below were fixed in `docs/DECISIONS.md` (D129) before anything was computed.

Input: data/processed/lincs/level5_bing_hiq.npy (90,988 signatures x 10,174 genes, float32, memory-mapped; the
first 978 columns are the landmarks, the remaining 9,196 best inferred), level5_bing_genes.tsv,
level5_bing_hiq_sigs.tsv.gz (all from scripts/76b_lincs_level5_bing.py).

Definitions, fixed in advance:
  (1) REPLICATE GROUP = signatures sharing (pert_id, cell_iname, pert_dose, pert_time). Level 5 signatures are
      already MODZ-consolidated over level 4 replicates, so two signatures in such a group are the same compound at
      the same dose and time in the same line on different plates / batches: independent repeats of one experiment.
      Groups with >= 2 signatures are used; larger groups are randomly subsampled to CAP_REPLICATES so that a few
      very heavily repeated compounds cannot dominate. Two statistics, both standard, both per gene:
        icc1        one-way random-effects ICC(1) over the groups: between-group variance / total variance,
                    unbalanced k0 = (N - sum n_g^2 / N) / (G - 1). Can be <= 0 and is reported as computed.
        splithalf_r Pearson r, across groups, between the mean of a random half of each group and the mean of the
                    other half (for n_g = 2, one signature each). Seeded once; the same split for every gene.
      Both conflate signal with noise on purpose: a gene that never responds reproducibly is not reproducible.
      The landmarks are the internal benchmark; an inferred gene below the landmark distribution looks like noise.
  (2) REDUNDANCY = held-out R^2 of each gene on the 978 landmark values of the SAME signature, OLS with intercept.
      Split is BY PERTURBATION (all signatures of a pert_id fall on one side, SHA1 of pert_id, TRAIN_FRAC), so no
      replicate structure can leak across the split. Inferred genes are regressed on all 978 landmarks; as the
      control, each landmark is regressed on the OTHER 977 (leave-one-out coefficients read off the inverse of the
      train scatter matrix, which is exact OLS). R^2 on the test rows uses the train fit and the test mean; it can
      be negative and is reported as computed.
  (3) EFFECTIVE INDEPENDENT INFORMATION, from the second-order structure of the gene x gene correlation matrix over
      all signatures:
        participation ratio   (sum lambda)^2 / sum(lambda^2)   -- the same measure scripts/86 and
                              scripts/diagnostics/drug_space_dimension.py use, so the numbers are comparable
        entropy effective rank exp(-sum p log p), p = lambda / sum lambda
        number of components needed for VAR_LEVELS of the variance
      reported for BING (10,174), the landmarks alone (978) and the best-inferred alone (9,196), on the correlation
      matrix (genes standardised) and on the covariance matrix (86's convention, center only). Assumptions: only
      linear / second-order dependence is seen (a gene that is a nonlinear function of the landmarks still counts as
      independent); signatures are treated as exchangeable draws, which they are not exactly (compounds repeat); the
      statistic is a property of this signature set, not of L1000 in general.
      Per gene, two independent-variance measures:
        indep_var_frac_landmark = 1 - r2_landmark_heldout   (what the landmarks do not explain)
        uniqueness_all_bing     = 1 / inv(C + RIDGE*I)_jj on the correlation matrix C, i.e. 1 - R^2 of the gene on
                                  ALL 10,173 other genes. The ridge is needed because C is near-singular; it biases
                                  uniqueness UPWARD, so this is a lower bound on redundancy.
  (4) VARIANCE COMPRESSION = per-gene sd over all signatures, IQR and the p01-p99 range (both on a seeded subsample
      of QUANTILE_N signatures), mean |z| and max |z|. Level 5 values are already z-like, so the informative split is
      sd_between (reproducible signal, from (1)) against sd_within (replicate noise): a shrunk inference model
      compresses the signal, not necessarily the total.
Missing values are reported as missing with a reason; nothing is imputed.

Outputs: results/lincs/gene_quality/bing_gene_quality.tsv          one row per gene, every statistic + its n
         results/lincs/gene_quality/bing_gene_quality_summary.tsv  landmark vs best-inferred distributions
         results/lincs/gene_quality/bing_gene_quality_run.json     run parameters and the independence numbers
"""
import os

# --- config block: every constant here and nowhere else -------------------------------------------------------
THREADS = 24                 # CPU rule: well under half of the 64 available
ROW_BLOCK = 4000             # signatures read from the memory-map at a time
TRAIN_FRAC = 0.70            # share of PERTURBATIONS (not signatures) in the regression training split
SPLIT_SALT = "d129-bing"     # salt for the SHA1 perturbation hash, so the split is reproducible and arbitrary
CAP_REPLICATES = 10          # at most this many signatures per replicate group (seeded subsample)
MIN_GROUP = 2                # a replicate group needs at least this many signatures
SEED = 0                     # replicate subsampling, split-half assignment, quantile subsample
QUANTILE_N = 20000           # signatures in the seeded subsample used for IQR / p01 / p99
RIDGE = 1e-6                 # ridge on the correlation matrix before inversion (uniqueness_all_bing)
VAR_LEVELS = (0.5, 0.8, 0.9, 0.95)
N_LANDMARK = 978             # the first N_LANDMARK columns of the matrix are the landmarks (asserted below)

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = str(THREADS)

import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/lincs/gene_quality"
GROUP_KEYS = ["pert_id", "cell_iname", "pert_dose", "pert_time"]

t0 = time.time()


def log(msg):
    print(f"[{time.time() - t0:7.1f}s] {msg}", flush=True)


def eff_dim(lam, n_obs):
    """Effective dimension of a spectrum: participation ratio, entropy rank, components per variance level.

    pr_null_iid is what the participation ratio would be for the same number of INDEPENDENT genes measured on the
    same number of signatures (Marchenko-Pastur: p / (1 + p / n)); it is the sampling-noise reference, not a null
    that anything is tested against."""
    lam = np.clip(np.asarray(lam, dtype=float), 0, None)
    tot = lam.sum()
    p = lam[lam > 0] / tot
    cum = np.cumsum(np.sort(lam)[::-1]) / tot
    n_g = len(lam)
    out = {"n_genes": int(n_g), "participation_ratio": float(tot ** 2 / (lam ** 2).sum()),
           "entropy_eff_rank": float(np.exp(-(p * np.log(p)).sum())),
           "pr_null_iid": float(n_g / (1 + n_g / n_obs)),
           "min_eigenvalue": float(lam.min()), "max_eigenvalue": float(lam.max()),
           "top1_var_share": float(lam.max() / tot)}
    for v in VAR_LEVELS:
        out[f"pcs_to_{int(v * 100)}"] = int(np.searchsorted(cum, v) + 1)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    mat = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")
    n_sig, n_gene = mat.shape
    assert (genes.feature_space.values[:N_LANDMARK] == "landmark").all()
    assert (genes.feature_space.values[N_LANDMARK:] == "best inferred").all()
    assert len(genes) == n_gene and len(sig) == n_sig
    is_lm = np.arange(n_gene) < N_LANDMARK
    log(f"{n_sig} signatures x {n_gene} genes ({N_LANDMARK} landmark, {n_gene - N_LANDMARK} best inferred)")

    # --- perturbation-level train / test split ------------------------------------------------------------
    def in_train(pid):
        h = hashlib.sha1((SPLIT_SALT + str(pid)).encode()).digest()
        return int.from_bytes(h[:8], "big") / 2 ** 64 < TRAIN_FRAC

    perts = sig.pert_id.unique()
    train_pert = {p for p in perts if in_train(p)}
    tr_mask = sig.pert_id.isin(train_pert).values
    log(f"split by pert_id: {len(train_pert)}/{len(perts)} perturbations train, "
        f"{tr_mask.sum()}/{n_sig} signatures train")

    # --- pass 1: streaming moments (train, test) and NaN count --------------------------------------------
    g_tr = np.zeros((n_gene, n_gene)); g_te = np.zeros((n_gene, n_gene))
    s_tr = np.zeros(n_gene); s_te = np.zeros(n_gene)
    n_tr = n_te = 0
    n_nan = 0
    gmin = np.full(n_gene, np.inf); gmax = np.full(n_gene, -np.inf)
    buf = np.empty((n_gene, n_gene))
    for i in range(0, n_sig, ROW_BLOCK):
        blk = np.asarray(mat[i:i + ROW_BLOCK], dtype=np.float64)
        n_nan += int(np.isnan(blk).sum())
        np.minimum(gmin, blk.min(0), out=gmin); np.maximum(gmax, blk.max(0), out=gmax)
        m = tr_mask[i:i + ROW_BLOCK]
        for sel, G, S, in_tr in ((m, g_tr, s_tr, True), (~m, g_te, s_te, False)):
            sub = blk[sel]
            if not len(sub):
                continue
            np.dot(sub.T, sub, out=buf)
            G += buf
            S += sub.sum(0)
            if in_tr:
                n_tr += len(sub)
            else:
                n_te += len(sub)
        if (i // ROW_BLOCK) % 5 == 0:
            log(f"  moments {i + ROW_BLOCK:>6} / {n_sig}")
    del buf
    assert n_nan == 0, f"{n_nan} NaN in the matrix: the moment accumulators assume none"
    log(f"moments done, train {n_tr} test {n_te}")

    mu_tr = s_tr / n_tr
    g_all = g_tr + g_te
    s_all = s_tr + s_te
    mean_all = s_all / n_sig
    var_all = (np.diag(g_all) - n_sig * mean_all ** 2) / (n_sig - 1)
    sd_all = np.sqrt(var_all)

    # --- (2) redundancy: held-out R^2 on the landmarks ----------------------------------------------------
    # train scatter about the train mean; only the 978 landmark rows are needed
    top_tr = g_tr[:N_LANDMARK] - np.outer(mu_tr[:N_LANDMARK], s_tr) - np.outer(s_tr[:N_LANDMARK], mu_tr) \
        + n_tr * np.outer(mu_tr[:N_LANDMARK], mu_tr)
    sxx = top_tr[:, :N_LANDMARK]
    cond = float(np.linalg.cond(sxx))
    log(f"train landmark scatter condition number {cond:.3e}")
    beta = np.empty((N_LANDMARK, n_gene))
    beta[:, N_LANDMARK:] = np.linalg.solve(sxx, top_tr[:, N_LANDMARK:])          # inferred on all 978
    pinv_lm = np.linalg.inv(sxx)                                                  # leave-one-out for landmarks
    d = np.diag(pinv_lm).copy()
    beta[:, :N_LANDMARK] = -pinv_lm / d[None, :]
    np.fill_diagonal(beta[:, :N_LANDMARK], 0.0)
    del pinv_lm

    # in-sample (train) R^2, from the same accumulators
    diag_tr = np.diag(g_tr) - n_tr * mu_tr ** 2
    sse_tr = diag_tr - 2 * np.einsum("ij,ij->j", beta, top_tr) \
        + np.einsum("ij,ij->j", beta, sxx @ beta)
    r2_train = 1 - sse_tr / diag_tr

    # held-out: test scatter about the TRAIN mean, plus the test's own total sum of squares
    top_te = g_te[:N_LANDMARK] - np.outer(mu_tr[:N_LANDMARK], s_te) - np.outer(s_te[:N_LANDMARK], mu_tr) \
        + n_te * np.outer(mu_tr[:N_LANDMARK], mu_tr)
    axx = top_te[:, :N_LANDMARK]
    diag_te = np.diag(g_te) - 2 * mu_tr * s_te + n_te * mu_tr ** 2
    sse_te = diag_te - 2 * np.einsum("ij,ij->j", beta, top_te) + np.einsum("ij,ij->j", beta, axx @ beta)
    sst_te = np.diag(g_te) - s_te ** 2 / n_te
    r2_heldout = 1 - sse_te / sst_te
    del top_tr, top_te, axx, sxx, beta, g_tr, g_te
    log(f"R2 done: landmark median {np.median(r2_heldout[is_lm]):.4f}, "
        f"best-inferred median {np.median(r2_heldout[~is_lm]):.4f}")

    # --- (3) effective independent information ------------------------------------------------------------
    cov = (g_all - np.outer(s_all, mean_all)) / (n_sig - 1)   # = (G - n*mu mu') / (n-1), symmetric
    cov = (cov + cov.T) / 2
    del g_all
    inv_sd = 1.0 / sd_all
    corr = cov * inv_sd[:, None] * inv_sd[None, :]
    np.fill_diagonal(corr, 1.0)
    indep = {}
    for name, sel in (("bing", slice(None)), ("landmark", slice(0, N_LANDMARK)),
                      ("best_inferred", slice(N_LANDMARK, n_gene))):
        indep[f"corr_{name}"] = eff_dim(np.linalg.eigvalsh(corr[sel, sel]), n_sig)
        indep[f"cov_{name}"] = eff_dim(np.linalg.eigvalsh(cov[sel, sel]), n_sig)
        log(f"  eff dim {name}: corr PR {indep[f'corr_{name}']['participation_ratio']:.1f}, "
            f"cov PR {indep[f'cov_{name}']['participation_ratio']:.1f}")
    del cov
    ridged = corr.copy()
    ridged[np.diag_indices(n_gene)] += RIDGE
    uniq = 1.0 / np.diag(np.linalg.inv(ridged))
    del ridged, corr
    indep["corr_bing"]["sum_uniqueness"] = float(uniq.sum())
    indep["corr_bing"]["ridge"] = RIDGE
    log(f"uniqueness (1 - R2 on all other BING genes, ridge {RIDGE}): sum {uniq.sum():.1f}, "
        f"landmark median {np.median(uniq[is_lm]):.4f}, inferred median {np.median(uniq[~is_lm]):.4f}")

    # --- (1) replicate reproducibility --------------------------------------------------------------------
    rng = np.random.default_rng(SEED)
    sig["_row"] = np.arange(n_sig)
    grp = sig.groupby(GROUP_KEYS, dropna=False)["_row"].apply(np.array)
    grp = grp[grp.map(len) >= MIN_GROUP]
    capped = [np.sort(rng.choice(r, CAP_REPLICATES, replace=False)) if len(r) > CAP_REPLICATES else r for r in grp]
    n_groups = len(capped)
    rows = np.concatenate(capped)
    uniq_rows, inv = np.unique(rows, return_inverse=True)
    log(f"{n_groups} replicate groups (>= {MIN_GROUP}, capped at {CAP_REPLICATES}), "
        f"{len(rows)} signature slots, {len(uniq_rows)} distinct signatures")
    rep = np.asarray(mat[uniq_rows], dtype=np.float32)   # ~0.6 GB

    sum_g = np.zeros((n_groups, n_gene)); ss_tot = np.zeros(n_gene); n_g = np.zeros(n_groups)
    a_sum = np.zeros((n_groups, n_gene)); b_sum = np.zeros((n_groups, n_gene))
    pos = 0
    for k, r in enumerate(capped):
        idx = inv[pos:pos + len(r)]; pos += len(r)
        x = rep[idx].astype(np.float64)
        sum_g[k] = x.sum(0); ss_tot += (x ** 2).sum(0); n_g[k] = len(r)
        perm = rng.permutation(len(r)); h = len(r) // 2
        a_sum[k] = x[perm[:h]].mean(0); b_sum[k] = x[perm[h:]].mean(0)
    del rep
    N = n_g.sum()
    ssb_terms = sum_g ** 2 / n_g[:, None]
    ssw = ss_tot - ssb_terms.sum(0)
    ssb = ssb_terms.sum(0) - sum_g.sum(0) ** 2 / N
    msb = ssb / (n_groups - 1); msw = ssw / (N - n_groups)
    k0 = (N - (n_g ** 2).sum() / N) / (n_groups - 1)
    icc1 = (msb - msw) / (msb + (k0 - 1) * msw)
    sd_between = np.sqrt(np.clip((msb - msw) / k0, 0, None))
    sd_within = np.sqrt(msw)

    def pearson(a, b):
        a = a - a.mean(0); b = b - b.mean(0)
        return (a * b).sum(0) / np.sqrt((a ** 2).sum(0) * (b ** 2).sum(0))

    splithalf_r = pearson(a_sum, b_sum)
    del sum_g, a_sum, b_sum, ssb_terms
    log(f"ICC(1) landmark median {np.median(icc1[is_lm]):.4f}, inferred median {np.median(icc1[~is_lm]):.4f}; "
        f"split-half r landmark {np.median(splithalf_r[is_lm]):.4f}, inferred {np.median(splithalf_r[~is_lm]):.4f}")

    # --- (4) variance compression: quantiles on a seeded subsample ----------------------------------------
    qrows = np.sort(np.random.default_rng(SEED).choice(n_sig, min(QUANTILE_N, n_sig), replace=False))
    qsub = np.asarray(mat[qrows], dtype=np.float32)
    q = np.percentile(qsub, [1, 25, 50, 75, 99], axis=0)
    mean_abs = np.abs(qsub).mean(0)
    del qsub
    iqr = q[3] - q[1]
    log(f"sd landmark median {np.median(sd_all[is_lm]):.4f}, inferred {np.median(sd_all[~is_lm]):.4f}")

    # --- per-gene table -----------------------------------------------------------------------------------
    tab = genes[["gene_id", "gene_symbol", "ensembl_id", "feature_space"]].copy()
    tab["n_signatures"] = n_sig
    tab["sd_all"] = sd_all
    tab["iqr"] = iqr
    tab["p01"] = q[0]; tab["p50"] = q[2]; tab["p99"] = q[4]
    tab["range_p01_p99"] = q[4] - q[0]
    tab["mean_abs"] = mean_abs
    tab["min_all"] = gmin; tab["max_all"] = gmax
    tab["n_quantile_subsample"] = len(qrows)
    tab["icc1"] = icc1
    tab["sd_between"] = sd_between
    tab["sd_within"] = sd_within
    tab["sd_between_over_sd_all"] = sd_between / sd_all
    tab["splithalf_r"] = splithalf_r
    tab["n_replicate_groups"] = n_groups
    tab["n_replicate_signatures"] = int(N)
    tab["r2_landmark_heldout"] = r2_heldout
    tab["r2_landmark_train"] = r2_train
    tab["indep_var_frac_landmark"] = 1 - r2_heldout
    tab["r2_predictors"] = np.where(is_lm, "other_977_landmarks", "978_landmarks")
    tab["n_train_signatures"] = n_tr
    tab["n_test_signatures"] = n_te
    tab["uniqueness_all_bing"] = uniq
    bad = ~np.isfinite(tab[["icc1", "splithalf_r", "r2_landmark_heldout", "uniqueness_all_bing"]].values).all(1)
    tab["missing_reason"] = np.where(bad, "non_finite_statistic_zero_variance_or_singular", "")
    tab.to_csv(OUT / "bing_gene_quality.tsv", sep="\t", index=False, float_format="%.6g")
    log(f"wrote {OUT / 'bing_gene_quality.tsv'} ({len(tab)} rows, {int(bad.sum())} with a missing statistic)")

    # --- summary ------------------------------------------------------------------------------------------
    stats = ["icc1", "splithalf_r", "sd_between", "sd_within", "sd_between_over_sd_all", "r2_landmark_heldout",
             "r2_landmark_train", "indep_var_frac_landmark", "uniqueness_all_bing", "sd_all", "iqr",
             "range_p01_p99", "mean_abs"]
    rows_out = []
    for st in stats:
        for fs in ["landmark", "best inferred"]:
            v = tab.loc[tab.feature_space == fs, st].astype(float).values
            v = v[np.isfinite(v)]
            rows_out.append(dict(statistic=st, feature_space=fs, n=len(v), mean=v.mean(), sd=v.std(ddof=1),
                                 **{f"p{p:02d}": np.percentile(v, p) for p in (1, 5, 10, 25, 50, 75, 90, 95, 99)}))
    summ = pd.DataFrame(rows_out)
    # how many inferred genes fall below the landmark distribution on each reproducibility statistic
    below = []
    for st in ("icc1", "splithalf_r", "sd_between"):
        lm = tab.loc[is_lm, st].astype(float).values
        inf = tab.loc[~is_lm, st].astype(float).values
        for p in (1, 5, 10, 25, 50):
            thr = np.percentile(lm, p)
            below.append(dict(statistic=st, landmark_percentile=p, threshold=thr,
                              n_best_inferred_below=int((inf < thr).sum()),
                              frac_best_inferred_below=float((inf < thr).mean()),
                              n_landmark_below=int((lm < thr).sum())))
    pd.concat([summ, pd.DataFrame(below)], ignore_index=True).to_csv(
        OUT / "bing_gene_quality_summary.tsv", sep="\t", index=False, float_format="%.6g")
    log(f"wrote {OUT / 'bing_gene_quality_summary.tsv'}")

    run = dict(script="scripts/96_bing_gene_quality.py", date="2026-11-16", decision="D129",
               n_signatures=n_sig, n_genes=n_gene, n_landmark=N_LANDMARK, n_best_inferred=n_gene - N_LANDMARK,
               threads=THREADS, seed=SEED, train_frac=TRAIN_FRAC, split_salt=SPLIT_SALT, split_unit="pert_id",
               n_pert=len(perts), n_pert_train=len(train_pert), n_train_signatures=n_tr, n_test_signatures=n_te,
               group_keys=GROUP_KEYS, min_group=MIN_GROUP, cap_replicates=CAP_REPLICATES, n_groups=n_groups,
               n_replicate_signatures=int(N), k0=float(k0), quantile_n=int(len(qrows)), ridge=RIDGE,
               landmark_scatter_cond=cond, n_nan=n_nan, effective_independence=indep,
               minutes=round((time.time() - t0) / 60, 2))
    (OUT / "bing_gene_quality_run.json").write_text(json.dumps(run, indent=2, default=float))
    print(summ.set_index(["statistic", "feature_space"])[["n", "p10", "p25", "p50", "p75", "p90"]].to_string())
    print("\neffective independence:")
    for k, v in indep.items():
        print(f"  {k:22s} n {v['n_genes']:>6} PR {v['participation_ratio']:8.2f} "
              f"entropy {v['entropy_eff_rank']:8.2f} pr_null_iid {v['pr_null_iid']:8.1f} "
              f"pcs_to_80 {v['pcs_to_80']:>5} pcs_to_95 {v['pcs_to_95']:>5} min_eig {v['min_eigenvalue']:.3e}")
    log("done")


if __name__ == "__main__":
    main()
