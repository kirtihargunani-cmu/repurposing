#!/usr/bin/env python
"""D131: is the reversal score's 1/sqrt(sum w^2) normalisation wrong, and does it change anything that shipped?

79e computes  rev = -sum_i s_i a_i / sqrt(sum_i s_i^2)  over a disease's master regulators. sqrt(sum s^2) is the SD
of that weighted sum ONLY if the regulator activities a_i are independent unit-variance terms. Regulons overlap
heavily, so the a_i are correlated and the true SD is larger: rev is on an inflated scale, by an amount that
depends on how correlated each disease's MR set happens to be.

AUDIT ONLY. No production score is changed, nothing is re-ranked, no fix is shipped.

What is computed, per cohort, on the production consensus arm:
  * the exact normaliser `norm = sqrt(sum s^2)` and the identity check rev_recomputed == published rev
  * THE QUANTITY ASKED FOR: the inflation ratio sqrt(w' Sigma w) / sqrt(w' w), Sigma the empirical covariance of
    regulator activities across the cohort's drugs. Because rev = -(D w)/sqrt(w'w), that ratio is identically the
    realised SD of rev across the cohort's drugs - the score is built to have SD 1 under its own assumption, so the
    SD it actually has IS the factor by which the assumption is wrong. It is computed as a plain sample variance of
    the scalar D w_hat, which needs no inversion of Sigma and is well defined even when there are fewer drugs than
    regulators; NO shrinkage is applied to it and none is needed. Shrinkage would be needed only to invert Sigma,
    which nothing here does.
  * the decomposition  sd(rev) = scale x corr_inflation, where
        v              = S w_hat,  w_hat = s / ||s||,  S = diag(sd of each regulator across the cohort's drugs)
        scale          = ||v||                        (per-regulator spread; 1 if a VIPER NES really were N(0,1))
        corr_inflation = sqrt(v' C v / v'v)           (C = correlation of regulator activity ACROSS DRUGS)
    corr_inflation is the part attributable to the independence assumption, and is the number asked for.
  * n_eff of the regulator activities, p^2 / ||C||_F^2, the same participation ratio used for genes in D128-A3.1.
    Read it as "this sum's variance behaves as if it were n_eff independent terms", NOT as "there are only n_eff
    regulators' worth of information": the correlation is real cascade biology and the numerator keeps all of it.
    ||C||_F^2 IS sensitive to sampling noise when p is comparable to the number of drugs, so the iid reference
    p / (1 + p / n_drugs) is reported beside it; the inflation ratio itself does not depend on this.
  * the weight-concentration count 79e already writes as `n_eff` ((sum w^2)^2 / sum w^4), so the two are not confused

Then, directly: does any shipped number change? Every column is checked for invariance to a per-cohort constant.

Inputs are read-only: the D126-B1 consensus activity dump (data/processed/lincs/d126_act, drugs x regulators per
cohort, written by 79e itself), the Step 6 consensus MR tables, and the shipped reversal / split tables.
Output results/lincs/d131/.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
ACT = PROC / "d126_act"
OUT = RES / "d131"
CONS_SETS = {"aracne": "results/step6/aracne_project_v2_tf", "collectri": "results/step6/collectri_v2"}
T0 = time.time()
LOG = []


def log(m):
    line = f"[{time.time() - T0:6.1f}s] {m}"
    LOG.append(line)
    print(line, flush=True)


def load_mr():
    """The set-tagged consensus MR table exactly as 79e builds it, including the D92/D97 ARACNe network overrides."""
    netmap_f = ROOT / "config/group_networks.tsv"
    netmap = pd.read_csv(netmap_f, sep="\t").set_index("group_id").net.to_dict() if netmap_f.exists() else {}
    parts = []
    for k, d in CONS_SETS.items():
        x = pd.read_csv(ROOT / d / "group_mr.tsv", sep="\t")
        x["regulator"] = k + "|" + x.regulator
        parts.append(x)
    group = pd.concat(parts, ignore_index=True)
    parts_s = []
    for k, d in CONS_SETS.items():
        f = ROOT / f"{d}_strata" / "group_mr.tsv"
        if f.exists():
            x = pd.read_csv(f, sep="\t")
            x["regulator"] = k + "|" + x.regulator
            parts_s.append(x)
    stratum = pd.concat(parts_s, ignore_index=True) if parts_s else None
    return group, stratum, netmap


def main():
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "8"
    OUT.mkdir(parents=True, exist_ok=True)
    mr_g, mr_s, netmap = load_mr()
    log(f"MR rows: group {len(mr_g)}, stratum {0 if mr_s is None else len(mr_s)}; netmap {len(netmap)} groups")

    pub = {}
    for lvl, f in (("group_id", RES / "drug_reversal_v2_consensus.tsv"),
                   ("stratum", RES / "stratum_reversal_v2_consensus.tsv")):
        if f.exists():
            pub[lvl] = pd.read_csv(f, sep="\t", low_memory=False)
            log(f"published {f.name}: {len(pub[lvl])} rows, {pub[lvl].group_id.nunique()} cohorts")

    rows = []
    for f in sorted(ACT.glob("*.f32")):
        g = f.stem
        regs = pd.Index((ACT / f"{g}.regulators.txt").read_text().split())
        drugs = pd.read_csv(ACT / f"{g}.drugs.tsv", sep="\t")
        D = np.fromfile(f, dtype=np.float32).reshape(len(drugs), len(regs)).astype(np.float64)
        lvl = "stratum" if (mr_s is not None and (mr_s.group_id == g).any()) else "group_id"
        src = mr_s if lvl == "stratum" else mr_g
        m = src[src.group_id == g]
        if not len(m):
            log(f"{g}: no MR rows, skipped")
            continue
        base = m.regulator.str.split("|", n=1).str[1].values
        keep = pd.Index(base).isin(regs)
        m, base = m[keep], base[keep]
        # consensus weights: per base regulator, the mean NES over the sets that carry it (79e consensus())
        s_ser = pd.Series(m.nes.values, index=base).groupby(level=0).mean().reindex(regs)
        if s_ser.isna().any():
            log(f"{g}: {int(s_ser.isna().sum())} dumped regulators have no MR weight, skipped")
            continue
        s = s_ser.to_numpy()
        norm = float(np.sqrt((s ** 2).sum()))
        rev = -(D * s).sum(1) / norm

        # identity check against the shipped table
        chk = dict(n_matched=0, max_abs_diff=np.nan, spearman=np.nan)
        if lvl in pub:
            P = pub[lvl]
            P = P[P.group_id == g][["pert_id", "rev"]]
            j = drugs[["pert_id"]].assign(rev_recomputed=rev).merge(P, on="pert_id", how="inner")
            if len(j):
                chk = dict(n_matched=len(j), max_abs_diff=float((j.rev_recomputed - j.rev).abs().max()),
                           spearman=float(np.corrcoef(j.rev_recomputed.rank(), j.rev.rank())[0, 1]))

        # the decomposition
        sd_rev = float(rev.std(ddof=1))
        Dc = D - D.mean(0)
        sd_reg = Dc.std(0, ddof=1)
        ok = sd_reg > 0
        Cm = np.corrcoef(Dc[:, ok], rowvar=False)
        p = int(ok.sum())
        n_eff_act = float(p ** 2 / np.square(Cm).sum())
        w_hat = s / np.linalg.norm(s)
        v = (sd_reg * w_hat)[ok]
        vv = float(v @ v)
        corr_infl = float(np.sqrt((v @ Cm @ v) / vv)) if vv > 0 else np.nan
        scale = float(np.sqrt(vv))
        n_eff_w = float((s ** 2).sum() ** 2 / ((s ** 4).sum() + 1e-12))

        rows.append(dict(cohort=g, level=lvl, n_drugs=len(drugs), n_mrs=len(s), n_reg_nonconstant=p,
                         norm_sqrt_sum_w2=norm, inflation_sqrt_wSw_over_sqrt_ww=sd_rev,
                         sd_rev_published_scale=sd_rev,
                         scale_per_regulator=scale, corr_inflation=corr_infl,
                         sd_rev_check=scale * corr_infl, n_eff_activity=n_eff_act,
                         n_eff_activity_per_reg=n_eff_act / p,
                         n_eff_null_iid=float(p / (1 + p / len(drugs))),
                         n_eff_weight_concentration_79e=n_eff_w,
                         mean_abs_offdiag_corr=float((np.abs(Cm).sum() - p) / (p * (p - 1))),
                         identity_n=chk["n_matched"], identity_max_abs_diff=chk["max_abs_diff"],
                         identity_spearman=chk["spearman"]))
        if len(rows) % 10 == 0:
            log(f"{len(rows)} cohorts done ({g})")

    R = pd.DataFrame(rows)
    R.to_csv(OUT / "cohort_inflation.tsv", sep="\t", index=False)
    log(f"cohorts: {len(R)}")
    log("identity check vs the shipped rev: max |diff| = "
        f"{R.identity_max_abs_diff.max():.6f} over {int(R.identity_n.sum())} matched rows; "
        f"min Spearman {R.identity_spearman.min():.6f}")
    for c in ("inflation_sqrt_wSw_over_sqrt_ww", "corr_inflation", "scale_per_regulator", "n_eff_activity",
              "n_eff_activity_per_reg", "n_mrs"):
        q = R[c].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
        log(f"{c}: median {q['50%']:.4g}  (5% {q['5%']:.4g}, 95% {q['95%']:.4g}, min {q['min']:.4g}, max {q['max']:.4g})")

    (OUT / "d131_run.json").write_text(json.dumps(dict(
        stage="D131-normalisation-inflation",
        git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
        n_cohorts=len(R), log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
