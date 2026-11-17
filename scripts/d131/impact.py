#!/usr/bin/env python
"""D131 part 2: does the wrong normaliser change anything that shipped?

Part 1 measured the inflation. This asks the only question that matters: which shipped numbers move when the
normaliser is corrected from sqrt(w'w) to sqrt(w' Sigma w).

  (1) THE EXACT INFLATION, with no reconstruction at all. rev = -(D w)/sqrt(w'w), so the SD of the SHIPPED rev
      across a cohort's drugs IS sqrt(w' Sigma w)/sqrt(w'w). Read straight off the published table.
  (2) INVARIANCE. The correction multiplies a cohort's whole rev vector by one positive constant. Every column
      built from rev by a within-cohort RANK (rev_pct, spec_pct, rank_rev, every tier) or by a within-cohort
      ROBUST Z (rev_z, and spec = rev_z - mean(reference rev_z)) is therefore unchanged exactly. Verified
      numerically by rescaling each cohort's rev by a random positive constant and recomputing.
  (3) prolif_share, the one place a part's magnitude is compared with another part's. 79i multiplies rev_part back
      by sqrt(sum w_part^2), so `a` and `b` are RAW summed reversals and prolif_share is ALSO invariant to the
      normaliser - but the comparison it makes is between two raw sums over gene sets of different size and
      different internal correlation. Under a covariance-aware comparison each part is divided by its own
      sqrt(w' Sigma w). This computes prolif_share both ways and counts how many rows cross 0.5 and how many
      change the shipped `prolif_driven` flag.

AUDIT ONLY: nothing is re-ranked, no production file is written, no fix is applied.
Output results/lincs/d131/.
"""
import json
import os
import subprocess
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
TOP_PCT = 10
T0 = time.time()
LOG = []


def log(m):
    line = f"[{time.time() - T0:6.1f}s] {m}"
    LOG.append(line)
    print(line, flush=True)


def net_of(cohort):
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def main():
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "8"
    OUT.mkdir(parents=True, exist_ok=True)
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    nmf = ROOT / "config/group_networks.tsv"
    netmap = pd.read_csv(nmf, sep="\t").set_index("group_id").net.to_dict() if nmf.exists() else {}
    PR = pd.read_csv(RES / "prolif_regulators_p2_v2_consensus.tsv", sep="\t")

    split = {"group_id": pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False),
             "stratum": pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)}
    # the subgroup table keys the cohort as `stratum`; rename so both levels share one key name
    split["stratum"] = split["stratum"].rename(columns={"stratum": "group_id"})
    par = split["stratum"].drop_duplicates("group_id").set_index("group_id").parent_group.to_dict()

    # ---- (1) exact inflation + (2) invariance, from the shipped tables only -------------------------------------
    rows, inv_fail = [], 0
    rng = np.random.default_rng(131)
    for lvl, S in split.items():
        for g, d in S.groupby("group_id"):
            rev = d.rev.to_numpy(float)
            if len(rev) < 3 or not np.isfinite(rev).all():
                continue
            c = float(rng.uniform(0.05, 20))                      # an arbitrary positive rescaling
            pct = pd.Series(rev).rank(method="max").values / len(rev) * 100
            pct_c = pd.Series(rev * c).rank(method="max").values / len(rev) * 100
            mad = 1.4826 * np.median(np.abs(rev - np.median(rev)))
            z = (rev - np.median(rev)) / mad
            revc = rev * c
            madc = 1.4826 * np.median(np.abs(revc - np.median(revc)))
            zc = (revc - np.median(revc)) / madc
            tier = (rev > 0) & (pct >= 100 - TOP_PCT)
            tierc = (revc > 0) & (pct_c >= 100 - TOP_PCT)
            ok = bool(np.array_equal(pct, pct_c) and np.allclose(z, zc, atol=1e-9) and np.array_equal(tier, tierc))
            inv_fail += (not ok)
            rows.append(dict(cohort=g, level=lvl, n_drugs=len(rev),
                             inflation_exact_sd_rev=float(rev.std(ddof=1)),
                             rank_invariant=ok, scale_tested=round(c, 4)))
    INV = pd.DataFrame(rows)
    INV.to_csv(OUT / "inflation_exact_and_invariance.tsv", sep="\t", index=False)
    q = INV.inflation_exact_sd_rev.describe(percentiles=[0.05, 0.5, 0.95])
    log(f"(1) EXACT inflation sqrt(w'Sigma w)/sqrt(w'w) = sd(shipped rev), {len(INV)} cohorts: "
        f"median {q['50%']:.2f}, 5-95% {q['5%']:.2f}-{q['95%']:.2f}, min {q['min']:.2f}, max {q['max']:.2f}")
    log(f"(2) within-cohort rank / robust-z invariance under an arbitrary positive rescaling: "
        f"{len(INV) - inv_fail}/{len(INV)} cohorts invariant in rev_pct, rev_z AND the tier mask "
        f"({inv_fail} failures)")

    # ---- (3) prolif_share, current vs covariance-aware ------------------------------------------------------------
    mr_parts = []
    for k, d in CONS_SETS.items():
        for suf in ("", "_strata"):
            f = ROOT / f"{d}{suf}" / "group_mr.tsv"
            if f.exists():
                x = pd.read_csv(f, sep="\t")
                x["set"] = k
                x["base"] = x.regulator
                mr_parts.append(x)
    MR = pd.concat(mr_parts, ignore_index=True)

    prows = []
    for f in sorted(ACT.glob("*.f32")):
        g = f.stem
        lvl = "stratum" if g in par else "group_id"
        parent = par.get(g, g)
        if parent not in grp.index:
            continue
        net = netmap.get(parent, net_of(grp.loc[parent, "cohort"]))
        regs = pd.Index((ACT / f"{g}.regulators.txt").read_text().split())
        drugs = pd.read_csv(ACT / f"{g}.drugs.tsv", sep="\t")
        D = np.fromfile(f, dtype=np.float32).reshape(len(drugs), len(regs)).astype(np.float64)
        m = MR[MR.group_id == g]
        if not len(m):
            continue
        s = m.groupby("base").nes.mean().reindex(regs)
        if s.isna().any():
            continue
        s = s.to_numpy()
        linked = set(PR.regulator[(PR.prolif_linked) & (((PR.set == "aracne") & (PR.net == net)) |
                                                        (PR.set == "collectri"))])
        isp = np.array([r in linked for r in regs])
        if isp.sum() < 5 or (~isp).sum() < 5:
            continue
        Dc = D - D.mean(0)
        out = {}
        for tag, mask in (("p", isp), ("n", ~isp)):
            w = s[mask]
            raw = np.abs(D[:, mask] @ w)                       # the raw summed reversal: normaliser-free
            sd_cov = float((Dc[:, mask] @ w).std(ddof=1))      # sqrt(w' Sigma w)
            out[tag] = (raw, float(np.sqrt((w ** 2).sum())), sd_cov, int(mask.sum()))
        a, b = out["p"][0], out["n"][0]
        share_cur = a / (a + b + 1e-12)
        # The covariance-aware share is applied to the SHIPPED prolif_share, not to the reconstruction, so the
        # reconstruction's small weight error cannot leak into the impact count. With r = b/a,
        # share_cur = 1/(1+r) and share_cov = 1/(1 + r*k), k = sd_cov_prolif / sd_cov_nonprolif - one constant per
        # cohort, a ratio of two SDs over ~1,000 regulators and insensitive to a small error in the weights.
        k = out["p"][2] / out["n"][2]
        S = split[lvl]
        S = S[S.group_id == g]
        j = drugs[["pert_id"]].assign(share_cur=share_cur).merge(
            S[["pert_id", "prolif_share", "prolif_driven", "rev_pct", "rev_mech_pct", "n_mrs_prolif"]],
            on="pert_id", how="inner")
        if not len(j):
            continue
        ps = j.prolif_share.to_numpy(float).clip(1e-9, 1 - 1e-9)
        r = (1 - ps) / ps
        j["share_cov"] = 1.0 / (1.0 + r * k)
        drv_cov = (j.share_cov > 0.5) & (j.rev_pct >= 100 - TOP_PCT) & (j.rev_mech_pct < 100 - TOP_PCT)
        prows.append(dict(cohort=g, level=lvl, n_drugs=len(j), n_prolif_regs=out["p"][3],
                          n_nonprolif_regs=out["n"][3],
                          published_n_mrs_prolif=float(j.n_mrs_prolif.iloc[0]),
                          recon_vs_published_share_maxdiff=float((j.share_cur - j.prolif_share).abs().max()),
                          recon_vs_published_share_meddiff=float((j.share_cur - j.prolif_share).abs().median()),
                          k_sd_ratio_prolif_over_nonprolif=float(k),
                          median_share_published=float(j.prolif_share.median()),
                          median_share_cov=float(j.share_cov.median()),
                          n_cross_half=int(((j.prolif_share > 0.5) != (j.share_cov > 0.5)).sum()),
                          n_prolif_driven_published=int(j.prolif_driven.fillna(False).astype(bool).sum()),
                          n_prolif_driven_cov=int(drv_cov.sum()),
                          n_prolif_driven_changed=int((j.prolif_driven.fillna(False).astype(bool) != drv_cov).sum()),
                          sd_cov_prolif=out["p"][2], sd_cov_nonprolif=out["n"][2],
                          norm_prolif=out["p"][1], norm_nonprolif=out["n"][1]))
    P = pd.DataFrame(prows)
    P.to_csv(OUT / "prolif_share_impact.tsv", sep="\t", index=False)
    if len(P):
        log(f"(3) prolif_share on {len(P)} cohorts, {int(P.n_drugs.sum())} cohort x drug rows. "
            f"Reconstruction check vs the shipped prolif_share: max |diff| "
            f"{P.recon_vs_published_share_maxdiff.max():.4f}, median per cohort "
            f"{P.recon_vs_published_share_meddiff.median():.4f}; regulator-count check "
            f"{int((P.n_prolif_regs == P.published_n_mrs_prolif).sum())}/{len(P)} cohorts match exactly")
        log(f"    k = sd_cov(prolif)/sd_cov(nonprolif): median {P.k_sd_ratio_prolif_over_nonprolif.median():.3f}, "
            f"range {P.k_sd_ratio_prolif_over_nonprolif.min():.3f}-{P.k_sd_ratio_prolif_over_nonprolif.max():.3f}")
        log(f"    median prolif_share published {P.median_share_published.median():.3f} -> "
            f"covariance-aware {P.median_share_cov.median():.3f}")
        log(f"    rows crossing the 0.5 line: {int(P.n_cross_half.sum())} of {int(P.n_drugs.sum())} "
            f"({100 * P.n_cross_half.sum() / P.n_drugs.sum():.1f}%)")
        log(f"    prolif_driven: {int(P.n_prolif_driven_published.sum())} published -> "
            f"{int(P.n_prolif_driven_cov.sum())} covariance-aware; "
            f"{int(P.n_prolif_driven_changed.sum())} rows change the flag")

    (OUT / "d131_impact_run.json").write_text(json.dumps(dict(
        stage="D131-impact",
        git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip(),
        n_cohorts_invariance=len(INV), n_cohorts_prolif=len(P), log=LOG,
        seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
