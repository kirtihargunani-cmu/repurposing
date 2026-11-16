#!/usr/bin/env python
"""D128-A3.1: why does M change sign between the all-BING and landmark-only gene arms? (diagnostic only)

The D128-A3 landmark sensitivity amendment found that 25.0% of units reverse the sign of M between `all_bing` and
`landmark_only` (rho 0.677; worst landmark-fraction quartile 37.3% flips). That number has several possible causes
with very different consequences, and the arm decision is the user's. This script separates the causes. It changes
no parameter, picks no arm and re-ranks nothing.

  Control A  is flipping benign - confined to units whose M is near zero? Flip rate per |M| bin, and the |M| above
             which flipping effectively stops (or the statement that it does not stop).
  Control B  is it IMPUTATION or just a small gene count? For each unit, seeded random subsets of the all-BING
             footprint of the SAME SIZE as its landmark set are scored through the same formula, and their
             sign-disagreement rate with the full all-BING M is compared with the landmark arm's.
  H1         hub concentration: flip rate against the largest single regulator's share of the footprint.
  H2         inferred-gene redundancy: the effective number of independent genes in a footprint,
             n_eff = p^2 / sum_ij corr(g_i, g_j)^2 over the LINCS profiles, landmark and inferred separately.
  H3         near-boundary genes: the share of genes, and of the absolute signed numerator, carried by small |z|.

Outputs results/lincs/d128/a3/landmark_flip_diagnostic{,_units,_neff,_summary}.tsv and _run.json.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import d128io as IO  # noqa: E402
import score as SC  # noqa: E402

A3 = C.OUT / "a3"
T0 = time.time()
LOG = []
N_DRAW = 200              # random same-size subsets per (unit, profile) row, seeded
N_PROFILES_PER_UNIT = 15  # seeded sample of profiles per unit for the per-row controls
N_CORR_PROFILES = 2000    # seeded sample of LINCS profiles for the correlation structure (H2)
N_NEFF_UNITS = 60         # seeded sample of units for H2, plus the named cases
NAMED = ("fulvestrant", "vemurafenib", "dabrafenib", "palbociclib", "lapatinib")


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def bins_of(x, edges):
    return pd.cut(x, edges, include_lowest=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--a2", default="a2")
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    A3.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(C.SEED)

    # --- inputs ---------------------------------------------------------------------------------------------------
    idx = IO.profile_index()
    genes = IO.profile_genes()
    col = {g: i for i, g in enumerate(genes.gene_id)}
    land = set(genes.gene_id[genes.feature_space == "landmark"])
    P = IO.profiles(C.PROFILES_SPECIFIC, len(idx), len(genes))
    F = pd.read_parquet(C.OUT / "a1/target_gene_footprints.parquet")
    amb = pd.read_csv(C.OUT / "a1/path_ambiguity.tsv", sep="\t")
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str}).fillna({"action_flag": ""})
    U["genes_str"] = U.genes
    U["unit_id"] = U.cmap_name + "::" + U.tid.fillna("annot:" + U.genes_str)
    dup = U.unit_id.duplicated(keep=False)
    if dup.any():
        U.loc[dup, "unit_id"] = U.unit_id[dup] + "#" + U[dup].groupby("unit_id").cumcount().add(1).astype(str)
    U_by_id = {r.unit_id: r for r in U.itertuples(index=False)}
    hub = dict(zip(amb.unit_id, amb.hub_fraction))
    hubtf = dict(zip(amb.unit_id, amb.hub_tf))

    FP = {}
    for uid, d in F.groupby("unit_id", sort=False):
        res = d[d.status == "resolved"]
        g = res.gene_id.to_numpy()
        lm = np.array([x in land for x in g], dtype=bool)
        FP[uid] = dict(cols=np.array([col[x] for x in g], dtype=np.int64),
                       expected=res.expected_sign_pre_action.to_numpy().astype(float), lm=lm,
                       subset=res.subset.to_numpy().astype(str), n=len(g), n_lm=int(lm.sum()))

    S = pd.read_parquet(C.OUT / a.a2 / "scores_by_line.parquet")
    # M does not depend on the cohort, so one row per (unit, pert_id, line) is the grain where M is computed
    W = (S[["gene_arm", "unit_id", "pert_id", "cmap_name", "ModelID", "M", "M_status", "M_n_genes"]]
         .drop_duplicates(["gene_arm", "unit_id", "pert_id", "ModelID"]))
    piv = W.pivot_table(index=["unit_id", "pert_id", "cmap_name", "ModelID"], columns="gene_arm",
                        values="M", aggfunc="first").reset_index()
    piv = piv.rename(columns={"all_bing": "M_all", "landmark_only": "M_lm"})
    both = piv[piv.M_all.notna() & piv.M_lm.notna()].copy()
    both["flip"] = np.sign(both.M_all) != np.sign(both.M_lm)
    both["absM"] = both.M_all.abs()
    both["absM_max"] = np.maximum(both.M_all.abs(), both.M_lm.abs())
    both["hub_fraction"] = both.unit_id.map(hub)
    both["hub_tf"] = both.unit_id.map(hubtf)
    both["n_footprint"] = both.unit_id.map(lambda u: FP[u]["n"] if u in FP else np.nan)
    both["n_landmark"] = both.unit_id.map(lambda u: FP[u]["n_lm"] if u in FP else np.nan)
    both["landmark_fraction"] = both.n_landmark / both.n_footprint
    log(f"(unit, pert_id, line) rows scored in both arms: {len(both)}; overall sign-flip rate "
        f"{both.flip.mean():.3f} over {both.unit_id.nunique()} units")

    # --- Control A: flip rate against |M| ---------------------------------------------------------------------------
    edges = [0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 1.0]
    A_rows = []
    for name, c in (("|M| (all_bing)", "absM"), ("max(|M_all|, |M_lm|)", "absM_max")):
        g = both.groupby(bins_of(both[c], edges), observed=True)
        for k, d in g:
            A_rows.append(dict(control="A", measure=name, bin=str(k), n=len(d), flip_rate=float(d.flip.mean()),
                               median_absM=float(d[c].median()), n_units=int(d.unit_id.nunique())))
    CA = pd.DataFrame(A_rows)
    log("Control A - flip rate by |M| bin:\n" + CA.to_string())
    # the |M| above which flipping effectively stops: the smallest bin edge beyond which the flip rate stays < 5%
    stop = np.nan
    for e in edges[1:-1]:
        d = both[both.absM >= e]
        if len(d) >= 50 and d.flip.mean() < 0.05:
            stop = e
            break

    # --- the per-row sample used by Control B and H3 ----------------------------------------------------------------
    rows = []
    for uid, d in both.groupby("unit_id", sort=False):
        take = d if len(d) <= N_PROFILES_PER_UNIT else d.iloc[rng.choice(len(d), N_PROFILES_PER_UNIT, replace=False)]
        rows.append(take)
    SAMP = pd.concat(rows, ignore_index=True)
    prow = {(p_, m_): i for i, (p_, m_) in enumerate(zip(idx.pert_id, idx.ModelID))}
    log(f"Control B / H3 sample: {len(SAMP)} (unit, pert_id, line) rows over {SAMP.unit_id.nunique()} units, "
        f"{N_DRAW} seeded draws each")

    B_rows = []
    for r in SAMP.itertuples(index=False):
        fp, i_row = FP[r.unit_id], prow.get((r.pert_id, r.ModelID))
        if i_row is None:
            continue
        u = U_by_id[r.unit_id]
        action = 1 if u.orient == "raise" else -1
        z = np.asarray(P[i_row], dtype=float)[fp["cols"]]
        c = fp["expected"] * action * SC.bounded(z)            # per-gene signed contribution
        k = int(fp["n_lm"])
        p = len(c)
        if k < 1 or k >= p:
            continue
        draw = rng.integers(0, p, size=(N_DRAW, k))            # with replacement is biased; use permutations
        draw = np.argsort(rng.random((N_DRAW, p)), axis=1)[:, :k]
        sub = c[draw].mean(axis=1)
        m_all = float(c.mean())
        rand_flip = float((np.sign(sub) != np.sign(m_all)).mean())
        small = np.abs(z) < 1.0
        B_rows.append(dict(unit_id=r.unit_id, cmap_name=r.cmap_name, pert_id=r.pert_id, ModelID=r.ModelID,
                           n_footprint=p, n_landmark=k, landmark_fraction=k / p,
                           M_all=r.M_all, M_lm=r.M_lm, landmark_flip=bool(r.flip),
                           random_subset_flip_rate=rand_flip,
                           random_subset_sd=float(sub.std(ddof=1)),
                           random_subset_median=float(np.median(sub)),
                           hub_fraction=r.hub_fraction,
                           # H3: near-boundary genes
                           share_genes_small_z=float(small.mean()),
                           share_numerator_small_z=float(np.abs(c[small]).sum() / np.abs(c).sum())
                           if np.abs(c).sum() else np.nan,
                           median_abs_z=float(np.median(np.abs(z)))))
    CB = pd.DataFrame(B_rows)
    CB.to_csv(A3 / "landmark_flip_diagnostic.tsv", sep="\t", index=False)
    per_unit = CB.groupby(["unit_id", "cmap_name"]).agg(
        n_rows=("M_all", "size"), n_footprint=("n_footprint", "first"), n_landmark=("n_landmark", "first"),
        landmark_fraction=("landmark_fraction", "first"), hub_fraction=("hub_fraction", "first"),
        median_absM=("M_all", lambda x: float(np.median(np.abs(x)))),
        landmark_flip_rate=("landmark_flip", "mean"),
        random_subset_flip_rate=("random_subset_flip_rate", "mean"),
        share_numerator_small_z=("share_numerator_small_z", "mean")).reset_index()
    per_unit["gap_landmark_minus_random"] = per_unit.landmark_flip_rate - per_unit.random_subset_flip_rate
    per_unit.to_csv(A3 / "landmark_flip_diagnostic_units.tsv", sep="\t", index=False)
    log(f"Control B - landmark flip rate {CB.landmark_flip.mean():.3f} vs random same-size subsets "
        f"{CB.random_subset_flip_rate.mean():.3f} (row level, n = {len(CB)}); per unit the gap has median "
        f"{per_unit.gap_landmark_minus_random.median():+.3f}, IQR "
        f"{per_unit.gap_landmark_minus_random.quantile(0.25):+.3f} to "
        f"{per_unit.gap_landmark_minus_random.quantile(0.75):+.3f}")

    # --- H1: flip rate against hub concentration --------------------------------------------------------------------
    hb = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 1.0]
    H1 = both.dropna(subset=["hub_fraction"]).groupby(bins_of(both.hub_fraction.dropna(), hb), observed=True).agg(
        n=("flip", "size"), flip_rate=("flip", "mean"), n_units=("unit_id", "nunique"),
        median_absM=("absM", "median")).reset_index().rename(columns={"hub_fraction": "bin"})
    log("H1 - flip rate by largest-regulator share:\n" + H1.to_string())
    # H1 is confounded with |M| (hub-dominated units have larger |M|), so the same cut is repeated WITHIN |M|
    # strata. Added after seeing the marginal table, as a confound control, not as a tuning step: the hypothesis and
    # the measure were both fixed in the decision row.
    both["absM_stratum"] = bins_of(both.absM, [0, 0.01, 0.03, 1.0]).astype(str)
    both["hub_stratum"] = bins_of(both.hub_fraction, [0, 0.3, 0.5, 1.0]).astype(str)
    H1b = both.dropna(subset=["hub_fraction"]).groupby(["absM_stratum", "hub_stratum"], observed=True).agg(
        n=("flip", "size"), flip_rate=("flip", "mean"), n_units=("unit_id", "nunique")).reset_index()
    log("H1 - flip rate by hub share WITHIN |M| strata (confound control):\n" + H1b.to_string())

    # --- H2: effective number of independent genes ------------------------------------------------------------------
    sel = sorted(set(rng.choice(sorted(FP), min(N_NEFF_UNITS, len(FP)), replace=False).tolist())
                 | {u for u in FP if any(u.startswith(n + "::") for n in NAMED)})
    pr = rng.choice(len(idx), min(N_CORR_PROFILES, len(idx)), replace=False)
    pr.sort()
    X = np.asarray(P[pr], dtype=np.float32)
    neff_rows = []
    for uid in sel:
        fp = FP[uid]
        for space, mask in (("all", np.ones(fp["n"], bool)), ("landmark", fp["lm"]), ("inferred", ~fp["lm"])):
            cols = fp["cols"][mask]
            p = len(cols)
            if p < 5:
                continue
            Z = X[:, cols].astype(np.float64)
            Z = Z - Z.mean(0)
            sd = Z.std(0, ddof=1)
            keep = sd > 0
            Z = Z[:, keep] / sd[keep]
            p = Z.shape[1]
            Cm = (Z.T @ Z) / (len(Z) - 1)
            # n_eff = (sum eigenvalues)^2 / sum eigenvalues^2 = p^2 / ||C||_F^2 (no eigendecomposition needed)
            neff = float(p ** 2 / np.square(Cm).sum())
            off = (np.abs(Cm).sum() - p) / (p * (p - 1)) if p > 1 else np.nan
            neff_rows.append(dict(unit_id=uid, cmap_name=U_by_id[uid].cmap_name, space=space, n_genes=p,
                                  n_eff=neff, n_eff_per_gene=neff / p, mean_abs_offdiag_corr=float(off)))
    # reference points for H2: the whole measured gene space, and a random footprint-sized gene set
    for label, cols in (("whole BING space", np.arange(len(genes))),
                        ("random 478 genes", rng.choice(len(genes), 478, replace=False)),
                        ("all 978 landmarks", np.flatnonzero(genes.feature_space.to_numpy() == "landmark"))):
        Z = X[:, np.sort(cols)].astype(np.float64)
        Z = Z - Z.mean(0)
        sd = Z.std(0, ddof=1)
        Z = Z[:, sd > 0] / sd[sd > 0]
        pp = Z.shape[1]
        Cm = (Z.T @ Z) / (len(Z) - 1)
        neff_rows.append(dict(unit_id=f"REFERENCE: {label}", cmap_name="(reference)", space="reference",
                              n_genes=pp, n_eff=float(pp ** 2 / np.square(Cm).sum()),
                              n_eff_per_gene=float(pp / np.square(Cm).sum()),
                              mean_abs_offdiag_corr=float((np.abs(Cm).sum() - pp) / (pp * (pp - 1)))))
    NE = pd.DataFrame(neff_rows)
    NE.to_csv(A3 / "landmark_flip_diagnostic_neff.tsv", sep="\t", index=False)
    log("H2 - effective independent genes (median over the sampled units):\n"
        + NE.groupby("space").agg(units=("unit_id", "nunique"), median_n_genes=("n_genes", "median"),
                                  median_n_eff=("n_eff", "median"),
                                  median_n_eff_per_gene=("n_eff_per_gene", "median"),
                                  median_abs_corr=("mean_abs_offdiag_corr", "median")).to_string())

    # --- the named case: fulvestrant / ESR1 -------------------------------------------------------------------------
    esr1 = [u for u in FP if u.startswith("fulvestrant::") and "ESR1" in str(U_by_id[u].genes_str).split("|")]
    named_rows = []
    for u in esr1 + [x for x in FP if any(x.startswith(n + "::") for n in NAMED[1:])]:
        d = CB[CB.unit_id == u]
        if not len(d):
            continue
        ne = NE[NE.unit_id == u].set_index("space")
        named_rows.append(dict(unit_id=u, cmap_name=U_by_id[u].cmap_name, genes=U_by_id[u].genes_str,
                               n_footprint=int(d.n_footprint.iloc[0]), n_landmark=int(d.n_landmark.iloc[0]),
                               hub_fraction=d.hub_fraction.iloc[0], hub_tf=hubtf.get(u, ""),
                               median_absM=float(np.median(np.abs(d.M_all))),
                               landmark_flip_rate=float(d.landmark_flip.mean()),
                               random_subset_flip_rate=float(d.random_subset_flip_rate.mean()),
                               share_numerator_small_z=float(d.share_numerator_small_z.mean()),
                               n_eff_all=float(ne.n_eff.get("all", np.nan)) if len(ne) else np.nan,
                               n_eff_landmark=float(ne.n_eff.get("landmark", np.nan)) if len(ne) else np.nan,
                               n_eff_inferred=float(ne.n_eff.get("inferred", np.nan)) if len(ne) else np.nan))
    NAMED_T = pd.DataFrame(named_rows)
    log("named cases:\n" + (NAMED_T.to_string() if len(NAMED_T) else "none"))

    SUMM = pd.concat([
        CA,
        H1.assign(control="H1", measure="largest single regulator's share").rename(columns={"bin": "bin"}),
        pd.DataFrame([dict(control="B", measure="landmark vs random same-size subsets", bin="row level",
                           n=len(CB), flip_rate=float(CB.landmark_flip.mean()),
                           median_absM=float(CB.random_subset_flip_rate.mean()),
                           n_units=int(CB.unit_id.nunique()))]),
        pd.DataFrame([dict(control="H3", measure="share of the signed numerator carried by |z| < 1",
                           bin="median over rows", n=len(CB),
                           flip_rate=float(CB.share_numerator_small_z.median()),
                           median_absM=float(CB.share_genes_small_z.median()),
                           n_units=int(CB.unit_id.nunique()))])], ignore_index=True)
    SUMM.to_csv(A3 / "landmark_flip_diagnostic_summary.tsv", sep="\t", index=False)

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A3 / "landmark_flip_diagnostic_run.json").write_text(json.dumps(dict(
        stage="D128-A3.1", git_sha=sha, frozen=C.FROZEN, seed=C.SEED, n_draw=N_DRAW,
        n_profiles_per_unit=N_PROFILES_PER_UNIT, n_corr_profiles=int(len(pr)), n_neff_units=int(NE.unit_id.nunique()),
        rows_both_arms=len(both), units_both_arms=int(both.unit_id.nunique()),
        overall_flip_rate=float(both.flip.mean()),
        flip_rate_by_absM=CA.to_dict("records"),
        absM_above_which_flipping_stops=stop,
        control_b=dict(row_level_landmark_flip=float(CB.landmark_flip.mean()),
                       row_level_random_flip=float(CB.random_subset_flip_rate.mean()),
                       per_unit_gap_median=float(per_unit.gap_landmark_minus_random.median()),
                       per_unit_gap_q25=float(per_unit.gap_landmark_minus_random.quantile(0.25)),
                       per_unit_gap_q75=float(per_unit.gap_landmark_minus_random.quantile(0.75)),
                       units_landmark_worse=int((per_unit.gap_landmark_minus_random > 0.05).sum()),
                       units_landmark_better=int((per_unit.gap_landmark_minus_random < -0.05).sum()),
                       units=int(len(per_unit))),
        h1_hub=H1.to_dict("records"), h1_hub_within_absM=H1b.to_dict("records"),
        h2_reference=NE[NE.space == "reference"].to_dict("records"),
        h2_neff=NE.groupby("space").agg(units=("unit_id", "nunique"), median_n_genes=("n_genes", "median"),
                                        median_n_eff=("n_eff", "median"),
                                        median_n_eff_per_gene=("n_eff_per_gene", "median"),
                                        median_abs_corr=("mean_abs_offdiag_corr", "median")).reset_index()
                                       .to_dict("records"),
        h3=dict(median_share_genes_small_z=float(CB.share_genes_small_z.median()),
                median_share_numerator_small_z=float(CB.share_numerator_small_z.median())),
        named_cases=NAMED_T.to_dict("records"), log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
