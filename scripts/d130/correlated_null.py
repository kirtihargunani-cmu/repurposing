#!/usr/bin/env python
"""D130: correlation-preserving nulls for the D128 mechanism score M.

Every permutation null in scripts/d128/nulls.py shuffles GENES as exchangeable units. D128-A3.1 measured that a
478-gene LINCS footprint carries n_eff ~ 46 independent measurements, so such a null can be too narrow and the
observed score can look more significant than it is. This script recomputes significance under a corrected
reference and reports old z beside new z. It changes NO score: M, R and every frozen parameter are imported from
the D128 engine untouched, and nothing is re-ranked or gated.

Method (fixed in the D130 decision row BEFORE running):

  CRP, the correlated-rank permutation. M = (1/n) sum_g e_g f(z_g) is a plain mean over footprint genes. A
  permutation null draws a UNIFORM random permutation of one of the two vectors. CRP replaces that uniform
  permutation with one drawn from the empirical gene-gene correlation, and changes nothing else:
      Z   = the m x n LINCS response matrix over the footprint's genes, column centred and scaled, m = 9,678
            profiles of the SAME profile arm, with the profile being scored excluded from the scatter
      u_k = Z^T g_k / sqrt(m - 2),  g_k ~ N(0, I_m) with the scored profile's coordinate set to 0
            => u_k ~ N(0, R_hat) EXACTLY, R_hat the leave-one-out gene-gene correlation (full-sample centring)
      pi_k = the rank order of u_k; the observed values are re-assigned to genes by that order
  so the null vector carries EXACTLY the observed multiset of values, arranged coherently with the measured
  covariance instead of at random. Applied to f(z) for the corrected `gene_label` null (within feature space,
  landmark and inferred separately, exactly as gene_label_permutation does) and to e for the corrected
  `sign_shuffle` null. With R_hat = I it reduces exactly to the uniform permutation.

Reported beside it, for every row:
  * the naive nulls recomputed at the same draw counts (so a change cannot be a draw-count artefact)
  * z_VIF = z_naive / sqrt(n_genes / n_eff), the cheap variance-inflation estimate
  * the correlated Gaussian null CGN (D130.1), which unlike any permutation is not invariant to a common shift
  * the analytic signed-statistic factor sqrt(e^T R e / n) and its ingredients (e_bar, 1^T R 1 / n), which is the
    correct inflation for a SIGN-weighted mean and is not in general equal to sqrt(n / n_eff)

Outputs results/lincs/d130/.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

D128 = Path(__file__).resolve().parent.parent / "d128"
sys.path.insert(0, str(D128))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a2_score as A2  # noqa: E402
import config as C  # noqa: E402
import nulls as NU  # noqa: E402
import score as SC  # noqa: E402

OUT = C.ROOT / "results/lincs/d130"
SEED = 130
T0 = time.time()
LOG = []


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line)
    print(line, flush=True)


# --------------------------------------------------------------------------------------------------- unit cache ----
class UnitCache:
    """Per (profile arm, unit, gene arm): the standardised profile matrix over the footprint, R-derived constants.

    R itself (n x n) is formed once, reduced to the three numbers that are reported, and dropped immediately, so
    only the m x n standardised matrix is retained. Total retained bytes are capped; on overflow the cache clears.
    """

    MAX_BYTES = 3_000_000_000

    def __init__(self, inp):
        self.inp = inp
        self.store = {}
        self.bytes = 0

    def get(self, unit_id, gene_arm, ev=None):
        key = (unit_id, gene_arm)
        e = self.store.get(key)
        if e is None:
            fp = self.inp["FP"][unit_id]
            cols = fp["lm_cols"] if gene_arm == "landmark_only" else fp["cols"]
            lm = fp["lm_mask"] if gene_arm == "landmark_only" else fp["lm"]
            n = len(cols)
            X = np.asarray(self.inp["M"][:, cols], dtype=np.float64)     # m x n, from the memmap
            mu = X.mean(0)
            sd = X.std(0, ddof=1)
            sd[sd == 0] = 1.0
            Zc = ((X - mu) / sd).astype(np.float32)   # float32 keeps the matmul cheap; R is exact to ~1e-6
            del X
            e = dict(cols=cols, lm=np.asarray(lm, bool), n=n, Zc=Zc, mu=mu, sd=sd)
            self._add_R(e, ev)
            if self.bytes + Zc.nbytes > self.MAX_BYTES:
                self.store.clear()
                self.bytes = 0
            self.store[key] = e
            self.bytes += Zc.nbytes
        elif ev is not None and "eRe_over_n" not in e:
            self._add_R(e, ev)
        return e

    def _add_R(self, e, ev):
        Zc = e["Zc"].astype(np.float64)
        m, n = Zc.shape
        R = (Zc.T @ Zc) / (m - 1)
        d = np.sqrt(np.diag(R))
        d[d == 0] = 1.0
        R /= np.outer(d, d)
        e["n_eff"] = float(n ** 2 / np.square(R).sum())
        e["one_R_one_over_n"] = float(R.sum() / n)                   # = n / n_eff_1, the all-+1 inflation
        if ev is not None:
            v = np.asarray(ev, dtype=float)
            e["eRe_over_n"] = float(v @ R @ v / len(v))
        del R, Zc

    def drop_R(self):
        self.store.clear()
        self.bytes = 0


# ------------------------------------------------------------------------------------------------------- nulls ----
def crp_ranks(Zc, i_excl, n_draws, rng, chunk=200):
    """n_draws x n matrix of correlated Gaussian draws u ~ N(0, R_hat), leave-one-out."""
    m, n = Zc.shape
    out = np.empty((n_draws, n), dtype=np.float32)
    k = 0
    while k < n_draws:
        b = min(chunk, n_draws - k)
        G = rng.standard_normal((b, m)).astype(np.float32)
        if i_excl is not None:
            G[:, i_excl] = 0.0
        out[k:k + b] = (G @ Zc) / np.sqrt(m - 2)
        k += b
    return out


def cgn_draws(e_signed, z_obs, ent, U):
    """D130.1: the correlated Gaussian null (unconditional).

    CRP is a permutation, so it is invariant to any shift common to every gene; the dominant eigenvector of the
    LINCS gene-gene correlation is close to such a common component, so a permutation null - naive OR corrected -
    conditions it away and cannot see it. M itself is NOT invariant to it: a common shift c moves M by e_bar * c.
    This null therefore draws the whole response instead of rearranging it:
        z* = mu_pop + sd_pop * u,   u ~ N(0, R_hat)   (so Cov(z*) = the empirical gene-gene covariance)
    rescaled so that the draw's spread across the footprint equals the observed profile's, which keeps the realised
    response magnitude fixed while letting the pattern - including the common component's share of it - be random.
    """
    mu, sd = ent["mu"], ent["sd"]
    s_obs = float(np.std(z_obs, ddof=1))
    vals = np.empty(len(U))
    for k in range(len(U)):
        zs = mu + sd * U[k]
        s = float(np.std(zs, ddof=1))
        if s > 0:
            zs = zs * (s_obs / s)
        vals[k] = float((e_signed * SC.bounded(zs)).mean())
    return vals


def _assign_by_rank(values, u):
    """Re-assign `values` to positions by the rank order of u (both 1-D, same length)."""
    order = np.argsort(u, kind="stable")
    out = np.empty_like(values)
    out[order] = np.sort(values)
    return out


def crp_gene_label(e_signed, fz, lm, U):
    """Corrected `gene_label` null: f(z) rank-mapped within feature space by each correlated draw."""
    L, I = np.flatnonzero(lm), np.flatnonzero(~lm)
    vals = np.empty(len(U))
    for k in range(len(U)):
        fz_star = fz.copy()
        for idx in (L, I):
            if len(idx) > 1:
                fz_star[idx] = _assign_by_rank(fz[idx], U[k][idx])
        vals[k] = float((e_signed * fz_star).mean())
    return vals


def crp_sign_shuffle(e_signed, fz, U):
    """Corrected `sign_shuffle` null: the expected signs rank-mapped by each correlated draw."""
    vals = np.empty(len(U))
    for k in range(len(U)):
        vals[k] = float((_assign_by_rank(e_signed, U[k]) * fz).mean())
    return vals


def naive_gene_label(e_signed, fz, lm, n_draws, rng):
    L, I = np.flatnonzero(lm), np.flatnonzero(~lm)
    vals = np.empty(n_draws)
    for k in range(n_draws):
        fz_star = fz.copy()
        for idx in (L, I):
            if len(idx) > 1:
                fz_star[idx] = fz[rng.permutation(idx)]
        vals[k] = float((e_signed * fz_star).mean())
    return vals


def naive_sign_shuffle(e_signed, fz, n_draws, rng):
    vals = np.empty(n_draws)
    for k in range(n_draws):
        vals[k] = float((rng.permutation(e_signed) * fz).mean())
    return vals


def summarise(v, observed, tag):
    v = np.asarray(v, dtype=float)
    sd = float(v.std(ddof=1))
    return {f"{tag}_null_mean": float(v.mean()), f"{tag}_null_sd": sd,
            f"{tag}_pct": float((v < observed).mean() * 100),
            f"{tag}_z": float((observed - v.mean()) / sd) if sd > 0 else np.nan}


# -------------------------------------------------------------------------------------------------------- rows ----
def score_row(inp, cache, unit_id, pert_id, model_id, gene_arm, n_small, n_big, rng):
    """Recompute M exactly as D128 does, then the four nulls (naive / corrected) x (sign / gene-label)."""
    i_row = inp["prow"].get((pert_id, model_id))
    if i_row is None:
        return None
    fp = inp["FP"].get(unit_id)
    if fp is None:
        return None
    zrow = np.asarray(inp["M"][i_row], dtype=float)
    if gene_arm == "landmark_only":
        gids, z, lm, expected, subset = fp["lm_genes"], zrow[fp["lm_cols"]], fp["lm_mask"], fp["lm_expected"], fp["lm_subset"]
    else:
        gids, z, lm, expected, subset = fp["genes"], zrow[fp["cols"]], fp["lm"], fp["expected"], fp["subset"]
    if not len(gids):
        return None
    u = inp["U_by_id"][unit_id]
    action = 1 if u.orient == "raise" else -1
    m = SC.score_mechanism(gids, expected, action, z, landmark_mask=lm, n_unresolved=fp["n_unresolved"],
                           action_flag=u.action_flag, subset=subset,
                           n_signaling_tfs=inp["n_sig_tfs"].get(unit_id, 0))
    if not m.scored:
        return None
    e_signed = expected * action
    if not NU.exact_shortcut_ok(e_signed, z):
        # the exact fast path does not apply (a zero expected sign or a missing response); skipped and recorded
        return dict(unit_id=unit_id, pert_id=pert_id, ModelID=model_id, gene_arm=gene_arm, M=m.value,
                    skipped="not_exact_shortcut")
    fz = SC.bounded(z)
    ent = cache.get(unit_id, gene_arm, ev=e_signed)
    n = len(e_signed)
    U = crp_ranks(ent["Zc"], i_row, n_big, rng)

    row = dict(unit_id=unit_id, cmap_name=u.cmap_name, unit_genes=u.genes_str, pert_id=pert_id, ModelID=model_id,
               gene_arm=gene_arm, profile_arm=inp["arm"], M=float(m.value), n_genes=n,
               n_landmark=int(lm.sum()), n_inferred=int((~lm).sum()),
               n_eff=ent["n_eff"], n_over_n_eff=n / ent["n_eff"],
               e_bar=float(e_signed.mean()), one_R_one_over_n=ent["one_R_one_over_n"],
               eRe_over_n=ent["eRe_over_n"])
    row["vif_factor_sqrt_n_over_neff"] = float(np.sqrt(row["n_over_n_eff"]))
    row["analytic_factor_sqrt_eRe_over_n"] = float(np.sqrt(max(row["eRe_over_n"], 0.0)))

    for tag, nd in (("200", n_small), ("1000", n_big)):
        row.update(summarise(naive_sign_shuffle(e_signed, fz, nd, np.random.default_rng(SEED + 1)),
                             m.value, f"naive_sign_{tag}"))
        row.update(summarise(naive_gene_label(e_signed, fz, lm, nd, np.random.default_rng(SEED + 2)),
                             m.value, f"naive_gene_{tag}"))
        row.update(summarise(crp_sign_shuffle(e_signed, fz, U[:nd]), m.value, f"crp_sign_{tag}"))
        row.update(summarise(crp_gene_label(e_signed, fz, lm, U[:nd]), m.value, f"crp_gene_{tag}"))
        row.update(summarise(cgn_draws(e_signed, z, ent, U[:nd]), m.value, f"cgn_{tag}"))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--nulls-small", type=int, default=200)
    ap.add_argument("--nulls-big", type=int, default=1000)
    ap.add_argument("--population-rows", type=int, default=1000)
    ap.add_argument("--only", default="", help="comma list of blocks to run: esr1,a4,population")
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen")
    OUT.mkdir(parents=True, exist_ok=True)
    blocks = set(a.only.split(",")) if a.only else {"esr1", "a4", "population"}
    A4 = C.OUT / "a4"

    inputs, caches = {}, {}

    def get(arm):
        if arm not in inputs:
            log(f"loading D128 inputs, profile arm {arm}")
            inputs[arm] = A2.load_inputs(arm)
            caches[arm] = UnitCache(inputs[arm])
        return inputs[arm], caches[arm]

    results = {}

    # --- block 1: the ESR1 / fulvestrant control ---------------------------------------------------------------
    if "esr1" in blocks:
        E = pd.read_csv(A4 / "esr1_control.tsv", sep="\t")
        E = E[E.M_status == "scored"]
        rows = []
        for r in E.itertuples(index=False):
            inp, cache = get(r.arm)
            out = score_row(inp, cache, r.unit_id, r.pert_id, "ACH-000019", "all_bing",
                            a.nulls_small, a.nulls_big, np.random.default_rng(SEED))
            if out is None:
                log(f"ESR1 {r.cmap_name} {r.pert_id} {r.arm}: not reproducible, skipped")
                continue
            out.update(set=r.set, old_M=r.M, old_sign_shuffle_z=r.sign_shuffle_z, old_gene_label_z=r.gene_label_z,
                       old_sign_shuffle_pct=r.sign_shuffle_pct, old_gene_label_pct=r.gene_label_pct)
            rows.append(out)
            log(f"ESR1 {r.arm} {r.cmap_name} {r.pert_id}: M={out['M']:.5f} old_sign_z={r.sign_shuffle_z} "
                f"crp_sign_z={out['crp_sign_1000_z']:.2f} old_gene_z={r.gene_label_z} "
                f"crp_gene_z={out['crp_gene_1000_z']:.2f}")
        D = pd.DataFrame(rows)
        D.to_csv(OUT / "esr1_control_corrected.tsv", sep="\t", index=False)
        results["esr1_rows"] = len(D)
        for c in caches.values():
            c.drop_R()

    # --- block 2: the A4 actions 2-5 controls ------------------------------------------------------------------
    if "a4" in blocks:
        C4 = pd.read_csv(A4 / "control_results.tsv", sep="\t")
        C4 = C4[(C4.M_status == "scored") & C4.sign_shuffle_z.notna()]
        rows = []
        inp, cache = get("specific")
        for r in C4.itertuples(index=False):
            out = score_row(inp, cache, r.unit_id, r.pert_id, r.ModelID, r.gene_arm,
                            a.nulls_small, a.nulls_big, np.random.default_rng(SEED))
            if out is None:
                continue
            out.update(wave=r.wave, context=r.context, anchor=r.anchor, line=r.line, lineage=r.lineage,
                       old_M=r.M, old_sign_shuffle_z=r.sign_shuffle_z, old_gene_label_z=r.gene_label_z,
                       old_sign_shuffle_pct=r.sign_shuffle_pct, old_gene_label_pct=r.gene_label_pct,
                       old_confident_positive=r.confident_positive)
            rows.append(out)
            if len(rows) % 20 == 0:
                log(f"A4 controls: {len(rows)} rows")
        D = pd.DataFrame(rows)
        D.to_csv(OUT / "a4_controls_corrected.tsv", sep="\t", index=False)
        results["a4_rows"] = len(D)
        cache.drop_R()

    # --- block 3: the population sample (the A3 population nulls were never run) --------------------------------
    if "population" in blocks:
        S = pd.read_parquet(C.OUT / "a2/scores_by_line.parquet")
        S = S[(S.gene_arm == "all_bing") & (S.M_status == "scored")]
        S = S.drop_duplicates(["unit_id", "pert_id", "ModelID"]).reset_index(drop=True)
        log(f"population: {len(S)} distinct (unit, pert_id, line) M-scored rows on the frozen arm")
        rng0 = np.random.default_rng(SEED)
        take = rng0.choice(len(S), min(a.population_rows, len(S)), replace=False)
        S = S.iloc[np.sort(take)]
        inp, cache = get("specific")
        rows = []
        for j, r in enumerate(S.itertuples(index=False)):
            out = score_row(inp, cache, r.unit_id, r.pert_id, r.ModelID, "all_bing",
                            a.nulls_small, a.nulls_big, np.random.default_rng(SEED))
            if out is None:
                continue
            out["old_M"] = r.M
            rows.append(out)
            if (j + 1) % 100 == 0:
                log(f"population: {j + 1}/{len(S)} rows, {len(rows)} scored, "
                    f"{len(cache.store)} unit matrices cached")
        D = pd.DataFrame(rows)
        D.to_csv(OUT / "population_sample_corrected.tsv", sep="\t", index=False)
        results["population_rows"] = len(D)

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (OUT / "d130_run.json").write_text(json.dumps(dict(
        stage="D130-correlated-null", git_sha=sha, seed=SEED, nulls_small=a.nulls_small, nulls_big=a.nulls_big,
        blocks=sorted(blocks), threads=a.threads, results=results, log=LOG,
        seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
