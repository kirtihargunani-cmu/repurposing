#!/usr/bin/env python
"""D132 step 1: build null profile matrices in the 79c layout, from REAL signatures only.

Nothing is shuffled at the gene level. Shuffling genes would reintroduce exactly the independence assumption D130
corrected; both nulls below are built from real data so the gene-gene correlation is carried by construction.

  replicate (PRIMARY)  the difference of two distinct level 5 signatures of the SAME (pert_id, cell_iname,
                       pert_dose, pert_time), divided by sqrt(2). The drug effect cancels exactly - same compound,
                       same dose, same time, same line - so what is left is real L1000 noise with real gene-gene
                       correlation and no regulator activity. This is the direct analogue of the disease side's
                       normal-split null (50_viper.R:131-136), which the D130 inventory cleared.
  rotation (SECONDARY) null_j = mu + (g_j . Pc)/sqrt(m-1), g ~ N(0, I) over the production `specific` profiles,
                       Pc gene-centred. Covariance = the empirical gene-gene covariance of the scored profiles
                       EXACTLY, so it carries the full across-drug structure including real mechanism-driven
                       covariance, and bounds the spread from above.

Both are written genes-contiguous-per-profile (the 79c .f32 layout) with a matching index file, so that the
PRODUCTION script 79d_drug_viper.R can be run on them unchanged via PROFILE_FILE / INDEX_FILE / RUN_TAG.

Outputs data/processed/lincs/d132/{replicate,rotation}_profiles.f32 + _index.tsv, and a build json.
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
OUT = PROC / "d132"
RES = ROOT / "results/lincs/d132"
SEED = 132
MAX_PER_GROUP = 4
T0 = time.time()
LOG = []


def log(m):
    line = f"[{time.time() - T0:6.1f}s] {m}"
    LOG.append(line)
    print(line, flush=True)


def gene_map():
    """level 5 BING row order -> the 79c profile gene order, matched by symbol (verified 1:1 and order-preserving)."""
    l5 = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
    p79 = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    m = l5.reset_index().merge(p79.reset_index(), left_on="gene_symbol", right_on="symbol",
                               suffixes=("_l5", "_79"))
    if len(m) != len(p79) or not (m["index_79"].to_numpy() == np.arange(len(p79))).all():
        raise SystemExit("gene mapping is not 1:1 onto the 79c gene order")
    return m["index_l5"].to_numpy()


def write_null(name, X, idx_names):
    """X is (n_profiles, n_genes) float32, C-order -> genes contiguous within each profile, the 79c layout."""
    f = OUT / f"{name}_profiles.f32"
    np.ascontiguousarray(X, dtype=np.float32).tofile(f)
    pd.DataFrame({"pert_id": idx_names, "ModelID": "NULL", "n_sig": 1}).to_csv(
        OUT / f"{name}_index.tsv", sep="\t", index=False)
    log(f"{name}: wrote {X.shape[0]} x {X.shape[1]} -> {f.name} ({f.stat().st_size / 1e6:.0f} MB)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-null", type=int, default=3000)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    OUT.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    cols = gene_map()
    n_genes = len(cols)
    log(f"gene map: {n_genes} level 5 columns in the 79c order")

    # --- primary: replicate differences --------------------------------------------------------------------------
    S = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    S["row"] = np.arange(len(S))
    key = ["pert_id", "cell_iname", "pert_dose", "pert_time"]
    groups = [g.row.to_numpy() for _, g in S.groupby(key, sort=True) if len(g) >= 2]
    log(f"replicate groups with >= 2 signatures: {len(groups)} "
        f"covering {sum(len(g) for g in groups)} signatures")
    pairs, used = [], np.zeros(len(groups), int)
    order = rng.permutation(len(groups))
    while len(pairs) < a.n_null:
        progressed = False
        for gi in order:
            if len(pairs) >= a.n_null:
                break
            if used[gi] >= MAX_PER_GROUP or len(groups[gi]) < 2:
                continue
            i, j = rng.choice(groups[gi], 2, replace=False)
            pairs.append((int(i), int(j)))
            used[gi] += 1
            progressed = True
        if not progressed:
            log(f"exhausted the pair pool at {len(pairs)} draws (cap {MAX_PER_GROUP} per group)")
            break
    pairs = np.array(pairs)
    log(f"replicate pairs drawn: {len(pairs)} from {int((used > 0).sum())} distinct groups")

    M = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")          # memory-mapped, never loaded whole
    log(f"level 5 matrix memory-mapped: {M.shape}")
    X = np.empty((len(pairs), n_genes), dtype=np.float32)
    CH = 250
    for s in range(0, len(pairs), CH):
        e = min(s + CH, len(pairs))
        A = np.asarray(M[pairs[s:e, 0]][:, cols], dtype=np.float32)
        B = np.asarray(M[pairs[s:e, 1]][:, cols], dtype=np.float32)
        X[s:e] = (A - B) / np.sqrt(2.0)
    log(f"replicate null built: per-gene sd median {np.median(X.std(0, ddof=1)):.3f}")
    write_null("replicate", X, [f"NULLREP{i:05d}" for i in range(len(X))])
    del X

    # --- secondary: rotation of the production `specific` profiles -------------------------------------------------
    pidx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    m = len(pidx)
    P = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype=np.float32, mode="r",
                  shape=(m, n_genes))
    mu = np.asarray(P, dtype=np.float64).mean(0)
    Pc = (np.asarray(P, dtype=np.float32) - mu.astype(np.float32))
    R = np.empty((a.n_null, n_genes), dtype=np.float32)
    for s in range(0, a.n_null, 500):
        e = min(s + 500, a.n_null)
        G = rng.standard_normal((e - s, m)).astype(np.float32)
        R[s:e] = (G @ Pc) / np.sqrt(m - 1)
    R += mu.astype(np.float32)
    log(f"rotation null built from {m} production profiles: per-gene sd median "
        f"{np.median(R.std(0, ddof=1)):.3f} (real profiles {np.median(Pc.std(0, ddof=1)):.3f})")
    write_null("rotation", R, [f"NULLROT{i:05d}" for i in range(len(R))])

    (RES / "d132_build.json").write_text(json.dumps(dict(
        stage="D132-null-profiles", seed=SEED, n_null=a.n_null, max_per_group=MAX_PER_GROUP,
        n_replicate_pairs=int(len(pairs)), n_groups_used=int((used > 0).sum()),
        n_genes=int(n_genes), n_production_profiles=int(m),
        git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                               text=True).stdout.strip(),
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
