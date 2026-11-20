#!/usr/bin/env python
"""O-3 / D136-S0: crossing-rate and censoring curves as a function of k.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python  (needs h5py)
SPEC: docs/SPEC_dose_response_bmd.md section 3.

The spec pre-registers a PROCEDURE, not a threshold. S0's only job is to characterise the
instrument: across a continuous range of k, what fraction of (drug, line, gene) dose curves ever
cross BMR_g = k * sd_within_g, and what fraction are right-censored at the top tested dose.

NO TARGET-SPECIFIC QUANTITY IS COMPUTED HERE. No target annotation, no regulon, no expected
direction, no null. That is what makes choosing k from this curve a calibration of measurement
range rather than fitting to an outcome (spec section 3).

Reads Level 5 directly - NOT 79c, which collapses each (drug, line) to one consensus profile at
TAS >= 0.4 and discards the dose structure this analysis is about.

Output: results/systems_validation/observability/o3_s0_crossing_curves.{tsv,json}
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd, h5py

ROOT = Path("/home/kirti/projects/repurposing")
RAW = ROOT / "data/raw/lincs2020"
OUT = ROOT / "results/systems_validation/observability"
GCTX = RAW / "level5_beta_trt_cp_n720216x12328.gctx"
SEED, N_CELLS, BLOCK = 20261118, 1500, 20000
KS = np.round(np.arange(0.25, 4.01, 0.25), 2)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "6"


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    q = pd.read_csv(ROOT / "results/lincs/gene_quality/bing_gene_quality.tsv", sep="\t",
                    dtype={"gene_id": str})
    sdw = q.set_index("gene_id").sd_within
    print(f"sd_within for {len(sdw):,} genes "
          f"(landmark median {q[q.feature_space=='landmark'].sd_within.median():.3f}, "
          f"best-inferred {q[q.feature_space!='landmark'].sd_within.median():.3f})")

    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname",
                               "pert_dose", "pert_time", "is_hiq", "tas"])
    sig = sig[(sig.pert_type == "trt_cp") & (sig.is_hiq == 1) & (sig.tas >= 0.4)]
    sig = sig[sig.pert_dose > 0].drop_duplicates("sig_id")
    # one signature per (cell, dose): the most common time, then highest tas - dose structure only
    sig = sig.sort_values("tas", ascending=False).drop_duplicates(
        ["pert_id", "cell_iname", "pert_dose"])
    nd = sig.groupby(["pert_id", "cell_iname"]).pert_dose.nunique()
    cells = nd[nd >= 3].index
    print(f"(compound, line) cells with >= 3 distinct doses at is_hiq & TAS>=0.4: {len(cells):,}")
    rng = np.random.default_rng(SEED)
    pick = cells[rng.choice(len(cells), min(N_CELLS, len(cells)), replace=False)] \
        if len(cells) > N_CELLS else cells
    sel = sig.set_index(["pert_id", "cell_iname"]).loc[pick].reset_index()
    print(f"sampled {len(pick):,} cells (seed {SEED}) -> {len(sel):,} signatures to read")

    with h5py.File(GCTX, "r") as f:
        mat = f["0/DATA/0/matrix"]
        col_ids = dec(f["0/META/COL/id"][:])
        row_ids = dec(f["0/META/ROW/id"][:])
        gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)
        keep_genes = [g for g in sdw.index if g in gpos.index]
        gidx = gpos[keep_genes].values
        want = set(sel.sig_id)
        keep = np.isin(col_ids, np.array(sorted(want), dtype=str))
        X = np.empty((int(keep.sum()), len(gidx)), dtype=np.float32)
        k = 0
        for i in range(0, len(col_ids), BLOCK):
            s = np.where(keep[i:i + BLOCK])[0]
            if len(s):
                X[k:k + len(s)] = mat[i:i + BLOCK, :][s][:, gidx]
                k += len(s)
            if (i // BLOCK) % 10 == 0:
                print(f"   {i+BLOCK:>7}/{len(col_ids)} scanned, {k:,} kept, {time.time()-t0:.0f}s",
                      flush=True)
        got = col_ids[keep]
    print(f"read {X.shape[0]:,} x {X.shape[1]:,} in {time.time()-t0:.0f}s")

    meta = sel.set_index("sig_id").loc[got].reset_index()
    sd = sdw[keep_genes].values.astype(np.float32)
    A = np.abs(X)
    # per (cell, gene): the MAXIMUM |z| over the dose series -> does the curve ever cross k*sd_g
    cell_key = meta.pert_id + "|" + meta.cell_iname
    order = np.argsort(cell_key.values, kind="stable")
    A, keys = A[order], cell_key.values[order]
    bounds = np.flatnonzero(np.r_[True, keys[1:] != keys[:-1], True])
    n_cells = len(bounds) - 1
    peak = np.empty((n_cells, A.shape[1]), dtype=np.float32)
    for i in range(n_cells):
        peak[i] = A[bounds[i]:bounds[i + 1]].max(axis=0)
    print(f"collapsed to {n_cells:,} (compound, line) cells x {peak.shape[1]:,} genes")

    rows = []
    for kk in KS:
        thr = kk * sd
        crossed = peak >= thr
        per_cell = crossed.mean(axis=1)
        rows.append(dict(k=float(kk),
                         pooled_crossing_rate=round(float(crossed.mean()), 5),
                         pooled_censored_rate=round(float(1 - crossed.mean()), 5),
                         median_cell_crossing_rate=round(float(np.median(per_cell)), 5),
                         p10_cell=round(float(np.percentile(per_cell, 10)), 5),
                         p90_cell=round(float(np.percentile(per_cell, 90)), 5),
                         cells_with_under_1pct_crossing=int((per_cell < 0.01).sum()),
                         cells_with_over_50pct_crossing=int((per_cell > 0.50).sum())))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "o3_s0_crossing_curves.tsv", sep="\t", index=False)
    (OUT / "o3_s0_crossing_curves.json").write_text(json.dumps(dict(
        step="O-3 / D136-S0 crossing and censoring curves",
        spec="docs/SPEC_dose_response_bmd.md section 3",
        seed=SEED, n_cells_available=int(len(cells)), n_cells_sampled=int(n_cells),
        n_signatures=int(X.shape[0]), n_genes=int(peak.shape[1]),
        source="level5_beta_trt_cp (NOT 79c consensus)",
        no_target_quantity_computed=True, curves=rows,
        runtime_s=round(time.time() - t0, 1)), indent=1))
    print("\nk    pooled_cross  median_cell  p10    p90    dead_cells(<1%)  saturated(>50%)")
    for r in rows:
        print(f"{r['k']:<4} {r['pooled_crossing_rate']:>11.4f} {r['median_cell_crossing_rate']:>12.4f} "
              f"{r['p10_cell']:>6.3f} {r['p90_cell']:>6.3f} {r['cells_with_under_1pct_crossing']:>14d} "
              f"{r['cells_with_over_50pct_crossing']:>15d}")
    print("\nwrote", OUT / "o3_s0_crossing_curves.tsv")


if __name__ == "__main__":
    sys.exit(main())
