#!/usr/bin/env python
"""LINCS 2020 level 5 compound signatures -> working subset for our cell lines (runs in the `celligner` env: h5py).

Input: data/raw/lincs2020/level5_beta_trt_cp_n720216x12328.gctx (scripts/08_download_lincs_level5.sh), MODZ
z-scores (replicate-consensus differential expression vs plate controls), 720,216 signatures x 12,328 genes.
Kept: signatures in LINCS lines that map to cancer DepMap models (results/lincs/cell_coverage.tsv), all of them with
the is_hiq flag carried along; genes: the 978 measured landmark genes only (the other 11,350 are inferred from them).
Read in blocks of signatures (the full matrix does not fit comfortably in one read).
Sanity check: correlation of high-quality signatures of the same compound across different lines vs random pairs.
Outputs: data/processed/lincs/level5_landmark.npy (float32, signatures x 978), level5_landmark_sigs.tsv.gz,
         level5_landmark_genes.tsv; results/lincs/level5_subset_info.tsv
"""
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/lincs2020"
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
GCTX = RAW / "level5_beta_trt_cp_n720216x12328.gctx"
BLOCK, SEED = 20000, 0


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def main():
    t0 = time.time()
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t")
    cov = cov.dropna(subset=["ModelID"])
    cov = cov[cov.OncotreePrimaryDisease != "Non-Cancerous"]
    lines = cov.set_index("cell_iname").ModelID
    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname", "pert_dose", "pert_dose_unit",
                               "pert_time", "is_hiq", "is_exemplar_sig", "tas", "nsample", "cc_q75"])
    sig = sig[(sig.pert_type == "trt_cp") & sig.cell_iname.isin(lines.index)].drop_duplicates("sig_id")
    genes = pd.read_csv(RAW / "geneinfo_beta.txt", sep="\t", dtype={"gene_id": str})
    lm = genes[genes.feature_space == "landmark"]

    with h5py.File(GCTX, "r") as f:
        mat = f["0/DATA/0/matrix"]
        col_ids = dec(f["0/META/COL/id"][:])   # signatures
        row_ids = dec(f["0/META/ROW/id"][:])   # genes (Entrez)
        print(f"matrix {mat.shape}, chunks {mat.chunks}, dtype {mat.dtype}; cols {len(col_ids)}, rows {len(row_ids)}")
        assert mat.shape == (len(col_ids), len(row_ids))
        gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)[lm.gene_id].values
        keep = np.isin(col_ids, sig.sig_id.values)
        out = np.empty((int(keep.sum()), len(gpos)), dtype=np.float32)
        k = 0
        for i in range(0, len(col_ids), BLOCK):
            sel = np.where(keep[i:i + BLOCK])[0]
            if len(sel):
                blk = mat[i:i + BLOCK, :]
                out[k:k + len(sel)] = blk[sel][:, gpos]
                k += len(sel)
            if (i // BLOCK) % 5 == 0:
                print(f"  {i + BLOCK:>7} / {len(col_ids)} read, {k} kept, {time.time() - t0:.0f}s", flush=True)
        kept_ids = col_ids[keep]

    PROC.mkdir(parents=True, exist_ok=True)
    np.save(PROC / "level5_landmark.npy", out)
    meta = sig.set_index("sig_id").loc[kept_ids].reset_index()
    meta.insert(meta.columns.get_loc("cell_iname") + 1, "ModelID", meta.cell_iname.map(lines).values)
    meta.to_csv(PROC / "level5_landmark_sigs.tsv.gz", sep="\t", index=False)
    lm[["gene_id", "gene_symbol", "ensembl_id"]].to_csv(PROC / "level5_landmark_genes.tsv", sep="\t", index=False)

    # sanity: same compound, different lines (high-quality) vs random pairs
    rng = np.random.default_rng(SEED)
    hq = meta.index[meta.is_hiq == 1].values
    by = meta.loc[hq].groupby("pert_id").apply(lambda d: d.drop_duplicates("cell_iname").index.values)
    pairs = [rng.choice(ix, 2, replace=False) for ix in by if len(ix) >= 2]
    pairs = [pairs[j] for j in rng.choice(len(pairs), min(20000, len(pairs)), replace=False)]
    same = np.array([np.corrcoef(out[a], out[b])[0, 1] for a, b in pairs])
    rnd = rng.choice(hq, (len(pairs), 2))
    rand = np.array([np.corrcoef(out[a], out[b])[0, 1] for a, b in rnd])
    info = {"n_signatures": len(kept_ids), "n_hiq": int((meta.is_hiq == 1).sum()), "n_lines": meta.cell_iname.nunique(),
            "n_compounds": meta.pert_id.nunique(), "n_genes_landmark": len(gpos),
            "value_median_abs": round(float(np.median(np.abs(out[rng.choice(len(out), 2000)]))), 3),
            "n_nan": int(np.isnan(out).sum()),
            "hiq_same_compound_diff_line_median_r": round(float(np.median(same)), 3),
            "hiq_random_pair_median_r": round(float(np.median(rand)), 3),
            "n_pairs": len(pairs), "minutes": round((time.time() - t0) / 60, 1)}
    pd.Series(info).to_csv(RES / "level5_subset_info.tsv", sep="\t", header=False)
    print(pd.Series(info).to_string())


if __name__ == "__main__":
    main()
