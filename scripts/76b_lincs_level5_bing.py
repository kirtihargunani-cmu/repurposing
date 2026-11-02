#!/usr/bin/env python
"""LINCS 2020 level 5 compound signatures on landmark + best-inferred genes, for regulator scoring (runs in `celligner`).

Why (user, 2026-10-12): scoring drug signatures with the project ARACNe regulons (the same regulators the disease side
uses) needs enough of each regulator's targets; on the 978 landmark genes only 11-20% of regulators have >= 10 targets,
on landmark + best inferred ("BING", 10,174 genes) 89-92%.
Input: data/raw/lincs2020/level5_beta_trt_cp_n720216x12328.gctx (scripts/08_download_lincs_level5.sh).
Kept: high-quality signatures (is_hiq == 1) in LINCS lines that map to cancer DepMap models, the same lines as
scripts/76_lincs_level5_subset.py (HA1E and other non-cancerous lines excluded: HA1E is never among any tumor's 10
nearest lines in Celligner, user 2026-10-12); genes: feature_space landmark or best inferred.
Check: the landmark columns equal data/processed/lincs/level5_landmark.npy for the same signatures.
Outputs: data/processed/lincs/level5_bing_hiq.npy (float32, signatures x genes), level5_bing_hiq_sigs.tsv.gz,
         level5_bing_genes.tsv (gene order, feature_space); results/lincs/level5_bing_info.tsv
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
    sig = sig[(sig.pert_type == "trt_cp") & sig.cell_iname.isin(lines.index) & (sig.is_hiq == 1)].drop_duplicates("sig_id")
    genes = pd.read_csv(RAW / "geneinfo_beta.txt", sep="\t", dtype={"gene_id": str})
    bing = genes[genes.feature_space.isin(["landmark", "best inferred"])]
    bing = pd.concat([bing[bing.feature_space == "landmark"], bing[bing.feature_space != "landmark"]])  # landmark first

    with h5py.File(GCTX, "r") as f:
        mat = f["0/DATA/0/matrix"]
        col_ids = dec(f["0/META/COL/id"][:])
        row_ids = dec(f["0/META/ROW/id"][:])
        assert mat.shape == (len(col_ids), len(row_ids))
        gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)[bing.gene_id].values
        keep = np.isin(col_ids, sig.sig_id.values.astype(str))  # str, not object: object arrays take numpy's slow loop
        out = np.empty((int(keep.sum()), len(gpos)), dtype=np.float32)
        k = 0
        for i in range(0, len(col_ids), BLOCK):
            sel = np.where(keep[i:i + BLOCK])[0]
            if len(sel):
                out[k:k + len(sel)] = mat[i:i + BLOCK, :][sel][:, gpos]
                k += len(sel)
            if (i // BLOCK) % 5 == 0:
                print(f"  {i + BLOCK:>7} / {len(col_ids)} read, {k} kept, {time.time() - t0:.0f}s", flush=True)
        kept_ids = col_ids[keep]

    meta = sig.set_index("sig_id").loc[kept_ids].reset_index()
    meta.insert(meta.columns.get_loc("cell_iname") + 1, "ModelID", meta.cell_iname.map(lines).values)

    # consistency with the landmark extraction (scripts/76)
    lm_sigs = pd.read_csv(PROC / "level5_landmark_sigs.tsv.gz", sep="\t", usecols=["sig_id"]).sig_id
    lm_pos = pd.Series(np.arange(len(lm_sigs)), index=lm_sigs.values)
    lm = np.load(PROC / "level5_landmark.npy", mmap_mode="r")
    rng = np.random.default_rng(SEED)
    chk = rng.choice(len(meta), min(2000, len(meta)), replace=False)
    n_lm = int((bing.feature_space == "landmark").sum())
    max_diff = float(np.abs(out[chk, :n_lm] - lm[lm_pos[meta.sig_id.values[chk]].values]).max())

    PROC.mkdir(parents=True, exist_ok=True)
    np.save(PROC / "level5_bing_hiq.npy", out)
    meta.to_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t", index=False)
    bing[["gene_id", "gene_symbol", "ensembl_id", "feature_space"]].to_csv(PROC / "level5_bing_genes.tsv", sep="\t",
                                                                          index=False)
    info = {"n_signatures_hiq": len(meta), "n_active_tas04": int((meta.tas >= 0.4).sum()),
            "n_lines": meta.cell_iname.nunique(), "n_compounds": meta.pert_id.nunique(),
            "n_genes": len(gpos), "n_landmark": n_lm, "n_best_inferred": len(gpos) - n_lm,
            "n_nan": int(np.isnan(out).sum()), "landmark_max_abs_diff_vs_76": max_diff,
            "gb": round(out.nbytes / 1e9, 2), "minutes": round((time.time() - t0) / 60, 1)}
    pd.Series(info).to_csv(RES / "level5_bing_info.tsv", sep="\t", header=False)
    print(pd.Series(info).to_string())


if __name__ == "__main__":
    main()
