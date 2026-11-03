#!/usr/bin/env python
"""LINCS 2020 TF overexpression signatures for the D84 perturbation benchmark (scripts/79n PERTURB=d84).

Why (user 2026-11-02, D84): ARACNe and CollecTRI read drugs differently (`79t`) but tie on knockdowns; overexpression
tests the opposite direction and was used to build neither set.
Runs in `celligner` (h5py). Signatures: trt_oe, wild-type ORF projects OEB / OEC only (one construct per gene; TA / MUT are allele-variant screens,
OFL mixes construct types and lncRNAs, HSF / HSFA are HEK293T), cancer DepMap lines (as scripts/76b / 79m),
is_hiq == 1, overexpressed gene = a regulator of CollecTRI, DoRothEA ABC or the v2 TF ARACNe set
(results/lincs/knockdown_regulators.txt, as 79m).
Input: data/raw/lincs2020/level5_beta_trt_oe_n34171x12328.gctx (scripts/08b)
Outputs: data/processed/lincs/overexpression_bing.f32 (float32, signatures x genes, gene order = level5_bing_genes.tsv),
         overexpression_bing_sigs.tsv.gz; results/lincs/overexpression_info.tsv
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
GCTX = RAW / "level5_beta_trt_oe_n34171x12328.gctx"
PROJECTS = ["OEB", "OEC"]
BLOCK = 20000


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def main():
    t0 = time.time()
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    cov = cov[cov.OncotreePrimaryDisease != "Non-Cancerous"]
    lines = cov.set_index("cell_iname").ModelID
    tfs = set(pd.read_csv(RES / "knockdown_regulators.txt", header=None)[0])
    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname", "pert_time", "is_hiq",
                               "nsample", "cc_q75", "project_code"])
    s = sig[(sig.pert_type == "trt_oe") & sig.project_code.isin(PROJECTS) & sig.cell_iname.isin(lines.index)
            & sig.cmap_name.isin(tfs) & (sig.is_hiq == 1)].drop_duplicates("sig_id")
    genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
    with h5py.File(GCTX, "r") as f:
        mat = f["0/DATA/0/matrix"]
        col_ids = dec(f["0/META/COL/id"][:]); row_ids = dec(f["0/META/ROW/id"][:])
        gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)[genes.gene_id].values
        keep = np.isin(col_ids, s.sig_id.values.astype(str))
        X = np.empty((int(keep.sum()), len(gpos)), dtype=np.float32); k = 0
        for i in range(0, len(col_ids), BLOCK):
            sel = np.where(keep[i:i + BLOCK])[0]
            if len(sel):
                X[k:k + len(sel)] = mat[i:i + BLOCK, :][sel][:, gpos]; k += len(sel)
        kept = col_ids[keep]
    assert k == len(kept) == len(s), (k, len(kept), len(s))
    meta = s.set_index("sig_id").loc[kept].reset_index()
    meta.insert(0, "technology", "ORF"); meta["ModelID"] = meta.cell_iname.map(lines).values
    X.tofile(PROC / "overexpression_bing.f32")   # profile-major, for R readBin (as 79m)
    meta.to_csv(PROC / "overexpression_bing_sigs.tsv.gz", sep="\t", index=False)
    info = meta.groupby("project_code").agg(signatures=("sig_id", "size"), tfs=("cmap_name", "nunique"),
                                            lines=("cell_iname", "nunique"))
    info["n_nan"] = int(np.isnan(X).sum())
    info.to_csv(RES / "overexpression_info.tsv", sep="\t")
    print(info.to_string(), f"\n({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
