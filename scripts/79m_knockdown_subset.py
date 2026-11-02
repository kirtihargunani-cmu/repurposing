#!/usr/bin/env python
"""LINCS 2020 TF knockdown signatures for the regulon benchmark (scripts/79n), on landmark + best-inferred genes (`celligner`).

Why (user 2026-10-28): the tumor-network ESR1 regulon reads anti-estrogens in MCF7 as ESR1 activation while the CollecTRI
regulon reads them as inhibition. A knockdown lowers the TF, so a regulon set that is right for that TF on perturbation
data should score it as inhibited; the benchmark picks, per TF, which regulon set the drug side uses.
Signatures (cancer DepMap lines as scripts/76b; knocked-down gene = a regulator of CollecTRI, DoRothEA ABC or the v2 TF
ARACNe set):
  shRNA   trt_sh.cgs  consensus gene signatures (several hairpins per gene combined, which damps seed off-target effects)
  CRISPR  trt_xpr     is_hiq == 1
Inputs: data/raw/lincs2020/level5_beta_trt_sh_n238351x12328.gctx, level5_beta_trt_xpr_n142901x12328.gctx
Outputs: data/processed/lincs/knockdown_bing.npy / .f32 (float32, signatures x genes, gene order = level5_bing_genes.tsv),
         knockdown_bing_sigs.tsv.gz; results/lincs/knockdown_info.tsv
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
FILES = {"shRNA": ("trt_sh.cgs", RAW / "level5_beta_trt_sh_n238351x12328.gctx"),
         "CRISPR": ("trt_xpr", RAW / "level5_beta_trt_xpr_n142901x12328.gctx")}
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
                               "nsample", "cc_q75"])
    genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
    blocks, metas = [], []
    for tech, (ptype, gctx) in FILES.items():
        s = sig[(sig.pert_type == ptype) & sig.cell_iname.isin(lines.index) & sig.cmap_name.isin(tfs)]
        if ptype == "trt_xpr":
            s = s[s.is_hiq == 1]
        s = s.drop_duplicates("sig_id")
        with h5py.File(gctx, "r") as f:
            mat = f["0/DATA/0/matrix"]
            col_ids = dec(f["0/META/COL/id"][:]); row_ids = dec(f["0/META/ROW/id"][:])
            gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)[genes.gene_id].values
            keep = np.isin(col_ids, s.sig_id.values.astype(str))
            out = np.empty((int(keep.sum()), len(gpos)), dtype=np.float32); k = 0
            for i in range(0, len(col_ids), BLOCK):
                sel = np.where(keep[i:i + BLOCK])[0]
                if len(sel):
                    out[k:k + len(sel)] = mat[i:i + BLOCK, :][sel][:, gpos]; k += len(sel)
            kept = col_ids[keep]
        m = s.set_index("sig_id").loc[kept].reset_index()
        m.insert(0, "technology", tech); m["ModelID"] = m.cell_iname.map(lines).values
        blocks.append(out); metas.append(m)
        print(f"{tech}: {len(m)} signatures, {m.cmap_name.nunique()} TFs, {m.cell_iname.nunique()} lines "
              f"({time.time() - t0:.0f}s)", flush=True)
    X = np.vstack(blocks); meta = pd.concat(metas, ignore_index=True)
    np.save(PROC / "knockdown_bing.npy", X)
    X.tofile(PROC / "knockdown_bing.f32")   # profile-major, for R readBin (as 79c / 79d)
    meta.to_csv(PROC / "knockdown_bing_sigs.tsv.gz", sep="\t", index=False)
    info = meta.groupby("technology").agg(signatures=("sig_id", "size"), tfs=("cmap_name", "nunique"),
                                          lines=("cell_iname", "nunique"))
    info["n_nan"] = [int(np.isnan(b).sum()) for b in blocks]
    info.to_csv(RES / "knockdown_info.tsv", sep="\t")
    print(info.to_string())


if __name__ == "__main__":
    main()
