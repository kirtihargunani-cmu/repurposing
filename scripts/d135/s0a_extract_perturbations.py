#!/usr/bin/env python
"""D135 stage S0, part A -- extract genetic-perturbation signatures of our target genes and measure them
PER SIGNATURE. Runs in the `celligner` env (h5py lives only there, as scripts/76b and 79m).

NO CONCORDANCE IS COMPUTED HERE. No drug profile is read. Nothing is correlated with anything. Every number is
a property of one perturbation signature on its own: how many genes it moved, and where the perturbed gene's
own mRNA sits. Definitions were fixed in docs/DECISIONS.md row D135-S0, committed and pushed before this ran.

Inputs
  results/lincs/target_role/d121/units.tsv            377 target genes (rows with non-empty tid)
  results/lincs/cell_coverage.tsv                     cell_iname -> ModelID (cancer lines only)
  data/processed/lincs/drug_line_profiles_index.tsv   the 73 drug-matrix lines
  data/processed/lincs/drug_line_profiles_genes.tsv   the 10,112-gene comparison space (ENSG)
  data/processed/lincs/level5_bing_genes.tsv          ENSG <-> entrez <-> symbol, feature_space
  data/raw/lincs2020/siginfo_beta.txt
  data/raw/lincs2020/level5_beta_trt_{xpr,sh,oe}_*.gctx   memory-mapped, read in blocks, never loaded whole

Output
  results/lincs/d135/perturbation_sigs.tsv.gz   one row per extracted signature: strength + target-mRNA QC
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "8"

import json
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RAW = ROOT / "data/raw/lincs2020"
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUT = RES / "d135"
OUT.mkdir(parents=True, exist_ok=True)

ARMS = {"trt_xpr": RAW / "level5_beta_trt_xpr_n142901x12328.gctx",
        "trt_sh.cgs": RAW / "level5_beta_trt_sh_n238351x12328.gctx",
        "trt_sh": RAW / "level5_beta_trt_sh_n238351x12328.gctx",
        "trt_oe": RAW / "level5_beta_trt_oe_n34171x12328.gctx"}
BLOCK = 10000


def dec(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a])


def main():
    t0 = time.time()
    prov = {}

    u = pd.read_csv(RES / "target_role/d121/units.tsv", sep="\t")
    u = u[u.tid.notna() & (u.tid.astype(str).str.strip() != "")]
    genes = sorted({g for s in u.genes.dropna() for g in str(s).split("|")})
    prov["target_genes"] = len(genes)

    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    cov = cov[cov.OncotreePrimaryDisease != "Non-Cancerous"]
    lmap = cov.drop_duplicates("cell_iname").set_index("cell_iname").ModelID

    didx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    drug_lines = set(didx.ModelID)
    prov["drug_matrix_lines"] = len(drug_lines)

    # gene space: the 10,112-gene drug-profile space, mapped to gctx entrez rows
    dg = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t", dtype={"gene_id": str})
    bg = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})
    gm = dg.merge(bg[["gene_id", "gene_symbol", "ensembl_id"]].rename(columns={"gene_id": "entrez"}),
                  left_on="gene_id", right_on="ensembl_id", how="left")
    unmapped = gm[gm.entrez.isna()]
    gm = gm[gm.entrez.notna()].reset_index(drop=True)
    prov["gene_space"] = dict(drug_space_genes=int(len(dg)), mapped_to_entrez=int(len(gm)),
                              unmapped=sorted(unmapped.symbol.astype(str)),
                              note="unmapped genes are dropped from the strength measures and reported, "
                                   "never imputed (D135-S0(4))")
    sym2pos = {s: i for i, s in enumerate(gm.symbol)}          # position within the extracted gene space
    fspace = dict(zip(gm.symbol, gm.feature_space))
    prov["target_genes_on_gene_space"] = int(sum(g in sym2pos for g in genes))
    prov["target_genes_landmark"] = int(sum(fspace.get(g) == "landmark" for g in genes))

    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname", "pert_time",
                               "pert_dose", "is_hiq", "qc_pass", "nsample", "cc_q75", "tas", "ss_ngene"])
    # one pass per FILE, not per arm: trt_sh and trt_sh.cgs share an 11.8 GB matrix and it is read once
    byfile = {}
    for a, p in ARMS.items():
        byfile.setdefault(p, []).append(a)
    out_frames = []
    for path, arms in byfile.items():
        s = sig[sig.pert_type.isin(arms) & sig.cmap_name.isin(genes)].copy()
        s["ModelID"] = s.cell_iname.map(lmap)
        s = s[s.ModelID.notna() & s.ModelID.isin(drug_lines)].drop_duplicates("sig_id")
        arm = "+".join(arms)
        if not len(s):
            continue
        with h5py.File(path, "r") as f:
            mat = f["0/DATA/0/matrix"]
            col_ids = dec(f["0/META/COL/id"][:])
            row_ids = dec(f["0/META/ROW/id"][:])
            gpos = pd.Series(np.arange(len(row_ids)), index=row_ids)[gm.entrez.values].values
            want = set(s.sig_id.astype(str))
            keep = np.isin(col_ids, np.array(sorted(want)))
            n = int(keep.sum())
            nz2 = np.empty(n, np.int32); nz4 = np.empty(n, np.int32)
            mab = np.empty(n, np.float32); nan = np.empty(n, np.int32)
            tz = np.full(n, np.nan, np.float32)
            kept_ids = col_ids[keep]
            tgt_sym = s.set_index("sig_id").cmap_name
            tgt_pos = np.array([sym2pos.get(tgt_sym[i], -1) for i in kept_ids])
            k = 0
            for i in range(0, len(col_ids), BLOCK):
                sel = np.where(keep[i:i + BLOCK])[0]
                if not len(sel):
                    continue
                blk = np.asarray(mat[i:i + BLOCK, :][sel][:, gpos], dtype=np.float32)
                a = np.abs(blk)
                nz2[k:k + len(sel)] = (a >= 2).sum(1)
                nz4[k:k + len(sel)] = (a >= 4).sum(1)
                mab[k:k + len(sel)] = np.nanmean(a, axis=1)
                nan[k:k + len(sel)] = np.isnan(blk).sum(1)
                tp = tgt_pos[k:k + len(sel)]
                ok = tp >= 0
                if ok.any():
                    tz[k:k + len(sel)][ok] = blk[np.where(ok)[0], tp[ok]]
                k += len(sel)
                print(f"  {arm}: {k}/{n} ({time.time() - t0:.0f}s)", flush=True)
        m = s.set_index("sig_id").loc[kept_ids].reset_index()
        m.insert(0, "arm", m.pert_type.values)
        m["nz2"], m["nz4"], m["mean_abs"], m["n_nan"] = nz2, nz4, mab, nan
        m["target_z"] = tz
        m["target_on_gene_space"] = [t >= 0 for t in tgt_pos]
        m["target_feature_space"] = [fspace.get(g, "absent") for g in m.cmap_name]
        out_frames.append(m)
        print(f"{arm}: {len(m)} signatures, {m.cmap_name.nunique()} genes, {m.ModelID.nunique()} lines "
              f"({time.time() - t0:.0f}s)", flush=True)
    res = pd.concat(out_frames, ignore_index=True)
    res.to_csv(OUT / "perturbation_sigs.tsv.gz", sep="\t", index=False)
    prov["extracted"] = {a: dict(signatures=int((res.arm == a).sum()),
                                 genes=int(res[res.arm == a].cmap_name.nunique()),
                                 lines=int(res[res.arm == a].ModelID.nunique())) for a in res.arm.unique()}
    prov["scope_note"] = ("signatures are extracted ONLY for the 73 drug-matrix lines, which is where any test "
                          "could live; coverage counts outside those lines are reported from siginfo in s0b")
    prov["constraints_held"] = ["no drug profile read", "no concordance", "no correlation of any kind",
                                "no positive control", "no threshold chosen", "nothing re-ranked"]
    prov["threads"] = 8
    prov["runtime_s"] = round(time.time() - t0, 1)
    json.dump(prov, open(OUT / "d135_s0a_extract.json", "w"), indent=1)
    print(json.dumps(prov, indent=1))


if __name__ == "__main__":
    main()
