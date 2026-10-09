#!/usr/bin/env python
"""Celligner: align our tumors to DepMap cell lines (Warren et al. 2021; broadinstitute/celligner d9c9246).

Runs in the `celligner` conda env (Python port pinned at commit d9c9246 + DeKegel/mnnpy bf1dd73, ~/tools/celligner).
  reference  DepMap 24Q4 OmicsExpressionProteinCodingGenesTPMLogp1 (all models in the file), log2(TPM+1)
  target     primary_per_patient tumors of every included group (Toil, all cohorts incl. hematologic),
             log2(TPM+0.001) -> log2(TPM+1). Uncorrected for purity: cPCA removes the tumor-vs-dish axes
             (largely microenvironment) itself; ESTIMATE / ABSOLUTE are used afterwards as diagnostics.
  genes      DepMap "SYMBOL (Entrez)" -> HGNC (Entrez) -> GENCODE v23 Ensembl via results/step5/id_map.tsv
             (primary_for_hgnc); protein-coding in both. Intersection count -> results/celligner/genes.tsv.
Full fit (fit on lines, transform tumors with cPCs computed): our tumors are TCGA/TARGET, i.e. the same kind
of target the published alignment uses (~10k tumors), not a small cohort projected into an existing space.
Parameters are the package defaults (PCA 70, cPCA 4, top 1,000 DE genes per dataset, MNN k1 = 5, k2 = 50).
Note: the package's limma wrapper casts the centered expression to int before the cluster DE step (DE-gene
selection only); left as published.
Outputs: data/processed/celligner/model.pkl (fitted object), aligned.npy (+ rows/genes), dist_tumor_model.tsv.gz
         results/celligner/{genes.tsv, umap_clusters.tsv, run_info.tsv}
"""
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
DEPMAP = RAW / "depmap_24q4"
SEED = 0


def reference(idmap):
    x = pd.read_csv(DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
    n_all = x.shape[1]
    entrez = pd.Series(x.columns, index=x.columns).str.extract(r"\((\d+)\)$")[0]
    m = idmap.dropna(subset=["entrez_id"]).copy()
    m["entrez_id"] = m.entrez_id.astype(int).astype(str)
    e2g = m.drop_duplicates("entrez_id").set_index("entrez_id").gene_id
    gene = entrez.map(e2g)
    keep = gene.notna() & ~gene.duplicated(keep=False)
    print(f"DepMap genes: {x.shape[1]}; mapped to v23 Ensembl (unique): {keep.sum()}")
    x = x.loc[:, keep.values]
    x.columns = gene[keep].values
    return x, n_all, int(keep.sum())


def target(genes):
    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
    inc = set(groups.loc[groups.status == "include", "group_id"])
    tum = cov[(cov.role == "tumor") & (cov.primary_per_patient == True) & cov.group_id.isin(inc)]
    tum = tum.drop_duplicates("sample")
    tpm_f = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
    use = [s for s in tum["sample"] if s in header]
    chunks = []
    for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + use, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        chunks.append(ch[ch.index.isin(genes)])
    x = pd.concat(chunks)
    x = np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0)  # log2(TPM+0.001) -> log2(TPM+1)
    print(f"target tumors: {x.shape[1]} of {len(tum)} selected")
    return x.T


def main():
    from celligner import Celligner

    PROC.mkdir(parents=True, exist_ok=True)
    RES.mkdir(parents=True, exist_ok=True)
    np.random.seed(SEED)  # cPCA uses randomized PCA without a fixed random_state
    t0 = time.time()

    idmap = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t")
    pc = set(pd.read_csv(ROOT / "results/step3/genes_protein_coding_v23.tsv", sep="\t").gene_id)
    idmap = idmap[idmap.primary_for_hgnc & idmap.gene_id.isin(pc)]

    ref, n_depmap, n_mapped = reference(idmap)
    tgt = target(set(ref.columns))
    common = [g for g in ref.columns if g in set(tgt.columns)]
    ref, tgt = ref[common], tgt[common]
    print(f"common genes: {len(common)}; reference {ref.shape}, target {tgt.shape}")
    pd.DataFrame({"gene_id": common}).merge(idmap[["gene_id", "symbol", "hgnc_id", "entrez_id"]], how="left") \
        .to_csv(RES / "genes.tsv", sep="\t", index=False)

    c = Celligner()
    c.fit(ref)
    c.transform(tgt)
    c.computeMetricsForOutput(model_ids=list(ref.index), tumor_ids=list(tgt.index))
    c.save(PROC / "model.pkl")

    out = c.combined_output
    np.save(PROC / "aligned.npy", out.values.astype(np.float32))
    pd.Series(out.index).to_csv(PROC / "aligned_rows.txt", index=False, header=False)
    pd.Series(out.columns).to_csv(PROC / "aligned_genes.txt", index=False, header=False)
    c.tumor_CL_dist.round(4).to_csv(PROC / "dist_tumor_model.tsv.gz", sep="\t", index_label="sample")

    um = c.umap_reduced.copy()
    um["cluster"] = pd.Series(c.output_clusters, index=out.index).reindex(um.index).values
    um["dataset"] = np.where(um.index.isin(ref.index), "depmap", "tumor")
    um.round(4).to_csv(RES / "umap_clusters.tsv", sep="\t", index_label="id")

    info = {
        "celligner_commit": "d9c9246", "mnnpy_commit": "bf1dd73", "depmap_release": "24Q4 (Figshare 27993248)",
        "n_models": ref.shape[0], "n_tumors": tgt.shape[0], "depmap_genes": n_depmap,
        "depmap_genes_mapped_v23": n_mapped, "common_genes": len(common), "de_genes_union": len(c.de_genes),
        "ref_clusters": len(set(c.ref_clusters)), "target_clusters": len(set(c.target_clusters)),
        "cpca_explained_var": ";".join(f"{v:.4g}" for v in c.cpca_explained_var),
        "mnn_k1": c.mnn_kwargs["k1"], "mnn_k2": c.mnn_kwargs["k2"], "pca_ncomp": c.pca_ncomp,
        "cpca_ncomp": c.cpca_ncomp, "top_k_genes": c.topKGenes, "seed": SEED,
        "minutes": round((time.time() - t0) / 60, 1),
    }
    pd.Series(info).to_csv(RES / "run_info.tsv", sep="\t", header=False)
    print(pd.Series(info).to_string())


if __name__ == "__main__":
    main()
