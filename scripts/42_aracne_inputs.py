#!/usr/bin/env python
"""Step 5c: inputs for project-built ARACNe networks (one per cohort).

Networks are built per cohort (subtypes pooled; ARACNe wants ~100+ samples), from the same Toil TPM
data as the Step 4 signatures, so regulons and signatures share pipeline and GENCODE v23 IDs.
  samples     primary_per_patient tumors of every included group in the cohort; COAD + READ -> COADREAD
  expression  log2(TPM + 1), protein-coding genes with TPM >= 1 in >= 10% of the cohort's tumors
  regulators  full set: Lambert 2018 TFs + every regulator in aracne.networks (TFs, co-factors,
              signaling proteins), mapped to v23 Ensembl; kinases flagged from HGNC gene groups
Only cohorts with at least one Step 4 signature are built.
Outputs: data/processed/step5/aracne_inputs/<NET>.exp.tsv (genes x samples, Ensembl IDs)
         results/step5/aracne_regulators.tsv, results/step5/aracne_cohorts.tsv
VERSION=v2 (rebuild spec, D46 / D51): genes and regulators restricted to the v2 gene set
(results/step3/genes_protein_coding_v23_v2.tsv: no mappability < 0.8, no readthrough loci); everything else as v1.
Outputs then go to aracne_inputs_v2/, aracne_regulators_v2.tsv, aracne_cohorts_v2.tsv. v2 also writes (D64, D65):
  <NET>.regulators_tf.txt   Lambert TFs only, for the TF-only build (43_run_aracne.sh REGSET=tf)
  BRCA_LUM.*                extra luminal-only breast network (BRCA-LumA + BRCA-LumB tumors)
"""
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
VERSION = os.environ.get("VERSION", "v1")
assert VERSION in ("v1", "v2")
SFX = "" if VERSION == "v1" else "_v2"
OUT = ROOT / f"data/processed/step5/aracne_inputs{SFX}"
MERGE = {"TCGA-COAD": "COADREAD", "TCGA-READ": "COADREAD"}
MIN_TPM, MIN_FRAC = 1.0, 0.10
LUMINAL = ["BRCA-LumA", "BRCA-LumB"]


def net_name(cohort):
    return MERGE.get(cohort, cohort.split("-", 1)[1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv(ROOT / "results/step4/signature_summary.tsv", sep="\t")
    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    signed = set(summ.loc[summ.status == "ok", "group_id"])
    cohorts = sorted(set(groups.loc[groups.group_id.isin(signed), "cohort"]))
    inc = set(groups.loc[groups.status == "include", "group_id"])
    tum = cov[(cov.role == "tumor") & cov.primary_per_patient & cov.group_id.isin(inc) & cov.cohort.isin(cohorts)]
    tum = tum.assign(net=tum.cohort.map(net_name)).drop_duplicates("sample")
    if VERSION == "v2":  # extra luminal-only breast network (D65); its tumors also stay in BRCA
        tum = pd.concat([tum, tum[tum.group_id.isin(LUMINAL)].assign(net="BRCA_LUM")], ignore_index=True)

    # regulators: Lambert TFs + aracne.networks regulators
    idmap = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t")
    gene_file = "genes_protein_coding_v23.tsv" if VERSION == "v1" else "genes_protein_coding_v23_v2.tsv"
    pc = set(pd.read_csv(ROOT / "results/step3" / gene_file, sep="\t").gene_id)
    tf = set(pd.read_csv(ROOT / "results/step5/tf_universe.tsv", sep="\t").gene_id)
    arac = set()
    for f in sorted((ROOT / "data/processed/step5").glob("regulons_aracne_*.tsv")):
        arac |= set(pd.read_csv(f, sep="\t", usecols=["source_gene_id"]).source_gene_id.dropna())
    hgnc = pd.read_csv(sorted((RAW / "hgnc").glob("hgnc_complete_set_*.txt"))[-1], sep="\t",
                       usecols=["hgnc_id", "gene_group", "locus_type"], dtype=str)
    kin = set(hgnc.loc[hgnc.gene_group.fillna("").str.contains(r"kinase", flags=re.I)
                       & ~hgnc.gene_group.fillna("").str.contains(r"kinase (?:regulatory|associated|anchor|inhibitor)", flags=re.I),
                       "hgnc_id"])
    regs = idmap[idmap.gene_id.isin((tf | arac) & pc)].copy()
    regs["is_lambert_tf"] = regs.gene_id.isin(tf)
    regs["in_aracne_networks"] = regs.gene_id.isin(arac)
    regs["is_kinase"] = regs.hgnc_id.isin(kin)
    regs = regs[["gene_id", "symbol", "is_lambert_tf", "in_aracne_networks", "is_kinase"]]

    tpm_f = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
    use = [s for s in tum["sample"].unique() if s in header]
    chunks = []
    for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + use, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        chunks.append(ch[ch.index.isin(pc)])
    x = pd.concat(chunks)
    x = np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0)  # log2(TPM+0.001) -> log2(TPM+1)

    rows = []
    for net, d in tum[tum["sample"].isin(x.columns)].groupby("net"):
        m = x[d["sample"].tolist()]
        keep = (m >= np.log2(MIN_TPM + 1)).mean(axis=1) >= MIN_FRAC
        m = m[keep & (m.std(axis=1) > 0)]
        m.round(4).to_csv(OUT / f"{net}.exp.tsv", sep="\t", index_label="gene")
        r = regs[regs.gene_id.isin(m.index)]
        r.gene_id.to_csv(OUT / f"{net}.regulators.txt", index=False, header=False)
        if VERSION == "v2":
            r.loc[r.is_lambert_tf, "gene_id"].to_csv(OUT / f"{net}.regulators_tf.txt", index=False, header=False)
        rows.append(dict(net=net, cohorts=";".join(sorted(d.group_id.unique() if net == "BRCA_LUM" else d.cohort.unique())), n_samples=m.shape[1],
                         n_genes=m.shape[0], n_regulators=len(r), n_regulators_tf=int(r.is_lambert_tf.sum()),
                         n_regulators_kinase=int(r.is_kinase.sum())))
    regs.to_csv(ROOT / f"results/step5/aracne_regulators{SFX}.tsv", sep="\t", index=False)
    rows = pd.DataFrame(rows)
    rows.to_csv(ROOT / f"results/step5/aracne_cohorts{SFX}.tsv", sep="\t", index=False)
    print(f"regulators: {len(regs)} (Lambert TF {regs.is_lambert_tf.sum()}, "
          f"in aracne.networks {regs.in_aracne_networks.sum()}, kinase {regs.is_kinase.sum()})")
    print(rows.to_string(index=False))


if __name__ == "__main__":
    main()
