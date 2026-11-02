#!/usr/bin/env python
"""Step 5 (D68 / D81): ARACNe inputs for the size-matched control networks of the BRCA_LUM test.

BRCA_LUM is built from the 769 luminal (LumA + LumB) tumors of the 1,040 in the pooled BRCA network. A network built
from fewer tumors differs from the pooled one even with identical wiring, so the test compares BRCA_LUM's change with
that of networks built from 769 tumors drawn at random from all BRCA subtypes.
  BRCA_SUB<k>, k = 1..N_CTRL (default 3): 769 tumors sampled without replacement from BRCA.exp.tsv's samples (seed k)
Built exactly like BRCA_LUM in 42_aracne_inputs.py (v2): log2(TPM + 1) from the Toil TPM file, v2 gene set, genes
with TPM >= 1 in >= 10% of the network's tumors and nonzero variance; regulators = the TF-only list (Lambert TFs)
restricted to the network's genes (also the full list, for completeness).
Outputs: data/processed/step5/aracne_inputs_v2/BRCA_SUB<k>.{exp.tsv,regulators.txt,regulators_tf.txt},
         results/step5/brca_size_controls.tsv (samples per control, subtype mix)
Then: VERSION=v2 REGSET=tf bash scripts/43_run_aracne.sh BRCA_SUB1 BRCA_SUB2 BRCA_SUB3
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
IN = ROOT / "data/processed/step5/aracne_inputs_v2"
N_CTRL = int(os.environ.get("N_CTRL", "3"))
MIN_TPM, MIN_FRAC = 1.0, 0.10  # as 42_aracne_inputs.py


def main():
    pool = list(pd.read_csv(IN / "BRCA.exp.tsv", sep="\t", nrows=0).columns[1:])
    n_lum = pd.read_csv(IN / "BRCA_LUM.exp.tsv", sep="\t", nrows=0).shape[1] - 1
    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    sub_of = cov.drop_duplicates("sample").set_index("sample").group_id
    pc = set(pd.read_csv(ROOT / "results/step3/genes_protein_coding_v23_v2.tsv", sep="\t").gene_id)
    regs = pd.read_csv(ROOT / "results/step5/aracne_regulators_v2.tsv", sep="\t")

    tpm_f = ROOT / "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
    chunks = []
    for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + pool, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        chunks.append(ch[ch.index.isin(pc)])
    x = pd.concat(chunks)
    x = np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0)  # log2(TPM+0.001) -> log2(TPM+1)

    rows = []
    for k in range(1, N_CTRL + 1):
        net = f"BRCA_SUB{k}"
        samp = list(np.random.default_rng(k).choice(pool, size=n_lum, replace=False))
        m = x[samp]
        keep = (m >= np.log2(MIN_TPM + 1)).mean(axis=1) >= MIN_FRAC
        m = m[keep & (m.std(axis=1) > 0)]
        m.round(4).to_csv(IN / f"{net}.exp.tsv", sep="\t", index_label="gene")
        r = regs[regs.gene_id.isin(m.index)]
        r.gene_id.to_csv(IN / f"{net}.regulators.txt", index=False, header=False)
        r.loc[r.is_lambert_tf, "gene_id"].to_csv(IN / f"{net}.regulators_tf.txt", index=False, header=False)
        mix = sub_of.reindex(samp).value_counts()
        rows.append(dict(net=net, seed=k, n_samples=m.shape[1], n_genes=m.shape[0], n_regulators_tf=int(r.is_lambert_tf.sum()),
                         **{f"n_{g.replace('BRCA-', '')}": int(mix.get(g, 0)) for g in ["BRCA-LumA", "BRCA-LumB", "BRCA-Basal", "BRCA-Her2"]}))
    rows = pd.DataFrame(rows)
    rows.to_csv(ROOT / "results/step5/brca_size_controls.tsv", sep="\t", index=False)
    print(rows.to_string(index=False))


if __name__ == "__main__":
    main()
