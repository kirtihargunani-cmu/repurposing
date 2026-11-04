#!/usr/bin/env python
"""D92 / D92a reported arms. NET / SUBTYPES env choose the network (default BRCA_HER2 = BRCA-Her2); D92a also BRCA_BAS = BRCA-Basal.
Original: D92 reported arm (user 2026-11-03: subtypes with > 50 tumors are large enough): a Her2-only breast network (BRCA_HER2,
81 tumors), built exactly like scripts/d92/d92_inputs.py. Not part of the D92 decision (Her2 is pooled with basal there)."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
IN = ROOT / "data/processed/step5/aracne_inputs_v2"
import os as _os
x0 = pd.read_csv(IN / f"{_os.environ.get('COHORT', 'BRCA')}.exp.tsv", sep="\t", nrows=0).columns[1:]   # COHORT (D96): source cohort input
cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
sub_of = cov.drop_duplicates("sample").set_index("sample").group_id.reindex(x0)
import os
NET = os.environ.get("NET", "BRCA_HER2"); GROUPS = os.environ.get("SUBTYPES", "BRCA-Her2").split(",")
samp = list(x0[sub_of.isin(GROUPS).values])
pc = set(pd.read_csv(ROOT / "results/step3/genes_protein_coding_v23_v2.tsv", sep="\t").gene_id)
regs = pd.read_csv(ROOT / "results/step5/aracne_regulators_v2.tsv", sep="\t")
tpm_f = ROOT / "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
ch = [c[c.index.isin(pc)] for c in (lambda it: (setattr(c, "index", c.index.str.split(".").str[0]) or c for c in it))(
      pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + samp, chunksize=5000, index_col=0))]
m = np.log2(np.power(2.0, pd.concat(ch)) - 0.001 + 1).clip(lower=0)
m = m[((m >= 1).mean(axis=1) >= 0.10) & (m.std(axis=1) > 0)]
m.round(4).to_csv(IN / f"{NET}.exp.tsv", sep="\t", index_label="gene")
r = regs[regs.gene_id.isin(m.index)]
r.gene_id.to_csv(IN / f"{NET}.regulators.txt", index=False, header=False)
r.loc[r.is_lambert_tf, "gene_id"].to_csv(IN / f"{NET}.regulators_tf.txt", index=False, header=False)
print(f"{NET}: {m.shape[1]} tumors x {m.shape[0]} genes; TF regulators {int(r.is_lambert_tf.sum())}")
