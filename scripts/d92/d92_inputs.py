#!/usr/bin/env python
"""D92 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): ARACNe inputs for cell-type breast networks.
Her2 (81 tumors, below the 100-tumor minimum) joins the cell-type group (luminal = LumA + LumB, basal = Basal) whose
tumor centroid (BRCA network input: log2 TPM + 1, genes centred over all BRCA tumors) correlates more with the Her2
centroid. Networks: Her2 -> luminal: BRCA_LUMX (LumA + LumB + Her2) and BRCA_BAS (Basal); Her2 -> basal: BRCA_LUM
(existing) and BRCA_BASX (Basal + Her2). Inputs built as scripts/47_brca_size_controls.py / 42 (v2 gene set, genes
with TPM >= 1 in >= 10% of the network's tumors and nonzero variance, TF-only regulator list).
Outputs: data/processed/step5/aracne_inputs_v2/<NET>.{exp.tsv,regulators.txt,regulators_tf.txt},
         results/step5/d92_breast_networks.tsv (assignment, correlations, sizes), results/step5/d92_group_net.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
IN = ROOT / "data/processed/step5/aracne_inputs_v2"
MIN_TPM, MIN_FRAC = 1.0, 0.10
x0 = pd.read_csv(IN / "BRCA.exp.tsv", sep="\t", index_col=0)
cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
sub_of = cov.drop_duplicates("sample").set_index("sample").group_id.reindex(x0.columns)
xc = x0.sub(x0.mean(axis=1), axis=0)
cen = {k: xc.loc[:, sub_of.isin(v).values].mean(axis=1) for k, v in
       {"luminal": ["BRCA-LumA", "BRCA-LumB"], "basal": ["BRCA-Basal"], "her2": ["BRCA-Her2"]}.items()}
r_lum, r_bas = np.corrcoef(cen["her2"], cen["luminal"])[0, 1], np.corrcoef(cen["her2"], cen["basal"])[0, 1]
her2_to = "luminal" if r_lum > r_bas else "basal"
print(f"Her2 centroid vs luminal r = {r_lum:.3f}, vs basal r = {r_bas:.3f} -> Her2 joins {her2_to}")
nets = ({"BRCA_LUMX": ["BRCA-LumA", "BRCA-LumB", "BRCA-Her2"], "BRCA_BAS": ["BRCA-Basal"]} if her2_to == "luminal"
        else {"BRCA_LUM": None, "BRCA_BASX": ["BRCA-Basal", "BRCA-Her2"]})
pc = set(pd.read_csv(ROOT / "results/step3/genes_protein_coding_v23_v2.tsv", sep="\t").gene_id)
regs = pd.read_csv(ROOT / "results/step5/aracne_regulators_v2.tsv", sep="\t")
build = {n: g for n, g in nets.items() if g is not None}
samp_all = [s for g in build.values() for s in x0.columns[sub_of.isin(g).values]]
tpm_f = ROOT / "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
chunks = []
for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + sorted(set(samp_all)), chunksize=5000, index_col=0):
    ch.index = ch.index.str.split(".").str[0]
    chunks.append(ch[ch.index.isin(pc)])
x = pd.concat(chunks)
x = np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0)
rows = []
for net, groups in build.items():
    samp = list(x0.columns[sub_of.isin(groups).values])
    m = x[samp]
    keep = (m >= np.log2(MIN_TPM + 1)).mean(axis=1) >= MIN_FRAC
    m = m[keep & (m.std(axis=1) > 0)]
    m.round(4).to_csv(IN / f"{net}.exp.tsv", sep="\t", index_label="gene")
    r = regs[regs.gene_id.isin(m.index)]
    r.gene_id.to_csv(IN / f"{net}.regulators.txt", index=False, header=False)
    r.loc[r.is_lambert_tf, "gene_id"].to_csv(IN / f"{net}.regulators_tf.txt", index=False, header=False)
    rows.append(dict(net=net, groups=";".join(groups), n_samples=m.shape[1], n_genes=m.shape[0], n_regulators_tf=int(r.is_lambert_tf.sum())))
gn = {"BRCA-LumA": [n for n in nets if n.startswith("BRCA_LUM")][0], "BRCA-LumB": [n for n in nets if n.startswith("BRCA_LUM")][0],
      "BRCA-Basal": [n for n in nets if n.startswith("BRCA_BAS")][0]}
gn["BRCA-Her2"] = gn["BRCA-LumA"] if her2_to == "luminal" else gn["BRCA-Basal"]
pd.DataFrame(rows).assign(her2_r_luminal=round(r_lum, 4), her2_r_basal=round(r_bas, 4), her2_to=her2_to) \
  .to_csv(ROOT / "results/step5/d92_breast_networks.tsv", sep="\t", index=False)
pd.DataFrame({"group_id": list(gn), "net": list(gn.values())}).to_csv(ROOT / "results/step5/d92_group_net.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False)); print(gn)
