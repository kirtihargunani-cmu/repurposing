#!/usr/bin/env python
"""Baseline cell-line networks for comparison with the D86 within-line drug networks (exploration 2026-11-03; runs in
`celligner`). DepMap 24Q4 log2(TPM + 1) of cancer lines, on the same 978 landmark genes (GENCODE v23 IDs, as
data/processed/lincs/pnet/*/all.exp.tsv) and the same 81 regulators. Sets: ALL (every cancer line) and BREAST
(OncotreeLineage Breast). Co-expression across lines is mostly lineage, so these networks stand for identity.
Outputs data/processed/lincs/pnet_depmap/<SET>/all.exp.tsv, regulators.txt; ARACNe with
PNET_DIR=data/processed/lincs/pnet_depmap scripts/pnet/pnet_run_aracne.sh
"""
from pathlib import Path
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
SRC = ROOT / "data/processed/lincs/pnet"
OUT = ROOT / "data/processed/lincs/pnet_depmap"

genes = pd.read_csv(SRC / "MCF7_specific/all.exp.tsv", sep="\t", usecols=["gene"]).gene
idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", dtype={"entrez_id": str})
idm = idm[idm.primary_for_hgnc == True].dropna(subset=["entrez_id"]).drop_duplicates("entrez_id")
e2g = dict(zip(idm.entrez_id, idm.gene_id))
ex = pd.read_csv(ROOT / "data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
ent = ex.columns.str.extract(r"\((\d+)\)$")[0]
ex.columns = [e2g.get(x) for x in ent]
ex = ex.loc[:, pd.Index(ex.columns).notna() & ~pd.Index(ex.columns).duplicated()]
miss = set(genes) - set(ex.columns)
print(f"landmark genes {len(genes)}, in DepMap {len(genes) - len(miss)}")
mod = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
can = mod.index[mod.OncotreePrimaryDisease != "Non-Cancerous"]
sets = {"ALL": can, "BREAST": mod.index[(mod.OncotreeLineage == "Breast") & mod.index.isin(can)]}
keep = [g for g in genes if g in ex.columns]
for name, ids in sets.items():
    x = ex.loc[ex.index.isin(ids), keep].T
    x = x[x.std(axis=1) > 0]
    d = OUT / name; d.mkdir(parents=True, exist_ok=True)
    x.index.name = "gene"
    x.to_csv(d / "all.exp.tsv", sep="\t", float_format="%.4f")
    print(f"{name}: {x.shape[1]} lines x {x.shape[0]} genes")
regs = [r for r in (SRC / "regulators.txt").read_text().split() if r in keep]
(OUT / "regulators.txt").write_text("\n".join(regs) + "\n")
print(f"regulators {len(regs)}")
