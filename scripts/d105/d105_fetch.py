#!/usr/bin/env python
"""D105: ChEMBL 37 mechanisms, molecule forms and human single-protein activities for every ranked drug (all cmap_names in
the production split table), with 79w's fetch code (read-only import) into its own cache:
data/raw/chembl/ChEMBL_37/d105_all/ (79w's receptor_gate cache is not touched). Refuses if the live API release differs
from config/datasets.yaml."""
import importlib.util
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
s = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(s); s.loader.exec_module(w)
w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
drugs = sorted(set(pd.read_csv(ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv", sep="\t", usecols=["cmap_name"]).cmap_name.dropna()))
print("ranked drugs:", len(drugs), flush=True)
d2m, mg, munit, best, n = w.chembl(drugs)
print(n, "| drugs with a ChEMBL id:", len(d2m), flush=True)
