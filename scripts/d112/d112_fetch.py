#!/usr/bin/env python
"""D112: ChEMBL 37 curated mechanisms for FDA drugs newly matched to Drugs@FDA (or otherwise missing from the D105a cache).

The D105a cache (data/raw/chembl/ChEMBL_37/d105_all/mechanisms*.json) was fetched for the FDA drugs of 2026-11-04; drugs
matched since (D112: lapatinib, rucaparib, cyclosporine, ...) have no mechanism records there, so 79i2 judged them on the
weaker tiers (e.g. lapatinib: EGFR only, ERBB2 lost). For every ranked FDA drug whose ChEMBL ids have no mechanism record
in the cache, mechanisms are fetched by parent_molecule_chembl_id and by molecule_chembl_id and appended (deduplicated by
mec_id) to mechanisms_by_parent.json; targets they reference that are not cached are appended to targets_by_parent.json.
The live API must report ChEMBL 37 (as 79w). Read by 79i2 unchanged."""
import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
s = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(s); s.loader.exec_module(w)
w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
live = w.get(f"{w.API}/status.json")["chembl_db_version"]
assert live == w.RELEASE, f"ChEMBL API is {live}, cache is {w.RELEASE}"

fm = pd.read_csv(ROOT / "results/lincs/fda_matches.tsv", sep="\t")
sp = pd.read_csv(ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv", sep="\t", usecols=["pert_id", "cmap_name"]).drop_duplicates()
drugs = sorted(set(sp[sp.pert_id.isin(fm.pert_id)].cmap_name))
d2m, _, _, _, _ = w.chembl(drugs)
fm_, fp_ = w.CACHE / "mechanisms_by_parent.json", w.CACHE / "targets_by_parent.json"
mech = json.load(open(w.CACHE / "mechanisms.json")) + json.load(open(fm_))
covered = {x for r in mech for x in (r["molecule_chembl_id"], r["parent_molecule_chembl_id"])}
missing = {d: ms for d, ms in d2m.items() if d in drugs and not (ms & covered)}
ids = sorted({m for ms in missing.values() for m in ms})
print(f"ranked FDA drugs {len(drugs)}; with a ChEMBL id {sum(d in d2m for d in drugs)}; without a cached mechanism record {len(missing)}")
new = []
for c in w.chunks(ids, 40):
    new += w.paged("mechanism", "mechanisms", {"parent_molecule_chembl_id__in": ",".join(c)})
    new += w.paged("mechanism", "mechanisms", {"molecule_chembl_id__in": ",".join(c)})
old_ids = {r["mec_id"] for r in mech}
add = list({r["mec_id"]: r for r in new if r["mec_id"] not in old_ids}.values())
json.dump(json.load(open(fm_)) + add, open(fm_, "w"))
found = {d for d, ms in missing.items() if any({r["molecule_chembl_id"], r["parent_molecule_chembl_id"]} & ms for r in add)}
print(f"mechanism records added {len(add)}; drugs now with a record {len(found)}: {sorted(found)}")
print(f"still without a record ({len(missing) - len(found)}): {sorted(set(missing) - found)}")
have = {t["target_chembl_id"] for t in json.load(open(w.CACHE / "targets.json")) + json.load(open(fp_))}
need = sorted({r["target_chembl_id"] for r in add if r["target_chembl_id"]} - have)
tx = sum([w.paged("target", "targets", {"target_chembl_id__in": ",".join(c),
          "only": "target_chembl_id,target_type,pref_name,target_components,organism"}) for c in w.chunks(need, 50)], [])
json.dump(json.load(open(fp_)) + tx, open(fp_, "w"))
print(f"targets added {len(tx)} of {len(need)}")
