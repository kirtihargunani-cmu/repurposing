#!/usr/bin/env python
"""D105a: ChEMBL 37 mechanisms looked up by parent_molecule_chembl_id (records stored on salt forms) for the FDA drugs
(results/lincs/fda_matches.tsv). 79w fetch helpers (read-only); cache data/raw/chembl/ChEMBL_37/d105_all/mechanisms_by_parent.json.
Prints mechanism coverage of the FDA drugs before / after."""
import importlib.util, json
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
s = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(s); s.loader.exec_module(w)
w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
fm = pd.read_csv(ROOT / "results/lincs/fda_matches.tsv", sep="\t")
sp = pd.read_csv(ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv", sep="\t", usecols=["pert_id", "cmap_name"]).drop_duplicates()
drugs = sorted(set(sp[sp.pert_id.isin(fm.pert_id)].cmap_name))
d2m, mg, _, best, _ = w.chembl(drugs)          # cache hit (d105_all)
parents = sorted({m for v in d2m.values() for m in v})
mp = w.cached("mechanisms_by_parent.json", lambda: sum([w.paged("mechanism", "mechanisms", {"parent_molecule_chembl_id__in": ",".join(c)}) for c in w.chunks(parents, 40)], []))
old = json.load(open(w.CACHE / "mechanisms.json"))
has = lambda recs: {d for d, ms in d2m.items() if d in drugs and any((r["molecule_chembl_id"] in ms or r["parent_molecule_chembl_id"] in ms) and r["direct_interaction"] == 1 for r in recs)}
a, b = has(old), has(old + mp)
print(f"FDA drugs ranked: {len(drugs)} | with a ChEMBL id: {len([d for d in drugs if d in d2m])} | curated mechanism: before {len(a)}, after parent lookup {len(a | b)} (+{len(b - a)})")
print("newly found:", sorted(b - a)[:40])
# targets referenced by the parent-lookup mechanisms that are not yet cached
have = {t["target_chembl_id"] for t in json.load(open(w.CACHE / "targets.json"))}
need = sorted({r["target_chembl_id"] for r in mp if r["target_chembl_id"]} - have)
tx = w.cached("targets_by_parent.json", lambda: sum([w.paged("target", "targets", {"target_chembl_id__in": ",".join(c),
              "only": "target_chembl_id,target_type,pref_name,target_components,organism"}) for c in w.chunks(need, 50)], []))
print("extra targets fetched:", len(tx), "of", len(need))
