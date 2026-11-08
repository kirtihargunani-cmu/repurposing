#!/usr/bin/env python
"""D113: ChEMBL ids of ranked FDA drugs resolved from the FDA ingredient name (docs/DECISIONS.md D113).

For every ranked FDA drug none of whose ChEMBL ids (79w.chembl: annotation + name lookup + parents) has a ChEMBL 37 mechanism
record: FDA ingredient minus salt words (fda_rows.drug_name) -> exact ChEMBL preferred name, else exact synonym -> parent.
Same InChIKey skeleton as the annotated entry -> added; different compound (or no annotated id) -> replaces. Zero hits or
more than one parent -> unresolved, listed. Writes d105_all/d113_ids.json {cmap_name: {"ids": [...], "mode": add|replace, ...}}
and appends mechanisms / targets / activities of the resolved parents to the *_by_parent.json caches. Amendment (D113):
resolutions whose parent has no mechanism record were then dropped from d113_ids.json (12 of 44; see DECISIONS).
DRY=1: print the mapping only (no fetch of mechanisms, nothing written)."""
import importlib.util
import json
import os
import sys
import urllib.parse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import drug_name  # noqa: E402

s = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(s); s.loader.exec_module(w)
w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
live = w.get(f"{w.API}/status.json")["chembl_db_version"]
assert live == w.RELEASE, f"ChEMBL API is {live}, cache is {w.RELEASE}"
DRY = os.environ.get("DRY") == "1"
OUT = w.CACHE / "d113_ids.json"
assert DRY or not OUT.exists(), f"{OUT} exists: D113 already applied"

fm = pd.read_csv(ROOT / "results/lincs/fda_matches.tsv", sep="\t")
sp = pd.read_csv(ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv", sep="\t", usecols=["pert_id", "cmap_name"]).drop_duplicates()
m = sp[sp.pert_id.isin(fm.pert_id)].merge(fm[["pert_id", "fda_ingredient"]], on="pert_id")
drugs = sorted(set(m.cmap_name))
d2m, _, _, _, _ = w.chembl(drugs)
mech = json.load(open(w.CACHE / "mechanisms.json")) + json.load(open(w.CACHE / "mechanisms_by_parent.json"))
covered = {x for r in mech for x in (r["molecule_chembl_id"], r["parent_molecule_chembl_id"])}
missing = [d for d in drugs if not (d2m.get(d, set()) & covered)]
print(f"ranked FDA drugs {len(drugs)}; without a mechanism record {len(missing)}")

ONLY = "molecule_chembl_id,pref_name,molecule_hierarchy,molecule_structures"


def mol(q):
    r = w.get(f"{w.API}/molecule.json?" + urllib.parse.urlencode({**q, "only": ONLY, "limit": 50}))["molecules"]
    return r


def skel(x):
    k = ((x or {}).get("molecule_structures") or {}).get("standard_inchi_key")
    return k[:14] if k else None


def parent(x):
    return ((x.get("molecule_hierarchy") or {}).get("parent_chembl_id")) or x["molecule_chembl_id"]


res, rows = {}, []
for d in missing:
    names = sorted({drug_name(i).upper() for i in m[m.cmap_name == d].fda_ingredient})
    hits, via = [], ""
    for n in names:
        h = mol({"pref_name__iexact": n})
        if h:
            hits += h; via = "pref_name"
        else:
            h = mol({"molecule_synonyms__molecule_synonym__iexact": n})
            hits += h; via = via or ("synonym" if h else "")
    pars = sorted({parent(x) for x in hits})
    ann_ids = sorted(d2m.get(d, set()))
    if len(pars) != 1:
        rows.append(dict(cmap_name=d, fda_name="|".join(names), ann_ids=";".join(ann_ids), resolved="", via=via,
                         mode="unresolved", note=f"{len(pars)} parents: {';'.join(pars)}"))
        continue
    p = pars[0]
    pm = mol({"molecule_chembl_id": p})
    ps = skel(pm[0]) if pm else None
    am = sum([mol({"molecule_chembl_id__in": ",".join(ann_ids)})], []) if ann_ids else []
    same = bool(ps) and any(skel(x) == ps for x in am)
    mode = "add" if same else "replace"
    an = ";".join(f"{x['molecule_chembl_id']}:{x.get('pref_name') or '-'}" for x in am)
    res[d] = {"ids": [p], "mode": mode, "fda_name": "|".join(names), "via": via, "annotated": ann_ids}
    rows.append(dict(cmap_name=d, fda_name="|".join(names), ann_ids=an, resolved=f"{p}:{pm[0].get('pref_name') if pm else ''}",
                     via=via, mode=mode, note=""))
tab = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60); pd.set_option("display.max_rows", 200)
print(tab.to_string(index=False))
print(tab["mode"].value_counts().to_string())
tab.to_csv(ROOT / "results/lincs/d113_chembl_ids.tsv", sep="\t", index=False)
if DRY:
    sys.exit(0)

ids = sorted({i for v in res.values() for i in v["ids"]})
fmp, ftp, fap = w.CACHE / "mechanisms_by_parent.json", w.CACHE / "targets_by_parent.json", w.CACHE / "activities_by_parent.json"
new = []
for c in w.chunks(ids, 40):
    new += w.paged("mechanism", "mechanisms", {"parent_molecule_chembl_id__in": ",".join(c)})
    new += w.paged("mechanism", "mechanisms", {"molecule_chembl_id__in": ",".join(c)})
old = {r["mec_id"] for r in mech}
add = list({r["mec_id"]: r for r in new if r["mec_id"] not in old}.values())
json.dump(json.load(open(fmp)) + add, open(fmp, "w"))
got = {d for d, v in res.items() if any({r["molecule_chembl_id"], r["parent_molecule_chembl_id"]} & set(v["ids"]) for r in add + mech)}
print(f"mechanism records added {len(add)}; resolved drugs with a record {len(got)} of {len(res)}; without: {sorted(set(res) - got)}")
act = sum([w.paged("activity", "activities", {
    "molecule_chembl_id__in": ",".join(c), "pchembl_value__isnull": "false", "target_type": "SINGLE PROTEIN",
    "target_organism": "Homo sapiens",
    "only": "molecule_chembl_id,target_chembl_id,pchembl_value,standard_type,data_validity_comment,potential_duplicate"})
    for c in w.chunks(ids, 20)], [])
json.dump((json.load(open(fap)) if fap.exists() else []) + act, open(fap, "w"))
print(f"activities added {len(act)}")
have = {t["target_chembl_id"] for t in json.load(open(w.CACHE / "targets.json")) + json.load(open(ftp))}
need = sorted(({r["target_chembl_id"] for r in add if r["target_chembl_id"]} | {a["target_chembl_id"] for a in act}) - have)
tx = sum([w.paged("target", "targets", {"target_chembl_id__in": ",".join(c),
          "only": "target_chembl_id,target_type,pref_name,target_components,organism"}) for c in w.chunks(need, 50)], [])
json.dump(json.load(open(ftp)) + tx, open(ftp, "w"))
print(f"targets added {len(tx)} of {len(need)}")
# Amendment (D113, before any pipeline step): keep a resolution only if its parent has a mechanism record
res = {d: v for d, v in res.items() if d in got}
tab["applied"] = tab.cmap_name.isin(res)
tab.loc[~tab.applied & tab["mode"].ne("unresolved"), "note"] = "no mechanism record for the resolved parent (amendment): unchanged"
tab.to_csv(ROOT / "results/lincs/d113_chembl_ids.tsv", sep="\t", index=False)
json.dump(res, open(OUT, "w"), indent=1)
print(f"wrote {OUT}: {len(res)} drugs")
