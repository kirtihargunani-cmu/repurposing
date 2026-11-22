#!/usr/bin/env python
"""Resolve ChEMBL ids by NAME for scored drugs that have none (D181b, 2026-11-22).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python  (run with -I)

WHY. D181 found that a drug with no ChEMBL id can be neither labelled nor marked `tried_elsewhere`, so it
lands silently in the "no record" bucket and reads as a negative. The 65 FDA-matched drugs still unresolved
after the name-pooling fix are mostly NON-oncology - acyclovir, citalopram, diclofenac, chloroquine - which
is exactly the repurposing population, so the bias runs against the thing we are looking for.

WHAT IT DOES. Exactly what D113 and 79w's `name_lookup` already do, applied to the remaining set: an exact
preferred-name lookup, then a synonym lookup, against the pinned ChEMBL release. Fetch only - it writes a
cache that `79aa_label_pairs.py` reads. No labelling, no scoring, no interpretation here.

TARGETS (env SCOPE): `fda` (default) = drugs flagged `fda_matched` in the D181 table with no resolved id;
`approved` = `approved_anywhere`; `named` = every named compound (skips `BRD-...` codes, which are unnamed
research compounds with no ChEMBL molecule record and are not worth querying).

Names tried per drug, first hit wins: cmap_name, its hyphens as spaces, the Drugs@FDA ingredient, the Hub
name, the ChEMBL preferred name. Then `molecule_synonyms__molecule_synonym__iexact` on the same list (the nested field; the flat
`molecule_synonyms__iexact` is a 400). A hit is recorded with
the name and endpoint that produced it, so a wrong match can be traced and removed by hand.

SAFETY. Generic User-Agent, no personal details, ever. Refuses if the live release differs from
config/datasets.yaml, as chembl_lib does.

Output: data/raw/chembl/ChEMBL_37/d105_all/d181_ids.json
        {cmap_name: {ids: [...], via: pref_name|synonym, matched_name: str, chembl_pref_name: str}}
        + d181_ids_run.json (counts, unresolved list)
"""
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
CACHE = ROOT / "data/raw/chembl/ChEMBL_37/d105_all"
OUT = CACHE / "d181_ids.json"
SCOPE = os.environ.get("SCOPE", "fda")
UA = "repurposing-pipeline/1.0 (academic research; python-urllib)"   # no personal details, ever
DS = yaml.safe_load(open(ROOT / "config/datasets.yaml"))
API = DS["chembl"]["url"]
RELEASE = DS["chembl"]["release"].split()[0]
norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def get(url, tries=5):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            print(f"  retry {i}: {e}", file=sys.stderr)
            time.sleep(3 * (i + 1))
    raise RuntimeError(url)


def query(field, name):
    q = urllib.parse.urlencode({field: name, "limit": 5, "only": "molecule_chembl_id,pref_name"})
    return get(f"{API}/molecule.json?{q}").get("molecules", [])


def main():
    live = get(f"{API}/status.json")["chembl_db_version"]
    assert live == RELEASE, f"ChEMBL API is {live}, config pins {RELEASE}"
    d = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    u = d.drop_duplicates("pert_id")
    col = {"fda": "fda_matched", "approved": "approved_anywhere"}.get(SCOPE)
    sel = u[u[col]] if col else u
    sel = sel[~sel.has_chembl_id & ~sel.cmap_name.astype(str).str.startswith("BRD-")]
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t", low_memory=False)
    fda_ing = fm.assign(k=fm.cmap_name.map(norm)).drop_duplicates("k").set_index("k").fda_ingredient.to_dict()

    todo = sorted(set(sel.cmap_name.dropna()))
    print(f"SCOPE={SCOPE}: {len(todo)} named drugs with no resolved ChEMBL id")
    out = json.load(open(OUT)) if OUT.exists() else {}
    todo = [n for n in todo if n not in out]
    print(f"{len(todo)} not already cached\n")

    unresolved = []
    for i, name in enumerate(todo, 1):
        p = sel[sel.cmap_name == name].pert_id.iloc[0]
        cands, seen = [], set()
        for c in (name, name.replace("-", " "), fda_ing.get(norm(name)),
                  ann.hub_name.get(p), ann.chembl_pref_name.get(p)):
            if isinstance(c, str) and c.strip() and c.lower() not in seen:
                seen.add(c.lower()); cands.append(c.strip())
        hit = None
        for field in ("pref_name__iexact", "molecule_synonyms__molecule_synonym__iexact"):
            for c in cands:
                m = query(field, c)
                if m:
                    hit = dict(ids=[m[0]["molecule_chembl_id"]], via="pref_name" if field.startswith("pref") else "synonym",
                               matched_name=c, chembl_pref_name=m[0].get("pref_name"))
                    break
            if hit:
                break
        if hit:
            out[name] = hit
            print(f"  [{i}/{len(todo)}] {name:<28} -> {hit['ids'][0]:<16} ({hit['via']}: {hit['matched_name']})")
        else:
            unresolved.append(name)
            print(f"  [{i}/{len(todo)}] {name:<28} -- no match")
        if i % 10 == 0:                      # save as we go: a 400 or a timeout must not lose the fetch
            json.dump(out, open(OUT, "w"), indent=1, sort_keys=True)
        time.sleep(0.15)

    json.dump(out, open(OUT, "w"), indent=1, sort_keys=True)
    run = dict(scope=SCOPE, release=live, queried=len(todo), resolved=len(todo) - len(unresolved),
               unresolved=unresolved, cache_total=len(out),
               fetched=time.strftime("%Y-%m-%dT%H:%M:%S"))
    json.dump(run, open(CACHE / "d181_ids_run.json", "w"), indent=1)
    print(f"\nresolved {run['resolved']} of {len(todo)}; cache now {len(out)} names -> {OUT.relative_to(ROOT)}")
    if unresolved:
        print(f"unresolved: {unresolved}")


if __name__ == "__main__":
    main()
