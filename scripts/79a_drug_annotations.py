#!/usr/bin/env python
"""Pharmacology annotations for LINCS drugs with an active signature in a cancer line (user, 2026-10-09: item 3).

Sources: Broad Drug Repurposing Hub 2020-03-24 (scripts/09_download_drug_annotations.sh; joined on the 13-character
Broad ID = LINCS pert_id, then InChIKey) and ChEMBL (REST API, release recorded; joined on InChIKey from LINCS
compoundinfo or the Hub). Raw ChEMBL JSON cached in data/raw/chembl/<release>/.
Flags (for review, not automatic exclusion; no open source records which organs a drug reaches):
  route_topical_only   ChEMBL topical and not oral / parenteral
  atc_local            every ATC code of the drug is in a locally acting group (LOCAL_ATC: mouth, gut-lumen,
                       skin except its oral subgroups D01B / D05B / D10B, gynecological anti-infective, nose /
                       throat, eye / ear, diagnostic / contrast); gut-acting groups count whatever the route
  organ_systems        ATC level-1 anatomical groups (e.g. A alimentary tract and metabolism, C cardiovascular)
  approved_indications ChEMBL indications at phase 4 (EFO terms)
Outputs: results/lincs/drug_annotations.tsv; data/raw/chembl/<release>/{molecules,indications,atc}.json
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
RES = ROOT / "results/lincs"
API = "https://www.ebi.ac.uk/chembl/api/data"
TAS_MIN = 0.4
LOCAL_ATC = ("A01", "A06", "A07", "D", "G01", "R01", "R02", "S01", "S02", "S03", "V04", "V08", "V09")
SYSTEMIC_IN_LOCAL = ("D01B", "D05B", "D10B")  # oral antifungals / antipsoriatics / anti-acne within group D
ATC1 = {"A": "alimentary tract and metabolism", "B": "blood and blood forming organs", "C": "cardiovascular",
        "D": "dermatologicals", "G": "genito-urinary and sex hormones", "H": "systemic hormones (excl. sex hormones)",
        "J": "anti-infectives for systemic use", "L": "antineoplastic and immunomodulating", "M": "musculo-skeletal",
        "N": "nervous system", "P": "antiparasitic", "R": "respiratory", "S": "sensory organs", "V": "various"}


def get(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return json.load(r)
        except Exception as e:  # transient API errors
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def paged(endpoint, params):
    out, url = [], f"{API}/{endpoint}.json?" + urllib.parse.urlencode(params)
    while url:
        page = get(url)
        key = [k for k in page if k != "page_meta"][0]
        out += page[key]
        nxt = page["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return out


def chembl(inchikeys, cache):
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / "molecules.json"
    if not f.exists():
        mols = []
        keys = sorted(inchikeys)
        for i in range(0, len(keys), 50):
            mols += paged("molecule", {"molecule_structures__standard_inchi_key__in": ",".join(keys[i:i + 50]),
                                       "only": "molecule_chembl_id,pref_name,max_phase,oral,parenteral,topical,"
                                               "atc_classifications,first_approval,withdrawn_flag,molecule_structures",
                                       "limit": 1000})
        json.dump(mols, open(f, "w"))
    mols = json.load(open(f))
    ids = sorted({m["molecule_chembl_id"] for m in mols})
    f = cache / "indications.json"
    if not f.exists():
        ind = []
        for i in range(0, len(ids), 50):
            ind += paged("drug_indication", {"molecule_chembl_id__in": ",".join(ids[i:i + 50]), "limit": 1000,
                                             "only": "molecule_chembl_id,efo_term,mesh_heading,max_phase_for_ind"})
        json.dump(ind, open(f, "w"))
    ind = json.load(open(f))
    return mols, ind


def main():
    status = get(f"{API}/status.json")
    rel = status["chembl_db_version"]
    cache = RAW / "chembl" / rel

    sig = pd.read_csv(RAW / "lincs2020/siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["pert_id", "pert_type", "cell_iname", "is_hiq", "tas"])
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    cancer = set(cov.loc[cov.OncotreePrimaryDisease != "Non-Cancerous", "cell_iname"])
    act = sig[(sig.pert_type == "trt_cp") & (sig.is_hiq == 1) & (sig.tas >= TAS_MIN) & sig.cell_iname.isin(cancer)]
    nlines = act.groupby("pert_id").cell_iname.nunique()
    drugs = pd.DataFrame({"pert_id": nlines.index, "n_active_lines": nlines.values})

    ci = pd.read_csv(RAW / "lincs2020/compoundinfo_beta.txt", sep="\t", na_values=['""'], keep_default_na=False)
    lin = ci.groupby("pert_id").agg(cmap_name=("cmap_name", "first"), inchi_key=("inchi_key", "first"),
                                    lincs_target=("target", lambda x: "|".join(sorted(set(x.dropna()) - {""}))),
                                    lincs_moa=("moa", lambda x: "|".join(sorted(set(x.dropna()) - {""}))))
    drugs = drugs.merge(lin, left_on="pert_id", right_index=True, how="left")

    hub_d = pd.read_csv(RAW / "repurposing_hub/repurposing_drugs_20200324.txt", sep="\t", comment="!")
    hub_s = pd.read_csv(RAW / "repurposing_hub/repurposing_samples_20200324.txt", sep="\t", comment="!")
    hub_s["pert_id"] = hub_s.broad_id.str[:13]
    hub = hub_s.drop_duplicates("pert_id")[["pert_id", "pert_iname", "InChIKey"]].merge(hub_d, on="pert_iname", how="left")
    hub_by_key = hub.dropna(subset=["InChIKey"]).drop_duplicates("InChIKey").set_index("InChIKey")
    h = drugs[["pert_id"]].merge(hub, on="pert_id", how="left")
    miss = h.pert_iname.isna().values & drugs.inchi_key.notna().values
    h.loc[miss, hub.columns[1:]] = hub_by_key.reindex(drugs.inchi_key[miss]).reset_index().reindex(columns=hub.columns[1:]).values
    drugs["hub_name"] = h.pert_iname.values
    for c in ["clinical_phase", "moa", "target", "disease_area", "indication"]:
        drugs[f"hub_{c}"] = h[c].values
    drugs["inchi_key"] = drugs.inchi_key.fillna(h.InChIKey)

    mols, ind = chembl(set(drugs.inchi_key.dropna()), cache)
    m = pd.DataFrame([{**{k: v for k, v in x.items() if k != "molecule_structures"},
                       "inchi_key": (x.get("molecule_structures") or {}).get("standard_inchi_key")} for x in mols])
    m = m.dropna(subset=["inchi_key"]).sort_values("max_phase", ascending=False).drop_duplicates("inchi_key")
    m["atc"] = m.atc_classifications.apply(lambda a: "|".join(sorted(a)) if a else "")
    m["atc_level1"] = m.atc_classifications.apply(lambda a: "|".join(sorted({c[0] for c in a})) if a else "")
    m["organ_systems"] = m.atc_level1.apply(lambda s: "|".join(ATC1[c] for c in s.split("|")) if s else "")
    m["atc_local"] = m.atc_classifications.apply(
        lambda a: bool(a) and all(c.startswith(LOCAL_ATC) and not c.startswith(SYSTEMIC_IN_LOCAL) for c in a))
    m["route_topical_only"] = m.topical.eq(True) & ~m.oral.eq(True) & ~m.parenteral.eq(True)
    ind = pd.DataFrame(ind)
    if len(ind):
        ind["max_phase_for_ind"] = pd.to_numeric(ind.max_phase_for_ind, errors="coerce")
        appr = ind[ind.max_phase_for_ind >= 4].groupby("molecule_chembl_id").efo_term.apply(
            lambda x: "|".join(sorted(set(x.dropna()))))
        m["approved_indications"] = m.molecule_chembl_id.map(appr).fillna("")
    keep = ["molecule_chembl_id", "pref_name", "max_phase", "first_approval", "withdrawn_flag", "oral", "parenteral",
            "topical", "route_topical_only", "atc", "organ_systems", "atc_local", "approved_indications"]
    drugs = drugs.merge(m.set_index("inchi_key")[keep].add_prefix("chembl_").rename(
        columns={"chembl_molecule_chembl_id": "chembl_id", "chembl_route_topical_only": "route_topical_only",
                 "chembl_atc_local": "atc_local", "chembl_organ_systems": "organ_systems",
                 "chembl_atc": "atc", "chembl_approved_indications": "approved_indications"}),
        left_on="inchi_key", right_index=True, how="left")
    drugs["pharmacology_flag"] = drugs.route_topical_only.eq(True) | drugs.atc_local.eq(True)
    drugs["chembl_release"] = rel
    drugs = drugs.sort_values("n_active_lines", ascending=False)
    drugs.to_csv(RES / "drug_annotations.tsv", sep="\t", index=False)

    n = len(drugs)
    print(f"{rel}; drugs active in a cancer line: {n}")
    print(f"  Repurposing Hub match: {drugs.hub_name.notna().sum()}; phase: "
          f"{drugs.hub_clinical_phase.value_counts().to_dict()}")
    print(f"  ChEMBL match: {drugs.chembl_id.notna().sum()}; with ATC: {(drugs.atc.fillna('') != '').sum()}; "
          f"max_phase 4: {(pd.to_numeric(drugs.chembl_max_phase, errors='coerce') >= 4).sum()}")
    print(f"  target annotation (LINCS or Hub): {((drugs.lincs_target.fillna('') != '') | drugs.hub_target.notna()).sum()}")
    print(f"  flagged locally acting: topical-only {drugs.route_topical_only.eq(True).sum()}, "
          f"ATC local {drugs.atc_local.eq(True).sum()}, either {drugs.pharmacology_flag.sum()}")
    print(drugs.loc[drugs.pharmacology_flag, ["cmap_name", "n_active_lines", "atc", "organ_systems"]].head(25).to_string(index=False))


if __name__ == "__main__":
    main()
