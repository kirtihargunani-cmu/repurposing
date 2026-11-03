"""ChEMBL fetching (cached) and the pipeline's binding-evidence rule, shared by 79y (stand-ins) and modelled on 79w.

Binding rule (config/drug_gate.yaml receptor_gate.binding, user 2026-11-02): a compound binds a gene if
  - a curated ChEMBL mechanism (direct_interaction = 1) on a single protein / complex / complex group / chimeric protein
    contains the gene (PROTEIN FAMILY / SELECTIVITY GROUP mechanisms give no member credit), or
  - its best valid pChEMBL on the gene is >= pchembl_min (1 uM);
  kinase genes (HGNC "Kinases") need a measured pChEMBL >= pchembl_min_kinase (100 nM) in both cases.
Valid activity = no data_validity_comment, not a potential duplicate, human single-protein target. Salt and parent
forms are pooled. Fetches refuse if the live API release differs from config/datasets.yaml; delete a cache folder to refetch.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
DS = yaml.safe_load(open(ROOT / "config/datasets.yaml"))
API = DS["chembl"]["url"]
RELEASE = DS["chembl"]["release"].split()[0]
BIND = yaml.safe_load(open(ROOT / "config/drug_gate.yaml"))["receptor_gate"]


def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception as e:  # ChEMBL returns sporadic 500s
            print("retry", i, e, file=sys.stderr)
            time.sleep(5 * (i + 1))
    raise RuntimeError(url)


def paged(ep, key, params):
    out, url = [], f"{API}/{ep}.json?" + urllib.parse.urlencode({**params, "limit": 1000})
    while url:
        p = get(url)
        out += p[key]
        nxt = p["page_meta"]["next"]
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else None
    return out


def chunks(x, n):
    return [x[i:i + n] for i in range(0, len(x), n)]


def kinase_genes():
    clo = pd.read_csv(ROOT / "data/raw/hgnc/genefamily_20261102/hierarchy_closure.csv")
    kf = set(clo[clo.parent_fam_id == BIND["kinase_family"]].child_fam_id)
    hg = pd.read_csv(sorted((ROOT / "data/raw/hgnc").glob("hgnc_complete_set_*.txt"))[-1], sep="\t", dtype=str,
                     usecols=["symbol", "gene_group_id"]).dropna()
    return {r.symbol for r in hg.itertuples() if set(map(int, r.gene_group_id.split("|"))) & kf}


class Binding:
    """Fetch (or read cached) ChEMBL data for LINCS compounds and answer binds(pert_id, gene)."""

    def __init__(self, perts, cache, log=print):
        self.cache = Path(cache)
        ann = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
        ann = ann[ann.pert_id.isin(set(perts))]
        p2c = ann.dropna(subset=["chembl_id"]).groupby("pert_id").chembl_id.apply(set).to_dict()
        ids = sorted({i for s in p2c.values() for i in s})
        forms = self._cached("forms.json", lambda: sum([paged("molecule_form", "molecule_forms",
                             {"molecule_chembl_id__in": ",".join(c)}) for c in chunks(ids, 40)], []))
        par = {f["molecule_chembl_id"]: f["parent_chembl_id"] for f in forms}
        parents = sorted(set(ids) | set(par.values()))
        log(f"ChEMBL: {len(p2c)} of {len(set(perts))} compounds have a ChEMBL id ({len(parents)} molecules incl. parents)")
        mech = self._cached("mechanisms.json", lambda: sum([paged("mechanism", "mechanisms",
                            {"molecule_chembl_id__in": ",".join(c)}) for c in chunks(parents, 40)], []))
        act = self._cached("activities.json", lambda: sum([paged("activity", "activities", {
            "molecule_chembl_id__in": ",".join(c), "pchembl_value__isnull": "false", "target_type": "SINGLE PROTEIN",
            "target_organism": "Homo sapiens",
            "only": "molecule_chembl_id,target_chembl_id,pchembl_value,data_validity_comment,potential_duplicate"})
            for c in chunks(parents, 20)], []))
        tids = sorted({m["target_chembl_id"] for m in mech if m["target_chembl_id"]} | {a["target_chembl_id"] for a in act})
        tg = self._cached("targets.json", lambda: sum([paged("target", "targets", {"target_chembl_id__in": ",".join(c),
                          "only": "target_chembl_id,target_type,pref_name,target_components"}) for c in chunks(tids, 50)], []))
        T = {}
        for t in tg:
            T[t["target_chembl_id"]] = (t["target_type"], frozenset(
                s["component_synonym"] for c in t.get("target_components") or [] for s in c.get("target_component_synonyms") or []
                if s["syn_type"] == "GENE_SYMBOL"))
        credit = set(BIND["binding"]["mechanism_target_types"])
        self.mech = {}
        for m in mech:
            if m["direct_interaction"] != 1 or not m["target_chembl_id"]:
                continue
            tt, gs = T.get(m["target_chembl_id"], ("", frozenset()))
            if tt in credit:
                for mol in {m["molecule_chembl_id"], m["parent_molecule_chembl_id"]}:
                    self.mech.setdefault(mol, set()).update(gs)
        a = pd.DataFrame(act)
        self.best = {}
        if len(a):
            a = a[a.data_validity_comment.isna() & (a.potential_duplicate.astype(str) != "1")].copy()
            a["v"] = a.pchembl_value.astype(float)
            a["g"] = a.target_chembl_id.map(lambda t: T.get(t, ("", frozenset()))[1])
            a = a[a.g.map(len) == 1].assign(g=lambda d: d.g.map(lambda s: next(iter(s))))
            for (mol, g), v in a.groupby(["molecule_chembl_id", "g"]).v.max().items():
                self.best.setdefault(mol, {})[g] = v
        self.p2m = {p: s | {par.get(m, m) for m in s} for p, s in p2c.items()}
        self.kin = kinase_genes()
        self.n = dict(compounds=len(set(perts)), with_chembl=len(p2c), molecules=len(parents), mechanisms=len(mech),
                      activities=len(act), targets=len(tg))

    def _cached(self, name, fn):
        f = self.cache / name
        if f.exists():
            return json.load(open(f))
        live = get(f"{API}/status.json")["chembl_db_version"]
        assert live == RELEASE, f"ChEMBL API is {live}, config pins {RELEASE}"
        self.cache.mkdir(parents=True, exist_ok=True)
        x = fn()
        json.dump(x, open(f, "w"))
        return x

    def binds(self, pert, gene):
        ms = self.p2m.get(pert, set())
        if not ms:
            return False
        pv = [self.best[m][gene] for m in ms if gene in self.best.get(m, {})]
        p = max(pv) if pv else None
        if gene in self.kin:
            return p is not None and p >= BIND["binding"]["pchembl_min_kinase"]
        return any(gene in self.mech.get(m, ()) for m in ms) or (p is not None and p >= BIND["binding"]["pchembl_min"])
