#!/usr/bin/env python
"""Receptor evidence gate: does a drug's expressed receptor target count as the drug's receptor in this cohort?
(user 2026-11-02; rules and thresholds in config/drug_gate.yaml: receptor_gate). Flags only: no drug leaves a ranking.

Rows: every disease (41) and v2 subgroup (14) x drug that passed the target gate (79b) x expressed target that is a
receptor. Receptor classes = HGNC gene-group roots in config (GPCR, nuclear receptor, RTK, receptor Ser/Thr kinase,
ligand-gated channel, cytokine / immune receptor, integrin, other ion channel), from the HGNC hierarchy closure.

1. Binding evidence (drug x gene, ChEMBL 37; parent and salt forms pooled):
     mechanism  curated mechanism (direct_interaction = 1) on a single protein / complex / complex group / chimeric
                protein containing the gene; PROTEIN FAMILY and SELECTIVITY GROUP mechanisms give no member credit
     potent     best valid pChEMBL >= pchembl_min (1 uM); kinase genes need pchembl_min_kinase (100 nM), also when
                mechanism-listed (kinase_mechanism_below_100nM / kinase_mechanism_unmeasured otherwise)
     weak / annotation_only / no_chembl_id -> no binding evidence
   Valid = no data_validity_comment, not a potential duplicate, human single-protein target.
2. Assembly: a gene whose mechanism target is a complex is judged as the complex (config assembly rules): tumors need
   each required subunit class expressed in >= 20% of tumors (aggregate approximation), lines need it per line.
3. Assay lines: the receptor (or assembly) is present (DepMap 24Q4 TPM >= 1) in >= half of the lines that 79e used for
   this cohort x drug, weighted as 79e weights them (direct: in-key lines, Celligner weights; borrowed: out-of-key
   lines, equal weights; all-zero weights -> equal). Same line rule as scripts/79e_drug_reversal.py (keep in sync).
receptor_gate: pass | binds_but_absent_in_matched_lines | assembly_incomplete_in_tumors | no_binding_evidence:<why>
single_line: pass resting on one matched line (kept, flagged; user 2026-11-02).
Selectivity: window_log10 = best pChEMBL over all of the drug's measured targets - pChEMBL on this gene (0 = the
drug's strongest target); drug_gap_1st_2nd = gap between the drug's two strongest distinct targets.

ChEMBL is fetched once into data/raw/chembl/ChEMBL_37/receptor_gate/ (refuses if the live API release differs from
config/datasets.yaml) and read from there afterwards; delete the folder to refetch.
Inputs: results/lincs/disease_drugs.tsv, stratum_drugs_v2.tsv, drug_annotations.tsv; data/processed/lincs/
        target_expression.tsv.gz, target_expression_strata_v2.tsv.gz, drug_line_profiles_index.tsv; results/celligner/
        disease_line_weights_long.tsv, stratum_line_weights_long_v2.tsv; config/celligner_tissue_keys.tsv;
        data/raw/depmap_24q4 (Model.csv, expression); data/raw/hgnc (complete set + genefamily tables)
Outputs: results/lincs/receptor_gate_v2.tsv (one row per cohort x drug x receptor), receptor_gate_summary_v2.tsv
         (per cohort), receptor_gate_selectivity_v2.tsv (window statistics), receptor_gate_run.json
"""
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
RAW, RES, PROC = ROOT / "data/raw", ROOT / "results/lincs", ROOT / "data/processed/lincs"
CFG = yaml.safe_load(open(ROOT / "config/drug_gate.yaml"))
RG, TPM_MIN, TUM_MIN = CFG["receptor_gate"], CFG["expression"]["tpm_min"], CFG["expression"]["min_fraction"]
DS = yaml.safe_load(open(ROOT / "config/datasets.yaml"))
RELEASE = DS["chembl"]["release"].split()[0]
CACHE = RAW / "chembl" / RELEASE / "receptor_gate"
GF = RAW / "hgnc/genefamily_20261102"
API = DS["chembl"]["url"]


# ---------------------------------------------------------------- HGNC classes
def hgnc():
    f = sorted((RAW / "hgnc").glob("hgnc_complete_set_*.txt"))[-1]
    return pd.read_csv(f, sep="\t", dtype=str, usecols=["hgnc_id", "symbol", "gene_group_id"])


def classes(hg):
    clo = pd.read_csv(GF / "hierarchy_closure.csv")
    g = hg.dropna(subset=["gene_group_id"]).assign(gid=lambda d: d.gene_group_id.str.split("|")).explode("gid")
    g["gid"] = g.gid.astype(int)
    cls = {}
    for c, roots in RG["receptor_families"].items():
        fams = set(clo[clo.parent_fam_id.isin(roots)].child_fam_id)
        for s in g[g.gid.isin(fams)].symbol:
            cls.setdefault(s, c)
    kin = set(g[g.gid.isin(set(clo[clo.parent_fam_id == RG["kinase_family"]].child_fam_id))].symbol)
    return cls, kin


# ---------------------------------------------------------------- rows: cohort x drug x expressed receptor
def receptor_rows(cls):
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    s2g = idm.set_index("symbol").gene_id
    out = []
    for f, ef, lvl in [("disease_drugs.tsv", "target_expression.tsv.gz", "disease"),
                       ("stratum_drugs_v2.tsv", "target_expression_strata_v2.tsv.gz", "subgroup")]:
        d = pd.read_csv(RES / f, sep="\t", dtype=str)
        d = d[d.passes_gate == "True"].dropna(subset=["targets_expressed"])
        t = d.assign(target=d.targets_expressed.str.split("|")).explode("target")
        t = t[t.target.isin(cls)].assign(receptor_class=lambda x: x.target.map(cls), level=lvl)
        t["gene_id"] = t.target.map(s2g)
        ex = pd.read_csv(PROC / ef, sep="\t")
        t = t.merge(ex, on=["gene_id", "group_id"], how="left")
        if "parent_group" not in t:
            t["parent_group"] = ""
        out.append(t[["level", "group_id", "parent_group", "pert_id", "cmap_name", "source", "n_targets", "hub_moa",
                      "target", "receptor_class", "median_tpm", "frac_expressed"]])
    r = pd.concat(out, ignore_index=True)
    r["n_targets"] = r.n_targets.astype(int)
    return r.drop_duplicates(["level", "group_id", "pert_id", "target"])


# ---------------------------------------------------------------- ChEMBL (cached)
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


def cached(name, fn):
    f = CACHE / name
    if f.exists():
        return json.load(open(f))
    live = get(f"{API}/status.json")["chembl_db_version"]
    assert live == RELEASE, f"ChEMBL API is {live}, config pins {RELEASE}: update config or keep the cache"
    CACHE.mkdir(parents=True, exist_ok=True)
    x = fn()
    json.dump(x, open(f, "w"))
    return x


def chembl(drugs):
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t")
    ann = ann[ann.cmap_name.isin(drugs)]
    missing = sorted(set(ann.cmap_name) - set(ann.dropna(subset=["chembl_id"]).cmap_name))

    def name_lookup():   # drugs without an InChIKey match in 79a: exact preferred-name lookup
        out = {}
        for n in missing:
            m = get(f"{API}/molecule.json?" + urllib.parse.urlencode({"pref_name__iexact": n.replace("-", " "), "limit": 5}))
            if m["molecules"]:
                out[n] = m["molecules"][0]["molecule_chembl_id"]
        return out
    byname = cached("name_lookup.json", name_lookup)
    ids = sorted(set(ann.chembl_id.dropna()) | set(byname.values()))
    forms = cached("forms.json", lambda: sum([paged("molecule_form", "molecule_forms",
                   {"molecule_chembl_id__in": ",".join(c)}) for c in chunks(ids, 40)], []))
    parents = sorted(set(ids) | {f["parent_chembl_id"] for f in forms})
    mech = cached("mechanisms.json", lambda: sum([paged("mechanism", "mechanisms",
                  {"molecule_chembl_id__in": ",".join(c)}) for c in chunks(parents, 40)], []))
    act = cached("activities.json", lambda: sum([paged("activity", "activities", {
        "molecule_chembl_id__in": ",".join(c), "pchembl_value__isnull": "false", "target_type": "SINGLE PROTEIN",
        "target_organism": "Homo sapiens",
        "only": "molecule_chembl_id,target_chembl_id,pchembl_value,standard_type,data_validity_comment,potential_duplicate"})
        for c in chunks(parents, 20)], []))
    tids = sorted({m["target_chembl_id"] for m in mech if m["target_chembl_id"]} | {a["target_chembl_id"] for a in act})
    tg = cached("targets.json", lambda: sum([paged("target", "targets", {"target_chembl_id__in": ",".join(c),
                "only": "target_chembl_id,target_type,pref_name,target_components,organism"}) for c in chunks(tids, 50)], []))
    # drug name -> molecule ids (salt + parent)
    par = {f["molecule_chembl_id"]: f["parent_chembl_id"] for f in forms}
    d2m = ann.dropna(subset=["chembl_id"]).groupby("cmap_name").chembl_id.apply(set).to_dict()
    for n, m in byname.items():
        d2m.setdefault(n, set()).add(m)
    d2m = {n: s | {par.get(m, m) for m in s} for n, s in d2m.items()}
    if (CACHE / "d113_ids.json").exists():   # D113: FDA drugs whose ChEMBL id was resolved from the FDA ingredient name
        for n, v in json.load(open(CACHE / "d113_ids.json")).items():
            d2m[n] = (d2m.get(n, set()) if v["mode"] == "add" else set()) | set(v["ids"])
    if (CACHE / "activities_by_parent.json").exists():   # D113: activities of those resolved ids
        act = act + json.load(open(CACHE / "activities_by_parent.json"))
    T = {}
    for t in tg:
        syms = {s["component_synonym"] for c in t.get("target_components") or [] for s in c.get("target_component_synonyms") or []
                if s["syn_type"] == "GENE_SYMBOL"}
        T[t["target_chembl_id"]] = (t["target_type"], t["pref_name"], frozenset(syms))
    credit = set(RG["binding"]["mechanism_target_types"])
    mg, munit = {}, {}
    for m in mech:
        if m["direct_interaction"] != 1 or not m["target_chembl_id"]:
            continue
        tt, nm, gs = T.get(m["target_chembl_id"], ("", "", frozenset()))
        if tt not in credit:
            continue
        for mol in {m["molecule_chembl_id"], m["parent_molecule_chembl_id"]}:
            mg.setdefault(mol, set()).update(gs)
            if tt in ("PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP"):
                for gene in gs:
                    munit[(mol, gene)] = (tt, nm, gs)
    a = pd.DataFrame(act)
    a = a[a.data_validity_comment.isna() & (a.potential_duplicate.astype(str) != "1")].copy()
    a["pchembl_value"] = a.pchembl_value.astype(float)
    a["gene"] = a.target_chembl_id.map(lambda t: T.get(t, ("", "", frozenset()))[2])
    a = a[a.gene.map(len) == 1].assign(gene=lambda d: d.gene.map(lambda s: next(iter(s))))
    best = {}
    for (mol, gene), v in a.groupby(["molecule_chembl_id", "gene"]).pchembl_value.max().items():
        best.setdefault(mol, {})[gene] = v
    n = dict(drugs=len(set(ann.cmap_name)), chembl_ids=len(ids), by_name=len(byname), mechanisms=len(mech),
             activities=len(act), activities_valid=len(a), targets=len(tg))
    return d2m, mg, munit, best, n


# ---------------------------------------------------------------- assembly rules
def unit_rule(unit):
    tt, nm, gs = unit
    if tt == "PROTEIN COMPLEX GROUP":
        for pre, req in RG["assembly"]["groups"].items():
            if nm.startswith(pre):
                if req == "any":
                    return [sorted(gs)]
                return [[g for g in gs if any(g.startswith(p) for p in cl)] for cl in req]
    if any(nm.startswith(p) for p in RG["assembly"]["any_component"]):
        return [sorted(gs)]
    return [[g] for g in sorted(gs)]


def main():
    hg = hgnc()
    cls, kinases = classes(hg)
    rows = receptor_rows(cls)
    d2m, mg, munit, best, nchem = chembl(sorted(rows.cmap_name.unique()))
    cut, cut_k = RG["binding"]["pchembl_min"], RG["binding"]["pchembl_min_kinase"]

    # ---- binding evidence per drug x gene
    ev = {}
    for dn, tgt in rows[["cmap_name", "target"]].drop_duplicates().itertuples(index=False):
        ms = d2m.get(dn, set())
        kin = tgt in kinases
        c = cut_k if kin else cut
        is_mech = any(tgt in mg.get(m, ()) for m in ms)
        pv = [best[m][tgt] for m in ms if tgt in best.get(m, {})]
        p = max(pv) if pv else np.nan
        allg = {}
        for m in ms:
            for gene, v in best.get(m, {}).items():
                allg[gene] = max(allg.get(gene, -1.0), v)
        vals = sorted(allg.values(), reverse=True)
        if is_mech and (not kin or p >= c):
            e = "mechanism"
        elif p >= c:
            e = "potent"
        elif is_mech and kin:
            e = "kinase_mechanism_below_100nM" if p == p else "kinase_mechanism_unmeasured"
        elif p == p:
            e = "weak"
        else:
            e = "annotation_only" if ms else "no_chembl_id"
        unit = next((munit[(m, tgt)] for m in ms if (m, tgt) in munit), None)
        ev[(dn, tgt)] = dict(kinase_target=kin, potency_cut=c, binding_evidence=e, pchembl_target=p,
                             pchembl_best_any=vals[0] if vals else np.nan,
                             window_log10=(vals[0] - p) if p == p else np.nan, n_chembl_targets=len(allg),
                             drug_gap_1st_2nd=(vals[0] - vals[1]) if len(vals) > 1 else np.nan,
                             assembly=unit[1] if unit else "", _unit=unit)

    # ---- line expression (DepMap) for receptors and assembly members
    genes = set(rows.target) | {g for e in ev.values() if e["_unit"] for g in e["_unit"][2]}
    f = RAW / "depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    hdr = pd.read_csv(f, nrows=0).columns
    X = pd.read_csv(f, usecols=[hdr[0]] + [c for c in hdr[1:] if c.split(" (")[0] in genes], index_col=0)
    X.columns = [c.split(" (")[0] for c in X.columns]
    X = X.loc[:, ~X.columns.duplicated()] >= np.log2(TPM_MIN + 1)
    # tumor fractions for assembly members
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    TE = pd.concat([pd.read_csv(PROC / "target_expression.tsv.gz", sep="\t"),
                    pd.read_csv(PROC / "target_expression_strata_v2.tsv.gz", sep="\t")])
    TE["symbol"] = TE.gene_id.map(dict(zip(idm.gene_id, idm.symbol)))
    tfrac = TE.dropna(subset=["symbol"]).drop_duplicates(["group_id", "symbol"]).set_index(["group_id", "symbol"]).frac_expressed

    # ---- lines used per cohort x drug (79e rule)
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    by_drug = {p: d.ModelID.values for p, d in idx.groupby("pert_id")}
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    W = pd.concat([pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")[["group_id", "ModelID", "weight"]],
                   pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t")
                   .rename(columns={"stratum": "group_id"})[["group_id", "ModelID", "weight"]]])
    wmap = {g: d.set_index("ModelID").weight for g, d in W.groupby("group_id")}
    keyset = {}

    def key_models(g, parent):
        P = parent or g
        if P not in keyset:
            k = keys.loc[P]
            km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
            if isinstance(k.diseases, str) and k.diseases:
                km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
            keyset[P] = set(km.index)
        return keyset[P]

    def present(unit, tgt, lines):
        if unit is None:
            return X.loc[lines, tgt] if tgt in X.columns else None
        ok = pd.Series(True, index=lines)
        for cl in unit_rule(unit):
            cl = [g for g in cl if g in X.columns]
            if not cl:
                return None
            ok &= X.loc[lines, cl].any(axis=1)
        return ok

    def tumor_ok(unit, g):
        if unit is None:
            return True
        for cl in unit_rule(unit):
            v = [tfrac.get((g, x)) for x in cl]
            v = [x for x in v if x is not None and x == x]
            if not v:
                return None
            if max(v) < TUM_MIN:
                return False
        return True

    out = []
    for r in rows.itertuples(index=False):
        e = ev[(r.cmap_name, r.target)]
        km = key_models(r.group_id, r.parent_group if isinstance(r.parent_group, str) else "")
        L = by_drug.get(r.pert_id, np.array([]))
        if r.source == "direct":
            L = [m for m in L if m in km]
            ww = wmap.get(r.group_id, pd.Series(dtype=float)).reindex(L).fillna(0).values
        else:
            L = [m for m in L if m not in km]
            ww = np.zeros(len(L))
        ww = ww if ww.sum() > 0 else np.ones(len(L))
        w = pd.Series(ww, index=L)
        w = w[w.index.isin(X.index)]
        w = w[~w.index.duplicated()]
        pr = present(e["_unit"], r.target, list(w.index)) if len(w) else None
        frac = float((w * pr.astype(float)).sum() / w.sum()) if pr is not None and w.sum() > 0 else np.nan
        tum = tumor_ok(e["_unit"], r.group_id)
        if e["binding_evidence"] not in ("mechanism", "potent"):
            gate = "no_binding_evidence:" + e["binding_evidence"]
        elif tum is False:
            gate = "assembly_incomplete_in_tumors"
        elif not frac >= RG["assay_lines"]["min_weighted_fraction"]:
            gate = "binds_but_absent_in_matched_lines"
        else:
            gate = "pass"
        out.append({**r._asdict(), **{k: v for k, v in e.items() if not k.startswith("_")},
                    "n_matched_lines": len(w), "frac_matched_lines_expr": round(frac, 3) if frac == frac else np.nan,
                    "assembly_in_tumors": tum, "receptor_gate": gate, "single_line": gate == "pass" and len(w) == 1})
    G = pd.DataFrame(out)
    for c in ["pchembl_target", "pchembl_best_any", "window_log10", "drug_gap_1st_2nd", "median_tpm"]:
        G[c] = G[c].round(2)
    G.to_csv(RES / "receptor_gate_v2.tsv", sep="\t", index=False)

    # ---- per cohort summary
    def summ(x):
        p = x[x.receptor_gate == "pass"]
        return pd.Series(dict(rows=len(x), drugs=x.pert_id.nunique(), receptors=x.target.nunique(), pass_rows=len(p),
                              pass_drugs=p.pert_id.nunique(), pass_receptors=p.target.nunique(),
                              single_line=int(p.single_line.sum()),
                              no_binding=int(x.receptor_gate.str.startswith("no_binding").sum()),
                              absent_in_lines=int((x.receptor_gate == "binds_but_absent_in_matched_lines").sum()),
                              assembly_incomplete=int((x.receptor_gate == "assembly_incomplete_in_tumors").sum()),
                              pass_receptor_list="|".join(sorted(p.target.unique()))))
    S = G.groupby(["level", "group_id"]).apply(summ).reset_index()
    S.to_csv(RES / "receptor_gate_summary_v2.tsv", sep="\t", index=False)

    # ---- selectivity window statistics (unique drug x receptor among passing rows)
    P = G[G.receptor_gate == "pass"].drop_duplicates(["cmap_name", "target"])
    def ws(x, set_, val):
        w = x.window_log10.dropna()
        return dict(set=set_, value=val, n=len(w), n_no_pchembl=int(x.window_log10.isna().sum()),
                    median=round(w.median(), 2) if len(w) else np.nan, q25=round(w.quantile(.25), 2) if len(w) else np.nan,
                    q75=round(w.quantile(.75), 2) if len(w) else np.nan, frac_best=round((w <= 0).mean(), 3) if len(w) else np.nan,
                    frac_lt10x=round((w < 1).mean(), 3) if len(w) else np.nan,
                    frac_10_100x=round(((w >= 1) & (w < 2)).mean(), 3) if len(w) else np.nan,
                    frac_ge100x=round((w >= 2).mean(), 3) if len(w) else np.nan)
    st = [ws(P, "all", "")]
    for col in ["binding_evidence", "kinase_target", "receptor_class"]:
        st += [ws(x, col, str(k)) for k, x in P.groupby(col)]
    pd.DataFrame(st).to_csv(RES / "receptor_gate_selectivity_v2.tsv", sep="\t", index=False)

    md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()
    run = dict(date=time.strftime("%Y-%m-%d"), config=RG, chembl_release=RELEASE, chembl=nchem,
               cache_md5={f.name: md5(f) for f in sorted(CACHE.glob("*.json"))},
               hgnc_genefamily_md5={f.name: md5(f) for f in sorted(GF.glob("*.csv"))},
               rows=len(G), gate=G.receptor_gate.value_counts().to_dict(),
               binding_evidence=G.drop_duplicates(["cmap_name", "target"]).binding_evidence.value_counts().to_dict(),
               single_line=int(G.single_line.sum()))
    json.dump(run, open(RES / "receptor_gate_run.json", "w"), indent=1, default=str)
    print(json.dumps({k: run[k] for k in ["chembl", "rows", "gate", "binding_evidence", "single_line"]}, indent=1))
    print(pd.DataFrame(st).to_string(index=False))


if __name__ == "__main__":
    main()
