#!/usr/bin/env python
"""D102 (rule fixed 2026-11-04, docs/DECISIONS.md): relative target-presence gate applied in the ranking.
STAGE=1 breast groups, STAGE=2 all 41. Reads production drug_reversal_split_v2_consensus.tsv; never edits 79b / 79w
(imports 79w helpers read-only; ChEMBL from its cache). Outputs results/lincs/d102/stage{N}_{flags,ranked,metrics,
controls,decision}.tsv"""
import importlib.util, os
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; RES = ROOT / "results/lincs"; PROC = ROOT / "data/processed/lincs"
OUT = RES / "d102"; OUT.mkdir(exist_ok=True); STAGE = int(os.environ.get("STAGE", "1")); AMEND = os.environ.get("AMEND", "")  # "a" = D102a
TAG = f"stage{STAGE}{AMEND}"
spec = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)
BREAST = ["BRCA-LumA", "BRCA-LumB", "BRCA-Her2", "BRCA-Basal"]

# ---- evidence targets per drug (79w binding rule, any target class)
hg = w.hgnc(); cls, kinases = w.classes(hg)
dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", dtype=str)
ann = {}
for r in dd[["cmap_name", "targets_expressed", "targets_not_expressed"]].itertuples(index=False):
    for c in r[1:]:
        if isinstance(c, str):
            ann.setdefault(r[0], set()).update(c.split("|"))
judged = {d for d, ts in ann.items() if ts & set(cls)}
d2m, mg, munit, best, _ = w.chembl(sorted(judged))
import json
MECH = json.load(open(w.CACHE / "mechanisms.json"))
TG = {t["target_chembl_id"]: {s_["component_synonym"] for c in t.get("target_components") or [] for s_ in c.get("target_component_synonyms") or [] if s_["syn_type"] == "GENE_SYMBOL"} for t in json.load(open(w.CACHE / "targets.json"))}
cut, cut_k = w.RG["binding"]["pchembl_min"], w.RG["binding"]["pchembl_min_kinase"]
E = {}
for d in judged:
    ms = d2m.get(d, set()); ev = set()
    for t in ann[d]:
        kin = t in kinases; c = cut_k if kin else cut
        mech = any(t in mg.get(m, ()) for m in ms)
        pv = [best[m][t] for m in ms if t in best.get(m, {})]; p = max(pv) if pv else np.nan
        if (mech and (not kin or p >= c)) or p >= c:
            ev.add(t)
    if AMEND == "a":   # D102a (2): curated direct mechanism genes (any target type) among the annotated targets
        mech_g = set()
        for m in MECH:
            if m["direct_interaction"] == 1 and m["target_chembl_id"] and ({m["molecule_chembl_id"], m["parent_molecule_chembl_id"]} & ms):
                mech_g |= TG.get(m["target_chembl_id"], set())
        mech_g &= ann[d]
        if mech_g:
            ev = mech_g
    if ev:
        E[d] = ev

# ---- presence: median TPM >= 1 and >= 10% of the gene's highest median over the 41 groups
idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
te = pd.read_csv(PROC / "target_expression.tsv.gz", sep="\t")
te["symbol"] = te.gene_id.map(dict(zip(idm.gene_id, idm.symbol))); te = te.dropna(subset=["symbol"]).drop_duplicates(["group_id", "symbol"])
te["gmax"] = te.groupby("symbol").median_tpm.transform(lambda v: v.quantile(0.9) if AMEND == "a" else v.max())
te["present"] = (te.median_tpm >= 1) & (te.median_tpm >= 0.1 * te.gmax)
P = te.set_index(["group_id", "symbol"])
pres = P.present.to_dict(); med = P.median_tpm.to_dict(); gmax = P.gmax.to_dict()

sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
groups = BREAST if STAGE == 1 else sorted(sp.group_id.unique())
sp = sp[sp.group_id.isin(groups)].copy()
def flag(g, d):
    if d not in E:
        return False, ""
    ok = [t for t in E[d] if pres.get((g, t), True)]
    return len(ok) == 0, ";".join(f"{t}:{med.get((g, t), np.nan):.2f}/{gmax.get((g, t), np.nan):.1f}" for t in sorted(E[d]))
f = [flag(g, d) for g, d in zip(sp.group_id, sp.cmap_name)]
sp["target_absent"] = [a for a, _ in f]; sp["evidence_targets_tpm"] = [b for _, b in f]
sp["spec_mech_pct_d102"] = np.nan
for g, x in sp.groupby("group_id"):
    key = (~x.target_absent).astype(float) * 1e6 + x.spec_mech.rank(method="max")
    sp.loc[x.index, "spec_mech_pct_d102"] = (key.rank(method="max", pct=True) * 100).round(1).values
sp[sp.target_absent][["group_id", "pert_id", "cmap_name", "source", "hub_moa", "spec_mech_pct", "spec_mech_pct_d102", "evidence_targets_tpm"]] \
    .sort_values(["group_id", "spec_mech_pct"], ascending=[True, False]).to_csv(OUT / f"{TAG}_flags.tsv", sep="\t", index=False)

# ---- metrics
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])
def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values; return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
PRAD = [("TCGA-PRAD", d) for d in ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]]
M, C, A = {}, [], []
auc_groups = ["BRCA-LumA", "BRCA-LumB", "BRCA-Her2"] if STAGE == 1 else groups
for arm, col in [("production", "spec_mech_pct"), ("d102", "spec_mech_pct_d102")]:
    best_ = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)][col].max()
    o = {k: float(np.nanmedian([best_(g, d) for g, d in v])) for k, v in [("ctl_pos", POS), ("ctl_neg", NEG)] + ([("ctl_prad", PRAD)] if STAGE == 2 else [])}
    C += [dict(arm=arm, group_id=g, drug=d, pct=best_(g, d)) for g, d in POS + NEG + (PRAD if STAGE == 2 else [])]
    s = sp.groupby(["group_id", "pert_id"])[col].max().reset_index()
    m = lab.merge(s, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        per = {c: auc(g[g[ls]][col].dropna().values, g[~g[ls]][col].dropna().values) for c, g in m.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20}
        A += [dict(arm=arm, label=ls, cohort=c, auc=v) for c, v in per.items()]
        o[f"auc_{ls}"] = float(np.median([per[c] for c in auc_groups if c in per]))
    M[arm] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / f"{TAG}_metrics.tsv", sep="\t", float_format="%.4f")
ct = pd.DataFrame(C).pivot_table(index=["group_id", "drug"], columns="arm", values="pct"); ct.to_csv(OUT / f"{TAG}_controls.tsv", sep="\t")
pd.DataFrame(A).pivot_table(index=["label", "cohort"], columns="arm", values="auc").to_csv(OUT / f"{TAG}_auc_by_cohort.tsv", sep="\t", float_format="%.4f")
labd = lab[lab.lab_hub | lab.lab_chembl].merge(sp[sp.target_absent][["group_id", "pert_id", "cmap_name", "evidence_targets_tpm"]],
                                                left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
labd.drop(columns="cohort").to_csv(OUT / f"{TAG}_labeled_demoted.tsv", sep="\t", index=False)
nd = sp.groupby("group_id").agg(drugs=("pert_id", "size"), demoted=("target_absent", "sum"))
nd.to_csv(OUT / f"{TAG}_demoted_per_group.tsv", sep="\t")
c, x = M["production"], M["d102"]
crit = [("B3 known-use AUC falls <= 0.01 (Hub, ChEMBL)", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}",
         x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01)]
if STAGE == 1:
    nl = labd[labd.group_id.isin(["BRCA-LumA", "BRCA-LumB"])]
    if AMEND == "a":   # D102a (3): only Hub (approved-use) labels stop; ChEMBL-only listed for review
        print("ChEMBL-only labeled demoted (review):", nl[~nl.lab_hub][["group_id", "cmap_name", "evidence_targets_tpm"]].to_string(index=False))
        nl = nl[nl.lab_hub]
    crit = [("B1 Basal negative controls fall >= 20", f"{x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}", x["ctl_neg"] <= c["ctl_neg"] - 20),
            ("B2 LumA / LumB positive controls not lower", f"{x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}", x["ctl_pos"] >= c["ctl_pos"])] + crit + \
           [("B4 no labeled drug demoted in LumA / LumB", f"{len(nl)} ({', '.join(sorted(set(nl.cmap_name)))})", len(nl) == 0)]
else:
    crit += [("PRAD controls not lower", f"{x['ctl_prad']:.1f} vs {c['ctl_prad']:.1f}", x["ctl_prad"] >= c["ctl_prad"])]
dec = pd.DataFrame(crit, columns=["criterion", "value", "met"])
dec.loc[len(dec)] = ("DECISION", ("go to stage 2" if STAGE == 1 else "adopt (needs the user's go-ahead)") if dec.met.all() else "stop", dec.met.all())
dec.to_csv(OUT / f"{TAG}_decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 90)
print("judged drugs", len(judged), "with evidence targets", len(E)); print(nd.to_string()); print(met.round(4).to_string()); print(ct.round(1).to_string())
print(dec.to_string(index=False))
