#!/usr/bin/env python
"""Relative target-presence step in the ranking (D102a rule, adopted as D102c by the user 2026-11-04).

Runs after 79i on the production split tables, in place:
  results/lincs/drug_reversal_split_v2_consensus.tsv     (41 groups; presence from the group's tumors)
  results/lincs/stratum_reversal_split_v2_consensus.tsv  (subgroups; presence from the subgroup's own tumors)
Rule (docs/DECISIONS.md D102, D102a):
  judged drugs   >= 1 annotated receptor target (the 79w drug set; ChEMBL 37 cache of 79w)
  evidence E     curated direct ChEMBL mechanism genes (any target type) among the drug's annotated non-ADME targets;
                 if none, annotated targets with pChEMBL >= 6 (kinases >= 7)
  present        median TPM >= 1 and >= 10% of R, R = 90th percentile of the gene's 41 group medians
                 (genes missing from the expression table count as present)
  target_absent  no gene of E present -> ranked below every other drug of the group (spec_mech_pct, and
                 spec_sib_mech_pct for subgroups, recomputed); tier_spec_mech (tier_spec_sib_mech) = 3
Reported: evidence_targets (gene:median/R), target_carrier_share (41 groups: share of tumors at TPM >= max(10, R),
from the D102b carrier table; shown on the page as "Minority target" when >= 3 tumors carry it).
Also refreshes the 79i summaries that use these percentiles: split_controls{SFX}.tsv (spec_mech, spec_sib_mech) and
split_overview{SFX}.tsv (drugs_top10_spec_mech_in_ge20_groups, controls_*_spec_mech_pct_mean).
Env SFX (default _v2_consensus) selects the tables, e.g. _v2_consensus_d100early for the D100 time-split rankings.
Idempotent: reruns start from the *_pregate columns. 79b / 79w are read, never edited (other agent's scripts).
"""
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, PROC = ROOT / "results/lincs", ROOT / "data/processed/lincs"
CARRIERS = RES / "d102/carriers_d102b.tsv.gz"
SFX = os.environ.get("SFX", "_v2_consensus")
TARGETS = os.environ.get("TARGETS", "")   # "d105": ChEMBL primary targets for every ranked drug (D105)
SOURCE, SECONDARY = {}, {}
TOP_PCT = 10   # as 79i
spec = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)


def evidence_targets_d105():
    """D105: primary targets = curated mechanism genes; else most potent target + those within 10-fold (cut 6, kinases 7);
    else annotated targets (unverified_binding). Secondary = other potent targets and other annotated targets."""
    import yaml
    hg = w.hgnc(); _, kinases = w.classes(hg)
    cfg = yaml.safe_load(open(ROOT / "config/drug_gate.yaml"))["adme"]
    gg = hg.dropna(subset=["gene_group_id"]).assign(gid=lambda d: d.gene_group_id.str.split("|")).explode("gid")
    adme = set(gg[gg.gid.isin({str(i) for i in cfg["hgnc_group_ids"]})].symbol) | set(cfg["genes"])
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", dtype=str)
    ann = {}
    for r in dd[["cmap_name", "targets_expressed", "targets_not_expressed"]].itertuples(index=False):
        for c in r[1:]:
            if isinstance(c, str):
                ann.setdefault(r[0], set()).update(c.split("|"))
    drugs = sorted(set(pd.read_csv(RES / f"drug_reversal_split{SFX}.tsv", sep="\t", usecols=["cmap_name"]).cmap_name.dropna()))
    w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
    d2m, _, _, best, _ = w.chembl(drugs)
    mech = json.load(open(w.CACHE / "mechanisms.json"))
    tg = {t["target_chembl_id"]: {x["component_synonym"] for c in t.get("target_components") or []
                                  for x in c.get("target_component_synonyms") or [] if x["syn_type"] == "GENE_SYMBOL"}
          for t in json.load(open(w.CACHE / "targets.json")) if t.get("organism") == "Homo sapiens"}
    m_by = {}
    for m in mech:
        if m["direct_interaction"] == 1 and m["target_chembl_id"] in tg:
            for mol in {m["molecule_chembl_id"], m["parent_molecule_chembl_id"]}:
                m_by.setdefault(mol, set()).update(tg[m["target_chembl_id"]])
    cut, cut_k = w.RG["binding"]["pchembl_min"], w.RG["binding"]["pchembl_min_kinase"]
    E = {}
    for d in drugs:
        ms = d2m.get(d, set())
        pv = {}
        for m in ms:
            for gene, v in best.get(m, {}).items():
                pv[gene] = max(pv.get(gene, -1.0), v)
        elig = {gene: v for gene, v in pv.items() if gene not in adme and v >= (cut_k if gene in kinases else cut)}
        mg = set().union(*(m_by.get(m, set()) for m in ms)) - adme if ms else set()
        an = ann.get(d, set()) - adme
        if mg:
            prim, src = mg, "chembl_mechanism"
        elif elig:
            top = max(elig.values()); prim, src = {gene for gene, v in elig.items() if v >= top - 1}, "chembl_potency"
        elif an:
            prim, src = an, "annotation"
        else:
            continue
        E[d] = prim; SOURCE[d] = src
        sec = sorted(((gene, v) for gene, v in elig.items() if gene not in prim), key=lambda x: -x[1])
        SECONDARY[d] = ";".join([f"{gene}:{v:.1f}" for gene, v in sec[:8]] + [f"{gene}:ann" for gene in sorted(an - prim - set(elig))][:8])
    return E


def evidence_targets():
    if TARGETS == "d105":
        return evidence_targets_d105()
    hg = w.hgnc(); cls, kinases = w.classes(hg)
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", dtype=str)
    ann = {}
    for r in dd[["cmap_name", "targets_expressed", "targets_not_expressed"]].itertuples(index=False):
        for c in r[1:]:
            if isinstance(c, str):
                ann.setdefault(r[0], set()).update(c.split("|"))
    judged = {d for d, ts in ann.items() if ts & set(cls)}
    d2m, _, _, best, _ = w.chembl(sorted(judged))
    mech = json.load(open(w.CACHE / "mechanisms.json"))
    tg = {t["target_chembl_id"]: {s["component_synonym"] for c in t.get("target_components") or []
                                  for s in c.get("target_component_synonyms") or [] if s["syn_type"] == "GENE_SYMBOL"}
          for t in json.load(open(w.CACHE / "targets.json"))}
    cut, cut_k = w.RG["binding"]["pchembl_min"], w.RG["binding"]["pchembl_min_kinase"]
    E = {}
    for d in judged:
        ms = d2m.get(d, set())
        mg = set()
        for m in mech:
            if m["direct_interaction"] == 1 and m["target_chembl_id"] and ({m["molecule_chembl_id"], m["parent_molecule_chembl_id"]} & ms):
                mg |= tg.get(m["target_chembl_id"], set())
        ev = mg & ann[d]
        if not ev:
            for t in ann[d]:
                pv = [best[m][t] for m in ms if t in best.get(m, {})]
                if pv and max(pv) >= (cut_k if t in kinases else cut):
                    ev.add(t)
        if ev:
            E[d] = ev
    return E


def extra_expression(sym, genes):
    """D105: group / subgroup medians (79b group_expression, read-only import) for target genes missing from the 79b
    tables; cached in data/processed/lincs/target_expression{,_strata_v2}_d105.tsv.gz"""
    s2g = {v: k for k, v in sym.items()}
    sp_ = importlib.util.spec_from_file_location("b79", ROOT / "scripts/79b_target_gate.py"); b = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(b)
    out = []
    for base, extra, strata in (("target_expression.tsv.gz", "target_expression_d105.tsv.gz", False),
                                ("target_expression_strata_v2.tsv.gz", "target_expression_strata_v2_d105.tsv.gz", True)):
        have = set(pd.read_csv(PROC / base, sep="\t", usecols=["gene_id"]).gene_id)
        fx = PROC / extra
        ex = pd.read_csv(fx, sep="\t") if fx.exists() else pd.DataFrame(columns=["gene_id", "group_id", "median_tpm", "frac_expressed"])
        need = {s2g[x] for x in genes if x in s2g} - have - set(ex.gene_id)
        if need:
            if strata:
                sp2 = pd.read_csv(ROOT / "results/step10_v2/stratum_patients.tsv", sep="\t")
                new = b.group_expression(need, None, 1.0, sp2.rename(columns={"group_id": "parent", "stratum": "group_id"})[["group_id", "sample"]])
            else:
                sig = pd.read_csv(ROOT / "results/step4/signature_summary.tsv", sep="\t")
                new = b.group_expression(need, sorted(set(sig.loc[sig.status == "ok", "group_id"])), 1.0)
            ex = pd.concat([ex, new], ignore_index=True); ex.to_csv(fx, sep="\t", index=False)
        out.append(ex)
    return out


def presence(genes=None):
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    sym = dict(zip(idm.gene_id, idm.symbol))
    out = {}
    g = pd.read_csv(PROC / "target_expression.tsv.gz", sep="\t")
    sx = pd.read_csv(PROC / "target_expression_strata_v2.tsv.gz", sep="\t")
    if TARGETS == "d105" and genes:
        gx, sxx = extra_expression(sym, genes)
        g, sx = pd.concat([g, gx], ignore_index=True), pd.concat([sx, sxx], ignore_index=True)
    g["symbol"] = g.gene_id.map(sym); g = g.dropna(subset=["symbol"]).drop_duplicates(["group_id", "symbol"])
    R = g.groupby("symbol").median_tpm.quantile(0.9)
    car = pd.read_csv(CARRIERS, sep="\t") if CARRIERS.exists() else pd.DataFrame(columns=["gene_id", "group_id", "n_tumors", "n_carriers"])
    g = g.merge(car, on=["gene_id", "group_id"], how="left")
    s = sx
    s["symbol"] = s.gene_id.map(sym); s = s.dropna(subset=["symbol"]).drop_duplicates(["group_id", "symbol"])
    for k, t in (("group", g), ("stratum", s)):
        t["R"] = t.symbol.map(R)
        t["present"] = (t.median_tpm >= 1) & (t.median_tpm >= 0.1 * t.R)
        out[k] = t.set_index(["group_id", "symbol"])
    return out


def apply(f, unit, P, E, pct_cols):
    df = pd.read_csv(f, sep="\t", low_memory=False)
    bak = f.with_name(f.stem + "_pre_d102.tsv")
    if not bak.exists():
        df.to_csv(bak, sep="\t", index=False)
    for c, _ in pct_cols:
        if c + "_pregate" not in df:
            df[c + "_pregate"] = df[c]
    pres, med, R = P.present.to_dict(), P.median_tpm.to_dict(), P.R.to_dict()
    nc = P.n_carriers.to_dict() if "n_carriers" in P else {}
    nt = P.n_tumors.to_dict() if "n_tumors" in P else {}
    absent, evs, share, ncar = [], [], [], []
    for g, d in zip(df[unit], df.cmap_name):
        e = E.get(d)
        if not e:
            absent.append(False); evs.append(""); share.append(np.nan); ncar.append(np.nan); continue
        absent.append(not any(pres.get((g, t), True) for t in e))
        evs.append(";".join(f"{t}:{med.get((g, t), np.nan):.2f}/{R.get((g, t), np.nan):.1f}" for t in sorted(e)))
        best = max(e, key=lambda t: nc.get((g, t), -1))
        n = nc.get((g, best), np.nan)
        ncar.append(n); share.append(n / nt[(g, best)] if n == n and nt.get((g, best)) else np.nan)
    df["target_absent"], df["evidence_targets"] = absent, evs
    if TARGETS == "d105":
        df["target_source"] = df.cmap_name.map(SOURCE).fillna("")
        df["unverified_binding"] = df.target_source.eq("annotation")
        df["secondary_targets"] = df.cmap_name.map(SECONDARY).fillna("")
    df["target_carriers"], df["target_carrier_share"] = ncar, np.round(share, 4)
    for c, tier in pct_cols:
        base = c.replace("_pct", "")
        df[c] = df[c + "_pregate"]
        for g, x in df.groupby(unit):
            key = (~x.target_absent).astype(float) * 1e6 + x[base].rank(method="max")
            df.loc[x.index, c] = (key.rank(method="max", pct=True) * 100).round(1).values
        df.loc[df.target_absent, tier] = 3
    df.to_csv(f, sep="\t", index=False)
    return df


def refresh_summaries(g, s):
    fc, fo = RES / f"split_controls{SFX}.tsv", RES / f"split_overview{SFX}.tsv"
    if not fc.exists():
        return
    ctl = pd.read_csv(fc, sep="\t")
    for i, r in ctl.iterrows():
        if r.unit in set(g.group_id):
            x = g[(g.group_id == r.unit) & (g.cmap_name == r.drug)]
            ctl.loc[i, "spec_mech"] = x.spec_mech_pct.max()
        else:
            x = s[(s.stratum == r.unit) & (s.cmap_name == r.drug)]
            ctl.loc[i, "spec_mech"] = x.spec_mech_pct.max()
            if "spec_sib_mech" in ctl:
                ctl.loc[i, "spec_sib_mech"] = x.spec_sib_mech_pct.max()
    ctl.to_csv(fc, sep="\t", index=False)
    if fo.exists():
        ov = pd.read_csv(fo, sep="\t").set_index("metric")
        ok = g[~g.low_mrs.astype(str).str.lower().eq("true")]
        top = ok[ok.spec_mech_pct >= 100 - TOP_PCT].groupby("pert_id").size()
        ov.loc["drugs_top10_spec_mech_in_ge20_groups", "value"] = int((top >= 20).sum())
        for exp in ["high", "low"]:
            ov.loc[f"controls_{exp}_spec_mech_pct_mean", "value"] = round(ctl[ctl.expected == exp].spec_mech.mean(), 1)
        ov.reset_index().to_csv(fo, sep="\t", index=False)


def main():
    E = evidence_targets(); P = presence(set().union(*E.values()))
    g = apply(RES / f"drug_reversal_split{SFX}.tsv", "group_id", P["group"], E, [("spec_mech_pct", "tier_spec_mech")])
    s = apply(RES / f"stratum_reversal_split{SFX}.tsv", "stratum", P["stratum"], E,
              [("spec_mech_pct", "tier_spec_mech"), ("spec_sib_mech_pct", "tier_spec_sib_mech")])
    refresh_summaries(g, s)
    run = dict(rule="D102a / D102c", tables=SFX, drugs_with_evidence_targets=len(E),
               groups_absent_rows=int(g.target_absent.sum()), subgroup_absent_rows=int(s.target_absent.sum()),
               groups_minority_rows=int((g.target_absent & (g.target_carriers >= 3)).sum()))
    json.dump(run, open(RES / f"target_presence_run{'' if SFX == '_v2_consensus' else SFX}.json", "w"), indent=1)
    print(run)


if __name__ == "__main__":
    main()
