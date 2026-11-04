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
Idempotent: reruns start from the *_pregate columns. 79b / 79w are read, never edited (other agent's scripts).
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, PROC = ROOT / "results/lincs", ROOT / "data/processed/lincs"
CARRIERS = RES / "d102/carriers_d102b.tsv.gz"
spec = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)


def evidence_targets():
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


def presence():
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    sym = dict(zip(idm.gene_id, idm.symbol))
    out = {}
    g = pd.read_csv(PROC / "target_expression.tsv.gz", sep="\t")
    g["symbol"] = g.gene_id.map(sym); g = g.dropna(subset=["symbol"]).drop_duplicates(["group_id", "symbol"])
    R = g.groupby("symbol").median_tpm.quantile(0.9)
    car = pd.read_csv(CARRIERS, sep="\t") if CARRIERS.exists() else pd.DataFrame(columns=["gene_id", "group_id", "n_tumors", "n_carriers"])
    g = g.merge(car, on=["gene_id", "group_id"], how="left")
    s = pd.read_csv(PROC / "target_expression_strata_v2.tsv.gz", sep="\t")
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
    df["target_carriers"], df["target_carrier_share"] = ncar, np.round(share, 4)
    for c, tier in pct_cols:
        base = c.replace("_pct", "")
        for g, x in df.groupby(unit):
            key = (~x.target_absent).astype(float) * 1e6 + x[base].rank(method="max")
            df.loc[x.index, c] = (key.rank(method="max", pct=True) * 100).round(1).values
        df.loc[df.target_absent, tier] = 3
    df.to_csv(f, sep="\t", index=False)
    return df


def main():
    E = evidence_targets(); P = presence()
    g = apply(RES / "drug_reversal_split_v2_consensus.tsv", "group_id", P["group"], E, [("spec_mech_pct", "tier_spec_mech")])
    s = apply(RES / "stratum_reversal_split_v2_consensus.tsv", "stratum", P["stratum"], E,
              [("spec_mech_pct", "tier_spec_mech"), ("spec_sib_mech_pct", "tier_spec_sib_mech")])
    run = dict(rule="D102a / D102c", drugs_with_evidence_targets=len(E),
               groups_absent_rows=int(g.target_absent.sum()), subgroup_absent_rows=int(s.target_absent.sum()),
               groups_minority_rows=int((g.target_absent & (g.target_carriers >= 3)).sum()))
    json.dump(run, open(RES / "target_presence_run.json", "w"), indent=1)
    print(run)


if __name__ == "__main__":
    main()
