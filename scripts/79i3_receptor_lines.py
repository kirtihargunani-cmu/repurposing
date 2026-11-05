#!/usr/bin/env python
"""Receptor present in the lines the score came from (D104 / D104b rule, adopted as D104c by the user 2026-11-04, cut 0.4).

Runs after 79i2, in place on the production split tables (41 groups and subgroups; env SFX as 79i2).
Judged: cohort x drug whose tumor-present evidence targets (79i2: E and presence rule) are all receptors (79w classes).
Scoring lines and weights as 79e / 79w (direct: in-key lines, Celligner weights; borrowed: out-of-key, equal weights;
subgroups: own weights, parent group's tissue key). Per receptor:
  a line supports it if DepMap 24Q4 TPM >= 1 or the drug induces it there (raw consensus profile z >= 2);
  support = weighted fraction of supporting lines
  matched            support >= 0.5
  tumor_environment  support < 0.5 and max(rho with ESTIMATE, rho with the endothelial score PECAM1 / VWF / CDH5)
                     >= 0.4 across the group's tumors (subgroups: parent group's rho)
  line_absent        support < 0.5 otherwise; not_tested = no DepMap data for the gene / lines
Drug level: any matched -> ok; else any tumor_environment -> acts_on_tumor_environment (flag); else any not_tested ->
receptor_not_tested (flag); else receptor_mismatch -> ranked below unflagged drugs and above target-absent drugs, tier 3.
Percentiles are recomputed from the *_pregate columns written by 79i2 (idempotent).
rho table: data/processed/lincs/receptor_environment_rho.tsv.gz (genes missing from it are computed and appended).
"""
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES, PROC, RAW = ROOT / "results/lincs", ROOT / "data/processed/lincs", ROOT / "data/raw"
SFX = os.environ.get("SFX", "_v2_consensus")
RHO_MIN, SUPPORT_MIN, INDUCE_Z = 0.4, 0.5, 2.0
ENDO = ["PECAM1", "VWF", "CDH5"]
RHO_FILE = PROC / "receptor_environment_rho.tsv.gz"
spec = importlib.util.spec_from_file_location("i2", ROOT / "scripts/79i2_target_presence.py")
i2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i2)
w = i2.w


def rho_table(genes, groups):
    have = pd.read_csv(RHO_FILE, sep="\t") if RHO_FILE.exists() else pd.DataFrame(columns=["group_id", "target", "rho", "n", "rho_endo"])
    need = set(genes) - set(have.target)
    if need:
        idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
        s2g = dict(zip(idm.symbol, idm.gene_id)); gids = {s2g[g]: g for g in need | set(ENDO) if g in s2g}
        cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
        tum = cov[(cov.role == "tumor") & (cov.primary_per_patient == True) & cov.group_id.isin(groups)]
        tum = tum.drop_duplicates(["group_id", "sample"]).dropna(subset=["estimate_score"])
        ft = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"; hd = pd.read_csv(ft, sep="\t", nrows=0).columns
        use = sorted(set(tum["sample"]) & set(hd)); parts = []
        for ch in pd.read_csv(ft, sep="\t", usecols=[hd[0]] + use, chunksize=5000, index_col=0):
            ch.index = ch.index.str.split(".").str[0]; parts.append(ch[ch.index.isin(set(gids))])
        T = pd.concat(parts); T.index = T.index.map(gids); T = T[~T.index.duplicated()]
        out = []
        for g, dg in tum[tum["sample"].isin(use)].groupby("group_id"):
            es = dg.set_index("sample").estimate_score.rank()
            lt = np.log2(np.power(2.0, T.loc[[x for x in ENDO if x in T.index], es.index]) - 0.001 + 1)
            en = lt.sub(lt.mean(axis=1), axis=0).div(lt.std(axis=1), axis=0).mean().rank()
            Tg = T.loc[[x for x in T.index if x in need], es.index].rank(axis=1)
            out.append(pd.DataFrame({"group_id": g, "target": Tg.index, "n": len(es),
                                     "rho": Tg.apply(lambda v: np.corrcoef(v.values, es.values)[0, 1], axis=1).values.round(3),
                                     "rho_endo": Tg.apply(lambda v: np.corrcoef(v.values, en.values)[0, 1], axis=1).values.round(3)}))
        have = pd.concat([have] + out, ignore_index=True); have.to_csv(RHO_FILE, sep="\t", index=False)
    have = have[~have.target.isin(ENDO)]
    return {(g, t): max(a, b) for g, t, a, b in zip(have.group_id, have.target, have.rho.fillna(-1), have.rho_endo.fillna(-1))}


def main():
    E = i2.evidence_targets(); P = i2.presence(set().union(*E.values()))
    cls, _ = w.classes(w.hgnc()); REC = set(cls)
    genes = set().union(*E.values()) & REC
    g = pd.read_csv(RES / f"drug_reversal_split{SFX}.tsv", sep="\t", low_memory=False)
    s = pd.read_csv(RES / f"stratum_reversal_split{SFX}.tsv", sep="\t", low_memory=False)
    RHO = rho_table(genes, set(g.group_id))
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"); idx["row"] = np.arange(len(idx))
    rowof = {(p, m): r for p, m, r in zip(idx.pert_id, idx.ModelID, idx.row)}
    by_drug = {p: d.ModelID.values for p, d in idx.groupby("pert_id")}
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0); models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    W = pd.concat([pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")[["group_id", "ModelID", "weight"]],
                   pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})[["group_id", "ModelID", "weight"]]])
    wmap = {u: d.set_index("ModelID").weight for u, d in W.groupby("group_id")}
    f = RAW / "depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"; hd = pd.read_csv(f, nrows=0).columns
    X = pd.read_csv(f, usecols=[hd[0]] + [c for c in hd[1:] if c.split(" (")[0] in genes], index_col=0)
    X.columns = [c.split(" (")[0] for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]; X = np.power(2.0, X) - 1
    G = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t"); gcol = {x: i for i, x in enumerate(G.symbol)}
    PR = np.memmap(PROC / "drug_line_profiles_bing.f32", dtype=np.float32, mode="r").reshape(len(idx), len(G))
    ks = {}

    def km(key_group):
        if key_group not in ks:
            k = keys.loc[key_group]; m = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
            if isinstance(k.diseases, str) and k.diseases:
                m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
            ks[key_group] = set(m.index)
        return ks[key_group]

    def judge(df, unit, parent, PP):
        pres = PP.present.to_dict(); pres_fda = PP.present_fda.to_dict() if "present_fda" in PP else pres
        status, detail = [], []
        for u, par, pert, d, src in zip(df[unit], df[parent], df.pert_id, df.cmap_name, df.source):
            e = E.get(d)
            pm = pres_fda if d in i2.FDA_NAMES else pres
            tp = [t for t in e if pm.get((u, t), True)] if e else []
            if not tp or not set(tp) <= REC:
                status.append(""); detail.append(""); continue
            L = by_drug.get(pert, np.array([]))
            if src == "direct":
                L = [m for m in L if m in km(par)]; ww = wmap.get(u, pd.Series(dtype=float)).reindex(L).fillna(0).values
            else:
                L = [m for m in L if m not in km(par)]; ww = np.zeros(len(L))
            ww = ww if ww.sum() > 0 else np.ones(len(L))
            sw = pd.Series(ww, index=L); sw = sw[~sw.index.duplicated()]
            lx = [m for m in sw.index if m in X.index]
            st = {}
            for t in tp:
                if t not in X.columns or not lx:
                    st[t] = "not_tested"; continue
                tpm = X.loc[lx, t]
                z = pd.Series([PR[rowof[(pert, m)], gcol[t]] if t in gcol and (pert, m) in rowof else np.nan for m in lx], index=lx)
                sup = float((sw[lx] * ((tpm >= 1) | (z >= INDUCE_Z))).sum() / sw[lx].sum())
                st[t] = "matched" if sup >= SUPPORT_MIN else ("tumor_environment" if RHO.get((par, t), -1) >= RHO_MIN else "line_absent")
            v = set(st.values())
            status.append("ok" if "matched" in v else "acts_on_tumor_environment" if "tumor_environment" in v
                          else "receptor_not_tested" if "not_tested" in v else "receptor_mismatch")
            detail.append(";".join(f"{t}:{x}" for t, x in sorted(st.items())))
        df["receptor_status"], df["receptor_detail"] = status, detail
        return df

    def rerank(df, unit, cols):
        mis = df.receptor_status.eq("receptor_mismatch"); ab = df.target_absent.astype(str).str.lower().eq("true")
        for c, tier in cols:
            for u, x in df.groupby(unit):
                key = (~ab[x.index]).astype(float) * 2e6 + (~mis[x.index]).astype(float) * 1e6 + x[c.replace("_pct", "")].rank(method="max")
                df.loc[x.index, c] = (key.rank(method="max", pct=True) * 100).round(1).values
            df.loc[mis | ab, tier] = 3
        return df

    g["_par"] = g.group_id
    g = rerank(judge(g, "group_id", "_par", P["group"]), "group_id", [("spec_mech_pct", "tier_spec_mech")]).drop(columns="_par")
    s = rerank(judge(s, "stratum", "parent_group", P["stratum"]), "stratum",
               [("spec_mech_pct", "tier_spec_mech"), ("spec_sib_mech_pct", "tier_spec_sib_mech")])
    g.to_csv(RES / f"drug_reversal_split{SFX}.tsv", sep="\t", index=False)
    s.to_csv(RES / f"stratum_reversal_split{SFX}.tsv", sep="\t", index=False)
    i2.refresh_summaries(g, s)
    run = dict(rule="D104b / D104c (rho 0.4)", tables=SFX, groups=g.receptor_status.value_counts().to_dict(),
               subgroups=s.receptor_status.value_counts().to_dict())
    json.dump(run, open(RES / f"receptor_lines_run{'' if SFX == '_v2_consensus' else SFX}.json", "w"), indent=1)
    print(run)


if __name__ == "__main__":
    main()
