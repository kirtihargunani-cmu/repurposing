#!/usr/bin/env python
"""D104 evaluation (rule fixed 2026-11-04, docs/DECISIONS.md): receptor level in tumor vs scoring lines, all pairs,
41 groups, on top of production (D102c). Reads 79i2 / 79w helpers; writes only results/lincs/d104/eval_*"""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import os
AMEND = os.environ.get("AMEND", ""); PFX = "evala_" if AMEND == "a" else "eval_"
ROOT = Path(__file__).resolve().parent.parent.parent; RES = ROOT / "results/lincs"; OUT = RES / "d104"; RAW = ROOT / "data/raw"; PROC = ROOT / "data/processed/lincs"
def load(name, f):
    s = importlib.util.spec_from_file_location(name, ROOT / f); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
i2 = load("i2", "scripts/79i2_target_presence.py"); w = i2.w
E = i2.evidence_targets(); P = i2.presence()["group"]
pres, med = P.present.to_dict(), P.median_tpm.to_dict()
cls, _ = w.classes(w.hgnc()); REC = set(cls)
sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
# lines, weights (79e / 79w rule)
idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"); idx["row"] = np.arange(len(idx))
rowof = {(p, m): r for p, m, r in zip(idx.pert_id, idx.ModelID, idx.row)}
by_drug = {p: d.ModelID.values for p, d in idx.groupby("pert_id")}
keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0); models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
W = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t"); wmap = {g: d.set_index("ModelID").weight for g, d in W.groupby("group_id")}
ks = {}
def km(g):
    if g not in ks:
        k = keys.loc[g]; m = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        ks[g] = set(m.index)
    return ks[g]
genes_all = set().union(*E.values()) & REC
f = RAW / "depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"; hdr = pd.read_csv(f, nrows=0).columns
X = pd.read_csv(f, usecols=[hdr[0]] + [c for c in hdr[1:] if c.split(" (")[0] in genes_all], index_col=0)
X.columns = [c.split(" (")[0] for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]; X = np.power(2.0, X) - 1
RHO = {}
if AMEND == "a":   # D104a (2): Spearman of receptor TPM with ESTIMATE score across each group's tumors
    rf = OUT / "estimate_rho.tsv.gz"
    if not rf.exists():
        idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
        s2g = dict(zip(idm.symbol, idm.gene_id)); gids = {s2g[g]: g for g in genes_all if g in s2g}
        cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
        tum = cov[(cov.role == "tumor") & (cov.primary_per_patient == True) & cov.group_id.isin(set(sp.group_id))].drop_duplicates(["group_id", "sample"]).dropna(subset=["estimate_score"])
        ft = ROOT / "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"; hd = pd.read_csv(ft, sep="\t", nrows=0).columns
        use = sorted(set(tum["sample"]) & set(hd)); parts = []
        for ch in pd.read_csv(ft, sep="\t", usecols=[hd[0]] + use, chunksize=5000, index_col=0):
            ch.index = ch.index.str.split(".").str[0]; parts.append(ch[ch.index.isin(set(gids))])
        T = pd.concat(parts); T.index = T.index.map(gids); T = T[~T.index.duplicated()]
        out = []
        for g, dg in tum[tum["sample"].isin(use)].groupby("group_id"):
            es = dg.set_index("sample").estimate_score.rank(); Tg = T[es.index].rank(axis=1)
            rho = Tg.apply(lambda v: np.corrcoef(v.values, es.values)[0, 1], axis=1)
            out.append(pd.DataFrame({"group_id": g, "target": rho.index, "rho": rho.values.round(3), "n": len(es)}))
        pd.concat(out).to_csv(rf, sep="\t", index=False)
    rr = pd.read_csv(rf, sep="\t"); RHO = {(g, t): v for g, t, v in zip(rr.group_id, rr.target, rr.rho)}
RHO_MIN = 0.3
G = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t"); gcol = {s: i for i, s in enumerate(G.symbol)}
PR = np.memmap(PROC / "drug_line_profiles_bing.f32", dtype=np.float32, mode="r").reshape(len(idx), len(G))

rows, recs = [], []
for x in sp[["group_id", "pert_id", "cmap_name", "source", "unmatched_lines"]].drop_duplicates(["group_id", "pert_id"]).itertuples(index=False):
    e = E.get(x.cmap_name)
    if not e:
        continue
    tp = [t for t in e if pres.get((x.group_id, t), True)]
    if not tp or not set(tp) <= REC:
        continue
    L = by_drug.get(x.pert_id, np.array([]))
    if x.source == "direct":
        L = [m for m in L if m in km(x.group_id)]; ww = wmap.get(x.group_id, pd.Series(dtype=float)).reindex(L).fillna(0).values
    else:
        L = [m for m in L if m not in km(x.group_id)]; ww = np.zeros(len(L))
    ww = ww if ww.sum() > 0 else np.ones(len(L))
    s = pd.Series(ww, index=L); s = s[~s.index.duplicated()]
    st = {}
    for t in tp:
        lx = [m for m in s.index if m in X.index]
        if t not in X.columns or not lx:
            st[t] = ("not_tested", np.nan, np.nan, False); continue
        ws = s[lx]; tpm = X.loc[lx, t]
        z = pd.Series([PR[rowof[(x.pert_id, m)], gcol[t]] if t in gcol and (x.pert_id, m) in rowof else np.nan for m in lx], index=lx)
        sup = (tpm >= 1) | (z >= 2)
        induced = bool(((tpm < 1) & (z >= 2)).any())
        frac = float((ws * sup).sum() / ws.sum()); lt = float((ws * tpm).sum() / ws.sum())
        tm = med.get((x.group_id, t), np.nan)
        rho = RHO.get((x.group_id, t), np.nan)
        if AMEND == "a":   # D104a (1) no line-relative abundance test, (2) tumor-environment exception
            stat = ("tumor_environment" if rho >= RHO_MIN else "line_absent") if frac < 0.5 else "matched"
        else:
            stat = "line_absent" if frac < 0.5 else ("tumor_far_lower" if tm < 0.1 * lt else "matched")
        st[t] = (stat, frac, lt, induced)
        recs.append(dict(group_id=x.group_id, pert_id=x.pert_id, cmap_name=x.cmap_name, target=t, tumor_tpm=tm, line_tpm=round(lt, 2),
                         support=round(frac, 3), rho_estimate=rho, induced=induced, in_lincs=t in gcol, status=stat, n_lines=len(lx)))
    S = [v[0] for v in st.values()]
    if "matched" in S:
        c = "ok"
    elif "tumor_environment" in S:
        c = "acts_on_tumor_environment"
    elif "not_tested" in S:
        c = "receptor_not_tested"
    else:
        c = "receptor_mismatch:" + ("line_absent" if all(v == "line_absent" for v in S) else ("tumor_far_lower" if all(v == "tumor_far_lower" for v in S) else "mixed"))
    poor = x.source == "borrowed" or str(x.unmatched_lines).lower() == "true"
    rows.append(dict(group_id=x.group_id, pert_id=x.pert_id, cmap_name=x.cmap_name, poor_match=poor, d104=c,
                     receptors=";".join(f"{t}:{v[0]}" for t, v in sorted(st.items())), any_induced=any(v[3] for v in st.values())))
R = pd.DataFrame(rows); RC = pd.DataFrame(recs)
RC.to_csv(OUT / (PFX + "receptors.tsv"), sep="\t", index=False, float_format="%.3f")
sp = sp.merge(R[["group_id", "pert_id", "d104", "poor_match", "receptors"]], on=["group_id", "pert_id"], how="left")
sp["mismatch"] = sp.d104.fillna("").str.startswith("receptor_mismatch")
sp["pct_d104"] = np.nan
for g, x in sp.groupby("group_id"):
    key = (~x.target_absent).astype(float) * 2e6 + (~x.mismatch).astype(float) * 1e6 + x.spec_mech.rank(method="max")
    sp.loc[x.index, "pct_d104"] = (key.rank(method="max", pct=True) * 100).round(1).values
sp[sp.mismatch | sp.d104.isin(["receptor_not_tested", "acts_on_tumor_environment"])][["group_id", "pert_id", "cmap_name", "source", "poor_match", "d104", "receptors", "spec_mech_pct", "pct_d104"]] \
    .sort_values("spec_mech_pct", ascending=False).to_csv(OUT / (PFX + "flags.tsv"), sep="\t", index=False)
# metrics (as D102 / D101)
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])
def auc(a, b):
    r = pd.Series(np.concatenate([a, b])).rank().values; return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
PRAD = [("TCGA-PRAD", d) for d in ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]]
M, CT = {}, []
for arm, col in [("production", "spec_mech_pct"), ("d104", "pct_d104")]:
    b = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)][col].max()
    o = {k: float(np.nanmedian([b(g, d) for g, d in v])) for k, v in [("ctl_pos", POS), ("ctl_neg", NEG), ("ctl_prad", PRAD)]}
    CT += [dict(arm=arm, group_id=g, drug=d, pct=b(g, d)) for g, d in POS + NEG + PRAD]
    m = lab.merge(sp.groupby(["group_id", "pert_id"])[col].max().reset_index(), left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        o[f"auc_{ls}"] = float(np.median([auc(g[g[ls]][col].dropna().values, g[~g[ls]][col].dropna().values) for c, g in m.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20]))
    M[arm] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / (PFX + "metrics.tsv"), sep="\t", float_format="%.4f")
pd.DataFrame(CT).pivot_table(index=["group_id", "drug"], columns="arm", values="pct").to_csv(OUT / (PFX + "controls.tsv"), sep="\t")
c, x = M["production"], M["d104"]
ok = [x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01,
      x["ctl_pos"] >= c["ctl_pos"] and x["ctl_prad"] >= c["ctl_prad"] and x["ctl_neg"] <= c["ctl_neg"]]
dec = pd.DataFrame([("1 known-use AUC falls <= 0.01", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}", ok[0]),
                    ("2 controls (breast pos, PRAD not lower; Basal not higher)", f"pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; PRAD {x['ctl_prad']:.1f} vs {c['ctl_prad']:.1f}; neg {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}", ok[1])], columns=["criterion", "value", "met"])
dec.loc[len(dec)] = ("DECISION", "adopt (needs the user's go-ahead)" if all(ok) else "do not adopt", all(ok)); dec.to_csv(OUT / (PFX + "decision.tsv"), sep="\t", index=False)
# reports
cnt = R.groupby(["poor_match", "d104"]).size().unstack(0).fillna(0).astype(int); cnt.to_csv(OUT / (PFX + "classes.tsv"), sep="\t")
top = sp[sp.mismatch & (sp.spec_mech_pct >= 90)]
labd = lab[lab.lab_hub | lab.lab_chembl].merge(sp[sp.mismatch][["group_id", "pert_id", "cmap_name", "receptors"]], left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
labd.to_csv(OUT / (PFX + "labeled_demoted.tsv"), sep="\t", index=False)
sens = {}
if AMEND == "a":
    for cut in [0.2, 0.3, 0.4]:
        st2 = RC.assign(s2=np.where(RC.support >= 0.5, "matched", np.where(RC.rho_estimate >= cut, "tumor_environment", "line_absent")))
        dl = st2.groupby(["group_id", "pert_id"]).s2.apply(lambda v: "ok" if "matched" in set(v) else ("env" if "tumor_environment" in set(v) else "mismatch"))
        sens[cut] = dl.value_counts().to_dict()
for cut in ([] if AMEND == "a" else [0.05, 0.1, 0.2]):
    st2 = np.where(RC.support < 0.5, "line_absent", np.where(RC.tumor_tpm < cut * RC.line_tpm, "tumor_far_lower", "matched"))
    sens[cut] = pd.Series(st2).value_counts().to_dict()
pd.DataFrame(sens).to_csv(OUT / (PFX + "sensitivity.tsv"), sep="\t")
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 70)
print("judged cohort x drug:", len(R)); print(cnt.to_string()); print("receptor status:", RC.status.value_counts().to_dict(), "| induced rescues:", int(RC.induced.sum()), "| not in LINCS genes:", int((~RC.in_lincs).sum()))
print("sensitivity (receptor rows):", sens)
print("demoted rows in top 10% before:", len(top), "| labeled demoted:", len(labd), "(hub", int(labd.lab_hub.sum()), ")")
print(labd[["group_id", "cmap_name", "lab_hub", "lab_chembl", "receptors"]].to_string(index=False))
print(met.round(4).to_string()); print(dec.to_string(index=False))
print(top.sort_values("spec_mech_pct", ascending=False)[["group_id", "cmap_name", "source", "receptors", "spec_mech_pct", "pct_d104"]].head(15).to_string(index=False))
