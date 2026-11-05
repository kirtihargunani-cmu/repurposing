#!/usr/bin/env python
"""D105 decision (rule fixed 2026-11-04, docs/DECISIONS.md): production (D104c) vs _v2_consensus_d105 (ChEMBL primary
targets for every drug). Outputs results/lincs/d105/{metrics,controls,decision,sources,changes,newly_demoted_top10,labeled_newly_demoted}.tsv"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; RES = ROOT / "results/lincs"; OUT = RES / "d105"; OUT.mkdir(exist_ok=True)
A = {"production": "_v2_consensus", "d105": "_v2_consensus_d105"}
T = {a: pd.read_csv(RES / f"drug_reversal_split{s}.tsv", sep="\t", low_memory=False) for a, s in A.items()}
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])
def auc(a, b):
    r = pd.Series(np.concatenate([a, b])).rank().values; return (r[:len(a)].sum() - len(a) * (len(a) + 1) / 2) / (len(a) * len(b))
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
PRAD = [("TCGA-PRAD", d) for d in ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]]
truthy = lambda s: s.astype(str).str.lower().eq("true")
M, C = {}, []
for a, sp in T.items():
    b = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {k: float(np.nanmedian([b(g, d) for g, d in v])) for k, v in [("ctl_pos", POS), ("ctl_neg", NEG), ("ctl_prad", PRAD)]}
    C += [dict(arm=a, group_id=g, drug=d, pct=b(g, d)) for g, d in POS + NEG + PRAD + [("BRCA-Basal", "toremifene"), ("BRCA-Basal", "bazedoxifene"), ("BRCA-Basal", "clomifene")]]
    m = lab.merge(sp.groupby(["group_id", "pert_id"]).spec_mech_pct.max().reset_index(), left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        o[f"auc_{ls}"] = float(np.median([auc(g[g[ls]].spec_mech_pct.dropna().values, g[~g[ls]].spec_mech_pct.dropna().values) for c, g in m.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20]))
    o["target_absent_rows"] = int(truthy(sp.target_absent).sum()); o["receptor_mismatch_rows"] = int(sp.receptor_status.eq("receptor_mismatch").sum())
    o["receptor_env_rows"] = int(sp.receptor_status.eq("acts_on_tumor_environment").sum())
    M[a] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
pd.DataFrame(C).pivot_table(index=["group_id", "drug"], columns="arm", values="pct").to_csv(OUT / "controls.tsv", sep="\t")
c, x = M["production"], M["d105"]
ok = [x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01,
      x["ctl_pos"] >= c["ctl_pos"] and x["ctl_prad"] >= c["ctl_prad"] and x["ctl_neg"] <= c["ctl_neg"]]
dec = pd.DataFrame([("1 known-use AUC falls <= 0.01", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}", ok[0]),
                    ("2 controls", f"pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; PRAD {x['ctl_prad']:.1f} vs {c['ctl_prad']:.1f}; Basal {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}", ok[1])], columns=["criterion", "value", "met"])
dec.loc[len(dec)] = ("DECISION", "adopt (needs the user's go-ahead)" if all(ok) else "do not adopt", all(ok)); dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
d = T["d105"]
src = d.drop_duplicates("cmap_name").target_source.replace("", "not judged").value_counts(); src.to_csv(OUT / "sources.tsv", sep="\t")
k = ["group_id", "pert_id"]
j = T["production"][k + ["cmap_name", "spec_mech_pct", "target_absent", "receptor_status", "evidence_targets"]].merge(
    d[k + ["spec_mech_pct", "target_absent", "receptor_status", "evidence_targets", "target_source", "secondary_targets"]], on=k, suffixes=("_prod", "_d105"))
dem = lambda s, r: truthy(s) | r.eq("receptor_mismatch")
j["dem_prod"], j["dem_d105"] = dem(j.target_absent_prod, j.receptor_status_prod), dem(j.target_absent_d105, j.receptor_status_d105)
ch = pd.crosstab(j.dem_prod, j.dem_d105); ch.to_csv(OUT / "changes.tsv", sep="\t")
nd = j[~j.dem_prod & j.dem_d105 & (j.spec_mech_pct_prod >= 90)].sort_values("spec_mech_pct_prod", ascending=False); nd.to_csv(OUT / "newly_demoted_top10.tsv", sep="\t", index=False)
ln = lab[lab.lab_hub | lab.lab_chembl].merge(j[~j.dem_prod & j.dem_d105], left_on=["cohort", "pert_id"], right_on=k); ln.to_csv(OUT / "labeled_newly_demoted.tsv", sep="\t", index=False)
lr = lab[lab.lab_hub | lab.lab_chembl].merge(j[j.dem_prod & ~j.dem_d105], left_on=["cohort", "pert_id"], right_on=k); lr.to_csv(OUT / "labeled_newly_restored.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 60)
print("sources:", src.to_dict()); print(met.round(4).to_string()); print(dec.to_string(index=False)); print("demoted prod x d105:\n", ch)
print(pd.read_csv(OUT / "controls.tsv", sep="\t").round(1).to_string(index=False))
print("newly demoted in top 10%:", len(nd)); print(nd[["group_id", "cmap_name", "target_source", "evidence_targets_d105", "spec_mech_pct_prod"]].head(15).to_string(index=False))
print("labeled newly demoted:", len(ln), "hub", int(ln.lab_hub.sum())); print(ln[["group_id", "cmap_name", "lab_hub", "lab_chembl", "target_source", "evidence_targets_d105"]].to_string(index=False))
print("labeled newly restored:", len(lr), "hub", int(lr.lab_hub.sum()))
