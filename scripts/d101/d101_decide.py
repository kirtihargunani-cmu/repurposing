#!/usr/bin/env python
"""D101 decision (rule fixed 2026-11-03, docs/DECISIONS.md). Production vs _d101 (ESR1 / AR read from PROGENy).
Adopt if: (1) breast positive controls rise >= 5 spec_mech percentile points, Basal negative controls not higher;
(2) PRAD controls (bicalutamide, enzalutamide, flutamide, nilutamide) not lower; (3) known-use AUC falls <= 0.01
(Hub, ChEMBL); (4) dominance rises <= 0.5 points. Outputs results/lincs/d101/{decision,metrics,controls}.tsv"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; RES = ROOT / "results/lincs"; OUT = RES / "d101"; OUT.mkdir(exist_ok=True)
ARMS = {"production": "_v2_consensus", "d101": "_v2_consensus_d101"}
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
PRAD = [("TCGA-PRAD", d) for d in ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]]
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])
def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values; return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))
M, C = {}, []
for arm, sfx in ARMS.items():
    sp = pd.read_csv(RES / f"drug_reversal_split{sfx}.tsv", sep="\t", low_memory=False)
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {k: float(np.nanmedian([best(g, d) for g, d in v])) for k, v in [("ctl_pos", POS), ("ctl_neg", NEG), ("ctl_prad", PRAD)]}
    C += [dict(arm=arm, group_id=g, drug=d, spec_mech_pct=best(g, d)) for g, d in POS + NEG + PRAD]
    s = sp.groupby(["group_id", "pert_id"]).spec_mech_pct.max().reset_index()
    dd = lab.merge(s, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        o[f"auc_{ls}"] = float(np.median([auc(g[g[ls]].spec_mech_pct.dropna().values, g[~g[ls]].spec_mech_pct.dropna().values)
                                          for c, g in dd.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20]))
    r = pd.read_csv(RES / f"drug_reversal{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean()); M[arm] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
ct = pd.DataFrame(C).pivot_table(index=["group_id", "drug"], columns="arm", values="spec_mech_pct"); ct.to_csv(OUT / "controls.tsv", sep="\t")
c, x = M["production"], M["d101"]
ok = [x["ctl_pos"] >= c["ctl_pos"] + 5 and x["ctl_neg"] <= c["ctl_neg"], x["ctl_prad"] >= c["ctl_prad"],
      x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01, x["dominance"] - c["dominance"] <= 0.005]
dec = pd.DataFrame([("1 breast: positive +5, Basal not higher", f"pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; neg {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}", ok[0]),
                    ("2 prostate controls not lower", f"{x['ctl_prad']:.1f} vs {c['ctl_prad']:.1f}", ok[1]),
                    ("3 known-use AUC falls <= 0.01", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}", ok[2]),
                    ("4 dominance rises <= 0.5", f"{100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%", ok[3])], columns=["criterion", "value", "met"])
dec.loc[len(dec)] = ("DECISION", "adopt (production change needs the user's go-ahead)" if all(ok) else "do not adopt", all(ok))
dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(met.round(4).to_string()); print(ct.round(1).to_string()); print(dec.to_string(index=False))
