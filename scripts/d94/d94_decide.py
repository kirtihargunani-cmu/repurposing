#!/usr/bin/env python
"""D94 decision (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner`.
Arms: CUR (_d85), LUM (_d94lum), SPLIT (_d94split), and GATED (_d94gated: the cohort gate applied, failing cohorts pooled).
Each arm's downstream outputs come from its worktree ../repurposing_d94_<arm>/results/lincs.
Adopt LUM or SPLIT if vs CUR: (1) Basal anti-estrogens fall >= 10 spec_mech percentile points and breast positive controls
do not fall; (2) C28 A1 lead-score AUC (all cohorts) falls <= 0.01 on Hub and ChEMBL; (3) C29 D1 (all) falls <= 0.02;
(4) dominance rises <= 0.5 points. LUM vs SPLIT: higher breast positive-control median, within 2 points SPLIT.
Cohort gate (results/lincs/d94/cohort_gate.tsv): cohorts failing stay pooled -> the production configuration is GATED,
whose checks (2)-(4) are reported alongside.
Outputs results/lincs/d94/{metrics.tsv, decision.tsv, top20_overlap.tsv, validation_<arm>.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results/lincs/d94"; OUT.mkdir(parents=True, exist_ok=True)
ARMS = {"cur": "_d85", "lum": "_d94lum", "split": "_d94split", "gated": "_d94gated"}
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
SPLITC = ["BRCA", "CESC", "ESCA", "LGG", "PCPG", "SARC", "TGCT", "UCEC"]
M, SP = {}, {}
for arm, tag in ARMS.items():
    R = ROOT.parent / f"repurposing_d94_{arm}" / "results/lincs"
    sp = pd.read_csv(R / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False); SP[arm] = sp
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {"ctl_pos": float(np.nanmedian([best(g, d) for g, d in POS])), "ctl_neg": float(np.nanmedian([best(g, d) for g, d in NEG]))}
    a1 = pd.read_csv(R / "known_use/a1_summary.tsv", sep="\t"); a1 = a1[(a1.subset == "all") & (a1.score == "spec_mech_pct")].set_index("source")
    o["auc_hub"], o["auc_chembl"] = float(a1.median_auc["hub"]), float(a1.median_auc["chembl"])
    a2 = pd.read_csv(R / "known_use/a2_summary.tsv", sep="\t").set_index("source")
    o["a2_hub_median_diff"], o["a2_chembl_median_diff"] = float(a2.median_diff["hub"]), float(a2.median_diff["chembl"])
    dp = pd.read_csv(R / "dependency/summary.tsv", sep="\t")
    o["d1_all"] = float(dp[(dp.test == "D1_lead_vs_dependency") & (dp.subset == "all")].estimate.iloc[0])
    o["d3_chembl"] = float(dp[(dp.test == "D3_known_use_dependency") & (dp.subset == "chembl")].estimate.iloc[0])
    gs = pd.read_csv(R / "gate_sensitivity/summary.tsv", sep="\t")
    for t, ds in [("B2_rho_matched", "prism_sec"), ("B2_rho_borrowed", "prism_sec"), ("B2_rho_matched", "prism_pri"), ("B2_rho_borrowed", "prism_pri")]:
        x = gs[(gs.test == t) & (gs.dataset == ds)]; o[f"{t}_{ds}"] = float(x.estimate.iloc[0]) if len(x) else np.nan
    r = pd.read_csv(R / "drug_reversal_v2_consensus.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean())
    M[arm] = o
    for f in ["known_use/a1_summary.tsv", "dependency/summary.tsv", "gate_sensitivity/summary.tsv"]:
        pd.read_csv(R / f, sep="\t").assign(arm=arm).to_csv(OUT / f"validation_{arm}_{f.split('/')[0]}.tsv", sep="\t", index=False)
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
# top-20 overlap with CUR, split-cohort groups
tw = []
for arm in ["lum", "split", "gated"]:
    for g in sorted(SP["cur"].group_id.unique()):
        if g.split("-")[0] not in SPLITC and not g.startswith(tuple(f"TCGA-{c}" for c in SPLITC)):
            continue
        a = SP["cur"][SP["cur"].group_id == g].groupby("pert_id").spec_mech.max().nlargest(20).index
        b = SP[arm][SP[arm].group_id == g].groupby("pert_id").spec_mech.max().nlargest(20).index
        tw.append(dict(arm=arm, group_id=g, top20_overlap=len(set(a) & set(b))))
tw = pd.DataFrame(tw); tw.to_csv(OUT / "top20_overlap.tsv", sep="\t", index=False)
gate = pd.read_csv(OUT / "cohort_gate.tsv", sep="\t")
c = M["cur"]
def checks(x):
    return [x["ctl_neg"] <= c["ctl_neg"] - 10 and x["ctl_pos"] >= c["ctl_pos"],
            x["auc_hub"] >= c["auc_hub"] - 0.01 and x["auc_chembl"] >= c["auc_chembl"] - 0.01,
            x["d1_all"] >= c["d1_all"] - 0.02, x["dominance"] - c["dominance"] <= 0.005]
rows = []
for arm in ["lum", "split", "gated"]:
    x = M[arm]; ok = checks(x)
    rows.append(dict(arm=arm, c1_breast=ok[0], c2_known_use=ok[1], c3_dependency=ok[2], c4_dominance=ok[3], all_met=all(ok),
                     detail=f"neg {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; AUC hub {x['auc_hub']:.3f} vs {c['auc_hub']:.3f}, "
                            f"chembl {x['auc_chembl']:.3f} vs {c['auc_chembl']:.3f}; D1 {x['d1_all']:+.4f} vs {c['d1_all']:+.4f}; dominance {100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%"))
dec = pd.DataFrame(rows)
cands = [a for a in ["lum", "split"] if dec.set_index("arm").all_met[a]]
if not cands:
    verdict = "keep CUR (pooled networks): neither LUM nor SPLIT meets (1)-(4)"
else:
    win = "split" if (len(cands) == 2 and M["lum"]["ctl_pos"] - M["split"]["ctl_pos"] <= 2) or cands == ["split"] else "lum"
    fail = gate[(gate.arm == win) & ~gate.passes].cohort.tolist()
    verdict = (f"adopt {win.upper()}" + (f"; cohort gate keeps {', '.join(fail)} pooled -> production configuration = GATED "
               f"(GATED checks (2)-(4): {'met' if all(dec.set_index('arm').loc['gated', ['c2_known_use', 'c3_dependency', 'c4_dominance']]) else 'NOT met'})" if fail else ""))
dec.loc[len(dec)] = dict(arm="DECISION", detail=verdict)
dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 200)
print(met.round(4).to_string()); print("\n", tw.pivot_table(index="group_id", columns="arm", values="top20_overlap").to_string())
print("\n", gate.to_string(index=False)); print("\n", dec.to_string(index=False))
