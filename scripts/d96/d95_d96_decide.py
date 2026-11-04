#!/usr/bin/env python
"""D95 / D96 decisions (rules fixed 2026-11-03, docs/DECISIONS.md). Runs in `celligner`.
Arms from worktrees ../repurposing_d94_<arm>: cur (CUR), lum / split (D94, reported), meta (D95), t85 (D96).
Adopt (each of meta, t85) if vs cur: (1) Basal anti-estrogens fall >= 10 spec_mech percentile points, breast positive
controls do not fall; (2) C28 A1 lead-score AUC falls <= 0.01 on Hub and ChEMBL; (3) C29 D1 falls <= 0.02; (4) dominance
rises <= 0.5 points. If both pass, the user chooses.
Outputs results/lincs/d96/{metrics.tsv, decision.tsv, top20_overlap.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
OUT = ROOT / "results/lincs/d96"; OUT.mkdir(parents=True, exist_ok=True)
ARMS = ["cur", "lum", "split", "meta", "t85"]
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
M, SP = {}, {}
for arm in ARMS:
    R = ROOT.parent / f"repurposing_d94_{arm}" / "results/lincs"
    sp = pd.read_csv(R / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False); SP[arm] = sp
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {"ctl_pos": float(np.nanmedian([best(g, d) for g, d in POS])), "ctl_neg": float(np.nanmedian([best(g, d) for g, d in NEG]))}
    a1 = pd.read_csv(R / "known_use/a1_summary.tsv", sep="\t"); a1 = a1[(a1.subset == "all") & (a1.score == "spec_mech_pct")].set_index("source")
    o["auc_hub"], o["auc_chembl"] = float(a1.median_auc["hub"]), float(a1.median_auc["chembl"])
    a2 = pd.read_csv(R / "known_use/a2_summary.tsv", sep="\t").set_index("source")
    o["a2_hub"], o["a2_chembl"] = float(a2.median_diff["hub"]), float(a2.median_diff["chembl"])
    dp = pd.read_csv(R / "dependency/summary.tsv", sep="\t")
    o["d1_all"] = float(dp[(dp.test == "D1_lead_vs_dependency") & (dp.subset == "all")].estimate.iloc[0])
    r = pd.read_csv(R / "drug_reversal_v2_consensus.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean()); o["dominant_drugs"] = int((t.pos >= 0.8).sum())
    M[arm] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
c = M["cur"]
rows = []
for arm in ["meta", "t85", "lum", "split"]:
    x = M[arm]
    ok = [x["ctl_neg"] <= c["ctl_neg"] - 10 and x["ctl_pos"] >= c["ctl_pos"],
          x["auc_hub"] >= c["auc_hub"] - 0.01 and x["auc_chembl"] >= c["auc_chembl"] - 0.01,
          x["d1_all"] >= c["d1_all"] - 0.02, x["dominance"] - c["dominance"] <= 0.005]
    rows.append(dict(arm=arm, rule="D95" if arm == "meta" else "D96" if arm == "t85" else "D94 (reported)", c1_breast=ok[0], c2_known_use=ok[1],
                     c3_dependency=ok[2], c4_dominance=ok[3], all_met=all(ok),
                     detail=f"neg {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; AUC hub {x['auc_hub']:.3f} vs {c['auc_hub']:.3f}, "
                            f"chembl {x['auc_chembl']:.3f} vs {c['auc_chembl']:.3f}; D1 {x['d1_all']:+.4f} vs {c['d1_all']:+.4f}; dominance {100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%"))
dec = pd.DataFrame(rows); dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
tw = []
for arm in ["meta", "t85"]:
    for g in sorted(SP["cur"].group_id.unique()):
        a = SP["cur"][SP["cur"].group_id == g].groupby("pert_id").spec_mech.max().nlargest(20).index
        b = SP[arm][SP[arm].group_id == g].groupby("pert_id").spec_mech.max().nlargest(20).index
        if len(set(a) ^ set(b)):
            tw.append(dict(arm=arm, group_id=g, top20_overlap=len(set(a) & set(b))))
tw = pd.DataFrame(tw); tw.to_csv(OUT / "top20_overlap.tsv", sep="\t", index=False)
ae = pd.concat([SP[a][SP[a].group_id.str.startswith("BRCA") & SP[a].cmap_name.isin(ANTI_ER)].groupby(["group_id", "cmap_name"]).spec_mech_pct.max()
                .rename(a) for a in ARMS], axis=1)
ae.to_csv(OUT / "antier.tsv", sep="\t")
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 220)
print(met.round(4).to_string()); print("\n", dec.to_string(index=False))
print("\n", tw.pivot_table(index="group_id", columns="arm", values="top20_overlap").to_string())
print("\n", ae.round(1).to_string())
