#!/usr/bin/env python
"""D98 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): selectivity flag for broad reversers at the tier step,
evaluated on the production consensus outputs (after D97). Runs in `celligner`. No production file is changed.
Broad reverser: scored in >= 33 of the 41 groups and rev > 0 in >= 80% of them. A broad reverser keeps tier 1-2
(tier_spec_mech) in a group only if its spec_mech percentile there is >= MARGIN points above its own median over the
groups where it is scored (subgroups: median over subgroups and groups where scored; reported).
Checks (groups): (1) broad reversers' tier 1-2 rows fall >= 50%; (2) C28 labels: labelled share of tier 1-2 FDA rows
(precision) does not fall (Hub, ChEMBL); labelled pairs reaching tier 1-2 (recall) fall <= 10%.
Outputs results/lincs/d98/{decision.tsv, margins.tsv, lost_tiers.tsv, broad_reversers.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results/lincs"; OUT = RES / "d98"; OUT.mkdir(exist_ok=True)
MARGIN, SENS = 25, [15, 35]
r = pd.read_csv(RES / "drug_reversal_v2_consensus.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean()))
broad = set(t[(t.n >= 33) & (t.pos >= 0.8)].index)
pd.DataFrame({"pert_id": sorted(broad)}).to_csv(OUT / "broad_reversers.tsv", sep="\t", index=False)
sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False,
                 usecols=["group_id", "pert_id", "cmap_name", "spec_mech_pct", "tier_spec_mech"])
st = pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
st = st.rename(columns={"stratum": "group_id"})[["group_id", "pert_id", "cmap_name", "spec_mech_pct", "tier_spec_mech"]] if "tier_spec_mech" in st else None
med = sp.groupby("pert_id").spec_mech_pct.median()
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])


def apply(d, margin, medians):
    d = d.copy()
    sel = (d.spec_mech_pct - d.pert_id.map(medians)) >= margin
    d["tier_new"] = np.where(d.pert_id.isin(broad) & (d.tier_spec_mech <= 2) & ~sel, 3, d.tier_spec_mech)
    return d


def metrics(d):
    b = d[d.pert_id.isin(broad)]
    o = {"broad_t12_before": int((b.tier_spec_mech <= 2).sum()), "broad_t12_after": int((b.tier_new <= 2).sum())}
    m = lab.merge(d, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        for k, col in [("before", "tier_spec_mech"), ("after", "tier_new")]:
            t12 = m[m[col] <= 2]
            o[f"precision_{ls}_{k}"] = float(t12[ls].mean()) if len(t12) else np.nan
            o[f"recall_{ls}_{k}"] = int((t12[ls]).sum())
    return o


res = {mg: metrics(apply(sp, mg, med)) for mg in [MARGIN] + SENS}
mt = pd.DataFrame(res).T; mt.index.name = "margin"; mt.to_csv(OUT / "margins.tsv", sep="\t", float_format="%.4f")
d = apply(sp, MARGIN, med)
lost = d[(d.tier_spec_mech <= 2) & (d.tier_new == 3)][["group_id", "cmap_name", "pert_id", "spec_mech_pct"]].assign(drug_median=lambda x: x.pert_id.map(med).round(1))
lost.sort_values(["group_id", "spec_mech_pct"], ascending=[True, False]).to_csv(OUT / "lost_tiers.tsv", sep="\t", index=False)
x = res[MARGIN]
c1 = x["broad_t12_after"] <= 0.5 * x["broad_t12_before"]
c2 = all(x[f"precision_{ls}_after"] >= x[f"precision_{ls}_before"] - 1e-12 and x[f"recall_{ls}_after"] >= 0.9 * x[f"recall_{ls}_before"]
         for ls in ["lab_hub", "lab_chembl"])
dec = pd.DataFrame([
    ("1 broad reversers' tier 1-2 rows fall >= 50%", f"{x['broad_t12_before']} -> {x['broad_t12_after']}", c1),
    ("2 known-use precision not lower; recall falls <= 10% (Hub, ChEMBL)",
     "; ".join(f"{ls}: precision {x[f'precision_{ls}_before']:.3f} -> {x[f'precision_{ls}_after']:.3f}, recall {x[f'recall_{ls}_before']} -> {x[f'recall_{ls}_after']}"
               for ls in ["lab_hub", "lab_chembl"]), c2)], columns=["criterion", "value", "met"])
dec.loc[len(dec)] = ("DECISION", "adopt the selectivity flag (production change needs the user's go-ahead)" if c1 and c2 else "do not adopt", bool(c1 and c2))
dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
if st is not None:
    medall = pd.concat([sp, st]).groupby("pert_id").spec_mech_pct.median()
    ds = apply(st, MARGIN, medall)
    print(f"subgroups (reported): broad tier 1-2 rows {int(((ds.tier_spec_mech <= 2) & ds.pert_id.isin(broad)).sum())} -> {int(((ds.tier_new <= 2) & ds.pert_id.isin(broad)).sum())}")
bc = d[d.group_id.str.startswith("BRCA") & d.cmap_name.isin(["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"])].groupby(["group_id", "cmap_name"])[["tier_spec_mech", "tier_new"]].min()
pd.set_option("display.width", 220); print(f"broad reversers: {len(broad)}"); print(mt.to_string()); print(bc.to_string()); print(lost.cmap_name.value_counts().head(15).to_string()); print(dec.to_string(index=False))
