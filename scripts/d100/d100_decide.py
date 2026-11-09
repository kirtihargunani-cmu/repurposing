#!/usr/bin/env python
"""D100 decision (rule fixed 2026-11-03, docs/DECISIONS.md). Runs in `celligner`.
AGG = production (post D97); EARLY / LATE = _d100early / _d100late consensus outputs.
Reported columns: per group x drug spec_mech_early_pct, spec_mech_late_pct, timing class (both: early and late rev_mech > 0
and both in the top 25%; early only; late only; neither; timing unknown = no timed pair among the drug's lines for the group:
in-key lines with Celligner weight > 0 for direct drugs, any line for borrowed drugs) -> results/lincs/d100/timing_scores.tsv
Lead test (EARLY as lead vs AGG): (1) known-use AUC (D118: 79l primary metric) falls <= 0.01 on Hub and ChEMBL and rises >= 0.01 on one; (2) breast
positive-control median not lower, Basal negative-control median not higher; (3) dominance rises <= 0.5 points.
Outputs results/lincs/d100/{decision.tsv, metrics.tsv, timing_scores.tsv, timing_classes.tsv, antier.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results/lincs"; OUT = RES / "d100"; OUT.mkdir(exist_ok=True); P = ROOT / "data/processed/lincs"
ARMS = {"agg": "_v2_consensus", "early": "_v2_consensus_d100early", "late": "_v2_consensus_d100late"}
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
# D118: known-use AUC = the primary metric of 79l (pooled within-cohort AUC, eligible rows of anticancer drugs, drug-level labels)
import importlib.util, sys  # noqa: E401
sys.path.insert(0, str(ROOT / "scripts"))
_s = importlib.util.spec_from_file_location("ku", ROOT / "scripts/79l_known_use.py"); KU = importlib.util.module_from_spec(_s); _s.loader.exec_module(KU)
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "fda_ingredient", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl", "anticancer"]]


def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values
    return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


M, SP = {}, {}
for arm, sfx in ARMS.items():
    sp = pd.read_csv(RES / f"drug_reversal_split{sfx}.tsv", sep="\t", low_memory=False); SP[arm] = sp
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {"ctl_pos": float(np.nanmedian([best(g, d) for g, d in POS])), "ctl_neg": float(np.nanmedian([best(g, d) for g, d in NEG]))}
    dd = KU.drug_table(sp).merge(lab, on=["cohort", "fda_ingredient"])
    for src in ["hub", "chembl"]:
        X = dd[dd.eligible_reversal & dd.anticancer & ~dd[f"excl_{src}"]]   # D122: d100 stays on the reversal-era rows
        o[f"auc_lab_{src}"] = float(KU.pooled_auc(X, f"lab_{src}", "spec_mech_pct")["auc"])
    r = pd.read_csv(RES / f"drug_reversal{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean())
    M[arm] = o
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
# timing classes
timed = pd.read_csv(P / "d100/timed_pairs.tsv", sep="\t"); tl = timed.groupby("pert_id").ModelID.apply(set).to_dict()
w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t"); wl = w[w.weight > 0].groupby("group_id").ModelID.apply(set).to_dict()
src = pd.read_csv(RES / "drug_reversal_v2_consensus.tsv", sep="\t", usecols=["group_id", "pert_id", "source"]).drop_duplicates(["group_id", "pert_id"])
k = ["group_id", "pert_id"]
sc = SP["agg"].groupby(k).agg(cmap_name=("cmap_name", "first"), spec_mech_pct=("spec_mech_pct", "max")).reset_index()
for arm in ["early", "late"]:
    a = SP[arm].groupby(k).agg(**{f"spec_mech_{arm}_pct": ("spec_mech_pct", "max"), f"rev_mech_{arm}": ("rev_mech", "max")}).reset_index()
    sc = sc.merge(a, on=k, how="left")
sc = sc.merge(src, on=k, how="left")
def is_timed(g, p, s):
    lines = tl.get(p)
    if not lines:
        return False
    return bool(lines & wl.get(g, set())) if s == "direct" else True
sc["timed"] = [is_timed(g, p, s) for g, p, s in zip(sc.group_id, sc.pert_id, sc.source)]
e_ok = (sc.rev_mech_early > 0) & (sc.spec_mech_early_pct >= 75); l_ok = (sc.rev_mech_late > 0) & (sc.spec_mech_late_pct >= 75)
sc["timing_class"] = np.where(~sc.timed, "timing unknown", np.where(e_ok & l_ok, "both", np.where(e_ok, "early only", np.where(l_ok, "late only", "neither"))))
sc.to_csv(OUT / "timing_scores.tsv", sep="\t", index=False, float_format="%.3f")
tc = sc.timing_class.value_counts(); tc.to_csv(OUT / "timing_classes.tsv", sep="\t")
ae = pd.concat([SP[a][SP[a].group_id.str.startswith("BRCA") & SP[a].cmap_name.isin(ANTI_ER)].groupby(["group_id", "cmap_name"]).spec_mech_pct.max().rename(a) for a in ARMS], axis=1)
ae.to_csv(OUT / "antier.tsv", sep="\t")
agree = sc[sc.timed].groupby("group_id").apply(lambda d: d.spec_mech_early_pct.corr(d.spec_mech_late_pct, method="spearman")).median()
c = M["agg"]
def lead(x):
    return [x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01 and
            max(x["auc_lab_hub"] - c["auc_lab_hub"], x["auc_lab_chembl"] - c["auc_lab_chembl"]) >= 0.01,
            x["ctl_pos"] >= c["ctl_pos"] and x["ctl_neg"] <= c["ctl_neg"], x["dominance"] - c["dominance"] <= 0.005]
rows = []
for arm in ["early", "late"]:
    x = M[arm]; ok = lead(x)
    rows.append(dict(arm=arm, role="lead test" if arm == "early" else "reported", c1_known_use=ok[0], c2_breast=ok[1], c3_dominance=ok[2], all_met=all(ok),
                     detail=f"AUC hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}, chembl {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}; "
                            f"pos {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; neg {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; dominance {100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%"))
dec = pd.DataFrame(rows)
dec.loc[len(dec)] = dict(arm="DECISION", detail=("early profiles become the lead for timed pairs" if dec.iloc[0].all_met else "aggregated stays the lead") +
                         "; early / late scores and timing class added as reported columns")
dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 200)
print(met.round(4).to_string()); print("\ntiming classes:\n", tc.to_string()); print(f"\nmedian per-group Spearman early vs late (timed rows): {agree:.3f}")
print("\n", ae.round(1).to_string()); print("\n", dec.to_string(index=False))
