#!/usr/bin/env python
"""D90 decision (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner`.
Arms (tumor ARACNe v2 TF only): CUR = 79e / 79i VERSION=v2 REGSET=tf RUN_TAG=_d85; MECH = the same with RUN_TAG=_d90mech
(profiles scripts/d90/d90_profiles.py). Reported: the consensus (REGULON=consensus) with _d85 and _d90mech.
Checks on the FDA rows with C28 known-use labels (results/lincs/known_use/labels.tsv), lead score spec_mech_pct:
  1 per cohort with >= 3 labelled and >= 20 unlabelled drugs: AUC labelled vs unlabelled; median over cohorts
  2 per labelled drug in >= 10 cohorts: mean percentile in labelled cohorts minus median elsewhere; mean over drugs
  3 controls (79i split_controls): positive / negative median spec_mech percentile
  4 growth-arrest dominance (D85): drugs scored in >= 33 of 41 groups, share with rev > 0 in >= 80% of them
  5 ESR1 in MCF7 (BRCA network) under anti-estrogens on MECH activity: share negative
Adopt MECH if: AUC +0.03 on Hub and not lower on ChEMBL; specificity gap not smaller on either; positive controls
fall <= 5 and negative rise <= 5 points; dominance not higher; ESR1 negative in most profiles.
Outputs results/lincs/d90/{decision.tsv, metrics.tsv, auc_per_cohort.tsv, agreement.tsv, esr1_mcf7.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; RES = ROOT / "results/lincs"; OUT = RES / "d90"; OUT.mkdir(exist_ok=True)
MIN_LAB, MIN_UNLAB, MIN_COH = 3, 20, 10
ARMS = {("aracne", "CUR"): "_v2_tf_d85", ("aracne", "MECH"): "_v2_tf_d90mech",
        ("consensus", "CUR"): "_v2_consensus_d85", ("consensus", "MECH"): "_v2_consensus_d90mech"}
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene", "ZK-164015", "ICI-182780"]
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "drug", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])


def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values
    return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


def metrics(sfx):
    sp = pd.read_csv(RES / f"drug_reversal_split{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "spec_mech_pct"], low_memory=False)
    sp = sp.groupby(["group_id", "pert_id"]).spec_mech_pct.max().reset_index()
    d = lab.merge(sp, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    out, per = {}, []
    for ls in ["lab_hub", "lab_chembl"]:
        a = []
        for c, g in d.groupby("cohort"):
            x, y = g[g[ls]].spec_mech_pct.dropna().values, g[~g[ls]].spec_mech_pct.dropna().values
            if len(x) >= MIN_LAB and len(y) >= MIN_UNLAB:
                a.append(auc(x, y)); per.append((ls, c, len(x), a[-1]))
        out[f"auc_{ls}"] = float(np.median(a)) if a else np.nan; out[f"cohorts_{ls}"] = len(a)
        gaps = []
        for dr, g in d.groupby("drug"):
            if g.cohort.nunique() >= MIN_COH and g[ls].any() and (~g[ls]).any():
                gaps.append(g[g[ls]].spec_mech_pct.mean() - g[~g[ls]].spec_mech_pct.median())
        out[f"spec_gap_{ls}"] = float(np.mean(gaps)) if gaps else np.nan; out[f"drugs_{ls}"] = len(gaps)
    ctl = pd.read_csv(RES / f"split_controls{sfx}.tsv", sep="\t")
    out["ctl_pos"] = float(ctl[ctl.expected == "high"].spec_mech.median()); out["ctl_neg"] = float(ctl[ctl.expected == "low"].spec_mech.median())
    r = pd.read_csv(RES / f"drug_reversal{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean()))
    t = t[t.n >= 33]; out["dominance"] = float((t.pos >= 0.8).mean()); out["dominance_drugs"] = len(t)
    return out, pd.DataFrame(per, columns=["labels", "cohort", "labelled", "auc"]), sp


M, PER, SP = {}, [], {}
for k, sfx in ARMS.items():
    M[k], per, SP[k] = metrics(sfx); PER.append(per.assign(set=k[0], arm=k[1]))
met = pd.DataFrame(M).T; met.index.names = ["set", "arm"]
met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
pd.concat(PER).to_csv(OUT / "auc_per_cohort.tsv", sep="\t", index=False, float_format="%.4f")

# ESR1 in MCF7 under anti-estrogens (BRCA network)
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); cmp = pd.read_csv(RES / "compounds.tsv", sep="\t")[["pert_id", "cmap_name"]]
rows = idx.reset_index().merge(cmp, on="pert_id")
rows = rows[(rows.ModelID == "ACH-000019") & rows.cmap_name.isin(ANTI_ER)]
src = pd.read_csv(OUT / "mech_source.tsv", sep="\t")
es = []
for arm, tag in [("CUR", "_d85"), ("MECH", "_d90mech")]:
    d = P / f"drug_viper_specific_v2_tf{tag}"; regs = (d / "BRCA.regulators.txt").read_text().split()
    A = np.fromfile(d / "BRCA.f32", dtype=np.float32).reshape(len(idx), len(regs))
    es.append(pd.DataFrame({"arm": arm, "drug": rows.cmap_name.values, "row": rows["index"].values, "ESR1": A[rows["index"].values, regs.index("ESR1")]}))
es = pd.concat(es).merge(src[["source"]].reset_index().rename(columns={"index": "row"}), on="row")
es.to_csv(OUT / "esr1_mcf7.tsv", sep="\t", index=False, float_format="%.3f")
esr_neg = float((es[es.arm == "MECH"].ESR1 < 0).mean())

# agreement between arms (aracne)
m = SP[("aracne", "CUR")].merge(SP[("aracne", "MECH")], on=["group_id", "pert_id"], suffixes=("_c", "_m"))
ag = []
for g, d in m.groupby("group_id"):
    ag.append((g, d.spec_mech_pct_c.corr(d.spec_mech_pct_m, method="spearman"),
               len(set(d.nlargest(20, "spec_mech_pct_c").pert_id) & set(d.nlargest(20, "spec_mech_pct_m").pert_id))))
ag = pd.DataFrame(ag, columns=["group_id", "spearman", "top20_overlap"]); ag.to_csv(OUT / "agreement.tsv", sep="\t", index=False, float_format="%.3f")

c, x = M[("aracne", "CUR")], M[("aracne", "MECH")]
crit = [("1 Hub AUC +0.03, ChEMBL not lower", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}",
         x["auc_lab_hub"] - c["auc_lab_hub"] >= 0.03 and x["auc_lab_chembl"] >= c["auc_lab_chembl"]),
        ("2 specificity gap not smaller", f"Hub {x['spec_gap_lab_hub']:+.1f} vs {c['spec_gap_lab_hub']:+.1f}; ChEMBL {x['spec_gap_lab_chembl']:+.1f} vs {c['spec_gap_lab_chembl']:+.1f}",
         x["spec_gap_lab_hub"] >= c["spec_gap_lab_hub"] and x["spec_gap_lab_chembl"] >= c["spec_gap_lab_chembl"]),
        ("3 controls", f"positive {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}; negative {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}",
         x["ctl_pos"] >= c["ctl_pos"] - 5 and x["ctl_neg"] <= c["ctl_neg"] + 5),
        ("4 dominance not higher", f"{100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%", x["dominance"] <= c["dominance"]),
        ("5 ESR1 negative in most MCF7 anti-estrogen profiles", f"{100*esr_neg:.0f}% negative (MECH, {int((es.arm == 'MECH').sum())} profiles)", esr_neg > 0.5)]
crit = pd.DataFrame(crit, columns=["criterion", "value", "met"])
crit.loc[len(crit)] = ("DECISION", "adopt MECH profiles" if crit.met.all() else "keep current profiles", bool(crit.met.all()))
crit.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(met.round(4).to_string()); print("\nsources:\n", src.source.value_counts().to_string())
print("\nagreement (aracne CUR vs MECH): median Spearman %.3f, median top-20 overlap %.0f" % (ag.spearman.median(), ag.top20_overlap.median()))
print("\nESR1 MCF7:\n", es.pivot_table(index=["drug", "row", "source"], columns="arm", values="ESR1").round(2).to_string())
print("\n", crit.to_string(index=False))
