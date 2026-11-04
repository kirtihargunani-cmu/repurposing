#!/usr/bin/env python
"""D91 decision (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner`.
Input: 79i with IDENTITY_FILE (results/lincs/drug_reversal_split_v2_tf_d85_id.tsv; consensus _v2_consensus_d85_id reported).
Compares the lead score spec_mech with spec_mech_id on the same rows:
  1 controls (79i CONTROLS: best LINCS ID per drug): negative median (anti-estrogens in BRCA-Basal) falls >= 10
    percentile points; positive median does not fall
  2 C28 known-use: median per-cohort AUC (>= 3 labelled, >= 20 unlabelled FDA drugs) falls <= 0.01 on Hub and on
    ChEMBL; specificity gap (labelled cohorts minus median elsewhere, drugs in >= 10 cohorts) shrinks <= 1 point on each
  3 growth-arrest dominance on the matching reversal score (rev_mech for spec_mech, rev_mech_id for spec_mech_id): share
    of drugs scored in >= 33 of 41 groups with score > 0 in >= 80% of them; does not rise
Outputs results/lincs/d91/{decision.tsv, metrics.tsv, antier.tsv, identity_summary.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results/lincs"; OUT = RES / "d91"
CONTROLS = [("BRCA-LumA", "high", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("BRCA-LumB", "high", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("TCGA-PRAD", "high", ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]),
            ("TCGA-COAD", "high", ["fluorouracil", "irinotecan", "oxaliplatin"]),
            ("TCGA-LUAD", "high", ["erlotinib", "gefitinib", "afatinib", "osimertinib"]),
            ("BRCA-Basal", "low", ["fulvestrant", "tamoxifen", "raloxifene"])]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "drug", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])


def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values
    return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


def metrics(sp, score, rev):
    s = sp.groupby(["group_id", "pert_id"])[score].max().reset_index()
    d = lab.merge(s, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    o = {}
    for ls in ["lab_hub", "lab_chembl"]:
        a = [auc(g[g[ls]][score].dropna().values, g[~g[ls]][score].dropna().values) for c, g in d.groupby("cohort")
             if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20]
        o[f"auc_{ls}"] = float(np.median(a)); o[f"cohorts_{ls}"] = len(a)
        gaps = [g[g[ls]][score].mean() - g[~g[ls]][score].median() for dr, g in d.groupby("drug")
                if g.cohort.nunique() >= 10 and g[ls].any() and (~g[ls]).any()]
        o[f"spec_gap_{ls}"] = float(np.mean(gaps)); o[f"drugs_{ls}"] = len(gaps)
    pos, neg = [], []
    for g, e, names in CONTROLS:
        x = sp[sp.group_id == g]
        for nm in names:
            r = x[x.cmap_name == nm]
            if len(r):
                (pos if e == "high" else neg).append(r[score].max())
    o["ctl_pos"], o["ctl_neg"], o["n_pos"], o["n_neg"] = float(np.median(pos)), float(np.median(neg)), len(pos), len(neg)
    t = sp.groupby("pert_id").agg(n=("group_id", "nunique"), pos=(rev, lambda x: (x > 0).mean()))
    t = t[t.n >= 33]; o["dominance"] = float((t.pos >= 0.8).mean()); o["dominance_drugs"] = len(t)
    return o


res = {}
for st, sfx in [("aracne", "_v2_tf_d85_id"), ("consensus", "_v2_consensus_d85_id")]:
    f = RES / f"drug_reversal_split{sfx}.tsv"
    if not f.exists():
        continue
    sp = pd.read_csv(f, sep="\t", low_memory=False)
    res[(st, "spec_mech")] = metrics(sp, "spec_mech_pct", "rev_mech")
    res[(st, "spec_mech_id")] = metrics(sp, "spec_mech_id_pct", "rev_mech_id")
    ae = sp[sp.group_id.isin(["BRCA-LumA", "BRCA-LumB", "BRCA-Basal"]) & sp.cmap_name.isin(ANTI_ER)]
    ae.groupby(["group_id", "cmap_name"])[["rev_pct", "spec_mech_pct", "spec_mech_id_pct", "rev_identity_pct"]].max().reset_index() \
      .assign(set=st).to_csv(OUT / f"antier_{st}.tsv", sep="\t", index=False)
    if st == "aracne":
        A = sp
met = pd.DataFrame(res).T; met.index.names = ["set", "score"]; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
it = pd.read_csv(OUT / "identity_regulators.tsv", sep="\t")
idsum = A.drop_duplicates("group_id")[["group_id", "n_mrs", "n_mrs_identity"]].assign(share=lambda x: (x.n_mrs_identity / x.n_mrs).round(3))
idsum.to_csv(OUT / "identity_summary.tsv", sep="\t", index=False)
c, x = res[("aracne", "spec_mech")], res[("aracne", "spec_mech_id")]
crit = pd.DataFrame([
    ("1 negative controls fall >= 10 points, positive do not fall",
     f"negative {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; positive {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}",
     x["ctl_neg"] <= c["ctl_neg"] - 10 and x["ctl_pos"] >= c["ctl_pos"]),
    ("2 known-use AUC falls <= 0.01 and specificity gap shrinks <= 1 point, Hub and ChEMBL",
     f"AUC Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}, ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}; "
     f"gap Hub {x['spec_gap_lab_hub']:+.1f} vs {c['spec_gap_lab_hub']:+.1f}, ChEMBL {x['spec_gap_lab_chembl']:+.1f} vs {c['spec_gap_lab_chembl']:+.1f}",
     x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01
     and x["spec_gap_lab_hub"] >= c["spec_gap_lab_hub"] - 1 and x["spec_gap_lab_chembl"] >= c["spec_gap_lab_chembl"] - 1),
    ("3 dominance does not rise", f"{100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%", x["dominance"] <= c["dominance"])],
    columns=["criterion", "value", "met"])
crit.loc[len(crit)] = ("DECISION", "adopt spec_mech_id as the lead score" if crit.met.all() else "keep spec_mech", bool(crit.met.all()))
crit.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(met.round(4).to_string())
print("\nidentity-linked share of tested regulators per group (median %.3f):\n" % idsum.share.median(), idsum.sort_values("share", ascending=False).head(12).to_string(index=False))
print("\nanti-estrogens (aracne):\n", pd.read_csv(OUT / "antier_aracne.tsv", sep="\t").pivot_table(index="cmap_name", columns="group_id", values=["spec_mech_pct", "spec_mech_id_pct"]).round(1).to_string())
print("\n", crit.to_string(index=False))
