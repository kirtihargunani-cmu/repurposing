#!/usr/bin/env python
"""D93 decision (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner`.
Pooled = _d85 (tumor ARACNe, current profiles); cell-type = _d93 (20 groups of the 8 split cohorts on their networks).
  1 breast benefit: anti-estrogens in BRCA-Basal (fulvestrant / tamoxifen / raloxifene, spec_mech percentile, best ID)
    fall >= 10 points; breast positive controls (fulvestrant / tamoxifen / raloxifene / palbociclib in LumA, LumB) do not fall
  2 no harm in the new cohorts: median per-group ChEMBL known-use AUC over the non-breast split groups (>= 3 labelled,
    >= 20 unlabelled FDA drugs) falls by <= 0.02
  3 C28 known-use AUC over the 41 groups (median per cohort) falls by <= 0.01 on Hub and ChEMBL
  4 growth-arrest dominance (D85, on rev) rises by <= 0.5 points
Reported: per-group AUC; disease side per group (weighted Spearman, sign flips with |NES| >= 1.96); anti-estrogens per
breast group; network sizes. Outputs results/lincs/d93/{decision,metrics,auc_per_group,disease_side,antier}.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
RES = ROOT / "results/lincs"; OUT = RES / "d93"; OUT.mkdir(exist_ok=True)
ARMS = {"pooled": "_v2_tf_d85", "cell_type": "_v2_tf_d93"}
SPLIT = pd.read_csv(ROOT / "results/step5/d93_group_net.tsv", sep="\t").group_id.tolist()
NONBRCA = [g for g in SPLIT if not g.startswith("BRCA")]
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])


def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values
    return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


M, PER, AE = {}, [], []
for arm, sfx in ARMS.items():
    sp = pd.read_csv(RES / f"drug_reversal_split{sfx}.tsv", sep="\t", low_memory=False)
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {"ctl_pos": float(np.nanmedian([best(g, d) for g, d in POS])), "ctl_neg": float(np.nanmedian([best(g, d) for g, d in NEG]))}
    s = sp.groupby(["group_id", "pert_id"]).spec_mech_pct.max().reset_index()
    d = lab.merge(s, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        a = {c: auc(g[g[ls]].spec_mech_pct.dropna().values, g[~g[ls]].spec_mech_pct.dropna().values)
             for c, g in d.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20}
        o[f"auc_{ls}"] = float(np.median(list(a.values())))
        PER += [(arm, ls, c, v) for c, v in a.items()]
    nb = [v for (ar, ls, c, v) in PER if ar == arm and ls == "lab_chembl" and c in NONBRCA]
    o["auc_chembl_nonbrca_split"] = float(np.median(nb)) if nb else np.nan; o["n_nonbrca_groups"] = len(nb)
    r = pd.read_csv(RES / f"drug_reversal{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean())
    M[arm] = o
    AE.append(sp[sp.group_id.str.startswith("BRCA") & sp.cmap_name.isin(ANTI_ER)].groupby(["group_id", "cmap_name"]).spec_mech_pct.max().reset_index().assign(arm=arm))
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
per = pd.DataFrame(PER, columns=["arm", "labels", "group_id", "auc"]).pivot_table(index=["labels", "group_id"], columns="arm", values="auc").reset_index()
per.to_csv(OUT / "auc_per_group.tsv", sep="\t", index=False, float_format="%.3f")
ae = pd.concat(AE); ae.to_csv(OUT / "antier.tsv", sep="\t", index=False)
a = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf_d93/group_mr.tsv", sep="\t")
b = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf_d93pooled/group_mr.tsv", sep="\t")
ds = []
for g in sorted(set(a.group_id)):
    x = a[a.group_id == g].merge(b[b.group_id == g], on="regulator", suffixes=("_c", "_p"))
    rc, rp = x.nes_c.rank(), x.nes_p.rank(); w = (x.nes_c.abs() + x.nes_p.abs()) / 2
    rw = np.cov(rc, rp, aweights=w)[0, 1] / np.sqrt(np.cov(rc, aweights=w) * np.cov(rp, aweights=w))
    fl = x[(np.sign(x.nes_c) != np.sign(x.nes_p)) & ((x.nes_c.abs() >= 1.96) | (x.nes_p.abs() >= 1.96))]
    ds.append(dict(group_id=g, shared_regulators=len(x), weighted_spearman=round(float(rw), 3), sign_flips=len(fl),
                   top_flips=", ".join(f"{r} {p:+.1f}->{c:+.1f}" for r, p, c in fl.reindex(fl.nes_p.abs().sort_values(ascending=False).index)[["regulator", "nes_p", "nes_c"]].head(5).values)))
ds = pd.DataFrame(ds); ds.to_csv(OUT / "disease_side.tsv", sep="\t", index=False)
c, x = M["pooled"], M["cell_type"]
crit = pd.DataFrame([
    ("1 breast: Basal anti-estrogens fall >= 10; positive controls do not fall",
     f"negative {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; positive {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}",
     x["ctl_neg"] <= c["ctl_neg"] - 10 and x["ctl_pos"] >= c["ctl_pos"]),
    ("2 non-breast split groups: ChEMBL AUC falls <= 0.02",
     f"{x['auc_chembl_nonbrca_split']:.3f} vs {c['auc_chembl_nonbrca_split']:.3f} ({int(x['n_nonbrca_groups'])} groups)",
     x["auc_chembl_nonbrca_split"] >= c["auc_chembl_nonbrca_split"] - 0.02),
    ("3 41 groups: AUC falls <= 0.01 (Hub, ChEMBL)", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}",
     x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01),
    ("4 dominance rises <= 0.5 points", f"{100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%", x["dominance"] - c["dominance"] <= 0.005)],
    columns=["criterion", "value", "met"])
crit.loc[len(crit)] = ("DECISION", "adopt cell-type networks for the 8 split cohorts" if crit.met.all() else "keep the pooled networks", bool(crit.met.all()))
crit.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30); pd.set_option("display.max_colwidth", 80)
print(met.round(4).to_string()); print("\n", per.round(3).to_string(index=False)); print("\n", ds.to_string(index=False))
print("\n", ae.pivot_table(index=["group_id", "cmap_name"], columns="arm", values="spec_mech_pct").round(1).to_string())
print("\n", crit.to_string(index=False))
