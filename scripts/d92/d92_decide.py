#!/usr/bin/env python
"""D92 decision (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner`.
Pooled = the _d85 run (tumor ARACNe, current profiles); cell-type = the _d92 run (breast groups on BRCA_LUM / BRCA_BASX).
  1 negative-control median (anti-estrogens fulvestrant / tamoxifen / raloxifene in BRCA-Basal, spec_mech percentile,
    best LINCS ID per drug) falls >= 10 points; positive-control median over the breast controls (fulvestrant,
    tamoxifen, raloxifene, palbociclib in LumA and LumB) does not fall
  2 C28 known-use AUC over the 41 groups (median per cohort, spec_mech_pct) falls <= 0.01 on Hub and ChEMBL
  3 growth-arrest dominance (D85, on rev) does not rise
Reported: anti-estrogens per breast group; ESR1 / FOXA1 / GATA3 / XBP1 / PGR in MCF7 under anti-estrogens per network;
disease side: weighted Spearman of group NES (cell-type vs pooled, same seeds) and sign of ESR1 / FOXA1 / GATA3 / XBP1.
Outputs results/lincs/d92/{decision.tsv, metrics.tsv, antier.tsv, mcf7_reading.tsv, disease_side.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
import os
TAG = os.environ.get("TAG", "_d92")   # D92a reported arm: TAG=_d92b MAP=results/step5/d92b_group_net.tsv
MAP = os.environ.get("MAP", "results/step5/d92_group_net.tsv")
RES = ROOT / "results/lincs"; OUT = RES / TAG.strip("_"); OUT.mkdir(exist_ok=True); P = ROOT / "data/processed/lincs"
ARMS = {"pooled": "_v2_tf_d85", "cell_type": f"_v2_tf{TAG}"}
POS = [(g, d) for g in ["BRCA-LumA", "BRCA-LumB"] for d in ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]]
NEG = [("BRCA-Basal", d) for d in ["fulvestrant", "tamoxifen", "raloxifene"]]
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "pert_id", "lab_hub", "lab_chembl"]].drop_duplicates(["cohort", "pert_id"])


def auc(x, y):
    r = pd.Series(np.concatenate([x, y])).rank().values
    return (r[:len(x)].sum() - len(x) * (len(x) + 1) / 2) / (len(x) * len(y))


M, AE = {}, []
for arm, sfx in ARMS.items():
    sp = pd.read_csv(RES / f"drug_reversal_split{sfx}.tsv", sep="\t", low_memory=False)
    best = lambda g, d: sp[(sp.group_id == g) & (sp.cmap_name == d)].spec_mech_pct.max()
    o = {"ctl_pos": float(np.nanmedian([best(g, d) for g, d in POS])), "ctl_neg": float(np.nanmedian([best(g, d) for g, d in NEG]))}
    s = sp.groupby(["group_id", "pert_id"]).spec_mech_pct.max().reset_index()
    d = lab.merge(s, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
    for ls in ["lab_hub", "lab_chembl"]:
        o[f"auc_{ls}"] = float(np.median([auc(g[g[ls]].spec_mech_pct.dropna().values, g[~g[ls]].spec_mech_pct.dropna().values)
                                          for c, g in d.groupby("cohort") if g[ls].sum() >= 3 and (~g[ls]).sum() >= 20]))
    r = pd.read_csv(RES / f"drug_reversal{sfx}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
    t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean())); t = t[t.n >= 33]
    o["dominance"] = float((t.pos >= 0.8).mean())
    M[arm] = o
    AE.append(sp[sp.group_id.str.startswith("BRCA") & sp.cmap_name.isin(ANTI_ER)].groupby(["group_id", "cmap_name"])[["rev_pct", "spec_mech_pct"]].max()
              .reset_index().assign(arm=arm))
met = pd.DataFrame(M).T; met.to_csv(OUT / "metrics.tsv", sep="\t", float_format="%.4f")
ae = pd.concat(AE); ae.to_csv(OUT / "antier.tsv", sep="\t", index=False)

# MCF7 reading per network
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); cmp = pd.read_csv(RES / "compounds.tsv", sep="\t")[["pert_id", "cmap_name"]]
rows = idx.reset_index().merge(cmp, on="pert_id"); rows = rows[(rows.ModelID == "ACH-000019") & rows.cmap_name.isin(ANTI_ER)]
mm = []
for net, d in [("BRCA", "drug_viper_specific_v2_tf_d85")] + [(n, f"drug_viper_specific_v2_tf{TAG}") for n in sorted(pd.read_csv(ROOT / MAP, sep="\t").net.unique())]:
    f = P / d / f"{net}.f32"
    if not f.exists():
        continue
    regs = (P / d / f"{net}.regulators.txt").read_text().split()
    A = np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs))[rows["index"].values]
    for r in ["ESR1", "FOXA1", "GATA3", "XBP1", "PGR"]:
        if r in regs:
            v = A[:, regs.index(r)]; mm.append((net, r, round(float(np.median(v)), 2), round(float((v < 0).mean()), 2)))
mm = pd.DataFrame(mm, columns=["net", "regulator", "median_nes", "share_negative"]); mm.to_csv(OUT / "mcf7_reading.tsv", sep="\t", index=False)

# disease side
a = pd.read_csv(ROOT / f"results/step6/aracne_project_v2_tf{TAG}/group_mr.tsv", sep="\t")
b = pd.read_csv(ROOT / f"results/step6/aracne_project_v2_tf{TAG}pooled/group_mr.tsv", sep="\t")
ds = []
for g in sorted(set(a.group_id)):
    x = a[a.group_id == g].merge(b[b.group_id == g], on="regulator", suffixes=("_c", "_p"))
    rc, rp = x.nes_c.rank(), x.nes_p.rank(); w = (x.nes_c.abs() + x.nes_p.abs()) / 2
    rw = np.cov(rc, rp, aweights=w)[0, 1] / np.sqrt(np.cov(rc, aweights=w) * np.cov(rp, aweights=w))
    kb = {r: (x[x.regulator == r].nes_p.round(2).tolist() or [np.nan])[0] for r in ["ESR1", "FOXA1", "GATA3", "XBP1"]}
    kc = {r: (x[x.regulator == r].nes_c.round(2).tolist() or [np.nan])[0] for r in ["ESR1", "FOXA1", "GATA3", "XBP1"]}
    ds.append(dict(group_id=g, shared_regulators=len(x), weighted_spearman=round(float(rw), 3),
                   **{f"{r}_pooled": kb[r] for r in kb}, **{f"{r}_cell_type": kc[r] for r in kc}))
ds = pd.DataFrame(ds); ds.to_csv(OUT / "disease_side.tsv", sep="\t", index=False)

c, x = M["pooled"], M["cell_type"]
crit = pd.DataFrame([
    ("1 negative controls fall >= 10 points; breast positive controls do not fall",
     f"negative {x['ctl_neg']:.1f} vs {c['ctl_neg']:.1f}; positive {x['ctl_pos']:.1f} vs {c['ctl_pos']:.1f}",
     x["ctl_neg"] <= c["ctl_neg"] - 10 and x["ctl_pos"] >= c["ctl_pos"]),
    ("2 known-use AUC falls <= 0.01 (Hub, ChEMBL)", f"Hub {x['auc_lab_hub']:.3f} vs {c['auc_lab_hub']:.3f}; ChEMBL {x['auc_lab_chembl']:.3f} vs {c['auc_lab_chembl']:.3f}",
     x["auc_lab_hub"] >= c["auc_lab_hub"] - 0.01 and x["auc_lab_chembl"] >= c["auc_lab_chembl"] - 0.01),
    ("3 dominance does not rise", f"{100*x['dominance']:.1f}% vs {100*c['dominance']:.1f}%", x["dominance"] <= c["dominance"])],
    columns=["criterion", "value", "met"])
crit.loc[len(crit)] = ("DECISION", "adopt cell-type networks for breast" if crit.met.all() else "keep the pooled breast network", bool(crit.met.all()))
crit.to_csv(OUT / "decision.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(met.round(4).to_string()); print("\n", mm.to_string(index=False)); print("\n", ds.to_string(index=False))
print("\n", ae.pivot_table(index=["group_id", "cmap_name"], columns="arm", values="spec_mech_pct").round(1).to_string())
print("\n", crit.to_string(index=False))
