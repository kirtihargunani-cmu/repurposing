#!/usr/bin/env python
"""D99 evaluation (rule fixed 2026-11-03 before running, docs/DECISIONS.md). Runs in `celligner` (scipy).
Inputs: scripts/d99/d99_scores.R outputs; production consensus outputs (post D97); D94 `cur` worktree (pre-switch).
V1 on-target pathway sign (pathways MAPK, EGFR, PI3K, JAK-STAT, NFkB, p53; VEGF excluded); V2 Estrogen in MCF7 under
anti-estrogens vs ARACNe ESR1; V3 disease contrasts; V4 ranking agreement pre vs post switch, tier 1-2 PROGENy percentile;
V5 DoRothEA A-B for the D84 ligands (reported). Validation only.
Outputs results/lincs/d99/{v1_on_target.tsv, v2_estrogen.tsv, v3_disease.tsv, v4_agreement.tsv, v5_dorothea.tsv, readings.tsv}
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; D = P / "d99"; RES = ROOT / "results/lincs"; OUT = RES / "d99"; OUT.mkdir(exist_ok=True)
PW_V1 = ["MAPK", "EGFR", "PI3K", "JAK-STAT", "NFkB", "p53"]
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); idx["condition"] = idx.pert_id + "|" + idx.ModelID
ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", dtype=str); ann["moa"] = ann.hub_moa.fillna(ann.lincs_moa).fillna("")
cfg = pd.read_csv(ROOT / "config/progeny_on_target.tsv", sep="\t", comment="#")
ex = pd.read_csv(ROOT / "data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0,
                 usecols=lambda c: c == "Unnamed: 0" or c.split(" (")[0] in ("ESR1", "AR"))
ex.columns = [c.split(" (")[0] for c in ex.columns]
ot = []
for r in cfg.itertuples():
    m = ann[ann.moa.str.contains(r.moa_regex, case=False, regex=True)]
    opp = cfg[(cfg.pathway == r.pathway) & (cfg.direction != r.direction)]
    if len(opp):
        m = m[~m.moa.str.contains("|".join(opp.moa_regex), case=False, regex=True)]
    p = idx[idx.pert_id.isin(m.pert_id)]
    if r.receptor != "-":
        p = p[p.ModelID.map(ex[r.receptor]).fillna(0) >= 1]
    ot.append(p.assign(pathway=r.pathway, direction=r.direction))
ot = pd.concat(ot, ignore_index=True)
cmp = pd.read_csv(RES / "compounds.tsv", sep="\t")[["pert_id", "cmap_name"]].drop_duplicates("pert_id")
ot = ot.merge(cmp, on="pert_id", how="left")


def load(name):
    return pd.read_csv(D / f"progeny_{name}.tsv.gz", sep="\t").set_index("condition")


# V1
rows = []
for name in ["specific", "raw", "T6", "T24", "DL"]:
    S = load(name)
    x = ot[ot.condition.isin(S.index) & ot.pathway.isin(PW_V1)].copy()
    x["score"] = [S.at[c, p] for c, p in zip(x.condition, x.pathway)]
    x["ok"] = np.sign(x.score) == x.direction
    for pw, g in list(x.groupby("pathway")) + [("pooled", x)]:
        rows.append(dict(profiles=name, pathway=pw, n=len(g), share_expected=round(g.ok.mean(), 4), p=binomtest(int(g.ok.sum()), len(g), 0.5).pvalue))
v1 = pd.DataFrame(rows); v1.to_csv(OUT / "v1_on_target.tsv", sep="\t", index=False)
# V2
S = load("specific")
e = ot[(ot.pathway == "Estrogen") & (ot.direction < 0) & (ot.ModelID == "ACH-000019")].copy()
e["estrogen"] = e.condition.map(S.Estrogen)
vd = P / "drug_viper_specific_v2_tf"
for net in ["BRCA", "BRCA_LUMA", "BRCA_LUMB", "BRCA_BAS"]:
    f = vd / f"{net}.f32"
    if f.exists():
        regs = (vd / f"{net}.regulators.txt").read_text().split()
        if "ESR1" in regs:
            A = np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs))
            e[f"ESR1_{net}"] = A[e.index.map(lambda i: idx.index[idx.condition == e.at[i, "condition"]][0]), regs.index("ESR1")]
e.to_csv(OUT / "v2_estrogen.tsv", sep="\t", index=False, float_format="%.3f")
# V3
dz = pd.read_csv(OUT / "progeny_disease.tsv", sep="\t").set_index("condition")
rk = lambda pw, g: int(dz[pw].rank(ascending=False)[g])
v3 = pd.DataFrame([("Estrogen LumA > Basal", dz.at["BRCA-LumA", "Estrogen"] > dz.at["BRCA-Basal", "Estrogen"]),
                   ("Estrogen LumB > Basal", dz.at["BRCA-LumB", "Estrogen"] > dz.at["BRCA-Basal", "Estrogen"]),
                   (f"Androgen PRAD top 5 (rank {rk('Androgen', 'TCGA-PRAD')})", rk("Androgen", "TCGA-PRAD") <= 5),
                   (f"Hypoxia KIRC top 5 (rank {rk('Hypoxia', 'TCGA-KIRC')})", rk("Hypoxia", "TCGA-KIRC") <= 5),
                   (f"MAPK THCA top 10 (rank {rk('MAPK', 'TCGA-THCA')})", rk("MAPK", "TCGA-THCA") <= 10)], columns=["contrast", "holds"])
v3.to_csv(OUT / "v3_disease.tsv", sep="\t", index=False)
# V4
w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
pws = list(dz.columns)
Sd = S[pws]
prof = idx.join(Sd, on="condition")
byd = {k: v for k, v in prof.groupby("pert_id")}
def prog_rev(arm_dir):
    rv = pd.read_csv(arm_dir / "drug_reversal_v2_consensus.tsv", sep="\t", usecols=["group_id", "pert_id", "source"])
    sp = pd.read_csv(arm_dir / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False, usecols=["group_id", "pert_id", "spec_mech", "tier_spec_mech"])
    out = []
    for g, d in rv.groupby("group_id"):
        if g not in dz.index:
            continue
        z = dz.loc[g, pws].values; lines = set(w[(w.group_id == g) & (w.weight > 0)].ModelID)
        vals = []
        for pid, src in zip(d.pert_id, d.source):
            q = byd.get(pid)
            if q is None:
                vals.append(np.nan); continue
            qq = q[q.ModelID.isin(lines)] if src == "direct" else q
            qq = qq if len(qq) else q
            vals.append(-(qq[pws].mean().values @ z) / np.linalg.norm(z))
        out.append(d.assign(prog_rev=vals))
    o = pd.concat(out).merge(sp.groupby(["group_id", "pert_id"]).agg(spec_mech=("spec_mech", "max"), tier=("tier_spec_mech", "min")).reset_index(), on=["group_id", "pert_id"])
    o["prog_pct"] = o.groupby("group_id").prog_rev.rank(pct=True) * 100
    return o
post, pre = prog_rev(RES), prog_rev(ROOT.parent / "repurposing_d94_cur/results/lincs")
agr = []
for name, o in [("post", post), ("pre", pre)]:
    for g, d in o.groupby("group_id"):
        agr.append(dict(arm=name, group_id=g, spearman=d.spec_mech.corr(d.prog_rev, method="spearman"), tier12_prog_pct_median=d[d.tier <= 2].prog_pct.median()))
v4 = pd.DataFrame(agr); v4.to_csv(OUT / "v4_agreement.tsv", sep="\t", index=False, float_format="%.4f")
m4 = v4.groupby("arm")[["spearman", "tier12_prog_pct_median"]].median()
# V5
lig = pd.read_csv(RES / "d84_ligands.tsv", sep="\t"); lig = lig[lig.profile == "specific"][["drug", "tf", "direction", "ModelID"]].drop_duplicates()
lig = lig.merge(cmp.rename(columns={"cmap_name": "drug"}), on="drug"); lig["condition"] = lig.pert_id + "|" + lig.ModelID
T = pd.read_csv(D / "dorothea_ab_specific.tsv.gz", sep="\t").set_index("condition")
lig = lig[lig.condition.isin(T.index) & lig.tf.isin(T.columns)]
lig["score"] = [T.at[c, t] for c, t in zip(lig.condition, lig.tf)]
v5 = pd.DataFrame([dict(rows=len(lig), tfs=lig.tf.nunique(), share_expected=round((np.sign(lig.score) == lig.direction).mean(), 4))])
v5.to_csv(OUT / "v5_dorothea.tsv", sep="\t", index=False)
pooled = v1[(v1.profiles == "specific") & (v1.pathway == "pooled")].iloc[0]
esr = e.estrogen
read = pd.DataFrame([
    ("V1 drug profiles carry on-target pathway signal", f"{100*pooled.share_expected:.1f}% expected sign over {pooled.n} profiles, p = {pooled.p:.2g}",
     "yes" if pooled.share_expected >= 0.55 and pooled.p < 0.05 else "no"),
    ("V2 ESR1 misread is network-side", f"PROGENy Estrogen negative in {(esr < 0).sum()} of {len(esr)}; ARACNe ESR1 (BRCA pooled) positive in "
     f"{(e.get('ESR1_BRCA', pd.Series(dtype=float)) > 0).sum()} of {len(esr)}", "yes" if (esr < 0).mean() > 0.5 and (e.get("ESR1_BRCA", pd.Series([0])) > 0).mean() > 0.5 else "no"),
    ("V3 disease contrasts hold", f"{int(v3.holds.sum())} of {len(v3)}", "all" if v3.holds.all() else "not all: V4 is weak evidence"),
    ("V4 D97 does not lose agreement", f"median Spearman post {m4.at['post', 'spearman']:.3f} vs pre {m4.at['pre', 'spearman']:.3f}",
     "ok" if m4.at["post", "spearman"] >= m4.at["pre", "spearman"] - 0.02 else "D97 flagged for review"),
    ("V4 tiers favour pathway-reversing drugs", f"tier 1-2 median PROGENy-reversal percentile post {m4.at['post', 'tier12_prog_pct_median']:.1f}",
     "yes" if m4.at["post", "tier12_prog_pct_median"] > 50 else "no"),
    ("V5 DoRothEA A-B ligand reading (reported)", f"{100*v5.share_expected[0]:.1f}% of {v5.rows[0]} ligand profiles", "-")], columns=["check", "value", "reading"])
read.to_csv(OUT / "readings.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 140)
print(v1[v1.profiles == "specific"].to_string(index=False)); print(v1[v1.pathway == "pooled"].to_string(index=False))
print(e[["cmap_name", "estrogen"] + [c for c in e.columns if c.startswith("ESR1_")]].round(2).to_string(index=False))
print(v3.to_string(index=False)); print(m4.round(3).to_string()); print(v5.to_string(index=False)); print(read.to_string(index=False))
