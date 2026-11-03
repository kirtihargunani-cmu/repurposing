#!/usr/bin/env python
"""D84 decisions: which regulon set reads TF perturbations correctly on LINCS profiles (runs in `celligner`, scipy).

Rule fixed 2026-11-02 before running (docs/DECISIONS.md D84). Per signature (79n PERTURB=d84): flipped percentile
fpct = pct for knockdowns, 1 - pct for overexpression (0 = read in the expected direction); correct = NES sign as
expected. Main analysis: overexpression signatures with own mRNA up where measured (own_ok); sensitivity: all.
  pass       per TF x set: >= 3 signatures, median fpct <= 0.25, >= 60% correct, and with both kinds no kind < 50% correct
  decision 1 set level, ARACNe v2 TF vs CollecTRI, TFs scored by both with >= 3 signatures: a set replaces the consensus
             only if its median per-TF correct-share difference is > 0 on overexpression and on knockdown separately,
             paired Wilcoxon p < 0.05 on the pooled shares, and the pooled median difference >= 0.05
  decision 2 per TF: >= 6 signatures in >= 2 lines, one set passes and the other < 50% correct -> that set
  decision 3 flags: 'not_supported' (>= 6 signatures, neither set passes), 'untested' (no usable signature),
             'insufficient' (1-5 signatures); with the 79t agreement class (networks where the TF disagrees)
  report 4   ligands: drugs annotated to act on a TF with a known direction, in lines expressing it (DepMap
             log2(TPM + 1) >= 1), ARACNe lineage-matched networks else all (as 79n); labelled as favouring CollecTRI
  report 5   the same ligand check on raw drug profiles (79d PROFILE=raw), and the anti-estrogen ESR1 reading in MCF7
Outputs: results/lincs/d84_tf.tsv (TF x set x analysis), d84_choice.tsv (per TF: choice, flag), d84_summary.tsv,
         d84_ligands.tsv
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
SETS = ["aracne_v2_tf", "collectri"]
MIN_SIGS, MAX_FPCT, MIN_CORRECT, MIN_KIND = 3, 0.25, 0.6, 0.5
TF_SIGS, TF_LINES, OTHER_MAX = 6, 2, 0.5
P_MAX, MIN_DIFF = 0.05, 0.05
MIN_EXPR = 1.0


def tf_table(sg):
    g = sg.groupby(["cmap_name", "set"])
    t = g.agg(n_sigs=("sig_id", "size"), n_lines=("cell_iname", "nunique"), n_kd=("direction", lambda d: int((d < 0).sum())),
              n_oe=("direction", lambda d: int((d > 0).sum())), median_fpct=("fpct", "median"),
              frac_correct=("correct", "mean")).reset_index()
    kind = sg.groupby(["cmap_name", "set", "direction"]).correct.mean().unstack("direction")
    kind.columns = ["correct_kd" if c < 0 else "correct_oe" for c in kind.columns]
    t = t.merge(kind.reset_index(), on=["cmap_name", "set"], how="left")
    both = (t.n_kd > 0) & (t.n_oe > 0)
    kinds_ok = ~both | ((t.correct_kd >= MIN_KIND) & (t.correct_oe >= MIN_KIND))
    t["pass"] = (t.n_sigs >= MIN_SIGS) & (t.median_fpct <= MAX_FPCT) & (t.frac_correct >= MIN_CORRECT) & kinds_ok
    return t


def set_level(sg, t):
    w = t[t.set.isin(SETS)].pivot(index="cmap_name", columns="set")
    elig = w.index[(w["n_sigs"][SETS] >= MIN_SIGS).all(axis=1)]
    a, c = (w["frac_correct"].loc[elig, s] for s in SETS)
    d = a - c
    out = {"tfs": len(elig), "median_diff_pooled": float(d.median()), "mean_diff_pooled": float(d.mean()),
           "aracne_better": int((d > 0).sum()), "collectri_better": int((d < 0).sum()), "ties": int((d == 0).sum())}
    out["wilcoxon_p"] = float(wilcoxon(a, c, zero_method="wilcox").pvalue) if (d != 0).any() else 1.0
    for kind, col in [("oe", "correct_oe"), ("kd", "correct_kd")]:
        x = w[col].loc[elig].dropna()
        dk = x[SETS[0]] - x[SETS[1]]
        out[f"tfs_{kind}"] = len(dk); out[f"median_diff_{kind}"] = float(dk.median()) if len(dk) else np.nan
        out[f"mean_diff_{kind}"] = float(dk.mean()) if len(dk) else np.nan
    win = None
    for s, sgn in [(SETS[0], 1), (SETS[1], -1)]:
        if (sgn * out["median_diff_oe"] > 0 and sgn * out["median_diff_kd"] > 0 and out["wilcoxon_p"] < P_MAX
                and sgn * out["median_diff_pooled"] >= MIN_DIFF):
            win = s
    out["decision_1"] = win or "keep_consensus"
    return out


def per_tf(t, agree):
    w = t[t.set.isin(SETS)].pivot(index="cmap_name", columns="set")
    rows = []
    for tf in w.index:
        n = int(w["n_sigs"].loc[tf].max()); nl = int(w["n_lines"].loc[tf].max())
        ps = {s: bool(w["pass"].loc[tf, s]) if pd.notna(w["pass"].loc[tf, s]) else False for s in SETS}
        fc = {s: w["frac_correct"].loc[tf, s] for s in SETS}
        choice = "consensus"
        if n >= TF_SIGS and nl >= TF_LINES:
            for s, o in [(SETS[0], SETS[1]), (SETS[1], SETS[0])]:
                if ps[s] and pd.notna(fc[o]) and fc[o] < OTHER_MAX:
                    choice = s
        flag = "tested" if n >= TF_SIGS else "insufficient"
        if n >= TF_SIGS and not any(ps.values()):
            flag = "not_supported"
        rows.append((tf, n, nl, ps[SETS[0]], ps[SETS[1]], fc[SETS[0]], fc[SETS[1]], choice, flag))
    ch = pd.DataFrame(rows, columns=["regulator", "n_sigs", "n_lines", "pass_aracne", "pass_collectri",
                                     "correct_aracne", "correct_collectri", "choice", "flag"])
    allreg = agree[["regulator", "n_nets", "n_disagree"]]
    ch = allreg.merge(ch, on="regulator", how="outer")
    ch["flag"] = ch.flag.fillna("untested"); ch["choice"] = ch.choice.fillna("consensus")
    return ch


def drug_nes(profile, regulon, net=None):
    d = PROC / (f"drug_viper_{profile}_v2_tf" if regulon == "aracne" else f"drug_viper_{profile}_collectri")
    net = net or "ALL"
    regs = (d / f"{net}.regulators.txt").read_text().split()
    n = len(pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"))
    return pd.DataFrame(np.fromfile(d / f"{net}.f32", dtype=np.float32).reshape(n, len(regs)), columns=regs)


def ligands(profile):
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", dtype=str)
    col = drug_nes(profile, "collectri")
    tfs = set(col.columns)
    rows = []
    for _, r in ann.iterrows():
        t = set(filter(None, "|".join(x for x in [r.lincs_target, r.hub_target] if isinstance(x, str)).split("|"))) & tfs
        moa = " ".join(x for x in [r.hub_moa, r.lincs_moa] if isinstance(x, str)).lower()
        if not t or not moa:
            continue
        ago = bool(re.search(r"agonist|activator", moa.replace("antagonist", "")))
        ant = bool(re.search(r"antagonist|inhibitor|degrader|blocker", moa))
        if ago == ant:   # unknown or both: no clean direction
            continue
        for g in t:
            rows.append((r.pert_id, r.cmap_name, g, 1 if ago else -1))
    lg = pd.DataFrame(rows, columns=["pert_id", "drug", "tf", "direction"]).drop_duplicates()
    lg = lg.merge(idx.reset_index().rename(columns={"index": "row"}), on="pert_id")
    ex = pd.read_csv(ROOT / "data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv", index_col=0)
    ex.columns = [c.split(" (")[0] for c in ex.columns]
    lg["expr"] = [ex.at[m, g] if (m in ex.index and g in ex.columns) else np.nan for m, g in zip(lg.ModelID, lg.tf)]
    lg = lg[lg.expr >= MIN_EXPR].copy()
    # ARACNe: lineage-matched networks (as 79n)
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#")
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
    net_of = dict(zip(grp.group_id, np.where(grp.cohort.isin(["TCGA-COAD", "TCGA-READ"]), "COADREAD",
                                              grp.cohort.str.replace(r"^[^-]+-", "", regex=True))))
    ad = PROC / f"drug_viper_{profile}_v2_tf"
    nets = sorted(p.stem for p in ad.glob("*.f32"))
    ln = {(l, net_of[gid]) for gid, ls in zip(keys.group_id, keys.lineages) if gid in net_of and isinstance(ls, str)
          for l in ls.split(";")}
    lin = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", usecols=["ModelID", "OncotreeLineage"]).set_index("ModelID").OncotreeLineage
    A = {n: drug_nes(profile, "aracne", n) for n in nets}
    out = []
    for _, r in lg.iterrows():
        c = col.at[r.row, r.tf]
        have = [n for n in nets if r.tf in A[n].columns]
        m = [n for n in have if (lin.get(r.ModelID), n) in ln]
        use = m or have
        a = float(np.mean([A[n].at[r.row, r.tf] for n in use])) if use else np.nan
        out.append((profile, r.drug, r.tf, r.direction, r.ModelID, a, bool(m), c))
    o = pd.DataFrame(out, columns=["profile", "drug", "tf", "direction", "ModelID", "nes_aracne", "aracne_matched", "nes_collectri"])
    o = o.dropna(subset=["nes_aracne", "nes_collectri"])   # compare only where both sets score the TF
    o["correct_aracne"] = np.sign(o.nes_aracne) == o.direction
    o["correct_collectri"] = np.sign(o.nes_collectri) == o.direction
    return o


def main():
    sg = pd.read_csv(RES / "perturbation_benchmark_d84_signatures.tsv.gz", sep="\t")
    sg["fpct"] = np.where(sg.direction < 0, sg.pct, 1 - sg.pct)
    sg["correct"] = np.sign(sg.nes) == sg.direction
    agree = pd.read_csv(RES / "regulon_agreement_v2_tf_summary.tsv", sep="\t")
    oe = sg[sg.direction > 0]
    print(f"overexpression signatures scored: {oe.sig_id.nunique()}; own mRNA measured "
          f"{oe.drop_duplicates('sig_id').own_z.notna().sum()}, own_ok {oe.drop_duplicates('sig_id').own_ok.sum()}")
    tabs, summ = [], []
    for analysis, d in [("main", sg[sg.own_ok.astype(bool)]), ("no_own_check", sg)]:
        t = tf_table(d); t.insert(0, "analysis", analysis); tabs.append(t)
        s = set_level(d, t); s["analysis"] = analysis
        for st in SETS + ["dorothea_ABC"]:
            x = d[d.set == st]
            s[f"{st}_signatures"] = len(x); s[f"{st}_correct_kd"] = round(x[x.direction < 0].correct.mean(), 3)
            s[f"{st}_correct_oe"] = round(x[x.direction > 0].correct.mean(), 3)
            s[f"{st}_auc_kd"] = round(1 - x[x.direction < 0].fpct.mean(), 3); s[f"{st}_auc_oe"] = round(1 - x[x.direction > 0].fpct.mean(), 3)
            tt = t[(t.set == st) & (t.n_sigs >= MIN_SIGS)]
            s[f"{st}_tfs_tested"] = len(tt); s[f"{st}_tfs_pass"] = int(tt["pass"].sum())
        summ.append(s)
        if analysis == "main":
            ch = per_tf(t, agree)
    tf = pd.concat(tabs, ignore_index=True)
    tf.to_csv(RES / "d84_tf.tsv", sep="\t", index=False, float_format="%.4g")
    ch.to_csv(RES / "d84_choice.tsv", sep="\t", index=False, float_format="%.4g")
    sm = pd.DataFrame(summ).set_index("analysis").T
    sm.to_csv(RES / "d84_summary.tsv", sep="\t")
    print(sm.to_string())
    print("\nper-TF choice (decision 2):", ch.choice.value_counts().to_dict())
    print("flags (decision 3):", ch.flag.value_counts().to_dict())
    print(ch[ch.choice != "consensus"].round(2).to_string(index=False))
    print("\nnot_supported x 79t disagreement (networks disagreeing, median):",
          ch.groupby("flag").n_disagree.median().round(1).to_dict())

    lig = []
    for prof in ["specific", "raw"]:
        if (PROC / f"drug_viper_{prof}_v2_tf").exists() and (PROC / f"drug_viper_{prof}_collectri").exists():
            lig.append(ligands(prof))
    lig = pd.concat(lig, ignore_index=True)
    lig.to_csv(RES / "d84_ligands.tsv", sep="\t", index=False, float_format="%.4g")
    ls = lig.groupby(["profile", "tf", "direction"]).agg(profiles=("drug", "size"), drugs=("drug", "nunique"),
                                                        correct_aracne=("correct_aracne", "mean"),
                                                        correct_collectri=("correct_collectri", "mean")).reset_index()
    print("\nligands (report only; favours CollecTRI), TFs with >= 5 profiles:\n",
          ls[ls.profiles >= 5].round(2).to_string(index=False))
    print("\nligands overall:\n", lig.groupby("profile")[["correct_aracne", "correct_collectri"]].mean().round(3).to_string())
    ae = lig[(lig.tf == "ESR1") & (lig.direction < 0) & (lig.ModelID == "ACH-000019")]
    print("\nanti-estrogens, ESR1 in MCF7 (median NES):\n",
          ae.groupby("profile")[["nes_aracne", "nes_collectri"]].median().round(2).to_string())


if __name__ == "__main__":
    main()
