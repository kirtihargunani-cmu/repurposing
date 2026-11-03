#!/usr/bin/env python
"""D85: raw vs line-specific drug profiles (rule fixed 2026-11-02 before running, docs/DECISIONS.md). Runs in `celligner`.

Both arms rerun with the current code and regulator filter (RUN_TAG=_d85): 79d (v2 TF + CollecTRI), 79e and 79i
(REGULON=consensus) with PROFILE=specific / raw; 79n PERTURB=d84 PROJECT=1 for the projected signatures.
Checks (arm raw vs arm specific):
  1 TF ground truth read as the ranking reads it (consensus NES = mean of ARACNe and CollecTRI where both, else one):
    pooled share in the expected direction over the D84 main-analysis signatures (own_ok; specific arm: line generic
    drug response projected out, raw arm: as measured; lines without drug profiles excluded from both) plus the D84
    ligand profiles (drug activity of each arm)
  2 known-drug controls (79i split_controls): median spec_mech percentile of the positive and of the negative pairs
  3 growth-arrest dominance: share of drugs with rev > 0 in >= 80% of the groups where they are scored, over drugs
    scored in >= 33 of the 41 groups (80%; interpretation fixed before the raw arm was looked at)
Adopt raw if: (1) rises >= 2 points; (2) positive median falls <= 5 points and negative median rises <= 5 points;
(3) dominance rises <= 10 points. Else keep specific.
Outputs: results/lincs/d85/{decision.tsv, check1.tsv, controls.tsv, antier.tsv, ranking_agreement.tsv, esr1_mcf7.tsv}
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUTD = RES / "d85"
TAG = "_d85"
MIN_GAIN, CTL_TOL, DOM_TOL, DOM_FRAC, DOM_MIN_GROUPS = 0.02, 5.0, 0.10, 0.8, 33
ARMS = {"specific": f"_v2_consensus{TAG}", "raw": f"_v2_consensus{TAG}_raw"}
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"]

spec = importlib.util.spec_from_file_location("d84", ROOT / "scripts/79v_d84_decisions.py")
d84 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d84)


def drug_nes(profile, regulon, net=None):   # 79v's reader, pointed at the D85 reruns
    d = PROC / (f"drug_viper_{profile}_v2_tf{TAG}" if regulon == "aracne" else f"drug_viper_{profile}_collectri{TAG}")
    net = net or "ALL"
    regs = (d / f"{net}.regulators.txt").read_text().split()
    n = len(pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"))
    return pd.DataFrame(np.fromfile(d / f"{net}.f32", dtype=np.float32).reshape(n, len(regs)), columns=regs)


d84.drug_nes = drug_nes


def cons(a, c):
    return np.where(np.isnan(a), c, np.where(np.isnan(c), a, (a + c) / 2))


def sig_table(path, ids=None):
    s = pd.read_csv(path, sep="\t")
    s = s[s.own_ok]
    w = s.pivot_table(index=["sig_id", "direction", "cmap_name", "cell_iname"], columns="set", values="nes").reset_index()
    if ids is not None:
        w = w[w.sig_id.isin(ids)]
    w["consensus"] = cons(w.get("aracne_v2_tf", np.nan * w.direction).values, w.get("collectri", np.nan * w.direction).values)
    return w


def main():
    OUTD.mkdir(exist_ok=True)
    # check 1: signatures
    proj = sig_table(RES / "perturbation_benchmark_d85_projected_signatures.tsv.gz")
    raw = sig_table(RES / "perturbation_benchmark_d84_signatures.tsv.gz", set(proj.sig_id))
    proj = proj[proj.sig_id.isin(set(raw.sig_id))]
    # check 1: ligands (79v logic on each arm's drug activity)
    lig = pd.concat([d84.ligands(p) for p in ARMS], ignore_index=True)
    lig["consensus"] = cons(lig.nes_aracne.values, lig.nes_collectri.values)
    keys = ["drug", "tf", "direction", "ModelID"]
    common = set.intersection(*(set(map(tuple, lig[lig.profile == p][keys].values)) for p in ARMS))
    lig = lig[lig[keys].apply(tuple, axis=1).isin(common)]
    lig.to_csv(OUTD / "ligands.tsv", sep="\t", index=False, float_format="%.4g")
    rows = []
    for arm, sg in [("specific", proj), ("raw", raw)]:
        lg = lig[lig.profile == arm]
        for col, lcol in [("consensus", "consensus"), ("aracne_v2_tf", "nes_aracne"), ("collectri", "nes_collectri")]:
            a = sg.dropna(subset=[col]); b = lg.dropna(subset=[lcol])
            ca = np.sign(a[col]) == a.direction; cb = np.sign(b[lcol]) == b.direction
            rows.append(dict(arm=arm, set=col, signatures=len(a), correct_signatures=ca.mean(), ligand_profiles=len(b),
                             correct_ligands=cb.mean(), pooled_n=len(a) + len(b), correct_pooled=(ca.sum() + cb.sum()) / (len(a) + len(b)),
                             correct_kd=ca[a.direction < 0].mean(), correct_oe=ca[a.direction > 0].mean()))
    c1 = pd.DataFrame(rows)
    c1.to_csv(OUTD / "check1.tsv", sep="\t", index=False, float_format="%.4f")
    g = c1[c1.set == "consensus"].set_index("arm").correct_pooled
    # check 2: controls
    ctl = pd.concat([pd.read_csv(RES / f"split_controls{s}.tsv", sep="\t").assign(arm=a) for a, s in ARMS.items()])
    ctl.to_csv(OUTD / "controls.tsv", sep="\t", index=False)
    cm = ctl.groupby(["arm", "expected"]).agg(pairs=("drug", "size"), median_spec_mech=("spec_mech", "median"),
                                              median_rev=("rev", "median")).reset_index()
    pos = cm[cm.expected == "high"].set_index("arm").median_spec_mech; neg = cm[cm.expected == "low"].set_index("arm").median_spec_mech
    # check 3: dominance
    dom = {}
    for a, s in ARMS.items():
        r = pd.read_csv(RES / f"drug_reversal{s}.tsv", sep="\t", usecols=["group_id", "pert_id", "rev"])
        t = r.groupby("pert_id").agg(n=("group_id", "nunique"), pos=("rev", lambda x: (x > 0).mean()))
        t = t[t.n >= DOM_MIN_GROUPS]
        dom[a] = ((t.pos >= DOM_FRAC).mean(), len(t))
    # reported: anti-estrogens, ranking agreement, split overview
    sp = {a: pd.read_csv(RES / f"drug_reversal_split{s}.tsv", sep="\t", low_memory=False) for a, s in ARMS.items()}
    ae = pd.concat([d[d.group_id.isin(["BRCA-LumA", "BRCA-LumB", "BRCA-Basal"]) & d.cmap_name.isin(ANTI_ER)]
                    .groupby(["group_id", "cmap_name"])[["rev_pct", "spec_pct", "spec_mech_pct", "prolif_share"]].max()
                    .reset_index().assign(arm=a) for a, d in sp.items()])
    ae.to_csv(OUTD / "antier.tsv", sep="\t", index=False)
    m = sp["specific"][["group_id", "pert_id", "spec_mech"]].merge(sp["raw"][["group_id", "pert_id", "spec_mech"]],
                                                                  on=["group_id", "pert_id"], suffixes=("_s", "_r"))
    ag = []
    for gid, d in m.groupby("group_id"):
        ts = set(d.nlargest(20, "spec_mech_s").pert_id); tr = set(d.nlargest(20, "spec_mech_r").pert_id)
        ag.append(dict(group_id=gid, drugs=len(d), spearman=d.spec_mech_s.corr(d.spec_mech_r, method="spearman"), top20_overlap=len(ts & tr)))
    ag = pd.DataFrame(ag); ag.to_csv(OUTD / "ranking_agreement.tsv", sep="\t", index=False, float_format="%.3f")
    ov = pd.concat([pd.read_csv(RES / f"split_overview{s}.tsv", sep="\t").set_index("metric").value.rename(a) for a, s in ARMS.items()], axis=1)
    ov.to_csv(OUTD / "split_overview.tsv", sep="\t")
    esr = lig[(lig.tf == "ESR1") & (lig.direction < 0) & (lig.ModelID == "ACH-000019")]
    esr.to_csv(OUTD / "esr1_mcf7.tsv", sep="\t", index=False, float_format="%.3f")

    crit = pd.DataFrame([
        ("1 correct share rises >= 2 points", f"raw {100*g['raw']:.1f}% vs specific {100*g['specific']:.1f}% ({100*(g['raw']-g['specific']):+.1f} points)",
         g["raw"] - g["specific"] >= MIN_GAIN),
        ("2a positive controls fall <= 5 points", f"median spec_mech pct raw {pos['raw']:.1f} vs specific {pos['specific']:.1f}", pos["raw"] >= pos["specific"] - CTL_TOL),
        ("2b negative controls rise <= 5 points", f"median spec_mech pct raw {neg['raw']:.1f} vs specific {neg['specific']:.1f}", neg["raw"] <= neg["specific"] + CTL_TOL),
        ("3 dominance rises <= 10 points", f"raw {100*dom['raw'][0]:.1f}% vs specific {100*dom['specific'][0]:.1f}% of {dom['raw'][1]} drugs",
         dom["raw"][0] - dom["specific"][0] <= DOM_TOL)], columns=["criterion", "value", "met"])
    crit.loc[len(crit)] = ("DECISION", "adopt raw profiles" if crit.met.all() else "keep line-specific profiles", crit.met.all())
    crit.to_csv(OUTD / "decision.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("== check 1\n", c1.round(4).to_string(index=False))
    print("\n== controls\n", cm.round(1).to_string(index=False))
    print("\n== dominance", {k: (round(v[0], 4), v[1]) for k, v in dom.items()})
    print("\n== anti-estrogens\n", ae.pivot_table(index=["group_id", "cmap_name"], columns="arm", values=["rev_pct", "spec_mech_pct"]).round(1).to_string())
    print("\n== ESR1 in MCF7, anti-estrogens\n", esr[["profile", "drug", "nes_aracne", "nes_collectri", "consensus"]].round(2).to_string(index=False))
    print("\n== ranking agreement (spec_mech, raw vs specific): median Spearman %.3f, median top-20 overlap %.0f, min Spearman %.3f"
          % (ag.spearman.median(), ag.top20_overlap.median(), ag.spearman.min()))
    print("\n== split overview\n", ov.to_string())
    print("\n== DECISION\n", crit.to_string(index=False))


if __name__ == "__main__":
    main()
