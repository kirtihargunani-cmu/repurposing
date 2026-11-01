#!/usr/bin/env python
"""D83: signaling-only target lists (sig1, sig2) vs the current lists rerun the same way (base). Rule fixed 2026-11-01.

Inputs (per arm a in base, sig1, sig2; ARACNE_SET=tf_<a> runs of 79e / 79i on the consensus):
  results/lincs/drug_reversal_split_v2_consensus_tf_<a>.tsv   ranking with spec_mech (diseases)
  results/lincs/knockdown_benchmark_compare_d83.tsv           79n COMPARE_SETS run
  main: results/step6/aracne_project_v2_tf_<a>/controls.tsv   51_compare_mr.R known-biology checks
  main: results/step5/identity_edges_v2_tf/signaling_lists.tsv regulators lost
Checks: (1) anti-estrogens (ANTI_ER) median spec_mech percentile and median rev percentile in LumA, LumB (Basal shown);
        (2) knockdown share read as lower activity (frac_neg); (3) controls pass / wrong-direction counts.
Adopt an arm if vs base: (1) anti-estrogen median spec_mech percentile rises in both LumA and LumB; (2) frac_neg not lower
by more than 0.02; (3) controls passing do not fall and wrong-direction ones do not rise.
Reported: drug ranking agreement per group (Spearman of spec_mech and of rev, arm vs base), regulators lost.
Also: base vs the production run (tf), to show what the shared-seed rerun alone changes.
Outputs: results/lincs/d83_compare_{antier,agreement,checks}.tsv
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
MAIN = Path.home() / "projects/repurposing"
ANTI_ER = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"]
ARMS = ["base", "sig1", "sig2"]
KD_TOL = 0.02


def load(arm):
    f = RES / ("drug_reversal_split_v2_consensus.tsv" if arm == "tf" else f"drug_reversal_split_v2_consensus_tf_{arm}.tsv")
    return pd.read_csv(f, sep="\t", low_memory=False)


def main():
    d = {a: load(a) for a in ["tf"] + ARMS}
    rows = []
    for a, x in d.items():
        for g in ["BRCA-LumA", "BRCA-LumB", "BRCA-Basal"]:
            y = x[(x.group_id == g) & x.cmap_name.isin(ANTI_ER)]
            rows.append(dict(arm=a, group_id=g, n_rows=len(y), median_spec_mech_pct=round(y.spec_mech_pct.median(), 1),
                             median_rev_pct=round(y.rev_pct.median(), 1), median_rev_mech_pct=round(y.rev_mech_pct.median(), 1)))
    ae = pd.DataFrame(rows)
    ae.to_csv(RES / "d83_compare_antier.tsv", sep="\t", index=False)

    ag = []
    for a in ["tf", "sig1", "sig2"]:
        b, x = d["base"], d[a]
        for g in sorted(set(b.group_id)):
            m = b[b.group_id == g].merge(x[x.group_id == g], on="pert_id", suffixes=("_b", "_x"))
            ag.append(dict(arm=a, group_id=g, n=len(m), rho_spec_mech=round(m.spec_mech_b.rank().corr(m.spec_mech_x.rank()), 3),
                           rho_rev=round(m.rev_b.rank().corr(m.rev_x.rank()), 3)))
    ag = pd.DataFrame(ag)
    ag.to_csv(RES / "d83_compare_agreement.tsv", sep="\t", index=False)

    kd = pd.read_csv(RES / "knockdown_benchmark_compare_d83.tsv", sep="\t").set_index("set")
    ctl = {}
    for a in ARMS:
        c = pd.read_csv(MAIN / f"results/step6/aracne_project_v2_tf_{a}/controls.tsv", sep="\t")
        ctl[a] = dict(pass_=int((c.verdict == "correct, passes").sum()), wrong=int((c.verdict == "wrong direction").sum()), n=len(c),
                      verdicts=c.verdict.value_counts().to_dict())
    sl = pd.read_csv(MAIN / "results/step5/identity_edges_v2_tf/signaling_lists.tsv", sep="\t")

    ck = []
    sb = ae[ae.arm == "base"].set_index("group_id")
    for a in ["sig1", "sig2"]:
        sa = ae[ae.arm == a].set_index("group_id")
        c1 = all(sa.loc[g, "median_spec_mech_pct"] > sb.loc[g, "median_spec_mech_pct"] for g in ["BRCA-LumA", "BRCA-LumB"])
        c2 = kd.loc[f"tf_{a}", "frac_neg"] >= kd.loc["tf_base", "frac_neg"] - KD_TOL
        c3 = ctl[a]["pass_"] >= ctl["base"]["pass_"] and ctl[a]["wrong"] <= ctl["base"]["wrong"]
        lost = sl[sl.arm == a]
        ck.append(dict(arm=a,
                       antier_LumA=f"{sb.loc['BRCA-LumA','median_spec_mech_pct']} -> {sa.loc['BRCA-LumA','median_spec_mech_pct']}",
                       antier_LumB=f"{sb.loc['BRCA-LumB','median_spec_mech_pct']} -> {sa.loc['BRCA-LumB','median_spec_mech_pct']}",
                       check1_antier=c1,
                       knockdown_frac_neg=f"{kd.loc['tf_base','frac_neg']} -> {kd.loc[f'tf_{a}','frac_neg']}",
                       knockdown_tfs_pass=f"{kd.loc['tf_base','tfs_pass']} -> {kd.loc[f'tf_{a}','tfs_pass']}",
                       check2_knockdown=bool(c2),
                       controls_pass=f"{ctl['base']['pass_']} -> {ctl[a]['pass_']}",
                       controls_wrong=f"{ctl['base']['wrong']} -> {ctl[a]['wrong']}",
                       check3_controls=c3,
                       adopt=bool(c1 and c2 and c3),
                       median_rho_spec_mech=ag[ag.arm == a].rho_spec_mech.median(),
                       min_rho_spec_mech=ag[ag.arm == a].rho_spec_mech.min(),
                       regulators_lost=int(lost.n_lost.sum()), edges_removed_pct=round(100 * (1 - lost.edges_after.sum() / lost.edges.sum()), 2)))
    ck = pd.DataFrame(ck)
    ck.to_csv(RES / "d83_compare_checks.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(ae.to_string(index=False)); print(ck.T.to_string())
    print("noise floor (production tf vs base rerun): median rho spec_mech",
          ag[ag.arm == "tf"].rho_spec_mech.median(), "min", ag[ag.arm == "tf"].rho_spec_mech.min())
    for a in ["sig1", "sig2"]:
        print(f"{a}: lowest agreement groups:", ag[ag.arm == a].nsmallest(6, "rho_spec_mech")[["group_id", "rho_spec_mech", "rho_rev"]].to_string(index=False))


if __name__ == "__main__":
    main()
