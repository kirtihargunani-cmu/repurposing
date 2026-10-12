#!/usr/bin/env python
"""Baseline summary of the first drug-reversal run: 41 Step 4 groups, tumor vs normal, no subgroups (runs in `celligner`).

Purpose (user, 2026-10-12): the reference to compare stratum results (Step 10 strata: COAD CMS-like groups, READ
clusters, breast axis ends) against. Strata signatures are 0.97-0.99 similar to their parents
(main:docs/notes/subgrouping_for_drug_side.md), so the key baseline is how similar the drug rankings of related groups
already are (breast subtypes with each other, COAD with READ), next to per-group statistics.
Input: results/lincs/drug_reversal.tsv, drug_reversal_summary.tsv, drug_reversal_controls.tsv (scripts/79e).
Outputs (results/lincs/reversal_baseline/):
  per_disease.tsv        per group: MRs, drugs, direct / borrowed, score spread, tiers by rev / spec / both,
                         approved and multi-target in tier 1, rank agreement rev vs spec, top tier-1 drugs
  drug_generality.tsv    per drug: number of groups where it is top 10% by rev / by spec (of groups it is scored in)
  pair_similarity.tsv    per pair of groups: Spearman of rev and of spec over shared drugs, Jaccard of top-10% lists,
                         relation (sibling = same parent family: BRCA subtypes, COAD-READ, other same tissue key; other)
  overview.tsv           headline numbers, incl. sibling vs other similarity and known-drug control percentiles
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "reversal_baseline"
TOP_PCT = 10
FAMILY = {"BRCA-Basal": "BRCA", "BRCA-Her2": "BRCA", "BRCA-LumA": "BRCA", "BRCA-LumB": "BRCA",
          "TCGA-COAD": "COADREAD", "TCGA-READ": "COADREAD"}   # the parents of the Step 10 strata


def main():
    OUT.mkdir(exist_ok=True)
    x = pd.read_csv(RES / "drug_reversal.tsv", sep="\t", low_memory=False)
    summ = pd.read_csv(RES / "drug_reversal_summary.tsv", sep="\t").set_index("group_id")
    ctl = pd.read_csv(RES / "drug_reversal_controls.tsv", sep="\t")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id").key
    x["top_rev"] = (x.rev > 0) & (x.rev_pct >= 100 - TOP_PCT)
    x["top_spec"] = (x.rev > 0) & (x.spec_pct >= 100 - TOP_PCT)

    # per disease
    rows = []
    for g, d in x.groupby("group_id"):
        t1r, t1s = d[d.tier_rev == 1], d[d.tier_spec == 1]
        both = d[(d.tier_rev == 1) & (d.tier_spec == 1)].sort_values("rev", ascending=False)
        rows.append(dict(
            group_id=g, key=keys[g], family=FAMILY.get(g, ""), n_mrs=int(d.n_mrs.iloc[0]), low_mrs=bool(d.low_mrs.iloc[0]),
            n_drugs=len(d), n_direct=int((d.source == "direct").sum()), n_borrowed=int((d.source == "borrowed").sum()),
            n_pass_gate=int(d.passes_gate.sum()), n_approved_pass=int((d.approved & d.passes_gate).sum()),
            rev_median=round(d.rev.median(), 2), rev_iqr=round(d.rev.quantile(.75) - d.rev.quantile(.25), 2),
            frac_rev_pos=round((d.rev > 0).mean(), 3),
            n_tier1_rev=len(t1r), n_tier2_rev=int((d.tier_rev == 2).sum()),
            n_tier1_spec=len(t1s), n_tier2_spec=int((d.tier_spec == 2).sum()), n_tier1_both=len(both),
            tier1_rev_multi_target=int(t1r.multi_target.sum()), tier1_spec_multi_target=int(t1s.multi_target.sum()),
            tier1_rev_direct=int((t1r.source == "direct").sum()), tier1_spec_direct=int((t1s.source == "direct").sum()),
            rho_rev_spec=round(d.rev.corr(d.spec, method="spearman"), 3),
            top_tier1_rev=", ".join(t1r.sort_values("rev", ascending=False).cmap_name.head(5)),
            top_tier1_spec=", ".join(t1s.sort_values("spec", ascending=False).cmap_name.head(5)),
            top_tier1_both=", ".join(both.cmap_name.head(5))))
    per = pd.DataFrame(rows)
    per.to_csv(OUT / "per_disease.tsv", sep="\t", index=False)

    # drug generality
    gen = x.groupby("pert_id").agg(cmap_name=("cmap_name", "first"), approved=("approved", "first"),
                                   hub_moa=("hub_moa", "first"), n_groups=("group_id", "nunique"),
                                   n_top_rev=("top_rev", "sum"), n_top_spec=("top_spec", "sum"))
    gen["frac_top_rev"] = (gen.n_top_rev / gen.n_groups).round(3)
    gen["frac_top_spec"] = (gen.n_top_spec / gen.n_groups).round(3)
    gen.sort_values(["n_top_rev", "n_top_spec"], ascending=False).to_csv(OUT / "drug_generality.tsv", sep="\t")

    # pairwise similarity of rankings
    R = x.pivot_table(index="pert_id", columns="group_id", values="rev")
    S = x.pivot_table(index="pert_id", columns="group_id", values="spec")
    top = {g: set(d.pert_id[d.top_rev]) for g, d in x.groupby("group_id")}
    tops = {g: set(d.pert_id[d.top_spec]) for g, d in x.groupby("group_id")}
    low = set(per.group_id[per.low_mrs])
    pr = []
    for a, b in combinations(sorted(R.columns), 2):
        fa, fb = FAMILY.get(a), FAMILY.get(b)
        rel = "parent_family" if fa and fa == fb else "same_tissue_key" if keys[a] == keys[b] else "other"
        jac = lambda s, t: len(s & t) / len(s | t) if s | t else np.nan
        pr.append(dict(group_a=a, group_b=b, relation=rel, low_mrs=a in low or b in low,
                       n_shared=int(R[[a, b]].dropna().shape[0]),
                       rho_rev=round(R[a].corr(R[b], method="spearman"), 3),
                       rho_spec=round(S[a].corr(S[b], method="spearman"), 3),
                       jaccard_top_rev=round(jac(top[a], top[b]), 3), jaccard_top_spec=round(jac(tops[a], tops[b]), 3)))
    pr = pd.DataFrame(pr)
    pr.to_csv(OUT / "pair_similarity.tsv", sep="\t", index=False)

    # overview
    ok = pr[~pr.low_mrs]
    ov = [("groups", len(per)), ("groups_low_mrs", int(per.low_mrs.sum())),
          ("drug_group_rows", len(x)), ("drugs", x.pert_id.nunique()),
          ("mrs_per_group_median", int(per.n_mrs.median())),
          ("tier1_rev_per_group_median", int(per.n_tier1_rev.median())),
          ("tier1_spec_per_group_median", int(per.n_tier1_spec.median())),
          ("tier1_both_per_group_median", int(per.n_tier1_both.median())),
          ("rho_rev_spec_within_group_median", round(per.rho_rev_spec.median(), 3)),
          ("drugs_top10_rev_in_ge20_groups", int((gen.n_top_rev >= 20).sum())),
          ("drugs_top10_spec_in_ge20_groups", int((gen.n_top_spec >= 20).sum()))]
    for rel, d in ok.groupby("relation"):
        ov += [(f"{rel}_pairs", len(d)), (f"{rel}_rho_rev_median", round(d.rho_rev.median(), 3)),
               (f"{rel}_rho_spec_median", round(d.rho_spec.median(), 3)),
               (f"{rel}_jaccard_top_rev_median", round(d.jaccard_top_rev.median(), 3)),
               (f"{rel}_jaccard_top_spec_median", round(d.jaccard_top_spec.median(), 3))]
    for fam in ["BRCA", "COADREAD"]:
        d = ok[(ok.relation == "parent_family") & ok.group_a.map(FAMILY).eq(fam)]
        ov += [(f"{fam}_siblings_rho_rev", ";".join(f"{a}~{b} {r}" for a, b, r in zip(d.group_a, d.group_b, d.rho_rev))),
               (f"{fam}_siblings_jaccard_top_rev",
                ";".join(f"{a}~{b} {r}" for a, b, r in zip(d.group_a, d.group_b, d.jaccard_top_rev)))]
    pos = ctl[(ctl.expected == "high") & ctl.rev_pct.notna()]
    neg = ctl[(ctl.expected == "low") & ctl.rev_pct.notna()]
    ov += [("controls_high_n", len(pos)), ("controls_high_rev_pct_mean", round(pos.rev_pct.mean(), 1)),
           ("controls_high_spec_pct_mean", round(pos.spec_pct.mean(), 1)),
           ("controls_low_rev_pct_mean", round(neg.rev_pct.mean(), 1)),
           ("controls_low_spec_pct_mean", round(neg.spec_pct.mean(), 1))]
    ov = pd.DataFrame(ov, columns=["metric", "value"])
    ov.to_csv(OUT / "overview.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(ov.to_string(index=False))


if __name__ == "__main__":
    main()
