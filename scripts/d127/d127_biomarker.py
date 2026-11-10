#!/usr/bin/env python
"""D127: biomarker-expectation test — the drug side's primary validation (docs/DECISIONS.md D127; D125.1).

The same FDA compound is compared against itself in two cohorts of a biomarker pair, so patents, sponsor, market size
and drug class cancel. One instance = (compound, high cohort, low cohort) where the compound is scored in both.

Score: `layered_pct`, the within-cohort percentile over FDA compounds that pass the expression gate and have a score
(79t3). That population is defined identically in both cohorts and is NOT affected by target_absent / local_route /
low_exposure, so the sides are like-for-like. Rows carrying those flags stay in: in the low cohort they are exactly the
expected-low cases, and dropping them deletes the negatives (only 1 of 7 anti-estrogens is eligible in Basal).

Primary: over the primary set, the share with layered_pct(high) > layered_pct(low), one-sided sign test.
Pass = share >= SHARE_MIN, p < P_MAX and median gap >= GAP_MIN. Thresholds are pre-registered; do not tune them.

Reported, not decision-bearing: per-pair share and median gap, and the per-pair median beside the pooled share; the
same test on spec_mech_pct_pregate (reversal) and target_role_pct (layer 1 alone, 41 groups only); the end-to-end
reading in which an exclusion flag in the low cohort counts as a correct demotion; the previously-seen pairs.

Outputs: results/lincs/d127/{instances,pairs,summary}.tsv, run.json.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import fda_rows as F  # noqa: E402

# --- config block: every threshold of the pre-registered rule lives here only (D127) ---------------------------------
SHARE_MIN = 0.70        # share of instances with layered_pct(high) > layered_pct(low)
P_MAX = 0.05            # one-sided sign test
GAP_MIN = 10.0          # median gap, percentile points
SCORE = "layered_pct"                       # the lead score; the test's primary quantity
ALT_SCORES = ["spec_mech_pct_pregate", "target_role_pct"]   # reported only
EXCLUSION_FLAGS = ["target_absent", "local_route", "low_exposure"]   # kept in the test; used by the end-to-end reading

CLASSES = {"BRAF/MEK": {"BRAF", "MAP2K1", "MAP2K2"}, "EGFR": {"EGFR"}, "HER2": {"ERBB2"},
           "PARP": {"PARP1", "PARP2"}, "anti-estrogen": {"ESR1", "ESR2", "CYP19A1"}, "ALK": {"ALK"},
           "CDK4/6": {"CDK4", "CDK6"}}
# (class, expected-high cohort, expected-low cohort, set). "primary" = the 49-instance decision set;
# "seen" = values already known or prior decisions taken on them (D122 controls, D92-D97, D102c, D110): reported only.
PAIRS = [
    ("BRAF/MEK", "THCA-BRAF-like-classical", "THCA-RAS-like-follicular", "primary"),
    ("BRAF/MEK", "TCGA-SKCM", "TCGA-LUSC", "primary"),
    ("EGFR", "TCGA-LUAD", "TCGA-LUSC", "primary"),
    ("EGFR", "TCGA-HNSC", "TCGA-LUSC", "primary"),
    ("HER2", "BRCA-Her2", "BRCA-Basal", "primary"),
    ("HER2", "BRCA-Her2", "BRCA-LumA", "primary"),
    ("PARP", "TCGA-OV", "TCGA-LUSC", "primary"),
    ("anti-estrogen", "BRCA-LumA", "BRCA-Basal", "seen"),
    ("anti-estrogen", "BRCA-LumB", "BRCA-Basal", "seen"),
    ("ALK", "TCGA-LUAD", "TCGA-LUSC", "seen"),
    ("CDK4/6", "BRCA-LumA", "BRCA-Basal", "seen"),
]
NOT_COMPUTABLE = [("AR antagonist", "TCGA-PRAD", "TCGA-BLCA", "no AR drug is scored in BLCA")]
OUT = ROOT / "results/lincs/d127"


def target_genes(x):
    """D105b evidence targets, written 'GENE:value/value' and separated by | or ;."""
    if not isinstance(x, str):
        return set()
    return {p.split(":")[0].strip() for p in re.split(r"[|;]", x) if p.strip()}


def load():
    g, s = F.fda_rows()
    a = pd.concat([g.assign(kind="group"), s.assign(kind="subgroup")], ignore_index=True)
    a["genes"] = [target_genes(x) for x in a.evidence_targets]
    for c in EXCLUSION_FLAGS:
        a[c] = a[c].map(F.truthy).fillna(False).astype(bool) if c in a else False
    a["excluded"] = a[EXCLUSION_FLAGS].any(axis=1)
    # target_role (D120) is per group only (no subgroups) and capitalises the drug name: join on a lowercase key
    tr = pd.read_csv(ROOT / "results/lincs/target_role/drug_target_role.tsv", sep="\t",
                     usecols=["cohort", "drug", "target_role_pct"])
    tr["_k"] = tr.drug.astype(str).str.lower()
    a["_k"] = a.cmap_name.astype(str).str.lower()
    a = a.merge(tr[["cohort", "_k", "target_role_pct"]].rename(columns={"cohort": "group_id"}),
                on=["group_id", "_k"], how="left").drop(columns="_k")
    return a


def instances(a):
    rows = []
    for cls, hi, lo, kind in PAIRS:
        gs = CLASSES[cls]
        m = a.genes.map(lambda t: bool(t & gs)) & a[SCORE].notna()
        H = a[(a.group_id == hi) & m].set_index("cmap_name")
        L = a[(a.group_id == lo) & m].set_index("cmap_name")
        for d in sorted(set(H.index) & set(L.index)):
            h, l = H.loc[d], L.loc[d]
            r = dict(pair_set=kind, drug_class=cls, high=hi, low=lo, drug=d,
                     hi_pct=float(h[SCORE]), lo_pct=float(l[SCORE]),
                     hi_excluded=bool(h.excluded), lo_excluded=bool(l.excluded),
                     hi_flags=";".join(c for c in EXCLUSION_FLAGS if h[c]),
                     lo_flags=";".join(c for c in EXCLUSION_FLAGS if l[c]))
            r["gap"] = round(r["hi_pct"] - r["lo_pct"], 1)
            r["correct"] = r["hi_pct"] > r["lo_pct"]
            # end-to-end (reported): an exclusion flag in the low cohort, absent in the high one, is a correct demotion
            r["correct_e2e"] = bool((l.excluded and not h.excluded) or
                                    (not (l.excluded and not h.excluded) and r["hi_pct"] > r["lo_pct"]))
            for sc in ALT_SCORES:
                hv, lv = h.get(sc, np.nan), l.get(sc, np.nan)
                r[f"{sc}_gap"] = round(float(hv) - float(lv), 1) if pd.notna(hv) and pd.notna(lv) else np.nan
            rows.append(r)
    return pd.DataFrame(rows)


def test(d, col="correct"):
    n = len(d)
    if not n:
        return dict(n=0, n_correct=0, share=np.nan, p=np.nan, median_gap=np.nan)
    k = int(d[col].sum())
    return dict(n=n, n_correct=k, share=round(k / n, 4),
                p=round(binomtest(k, n, 0.5, alternative="greater").pvalue, 5),
                median_gap=round(float(d.gap.median()), 1))


def alt_test(d, sc):
    x = d[d[f"{sc}_gap"].notna()]
    if not len(x):
        return dict(n=0, share=np.nan, p=np.nan, median_gap=np.nan)
    k = int((x[f"{sc}_gap"] > 0).sum())
    return dict(n=len(x), share=round(k / len(x), 4),
                p=round(binomtest(k, len(x), 0.5, alternative="greater").pvalue, 5),
                median_gap=round(float(x[f"{sc}_gap"].median()), 1))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = load()
    inst = instances(a)
    inst.to_csv(OUT / "instances.tsv", sep="\t", index=False)

    per = []
    for (ps, cls, hi, lo), d in inst.groupby(["pair_set", "drug_class", "high", "low"], sort=False):
        per.append(dict(pair_set=ps, drug_class=cls, high=hi, low=lo, **test(d),
                        n_lo_excluded=int(d.lo_excluded.sum()), n_hi_excluded=int(d.hi_excluded.sum()),
                        share_e2e=round(float(d.correct_e2e.mean()), 4)))
    per = pd.DataFrame(per)
    per.to_csv(OUT / "pairs.tsv", sep="\t", index=False)

    P = inst[inst.pair_set == "primary"]
    pri = test(P)
    ties = int((P.gap == 0).sum())
    nt = P[P.gap != 0]
    pass_share = pri["share"] >= SHARE_MIN
    pass_p = pri["p"] < P_MAX
    pass_gap = pri["median_gap"] >= GAP_MIN
    verdict = bool(pass_share and pass_p and pass_gap)

    rows = [dict(item="PRIMARY: instances", value=pri["n"]),
            dict(item="PRIMARY: correct (high > low)", value=pri["n_correct"]),
            dict(item=f"PRIMARY: share (pass >= {SHARE_MIN})", value=pri["share"]),
            dict(item=f"PRIMARY: sign-test p (pass < {P_MAX})", value=pri["p"]),
            dict(item=f"PRIMARY: median gap (pass >= {GAP_MIN})", value=pri["median_gap"]),
            dict(item="PRIMARY VERDICT", value="PASS" if verdict else "NOT MET"),
            dict(item="  criteria met", value=f"share {pass_share}, p {pass_p}, gap {pass_gap}"),
            dict(item="ties (gap == 0), counted as not-correct", value=ties),
            dict(item="sign test excluding ties", value=f"{test(nt)['share']} (n={len(nt)}, p={test(nt)['p']})"),
            dict(item="per-pair median of the share", value=round(float(per[per.pair_set == "primary"].share.median()), 4)),
            dict(item="end-to-end share (reported)", value=round(float(P.correct_e2e.mean()), 4))]
    for sc in ALT_SCORES:
        t = alt_test(P, sc)
        rows.append(dict(item=f"reported: {sc}", value=f"share {t['share']} (n={t['n']}, p={t['p']}, median gap {t['median_gap']})"))
    S = inst[inst.pair_set == "seen"]
    if len(S):
        ts = test(S)
        rows.append(dict(item="previously seen pairs (reported, not primary)",
                         value=f"share {ts['share']} (n={ts['n']}, median gap {ts['median_gap']})"))
    for cls, hi, lo, why in NOT_COMPUTABLE:
        rows.append(dict(item=f"not computable: {cls} {hi} vs {lo}", value=why))
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    run = dict(score=SCORE, share_min=SHARE_MIN, p_max=P_MAX, gap_min=GAP_MIN, primary=pri, verdict=verdict,
               ties=ties, n_pairs_primary=int((per.pair_set == "primary").sum()),
               per_pair_median_share=float(per[per.pair_set == "primary"].share.median()))
    (OUT / "run.json").write_text(json.dumps(run, indent=1, default=float))

    pd.set_option("display.width", 250)
    print(summ.to_string(index=False))
    print("\nPer pair:")
    print(per.to_string(index=False))


if __name__ == "__main__":
    main()
