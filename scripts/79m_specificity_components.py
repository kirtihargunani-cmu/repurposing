#!/usr/bin/env python
"""Two specificity components the ranking does not currently use, as REPORTED columns (D173).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

Measured gap (D173): target specificity carries **no weight** in the current ranking.
Spearman(target tau, `layered_pct`) = **-0.017**; mean target tau is flat across all ten lead-score
deciles (0.690-0.751). A drug targeting ESR1 (tau 0.961) and one targeting TUBB (0.700) rank alike,
which is why ubiquitously-expressed targets - NR3C1 at 0.634 - let corticosteroids clear the gate in
all 41 cohorts and crowd the top of several.

COMPONENT 1 - target specificity across cohorts (`tau`).
    tau(gene) = sum_i (1 - x_i / max_i x) / (n_cohorts - 1)   over the 41 cohort median TPMs
    0 = expressed evenly everywhere (carries no disease-specific information)
    1 = expressed in a single cohort
Per drug x cohort we report the mean and max tau of the targets that passed the gate THERE, so the
statistic describes the targets actually supporting the drug in that disease.

COMPONENT 2 - cross-disease promiscuity. `split_overview` already counts 33 drugs in the top 10 of
>= 20 cohorts by `rev` (100 by `rev_prolif`) and that count feeds nothing. A drug that tops half the
atlas is making a statement about cytotoxicity, not about any disease.

REPORTED ONLY. No gate, no re-scoring, no change to any existing column. Written in place and
idempotently, as 79i2 / 79j do.

LIMITS, because these columns invite over-reading:
  * **tau is a property of the cohort PANEL, not of the gene.** It is computed over these 41
    cohorts; a different panel gives a different number. It is not an absolute tissue-specificity.
  * **A ubiquitous target is not a bad drug target.** Topoisomerase and tubulin inhibitors work.
    This measures *discoverability*, not efficacy, and demoting low-tau targets deliberately trades
    away real drugs for repurposing-relevant ones ([[discovery-is-the-goal]]). Known-use metrics
    should be expected to get WORSE if these are ever promoted into the score.
  * **Promiscuity is confounded with coverage.** A drug scored in 41 cohorts can top 20; one scored
    in 5 cannot. `promiscuity` normalises by cohorts scored, but the raw counts do not.
  * tau uses the targets EXPRESSED in that cohort (79b `targets_expressed`). A drug whose gate pass
    rests on an off-mechanism target (D168: 1,192 such rows) gets that target's tau, not its
    mechanism target's - read it beside `target_absent`.

Env: SPLIT_SFX (default _v2_tf). Output: columns on results/lincs/drug_reversal_split<SFX>.tsv
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
RES = ROOT / "results/lincs"
PROC = ROOT / "data/processed/lincs"
SFX = os.environ.get("SPLIT_SFX", "_v2_tf")
TOP_PCT = 90.0          # "tops a cohort" = lead percentile at or above this
COLS = ["target_tau", "target_tau_max", "target_tau_pct", "n_cohorts_scored",
        "n_cohorts_lead_top", "promiscuity", "n_cohorts_rev_top10", "n_cohorts_rev_topdecile"]


def gene_tau():
    """tau per gene over the cohort median-TPM panel, keyed by gene symbol."""
    e = pd.read_csv(PROC / "target_expression.tsv.gz", sep="\t")
    x = (e.pivot_table(index="gene_id", columns="group_id", values="median_tpm", aggfunc="first")
         .fillna(0).clip(lower=0))
    mx = x.max(axis=1)
    tau = (1 - x.div(mx.replace(0, np.nan), axis=0)).sum(axis=1) / (x.shape[1] - 1)
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", low_memory=False)
    sc = [c for c in idm.columns if "symbol" in c.lower()][0]
    gc = [c for c in idm.columns if "gene_id" in c.lower()][0]
    out = {}
    for s, g in zip(idm[sc].astype(str), idm[gc].astype(str)):
        if g in tau.index and np.isfinite(tau[g]):
            out[s] = float(tau[g])
    print(f"  tau over {x.shape[0]:,} genes x {x.shape[1]} cohorts; "
          f"median {tau.median():.3f}; mapped to {len(out):,} symbols")
    return out


def main():
    sp = RES / f"drug_reversal_split{SFX}.tsv"
    if not sp.exists():
        print(f"missing {sp}", file=sys.stderr)
        return 1
    d = pd.read_csv(sp, sep="\t", low_memory=False)
    d = d.drop(columns=[c for c in COLS if c in d.columns])          # idempotent
    tau = gene_tau()

    gate = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False)
    gate = gate[["group_id", "pert_id", "targets_expressed"]].drop_duplicates(["group_id", "pert_id"])

    def agg(s, f):
        if not isinstance(s, str) or not s:
            return np.nan
        v = [tau[t] for t in s.split("|") if t in tau]
        return float(f(v)) if v else np.nan

    gate["target_tau"] = [agg(s, np.mean) for s in gate.targets_expressed]
    gate["target_tau_max"] = [agg(s, np.max) for s in gate.targets_expressed]
    d = d.merge(gate.drop(columns="targets_expressed"), on=["group_id", "pert_id"], how="left")
    d["target_tau_pct"] = (d.groupby("group_id").target_tau
                           .transform(lambda v: v.rank(method="max") / v.notna().sum() * 100
                                      if v.notna().any() else v).round(1))

    # component 2: how much of the atlas does this drug top?
    lead = pd.to_numeric(d.layered_pct, errors="coerce")
    revp = pd.to_numeric(d.get("rev_pct"), errors="coerce")
    sc = d.groupby("pert_id").group_id.nunique().rename("n_cohorts_scored")
    tp = d[lead >= TOP_PCT].groupby("pert_id").group_id.nunique().rename("n_cohorts_lead_top")
    # Two different "tops a cohort" definitions, kept apart because 79i conflates them:
    #   n_cohorts_rev_top10      literally the 10 highest-rev drugs in the cohort (strict)
    #   n_cohorts_rev_topdecile  rev_pct >= 90, i.e. ~150 drugs of ~1,500 (what 79i:550 actually does)
    # NOTE 79i:551 names its metric `drugs_top10_rev_in_ge20_groups` but computes the DECILE
    # (`ok[c+"_pct"] >= 100 - TOP_PCT`), which is why it reports 33 where a true top-10 gives 0.
    r10 = (d.assign(_r=revp).sort_values("_r", ascending=False)
           .groupby("group_id").head(10).groupby("pert_id").group_id.nunique()
           .rename("n_cohorts_rev_top10"))
    rdec = d[revp >= 90].groupby("pert_id").group_id.nunique().rename("n_cohorts_rev_topdecile")
    d = d.join(sc, on="pert_id").join(tp, on="pert_id").join(r10, on="pert_id").join(rdec, on="pert_id")
    for c in ("n_cohorts_lead_top", "n_cohorts_rev_top10", "n_cohorts_rev_topdecile"):
        d[c] = d[c].fillna(0).astype(int)
    d["promiscuity"] = (d.n_cohorts_lead_top / d.n_cohorts_scored.replace(0, np.nan)).round(3)
    d.to_csv(sp, sep="\t", index=False)

    tb = lambda v: v.astype(str).isin(["True", "TRUE", "1"])
    g = d[tb(d.passes_gate)]
    print(f"\n  wrote {len(COLS)} columns to {sp.name} ({len(d):,} rows)")
    print(f"  target_tau present on {g.target_tau.notna().sum():,} of {len(g):,} gate-passing rows "
          f"(mean {g.target_tau.mean():.3f})")
    pr = g.drop_duplicates("pert_id")
    print(f"  drugs in >= 20 cohorts' top DECILE by lead: {int((pr.n_cohorts_lead_top >= 20).sum()):,}; "
          f"by rev: {int((pr.n_cohorts_rev_topdecile >= 20).sum()):,} "
          f"(reconciles with split_overview's mis-named drugs_top10_rev_in_ge20_groups)")
    print(f"  drugs in >= 20 cohorts' true TOP 10 by rev: "
          f"{int((pr.n_cohorts_rev_top10 >= 20).sum()):,} (max {int(pr.n_cohorts_rev_top10.max())} cohorts)")
    print("\n  mean target_tau by lead-score decile (flat = the score ignores target specificity):")
    dec = pd.qcut(pd.to_numeric(g.layered_pct, errors="coerce"), 10, labels=False, duplicates="drop")
    print(g.groupby(dec).target_tau.agg(["mean", "size"]).round(3).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
