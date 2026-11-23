#!/usr/bin/env python
"""Two-stage scorer: gate on the lead score, then rank survivors on gene-level activity (D204).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-22: "a scorer that performed the initial filtering with the lead score which we saw
is effective and then calculated a score among the drugs kept with this gene-level master
regulator activity profiles ... timepoint separated, drop the plate and non-biological
confounders".

EXPLORATORY. Not a registered test - it measures a design before deciding whether to register it.
No acceptance criterion is claimed here, and nothing is adopted on the strength of it.

THE PREMISE, MEASURED FIRST. On the rows where both exist, the lead score beats `rev`:
median AUC **0.603 vs 0.545**, +0.035, better in 16/25 cells. `layered_pct` is populated for
**15.1%** of rows (349 of 3,464 drugs) because it needs an identified target unit - so the gate is
a precision filter that discards 85% of rows by construction, which is the point of a gate but
must be stated rather than discovered later (the D184 first-run failure).

STAGE 1 - gate.  Keep rows that pass the production gate, are not `target_absent`, and clear a
`layered_pct` percentile. Window-agnostic, because the lead score has no window-resolved form.

STAGE 2 - rank.  Among survivors, rank by the gene-level score, WINDOW-SEPARATED.

CONFOUNDER REMOVAL. Regulators belonging to a program flagged `artefact` (plate / RNA quality) or
`proliferation` by D199 are dropped from the sum. Only FLAGGED regulators are dropped - a
regulator in no program is not evidence of a confound, so it is kept. Dropping everything outside
a clean program would cut the regulator space to ~14% and repeat D201's coverage collapse.

Judged on the anchors that are independent of approval (D120 dependency / alteration /
over-activity) as well as on approval recovery, because D192 showed the latter is bimodal and
rewards cytotoxics.

Output: results/lincs/interpretation/d204_two_stage.tsv
"""
import json
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "interpretation"
ACT = RES / "regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz"
GATES = [0, 25, 50, 75]          # layered_pct percentile floors to try


def auc(y, s):
    m = ~np.isnan(s)
    y, s = np.asarray(y)[m], np.asarray(s)[m]
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 < 3 or n0 < 3:
        return None
    r = pd.Series(s).rank().to_numpy()
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # regulators to drop: those in an artefact or proliferation program (D199)
    P = pd.read_csv(RES / "programs/programs.tsv", sep="\t")
    F = pd.read_csv(RES / "programs/program_filter.tsv", sep="\t")
    bad = F[(F.artefact) | (F.fail_proliferation)][["cohort", "program"]]
    Pb = P.merge(bad, on=["cohort", "program"])
    drop = {}
    for r in Pb.itertuples():
        drop.setdefault(r.cohort, set()).update(str(r.members).split(","))
    print(f"confounded regulators dropped: "
          f"{sum(len(v) for v in drop.values()):,} over {len(drop)} cohorts "
          f"(median {int(np.median([len(v) for v in drop.values()]))} per cohort)")

    acc = []
    for ch in pd.read_csv(ACT, sep="\t", chunksize=3_000_000, low_memory=False,
                          usecols=["cohort", "pert_id", "cmap_name", "window", "regulator",
                                   "rev_contrib", "rev_gene_contrib"]):
        ch["bad"] = [reg in drop.get(c, ()) for c, reg in zip(ch.cohort, ch.regulator)]
        k = ["cohort", "pert_id", "cmap_name", "window"]
        a = ch.groupby(k, sort=False).agg(rev=("rev_contrib", "sum"),
                                          rev_gene=("rev_gene_contrib", "sum")).reset_index()
        b = (ch[~ch.bad].groupby(k, sort=False)
             .agg(rev_clean=("rev_contrib", "sum"),
                  rev_gene_clean=("rev_gene_contrib", "sum"),
                  n_reg_clean=("regulator", "size")).reset_index())
        acc.append(a.merge(b, on=k, how="left"))
    S = (pd.concat(acc, ignore_index=True)
         .groupby(["cohort", "pert_id", "cmap_name", "window"], as_index=False)
         [["rev", "rev_gene", "rev_clean", "rev_gene_clean", "n_reg_clean"]].sum())
    del acc

    sp = pd.read_csv(RES / "drug_reversal_split_v2_tf.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "layered_pct", "passes_gate", "target_absent"])
    lab = pd.read_csv(RES / "labels/pair_labels_v2_tf.tsv.gz", sep="\t", low_memory=False)
    lab = lab[lab.unit_kind == "group"][["unit_id", "pert_id", "label_phase"]]
    TR = pd.read_csv(RES / "target_role/drug_target_role.tsv", sep="\t", low_memory=False)
    TR["k"] = TR.drug.str.lower().str.strip()

    M = (S.merge(sp, left_on=["cohort", "pert_id"], right_on=["group_id", "pert_id"])
         .merge(lab, left_on=["cohort", "pert_id"], right_on=["unit_id", "pert_id"]))
    M["y"] = (M.label_phase >= 3).astype(int)
    M["k"] = M.cmap_name.str.lower().str.strip()
    M = M.merge(TR[["cohort", "k", "T", "D", "A"]], on=["cohort", "k"], how="left")
    truthy = lambda s: s.astype(str).str.lower().isin(("true", "1", "yes"))        # noqa: E731
    M["ok_gate"] = truthy(M.passes_gate) & ~truthy(M.target_absent)
    M.to_csv(OUT / "d204_scores.tsv.gz", sep="\t", index=False, float_format="%.6g")

    rows = []
    for g0 in GATES:
        sub = M[M.ok_gate & M.layered_pct.notna()]
        if g0:
            sub = sub[sub.groupby(["cohort", "window"]).layered_pct
                      .transform(lambda s: s >= np.percentile(s, g0))]
        for (c, w), g in sub.groupby(["cohort", "window"]):
            if len(g) < 25 or g.y.sum() < 3:
                continue
            r = dict(gate=g0, cohort=c, window=w, n=len(g), n_pos=int(g.y.sum()))
            for nm, col in (("lead", "layered_pct"), ("rev", "rev"), ("gene", "rev_gene"),
                            ("gene_clean", "rev_gene_clean"), ("rev_clean", "rev_clean")):
                r[f"auc_{nm}"] = auc(g.y.values, g[col].values.astype(float))
                gg = g.dropna(subset=["T"])
                r[f"T_{nm}"] = (gg[col].corr(gg["T"], method="spearman")
                                if len(gg) >= 20 and gg[col].nunique() > 4 else np.nan)
            rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "d204_two_stage.tsv", sep="\t", index=False)

    print(f"\n{'gate':>5s} {'cells':>5s} | " + " ".join(f"{n:>11s}" for n in
          ("AUC lead", "AUC gene_cl", "AUC rev", "T lead", "T gene_cl", "T rev")))
    for g0, g in R.groupby("gate"):
        print(f"{g0:>5d} {len(g):>5d} | "
              f"{g.auc_lead.median():11.3f} {g.auc_gene_clean.median():11.3f} {g.auc_rev.median():11.3f} "
              f"{g.T_lead.median():11.3f} {g.T_gene_clean.median():11.3f} {g.T_rev.median():11.3f}")
    for g0, g in R.groupby("gate"):
        d = (g.auc_gene_clean - g.auc_lead).dropna()
        dt = (g.T_gene_clean - g.T_lead).dropna()
        print(f"  gate {g0:>2d}: gene_clean vs lead — AUC {d.median():+.3f} ({int((d>0).sum())}/{len(d)}), "
              f"target-role {dt.median():+.3f} ({int((dt>0).sum())}/{len(dt)})")
    (OUT / "d204_run.json").write_text(json.dumps(dict(
        script="79ao_two_stage.py", exploratory=True, gates=GATES,
        n_cells=int(len(R)), note="not a registered test; measures a design before registering"),
        indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
