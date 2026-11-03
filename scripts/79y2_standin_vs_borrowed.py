#!/usr/bin/env python
"""Stand-in vs the drug's own borrowed score (user 2026-11-02; follows 79y, C27).

Question: for an FDA drug never profiled in a cohort's tissue, is a stand-in (79y) a better estimate of the drug's
matched-line score than the drug's own borrowed score (its average over other tissues' lines, as 79e scores borrowed
drugs)?
Counterfactual run (scripts/79y2_borrowed_counterfactual.sh): a copy of disease_drugs.tsv in which every FDA drug that
was tested in its cohort's lines (direct) and has >= 3 lines elsewhere (79b's borrowing rule) is switched to borrowed;
79e and 79i are rerun on it in a sandbox (production results untouched). Its borrowed spec_mech is placed in the
cohort's distribution in that run (percentile among the other drugs' spec_mech; specificity is scaled over all drugs of
a cohort, so the counterfactual's scale differs slightly from production's).
Set: the 79y validation pairs (FDA drugs tested in their cohort's lines with an accepted stand-in, similarity measured
outside those lines) that also have a borrowed score. Truth = the drug's own production spec_mech percentile.
Decision rule (fixed before looking): Spearman with the truth for stand-in and for borrowed, paired bootstrap (1,000)
of the difference; "stand-in better" if the 95% CI of (stand-in - borrowed) is above 0, "borrowed better" if below,
otherwise "no clear difference". The average of the two is reported as a third estimate under the same test vs borrowed.
Usage: SBX=<sandbox root> python scripts/79y2_standin_vs_borrowed.py
Outputs: results/lincs/standin_vs_borrowed_v2.tsv (per drug x cohort), standin_vs_borrowed_summary_v2.tsv
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SBX = Path(os.environ["SBX"])
N_BOOT, SEED = 1000, 0


def main():
    orig = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    cf = pd.read_csv(SBX / "results/lincs/drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    V = pd.read_csv(RES / "standin_validation_pairs_v2.tsv", sep="\t")
    flipped = cf[cf.source == "borrowed"].merge(orig[["group_id", "pert_id", "source"]], on=["group_id", "pert_id"],
                                                suffixes=("", "_orig"))
    nl = pd.read_csv(SBX / "results/lincs/drug_reversal_v2_consensus.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_lines_used"])
    flipped = flipped[flipped.source_orig == "direct"].merge(nl, on=["group_id", "pert_id"], how="left")
    flipped = flipped[["group_id", "pert_id", "spec_mech", "n_lines_used", "transfer_conf"]]
    # specificity is scaled over all drugs of a cohort, so switching drugs shifts every drug's spec slightly (max 0.24
    # spec_mech units, rev unchanged): the borrowed score is placed in the counterfactual run's own distribution
    dist = {g: x.set_index("pert_id").spec_mech for g, x in cf.groupby("group_id")}

    def pct(g, p, v):
        o = dist[g].drop(p, errors="ignore").dropna().values
        return 100.0 * (o <= v).mean() if len(o) and v == v else np.nan
    flipped["borrowed_pct"] = [pct(g, p, v) for g, p, v in zip(flipped.group_id, flipped.pert_id, flipped.spec_mech)]
    T = V.merge(flipped.rename(columns={"pert_id": "fda_pert", "spec_mech": "borrowed_spec_mech",
                                        "n_lines_used": "borrowed_n_lines"}), on=["group_id", "fda_pert"], how="inner")
    T = T.dropna(subset=["own", "standin", "borrowed_pct"])
    T["average"] = (T.standin + T.borrowed_pct) / 2
    T.to_csv(RES / "standin_vs_borrowed_v2.tsv", sep="\t", index=False)

    rng = np.random.default_rng(SEED)
    own = T.own.values
    est = {k: T[k].values for k in ("standin", "borrowed_pct", "average")}
    rho = {k: spearmanr(own, v)[0] for k, v in est.items()}
    mae = {k: float(np.abs(own - v).mean()) for k, v in est.items()}
    boots = {k: [] for k in est}
    for _ in range(N_BOOT):
        i = rng.integers(0, len(T), len(T))
        for k, v in est.items():
            boots[k].append(spearmanr(own[i], v[i])[0])
    boots = {k: np.array(v) for k, v in boots.items()}
    rows = []
    for k in est:
        lo, hi = np.nanpercentile(boots[k], [2.5, 97.5])
        rows.append(dict(estimate=k, n=len(T), spearman=round(rho[k], 3), ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                         mae_pct=round(mae[k], 1), diff_vs_borrowed=np.nan, diff_ci_lo=np.nan, diff_ci_hi=np.nan, verdict=""))
    for r in rows:
        if r["estimate"] == "borrowed_pct":
            continue
        d = boots[r["estimate"]] - boots["borrowed_pct"]
        lo, hi = np.nanpercentile(d, [2.5, 97.5])
        r.update(diff_vs_borrowed=round(rho[r["estimate"]] - rho["borrowed_pct"], 3), diff_ci_lo=round(lo, 3), diff_ci_hi=round(hi, 3),
                 verdict=f"{r['estimate']} better" if lo > 0 else ("borrowed better" if hi < 0 else "no clear difference"))
    # context: all switched drugs (with or without a stand-in)
    allf = flipped.merge(orig[["group_id", "pert_id", "spec_mech_pct"]], on=["group_id", "pert_id"]).dropna(subset=["borrowed_pct"])
    rows.append(dict(estimate="borrowed_pct, all switched drugs", n=len(allf),
                     spearman=round(spearmanr(allf.spec_mech_pct, allf.borrowed_pct)[0], 3), ci_lo=np.nan, ci_hi=np.nan,
                     mae_pct=round(float((allf.spec_mech_pct - allf.borrowed_pct).abs().mean()), 1),
                     diff_vs_borrowed=np.nan, diff_ci_lo=np.nan, diff_ci_hi=np.nan, verdict=""))
    # by similarity / shared lines
    for lab, m in [("standin, similarity >= 0.6", T.similarity >= 0.6), ("standin, similarity < 0.6", T.similarity < 0.6),
                   ("standin, shared lines >= 3", T.n_shared_lines >= 3)]:
        x = T[m]
        rows.append(dict(estimate=lab, n=len(x), spearman=round(spearmanr(x.own, x.standin)[0], 3),
                         ci_lo=np.nan, ci_hi=np.nan, mae_pct=round(float((x.own - x.standin).abs().mean()), 1),
                         diff_vs_borrowed=round(spearmanr(x.own, x.standin)[0] - spearmanr(x.own, x.borrowed_pct)[0], 3),
                         diff_ci_lo=np.nan, diff_ci_hi=np.nan, verdict=f"borrowed on same set {spearmanr(x.own, x.borrowed_pct)[0]:.3f}"))
    S = pd.DataFrame(rows)
    S.to_csv(RES / "standin_vs_borrowed_summary_v2.tsv", sep="\t", index=False)
    print(S.to_string(index=False))


if __name__ == "__main__":
    main()
