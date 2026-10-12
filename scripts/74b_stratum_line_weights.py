#!/usr/bin/env python
"""Celligner line weights for the Step 10 strata (runs in `celligner`).

User, 2026-10-12: a stratum's line weights are the mean of its own patients' Celligner weights (not the parent
group's); the tissue key is the parent's. Same rule as the disease weights of scripts/74 (hybrid per-patient softmax
weights of scripts/73, mean over patients, absent = 0, renormalized to 1).
Strata (disease side, results/step8/step10_strata.tsv):
  subgroups  COAD-CMS1/3-like, COAD-CMS2/4-like, READ-C1..C3: patients from results/step8/step10_patient_strata.tsv
  breast     each kept axis cut into thirds of the axis score within the subtype (results/step8/breast_axes_scores.tsv);
             <axis>-high = top third (score >= 2/3 quantile), <axis>-low = bottom third (<= 1/3 quantile), as in
             scripts/92_step10_signatures.R
Outputs: results/celligner/stratum_line_weights_long.tsv (stratum, parent, ModelID, weight, in_key, line info,
         cum_weight), stratum_line_summary.tsv (per stratum: patients, effective lines, top lines, overlap with the
         parent's weights and with the sibling strata's)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
S8 = ROOT / "results/step8"


def strata_members():
    strata = pd.read_csv(S8 / "step10_strata.tsv", sep="\t")
    pat = pd.read_csv(S8 / "step10_patient_strata.tsv", sep="\t")
    ax = pd.read_csv(S8 / "breast_axes_scores.tsv", sep="\t")
    rows = [pat.rename(columns={"group_id": "parent"})[["stratum", "parent", "sample"]].assign(sibling_set=pat.group_id)]
    for _, r in strata[strata.kind == "continuous axis"].iterrows():
        a = r.stratum.rsplit("-", 1)[1]
        sc = ax[ax.group_id == r.group_id][["sample", a]]
        lo, hi = sc[a].quantile([1 / 3, 2 / 3])
        for end, sel in [("high", sc[a] >= hi), ("low", sc[a] <= lo)]:
            rows.append(pd.DataFrame({"stratum": f"{r.stratum}-{end}", "parent": r.group_id,
                                      "sample": sc["sample"][sel].values, "sibling_set": r.stratum}))
    return pd.concat(rows, ignore_index=True)


def main():
    mem = strata_members()
    w = pd.read_csv(PROC / "weights_hybrid.tsv.gz", sep="\t")
    parent_w = pd.read_csv(RES / "disease_line_weights_long.tsv", sep="\t")
    line_info = parent_w.drop_duplicates("ModelID").set_index("ModelID")[["CellLineName", "OncotreeLineage",
                                                                          "OncotreePrimaryDisease"]]
    out = []
    for s, m in mem.groupby("stratum"):
        parent = m.parent.iloc[0]
        ws = w[w["sample"].isin(m["sample"]) & (w.group_id == parent)]
        n = ws["sample"].nunique()
        g = ws.groupby("ModelID").weight.sum() / n
        g = (g / g.sum()).rename("weight").reset_index()
        g.insert(0, "stratum", s)
        g.insert(1, "parent", parent)
        g.insert(2, "sibling_set", m.sibling_set.iloc[0])
        g["n_patients"] = n
        g["n_members"] = len(m)
        out.append(g)
    out = pd.concat(out, ignore_index=True)
    pk = parent_w.drop_duplicates(["group_id", "ModelID"]).set_index(["group_id", "ModelID"]).in_key
    out["in_key"] = [bool(pk.get((p, mid), False)) for p, mid in zip(out.parent, out.ModelID)]
    out = out.join(line_info, on="ModelID").sort_values(["stratum", "weight"], ascending=[True, False])
    out["cum_weight"] = out.groupby("stratum").weight.cumsum()
    out.round(6).to_csv(RES / "stratum_line_weights_long.tsv", sep="\t", index=False)

    # summary: concentration, top lines, similarity to parent and to siblings (1 - total variation distance)
    W = out.pivot_table(index="ModelID", columns="stratum", values="weight", fill_value=0)
    P = parent_w.pivot_table(index="ModelID", columns="group_id", values="weight", fill_value=0)
    sim = lambda a, b: 1 - 0.5 * np.abs(a.reindex(a.index.union(b.index), fill_value=0)
                                        - b.reindex(a.index.union(b.index), fill_value=0)).sum()
    summ = []
    for s, d in out.groupby("stratum"):
        sibs = [t for t in out.stratum[out.sibling_set == d.sibling_set.iloc[0]].unique() if t != s]
        summ.append(dict(stratum=s, parent=d.parent.iloc[0], n_members=int(d.n_members.iloc[0]),
                         n_patients_weighted=int(d.n_patients.iloc[0]),
                         n_eff_lines=round(1 / (d.weight ** 2).sum(), 1), weight_in_key=round(d.weight[d.in_key].sum(), 3),
                         overlap_parent=round(sim(W[s], P[d.parent.iloc[0]]), 3),
                         overlap_siblings=round(np.mean([sim(W[s], W[t]) for t in sibs]), 3) if sibs else np.nan,
                         top_lines=", ".join(f"{a} {b:.3f}" for a, b in zip(d.CellLineName.head(5), d.weight.head(5)))))
    summ = pd.DataFrame(summ)
    summ.to_csv(RES / "stratum_line_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
