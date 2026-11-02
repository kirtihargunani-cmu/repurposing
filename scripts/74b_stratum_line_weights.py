#!/usr/bin/env python
"""Celligner line weights for the Step 10 strata (runs in `celligner`).

User, 2026-10-12: a stratum's line weights are the mean of its own patients' Celligner weights (not the parent
group's); the tissue key is the parent's. Same rule as the disease weights of scripts/74 (hybrid per-patient softmax
weights of scripts/73, mean over patients, absent = 0, renormalized to 1).
Strata and their patients: results/step10/stratum_patients.tsv (disease side, scripts/95; exactly the tumors of each
stratum signature): COAD-CMS1/3-like, COAD-CMS2/4-like, READ-C1..C3, and the high / low thirds of the 13 breast axes.
Siblings: subgroups of the same parent; the two ends of the same breast axis.
Outputs: results/celligner/stratum_line_weights_long.tsv (stratum, parent, ModelID, weight, in_key, line info,
         cum_weight), stratum_line_summary.tsv (per stratum: patients, effective lines, top lines, overlap with the
         parent's weights and with the sibling strata's)
VERSION=v2 (D72, 2026-10-28): the v2 subgroups (results/step10_v2/stratum_patients.tsv, copied from main: 10 subgroups
of BLCA, HNSC, COAD, ACC, LGG-IDHmut-noncodel; siblings = the other subgroup of the same parent); outputs suffixed _v2.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
V2 = os.environ.get("VERSION", "v1") == "v2"
SFX = "_v2" if V2 else ""


def strata_members():
    sp = pd.read_csv(ROOT / f"results/step10{SFX}/stratum_patients.tsv", sep="\t")
    if "weak" not in sp:
        sp["weak"] = False
    sp = sp.rename(columns={"group_id": "parent"})
    # siblings: subgroups of the same parent; the two ends of the same breast axis
    sp["sibling_set"] = np.where(sp.stratum.str.endswith(("-high", "-low")), sp.stratum.str.rsplit("-", n=1).str[0],
                                 sp.parent)
    return sp[["stratum", "parent", "sample", "sibling_set", "weak"]]


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
    out.round(6).to_csv(RES / f"stratum_line_weights_long{SFX}.tsv", sep="\t", index=False)

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
    summ.to_csv(RES / f"stratum_line_summary{SFX}.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
