#!/usr/bin/env python
"""Celligner matrices: disease x cell line weights, and cell line x cell line distances (for borrowing drug data).

Runs in the `celligner` env (sklearn).
1. Disease x line weights: for each group, the mean over its patients of their softmax weights from
   scripts/73_celligner_tissue_keys.py (WEIGHTS=hybrid default: organ key + closer non-key lines; or strict),
   renormalized to sum 1. A line many patients resemble gets more weight.
   Bulk tumors give similarity, not cell-type fractions, so a weight means "share of the disease's patients that
   resemble this line", not "share of tumor cells like this line".
2. Line x line distances: Euclidean distance between cancer lines in the Celligner aligned space (70 PCs of the
   combined aligned matrix, as for tumor-line distances), the basis for extrapolating LINCS drug data from
   profiled lines to unprofiled ones. The package's PCA is unseeded, so the 70-PC space is recomputed here with a
   fixed seed and saved; vs the main run's tumor-line distances: r = 0.99987, same nearest line 93%.
   Assumption to test once LINCS is in: lines close in baseline expression respond alike to drugs.
Outputs: results/celligner/disease_line_weights.tsv (groups x lines, wide), disease_line_weights_long.tsv,
         line_neighbors.tsv (10 nearest lines per line), line_distance_summary.tsv (per lineage)
         data/processed/celligner/{pcs_70.npy, pcs_70_rows.txt, dist_line_line.tsv.gz}
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
DEPMAP = ROOT / "data/raw/depmap_24q4"
NEIGHBORS, SEED = 10, 0
PFX = {"hybrid": "hybrid", "strict": "keyed"}[os.environ.get("WEIGHTS", "hybrid")]


def disease_weights(models):
    w = pd.read_csv(PROC / f"weights_{PFX}.tsv.gz", sep="\t")
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    n_pat = w.groupby("group_id")["sample"].nunique()
    g = w.groupby(["group_id", "ModelID"]).weight.sum().rename("weight").reset_index()
    g["weight"] = g.weight / g.group_id.map(n_pat)          # mean over patients (absent = 0)
    g["weight"] = g.weight / g.groupby("group_id").weight.transform("sum")
    top = pd.read_csv(RES / f"{PFX}_top_lines.tsv", sep="\t")
    g = g.merge(top.groupby(["group_id", "ModelID"])["sample"].nunique().rename("n_patients_top10").reset_index(),
                how="left").fillna({"n_patients_top10": 0})
    g = g.merge(top[top["rank"] == 1].groupby(["group_id", "ModelID"]).size().rename("n_patients_nearest").reset_index(),
                how="left").fillna({"n_patients_nearest": 0})
    g["key"] = g.group_id.map(keys.key)
    def key_set(k):  # same rule as scripts/73_celligner_tissue_keys.py
        m = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        return set(m.index)
    ks = {gid: key_set(k) for gid, k in keys.iterrows()}
    g["in_key"] = [m in ks[gid] for gid, m in zip(g.group_id, g.ModelID)]
    g = g.join(models[["CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]], on="ModelID")
    g = g.sort_values(["group_id", "weight"], ascending=[True, False])
    g["cum_weight"] = g.groupby("group_id").weight.cumsum()
    g.round(6).to_csv(RES / "disease_line_weights_long.tsv", sep="\t", index=False)
    wide = g.pivot(index="group_id", columns="ModelID", values="weight").fillna(0)
    wide = wide[wide.sum().sort_values(ascending=False).index]
    wide.round(6).to_csv(RES / "disease_line_weights.tsv", sep="\t")
    # lines needed to cover 80% of each disease's weight
    n80 = g.groupby("group_id").apply(lambda d: int((d.cum_weight < 0.8).sum() + 1)).rename("n_lines_80pct")
    print(f"disease x line matrix: {wide.shape[0]} groups x {wide.shape[1]} lines with weight >= 1e-4 for some patient")
    print(n80.describe().round(1).to_string())
    return n80


def main():
    models = pd.read_csv(DEPMAP / "Model.csv", index_col=0)
    rows = pd.read_csv(PROC / "aligned_rows.txt", header=None)[0].values
    cancer = models.index[(models.OncotreePrimaryDisease != "Non-Cancerous") & models.index.isin(rows)]
    models = models.loc[cancer]

    n80 = disease_weights(models)

    x = np.load(PROC / "aligned.npy").astype(np.float64)
    pcs = PCA(70, random_state=SEED).fit_transform(x)
    np.save(PROC / "pcs_70.npy", pcs.astype(np.float32))
    pd.Series(rows).to_csv(PROC / "pcs_70_rows.txt", index=False, header=False)
    idx = pd.Series(np.arange(len(rows)), index=rows)

    d = pd.DataFrame(pairwise_distances(pcs[idx[cancer]]), index=cancer, columns=cancer)
    d.round(3).to_csv(PROC / "dist_line_line.tsv.gz", sep="\t", index_label="ModelID")

    v = d.values.copy()
    np.fill_diagonal(v, np.inf)
    o = np.argsort(v, axis=1)[:, :NEIGHBORS]
    lin = models.OncotreeLineage
    nb = []
    for k in range(NEIGHBORS):
        mid = cancer[o[:, k]]
        nb.append(pd.DataFrame({"ModelID": cancer, "rank": k + 1, "neighbor": mid,
                                "distance": np.take_along_axis(v, o[:, [k]], 1)[:, 0].round(3),
                                "same_lineage": lin[cancer].values == lin[mid].values}))
    nb = pd.concat(nb).sort_values(["ModelID", "rank"])
    nb.insert(1, "CellLineName", nb.ModelID.map(models.CellLineName))
    nb.insert(2, "OncotreeLineage", nb.ModelID.map(lin))
    nb.insert(5, "neighbor_name", nb.neighbor.map(models.CellLineName))
    nb.insert(6, "neighbor_lineage", nb.neighbor.map(lin))
    nb.to_csv(RES / "line_neighbors.tsv", sep="\t", index=False)

    # how far lines sit from each other, within vs across lineage, against the tumor-line scale
    same = lin.values[:, None] == lin.values[None, :]
    iu = np.triu_indices(len(cancer), 1)
    nn = v.min(axis=1)
    rows_s = []
    for L, ix in pd.Series(np.arange(len(cancer)), index=lin.values).groupby(level=0):
        ix = ix.values
        sub = d.values[np.ix_(ix, ix)][np.triu_indices(len(ix), 1)]
        rows_s.append(dict(lineage=L, n_lines=len(ix), median_nn_distance=round(np.median(nn[ix]), 2),
                           median_within_lineage=round(np.median(sub), 2) if len(sub) else np.nan,
                           median_to_other_lineages=round(np.median(d.values[np.ix_(ix, np.setdiff1d(np.arange(len(cancer)), ix))]), 2),
                           pct_nn_same_lineage=round(100 * nb[(nb["rank"] == 1) & nb.ModelID.isin(cancer[ix])].same_lineage.mean(), 1)))
    s = pd.DataFrame(rows_s).sort_values("n_lines", ascending=False)
    s.to_csv(RES / "line_distance_summary.tsv", sep="\t", index=False)
    tl = pd.read_csv(RES / "patient_qc.tsv", sep="\t", usecols=["d1"]).d1
    print(f"line-line: {len(cancer)} lines; median nearest-neighbor distance {np.median(nn):.1f}; "
          f"within-lineage median {np.median(d.values[iu][same[iu]]):.1f}, across {np.median(d.values[iu][~same[iu]]):.1f}; "
          f"nearest neighbor same lineage {100 * nb[nb['rank'] == 1].same_lineage.mean():.1f}%; "
          f"tumor -> nearest line median {tl.median():.1f}")
    pd.set_option("display.width", 200)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
