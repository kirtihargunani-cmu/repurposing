#!/usr/bin/env python
"""Celligner matches restricted to each group's tissue key (config/celligner_tissue_keys.tsv).

Goal (user, 2026-10-08): matching is for drug applicability to the disease, so each patient is matched only to
DepMap lines of the group's key: by default the organ (DepMap OncotreeLineage = OncoTree tissue of the group's
code); a substitute where the organ has no usable lines; Myeloid or Lymphoid lines for hematologic groups
(user 2026-10-12; was one "Blood" key, myeloid + lymphoid).
MODE=strict: only the key's lines. MODE=hybrid (user 2026-10-08): the key's lines plus, per patient, any other
line closer than that patient's nearest in-key line; patients whose nearest line is already in the key are unchanged.
Distances are the main run's (scripts/70_celligner.py), cancer models only; weights use the same softmax and tau
rule as scripts/71_celligner_outputs.py, computed over the key's lines (tau from the unrestricted distances).
Outputs (<p> = keyed for strict, hybrid for hybrid): results/celligner/<p>_top_lines.tsv (per patient, top 10),
         <p>_group_summary.tsv, <p>_line_coverage.tsv; data/processed/celligner/weights_<p>.tsv.gz (long form)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
DEPMAP = ROOT / "data/raw/depmap_24q4"
TOPK = 10
MODE = os.environ.get("MODE", "strict")
PFX = {"strict": "keyed", "hybrid": "hybrid"}[MODE]


def key_lines(row, models):
    m = models[models.OncotreeLineage.isin(str(row.lineages).split(";"))]
    if isinstance(row.diseases, str) and row.diseases:
        m = m[m.OncotreePrimaryDisease.isin(row.diseases.split(";"))]
    return list(m.index)


def main():
    dist = pd.read_csv(PROC / "dist_tumor_model.tsv.gz", sep="\t", index_col=0)
    models = pd.read_csv(DEPMAP / "Model.csv", index_col=0)
    models = models[models.index.isin(dist.columns) & (models.OncotreePrimaryDisease != "Non-Cancerous")]
    dist = dist[models.index]
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    qc = pd.read_csv(RES / "patient_qc.tsv", sep="\t", index_col=0)  # sample -> group, unrestricted nearest

    ds = np.sort(dist.values, axis=1)
    tau = float(os.environ.get("TAU", np.median(ds[:, TOPK - 1] - ds[:, 0]) / np.log(10)))
    mcols = ["CellLineName", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeCode"]

    top, wlong, summ = [], [], []
    for gid, g in qc.groupby("group_id"):
        k = keys.loc[gid]
        lines = key_lines(k, models)
        if not lines:
            raise ValueError(f"{gid}: key {k.key} has no lines")
        if MODE == "strict":
            d = dist.loc[g.index, lines]
        else:  # key lines + any line closer than the patient's nearest in-key line (inf = not a candidate)
            full = dist.loc[g.index]
            near_key = full[lines].min(axis=1).values
            keep = full.columns.isin(lines)[None, :] | (full.values < near_key[:, None])
            d = full.where(keep, np.inf)
            d = d.loc[:, np.isfinite(d.values).any(axis=0)]
        o = np.argsort(d.values, axis=1)
        dsort = np.take_along_axis(d.values, o, axis=1)
        w = np.exp(-(d.values - dsort[:, [0]]) / tau)
        w /= w.sum(axis=1, keepdims=True)
        out_key = ~d.columns.isin(lines)
        wsort = np.take_along_axis(w, o, axis=1)
        for r in range(min(TOPK, len(lines))):
            top.append(pd.DataFrame({"sample": d.index, "patient": g.patient.values, "group_id": gid, "key": k.key,
                                     "rank": r + 1, "ModelID": d.columns[o[:, r]], "distance": dsort[:, r].round(3),
                                     "weight": wsort[:, r].round(4)}))
        wl = pd.DataFrame(w, index=d.index, columns=d.columns).stack().rename("weight").reset_index()
        wl.columns = ["sample", "ModelID", "weight"]
        wlong.append(wl[wl.weight >= 1e-4].assign(group_id=gid))
        wm = pd.Series(w.mean(axis=0), index=d.columns).sort_values(ascending=False)
        n_eff = 1 / (w ** 2).sum(axis=1)
        in_key = g.nearest_line.isin(models.loc[lines, "CellLineName"])
        summ.append(dict(
            group_id=gid, key=k.key, key_type=k.key_type, n_patients=len(g), n_key_lines=len(lines),
            pct_unrestricted_nearest_in_key=round(100 * in_key.mean(), 1),
            median_d1_unrestricted=round(g.d1.median(), 2), median_d1_in_key=round(np.median(dist.loc[g.index, lines].min(axis=1)), 2),
            median_extra_distance=round(np.median(dist.loc[g.index, lines].min(axis=1).values - g.d1.values), 2),
            median_n_eff_lines=round(np.median(n_eff), 1),
            top_weight_share=round(wm.iloc[0], 3),
            pct_patients_with_out_of_key=round(100 * (w[:, out_key].sum(axis=1) > 0).mean(), 1) if out_key.any() else 0.0,
            mean_out_of_key_weight=round(w[:, out_key].sum(axis=1).mean(), 3) if out_key.any() else 0.0,
            top_lines_by_mean_weight=", ".join(f"{models.CellLineName[m]} {v:.3f}" for m, v in wm.head(5).items())))

    top = pd.concat(top)
    top = top[np.isfinite(top.distance)].join(models[mcols], on="ModelID")  # hybrid: fewer than 10 candidates.sort_values(["group_id", "sample", "rank"])
    top.to_csv(RES / f"{PFX}_top_lines.tsv", sep="\t", index=False)
    pd.concat(wlong).round(6).to_csv(PROC / f"weights_{PFX}.tsv.gz", sep="\t", index=False)
    summ = pd.DataFrame(summ)
    summ.to_csv(RES / f"{PFX}_group_summary.tsv", sep="\t", index=False)
    cov = (top.groupby(["key", "ModelID"]).agg(n_patients_top10=("sample", "nunique"), n_rank1=("rank", lambda r: int((r == 1).sum())),
                                              groups=("group_id", lambda x: ",".join(sorted(set(x)))))
           .reset_index().join(models[mcols], on="ModelID").sort_values(["key", "n_patients_top10"], ascending=[True, False]))
    cov.to_csv(RES / f"{PFX}_line_coverage.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(f"mode = {MODE}; tau = {tau:.3f}")
    print(summ.drop(columns="top_lines_by_mean_weight").to_string(index=False))
    print(summ[["group_id", "top_lines_by_mean_weight"]].to_string(index=False))


if __name__ == "__main__":
    main()
