#!/usr/bin/env python
"""Celligner negative controls: where do normal tissues land relative to the cell lines?

Question: is "tumor X is close to line Y" cancer-specific, or would the same tissue without cancer be just as close?
Normals are the baseline: a tumor's match is informative only to the extent the tumor sits closer to lines than
normal tissue of the same organ (ideally the same patient's).

Two arms (CONTROL env var; one run each, so neither normal set shifts the other):
  adjacent  all 737 adjacent normals (TCGA + TARGET); paired with the same patient's tumor. Same organ and
            person, but field effect / tumor contamination can make it tumor-like (conservative control).
  gtex      GTEx healthy donors of each group's Step 4 tissues (primary_per_patient rows), at most GTEX_CAP donors
            per group (seeded sample); unpaired. Truly healthy tissue, but post-mortem and a different collection.

Design (runs in the `celligner` env, reuses data/processed/celligner/model.pkl from scripts/70_celligner.py):
  target   all main-run tumors + the arm's normals, as ONE batch. One batch matters: Celligner mean-centers each
           dataset on its own mean, so normals projected alone would be re-centered onto the tumor mean and the
           control would be void. Normals are kept a small share of the batch so centering is ~ unchanged.
  method   transform(compute_cPCs=False): the main run's cPCs, DE genes and MNN settings (k1 5, k2 50) against the
           same DepMap reference; tumors and normals get identical treatment. Tumor positions should reproduce the
           main run (checked: nearest-line distance and identity).
  metrics  per group with >= 8 normals (adjacent: >= 8 same-patient pairs): nearest-line distance d1 of normals vs
           tumors (AUC = P(normal farther); adjacent also paired Wilcoxon on own tumor), fraction of the group's
           tumors closer to a line than the median normal, nearest-line tissue match, and whether normals pick the
           tumors' lines (adjacent: same nearest line / top-10 overlap with own tumor; both: share of normals whose
           nearest line is one of the group tumors' 5 most frequent nearest lines).
CONTROL=none re-projects the tumors alone (no normals) and only reports reproducibility vs the main run: the
baseline for how much of the arms' tumor drift is due to the normals rather than the projection itself.
Outputs (<arm> = adjacent | gtex): results/celligner/normal_control_<arm>.tsv (group), normal_control_<arm>_samples.tsv,
         normal_control_<arm>_info.tsv; data/processed/celligner/normal_control_<arm>_dist.tsv.gz
"""
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
DEPMAP = RAW / "depmap_24q4"
TOPK, MIN_N, SEED, GTEX_CAP = 10, 8, 0, 50
ARM = os.environ.get("CONTROL", "adjacent")


def load_toil(samples, genes):
    tpm_f = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
    use = [s for s in samples if s in header]
    chunks = []
    for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + use, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        chunks.append(ch[ch.index.isin(genes)])
    x = pd.concat(chunks)
    x = np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0)  # log2(TPM+0.001) -> log2(TPM+1)
    return x.T


def select_normals(cov, tumors):
    """Normal samples to project and their (sample, group_id[, patient]) rows for the analysis."""
    if ARM == "none":
        return [], pd.DataFrame(columns=["sample", "patient"])
    if ARM == "adjacent":
        nrm = cov[cov.role == "adjacent_normal"].drop_duplicates("sample")
        return list(nrm["sample"]), nrm[["sample", "patient"]]
    tcov = cov[cov["sample"].isin(tumors) & (cov.role == "tumor")]
    g = cov[(cov.role == "gtex_normal") & (cov.primary_per_patient == True) & cov.group_id.isin(set(tcov.group_id))]
    g = g.drop_duplicates(["group_id", "patient"])
    g = g.groupby("group_id", group_keys=False).apply(lambda d: d.sample(min(len(d), GTEX_CAP), random_state=SEED))
    return sorted(set(g["sample"])), g[["sample", "patient", "group_id"]]


def main():
    from celligner import Celligner

    assert ARM in ("adjacent", "gtex", "none"), ARM
    np.random.seed(SEED)
    t0 = time.time()
    c = Celligner().load(PROC / "model.pkl")
    genes = list(c.ref_input.columns)
    tumors = list(c.target_input.index)
    lines = list(c.ref_input.index)
    # transform(compute_cPCs=False) uses ref_input as-is (it expects makeNewReference() output), so give it the
    # main run's cPC-regressed lines (their rows of combined_output; MNN corrects only the target). With the raw
    # centered lines instead, tumors failed to reproduce the main run (same nearest line 51%).
    c.ref_input = c.combined_output.loc[lines, genes]
    main_dist = pd.read_csv(PROC / "dist_tumor_model.tsv.gz", sep="\t", index_col=0)

    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    normals, nrows = select_normals(cov, tumors)
    x = load_toil(tumors + normals, set(genes))[genes]
    normals = [s for s in normals if s in x.index]
    print(f"[{ARM}] projecting {len(tumors)} tumors + {len(normals)} normals onto {len(lines)} lines "
          f"(normals {100 * len(normals) / (len(tumors) + len(normals)):.1f}% of the batch)")

    c.transform(x.loc[tumors + normals], compute_cPCs=False)
    c.computeMetricsForOutput(model_ids=lines, tumor_ids=tumors + normals)
    dist = c.tumor_CL_dist
    dist.round(4).to_csv(PROC / f"normal_control_{ARM}_dist.tsv.gz", sep="\t", index_label="sample")

    # candidates as in scripts/71_celligner_outputs.py: cancer models only
    models = pd.read_csv(DEPMAP / "Model.csv", index_col=0)
    cancer = [m for m in dist.columns if m in set(models.index[models.OncotreePrimaryDisease != "Non-Cancerous"])]
    dist, main_dist = dist[cancer], main_dist[cancer]

    def nearest(d):
        o = np.argsort(d.values, axis=1)
        return (pd.Series(np.take_along_axis(d.values, o[:, :1], 1)[:, 0], index=d.index),
                pd.Series(d.columns[o[:, 0]], index=d.index),
                pd.Series([set(d.columns[r[:TOPK]]) for r in o], index=d.index))

    d1, l1, top = nearest(dist)
    md1, ml1, _ = nearest(main_dist.loc[tumors])
    repro = dict(tumor_d1_spearman_vs_main=round(d1[tumors].rank().corr(md1.rank()), 4),
                 tumor_same_nearest_line_vs_main=round((l1[tumors] == ml1).mean(), 4))
    md1_top3 = pd.Series([set(main_dist.columns[r[:3]]) for r in np.argsort(main_dist.loc[tumors].values, axis=1)],
                         index=tumors)
    repro["tumor_nearest_in_main_top3"] = round(np.mean([l in t for l, t in zip(l1[tumors], md1_top3)]), 4)
    print(repro)
    if ARM == "none":
        pd.Series({**repro, "arm": ARM, "minutes": round((time.time() - t0) / 60, 1)}).to_csv(
            RES / "normal_control_none_info.tsv", sep="\t", header=False)
        return

    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    ot = {t["code"]: t["tissue"] for t in json.load(open(RAW / "oncotree/tumorTypes.json"))}
    gtissue = groups.oncotree.fillna("").str.split(";").str[0].map(ot)
    lin = models.OncotreeLineage
    name = models.CellLineName

    tcov = cov[cov["sample"].isin(tumors) & (cov.role == "tumor")].drop_duplicates("sample")
    t_by_pt = tcov.drop_duplicates("patient").set_index("patient")
    nrows = nrows[nrows["sample"].isin(normals)].copy()
    if ARM == "adjacent":  # analyse same-patient pairs only; group = the patient's tumor group
        nrows = nrows.drop_duplicates("patient")
        nrows = nrows[nrows.patient.isin(t_by_pt.index)]
        nrows["tumor"] = t_by_pt.loc[nrows.patient, "sample"].values
        nrows["group_id"] = t_by_pt.loc[nrows.patient, "group_id"].values
    smp = nrows.copy()
    smp["d1_normal"] = d1[smp["sample"]].values.round(3)
    smp["nearest_normal"] = name.reindex(l1[smp["sample"]]).values
    smp["normal_tissue_match"] = lin.reindex(l1[smp["sample"]]).values == smp.group_id.map(gtissue).values
    if ARM == "adjacent":
        smp["d1_tumor"] = d1[smp.tumor].values.round(3)
        smp["normal_farther"] = smp.d1_normal > smp.d1_tumor
        smp["nearest_tumor"] = name.reindex(l1[smp.tumor]).values
        smp["same_nearest"] = l1[smp.tumor].values == l1[smp["sample"]].values
        smp["top10_overlap"] = [len(top[a] & top[b]) / TOPK for a, b in zip(smp.tumor, smp["sample"])]
    smp.to_csv(RES / f"normal_control_{ARM}_samples.tsv", sep="\t", index=False)

    rows = []
    for gid, p in smp.groupby("group_id"):
        if len(p) < MIN_N:
            continue
        gt = tcov.loc[tcov.group_id == gid, "sample"]
        gn = p["sample"]
        auc = mannwhitneyu(d1[gn], d1[gt]).statistic / (len(gn) * len(gt))  # P(normal d1 > tumor d1)
        top5 = set(l1[gt].value_counts().index[:5])
        r = dict(group_id=gid, n_normals=len(p), n_tumors=len(gt), median_d1_tumor=round(d1[gt].median(), 2),
                 median_d1_normal=round(d1[gn].median(), 2), auc_normal_farther=round(auc, 3),
                 pct_tumors_closer_than_median_normal=round(100 * (d1[gt] < d1[gn].median()).mean(), 1),
                 pct_tumor_tissue_match=round(100 * (lin.reindex(l1[gt]).values == gtissue.get(gid)).mean(), 1),
                 pct_normal_tissue_match=round(100 * p.normal_tissue_match.mean(), 1),
                 pct_normal_nearest_in_tumor_top5=round(100 * l1[gn].isin(top5).mean(), 1))
        if ARM == "adjacent":
            r.update(median_paired_diff=round((p.d1_normal - p.d1_tumor).median(), 2),
                     pct_normal_farther_than_own_tumor=round(100 * p.normal_farther.mean(), 1),
                     wilcoxon_p=wilcoxon(p.d1_normal, p.d1_tumor).pvalue,
                     pct_same_nearest_line=round(100 * p.same_nearest.mean(), 1),
                     mean_top10_overlap=round(p.top10_overlap.mean(), 2))
        r["top_nearest_normal"] = ", ".join(f"{k} {v}" for k, v in p.nearest_normal.value_counts().head(3).items())
        rows.append(r)
    out = pd.DataFrame(rows).sort_values("auc_normal_farther")
    out.to_csv(RES / f"normal_control_{ARM}.tsv", sep="\t", index=False, float_format="%.4g")

    info = {**repro, "arm": ARM, "n_normals_projected": len(normals), "n_normals_analysed": len(smp),
            "n_groups_analysed": len(out), "median_auc": round(out.auc_normal_farther.median(), 3),
            "minutes": round((time.time() - t0) / 60, 1)}
    if ARM == "adjacent":
        info.update(all_pct_normal_farther=round(100 * smp.normal_farther.mean(), 1),
                    all_wilcoxon_p=wilcoxon(smp.d1_normal, smp.d1_tumor).pvalue,
                    all_pct_same_nearest_line=round(100 * smp.same_nearest.mean(), 1))
    pd.Series(info).to_csv(RES / f"normal_control_{ARM}_info.tsv", sep="\t", header=False)
    pd.set_option("display.width", 250)
    print(pd.Series(info).to_string())
    print(out.drop(columns="top_nearest_normal").to_string(index=False))


if __name__ == "__main__":
    main()
