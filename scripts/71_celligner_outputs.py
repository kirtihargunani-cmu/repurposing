#!/usr/bin/env python
"""Celligner outputs and checks: patient -> cell-line matches from the scripts/70_celligner.py alignment.

Distances: Euclidean in the 70-PC space of the aligned matrix (Celligner tumor_CL_dist). Non-cancerous DepMap
models (OncotreePrimaryDisease == "Non-Cancerous": fibroblasts, normal / immortalized lines) stay in the
alignment but are not candidate matches.
  weights   per patient softmax over negative distance, w_l = exp(-d_l / tau) / sum_j exp(-d_j / tau): a similarity
            weighting, not a mixture fraction (aligned values are MNN-shifted, so NNLS deconvolution is not meaningful).
            tau (provisional, env TAU overrides): median over patients of (d_10 - d_1) / ln 10, i.e. for the median
            patient the 10th-nearest line has 1/10 the weight of the nearest.
  lineage   group tissue = OncoTree tissue of the group's code; line tissue = DepMap OncotreeLineage (same vocabulary).
  purity    nearest-line distance vs ESTIMATE score, ABSOLUTE purity, blast % (AML) per group; if distance rises as
            purity falls, the alignment under-corrected low-purity samples (flagged, not corrected).
  lesions   TCGA patients only (MC3 PASS + GISTIC): oncogene hotspots (same positions as scripts/60_driver_classes.py,
            plus IDH1 R132 / IDH2 R140, R172) and high-level amplification (GISTIC 2; DepMap relative copy ratio >= 3).
            For carriers: fraction of their top-10 lines with the same lesion vs the same-tissue line background.
Outputs: results/celligner/{patient_top_lines, patient_qc, group_summary, line_coverage, purity_diagnostic,
         lesion_check}.tsv; data/processed/celligner/weights_softmax.tsv.gz
"""
import gzip
import importlib.util
import json
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
DEPMAP = RAW / "depmap_24q4"
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
TOPK = 10
AMP_RATIO = 3.0
AMP_GENES = ["ERBB2", "EGFR", "MYC", "MYCN", "CCND1", "CDK4", "MDM2", "MET", "FGFR1", "KRAS"]

spec = importlib.util.spec_from_file_location("dc", ROOT / "scripts/60_driver_classes.py")
dc = importlib.util.module_from_spec(spec); spec.loader.exec_module(dc)
HOTSPOTS = {g: h for g, (_, role, h) in dc.DRIVERS.items() if role == "onc" and h}
HOTSPOTS.update({"IDH1": {132}, "IDH2": {140, 172}})


def tumor_lesions(patients):
    maf = pd.read_csv(RAW / "tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=[0, 8, 15, 53, 108],
                      dtype=str, low_memory=False)
    maf.columns = ["gene", "vclass", "barcode", "pos", "filter"]
    maf["patient"] = maf.barcode.str[:12]
    profiled = set(maf.patient)
    maf = maf[(maf["filter"] == "PASS") & maf.gene.isin(HOTSPOTS) & maf.vclass.isin(dc.HOT)]
    maf["p1"] = pd.to_numeric(maf.pos.str.split("/").str[0].str.split("-").str[0], errors="coerce")
    mut = maf[[p in HOTSPOTS[g] for g, p in zip(maf.gene, maf.p1)]]
    mut = mut.assign(lesion=mut.gene + " hotspot")[["patient", "lesion"]]
    with gzip.open(RAW / "tcga_pancan/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", "rt") as fh:
        cn = pd.read_csv(fh, sep="\t", index_col=0)
    cn = cn.loc[cn.index.intersection(AMP_GENES)]
    cn.columns = cn.columns.str[:12]
    cn = cn.loc[:, ~cn.columns.duplicated()]
    amp = [(p, f"{g} amp") for g in cn.index for p in cn.columns[cn.loc[g] == 2]]
    les = pd.concat([mut, pd.DataFrame(amp, columns=["patient", "lesion"])]).drop_duplicates()
    covered = set(patients) & profiled & set(cn.columns)
    return les[les.patient.isin(covered)], covered


def line_lesions():
    m = pd.read_csv(DEPMAP / "OmicsSomaticMutations.csv", usecols=["ModelID", "HugoSymbol", "ProteinChange",
                    "VariantInfo"], low_memory=False)
    m = m[m.HugoSymbol.isin(HOTSPOTS) & m.VariantInfo.str.contains("missense_variant|inframe_", na=False)]
    m["p1"] = pd.to_numeric(m.ProteinChange.str.extract(r"^p\.[A-Za-z*]+?(\d+)")[0], errors="coerce")
    mut = m[[p in HOTSPOTS[g] for g, p in zip(m.HugoSymbol, m.p1)]]
    mut = mut.assign(lesion=mut.HugoSymbol + " hotspot")[["ModelID", "lesion"]]
    h = pd.read_csv(DEPMAP / "OmicsCNGene.csv", nrows=0).columns
    cols = [c for c in h[1:] if c.split(" ")[0] in AMP_GENES]
    cn = pd.read_csv(DEPMAP / "OmicsCNGene.csv", usecols=[h[0]] + cols, index_col=0)
    cn.columns = cn.columns.str.split(" ").str[0]
    amp = [(mid, f"{g} amp") for g in cn.columns for mid in cn.index[cn[g] >= AMP_RATIO]]
    return pd.concat([mut, pd.DataFrame(amp, columns=["ModelID", "lesion"])]).drop_duplicates()


def rho(a, b):
    ok = a.notna() & b.notna()
    if ok.sum() < 10:
        return np.nan, np.nan, int(ok.sum())
    n = int(ok.sum())
    r = a[ok].rank().corr(b[ok].rank())
    z = abs(r) * math.sqrt((n - 2) / max(1e-12, 1 - r * r))  # t statistic; normal tail (n >= 10)
    return round(r, 3), math.erfc(z / math.sqrt(2)), n


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1); return out


def main():
    dist = pd.read_csv(PROC / "dist_tumor_model.tsv.gz", sep="\t", index_col=0)
    models = pd.read_csv(DEPMAP / "Model.csv", index_col=0)
    cancer = models.index[models.OncotreePrimaryDisease != "Non-Cancerous"]
    dist = dist.loc[:, dist.columns.isin(cancer)]
    print(f"patients {dist.shape[0]}, candidate cancer models {dist.shape[1]}")

    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    cov = cov[cov["sample"].isin(dist.index) & (cov.role == "tumor") & (cov.primary_per_patient == True)]
    cov = cov.drop_duplicates("sample").set_index("sample").loc[dist.index]
    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    ot = {t["code"]: t["tissue"] for t in json.load(open(RAW / "oncotree/tumorTypes.json"))}
    gtissue = groups.oncotree.fillna("").str.split(";").str[0].map(ot)

    d = dist.values
    order = np.argsort(d, axis=1)
    ds = np.take_along_axis(d, order, axis=1)
    tau = float(os.environ.get("TAU", np.median(ds[:, TOPK - 1] - ds[:, 0]) / np.log(10)))
    w = np.exp(-(d - ds[:, [0]]) / tau)
    w /= w.sum(axis=1, keepdims=True)
    weights = pd.DataFrame(w, index=dist.index, columns=dist.columns)
    weights.round(6).to_csv(PROC / "weights_softmax.tsv.gz", sep="\t", index_label="sample")
    ws = np.take_along_axis(w, order, axis=1)
    # effective number of lines per patient (inverse Simpson of the weights)
    n_eff = 1 / (w ** 2).sum(axis=1)

    mcols = ["CellLineName", "OncotreeLineage", "OncotreeCode", "OncotreePrimaryDisease"]
    lin = models.OncotreeLineage
    top = []
    for k in range(TOPK):
        mid = dist.columns[order[:, k]]
        top.append(pd.DataFrame({"sample": dist.index, "patient": cov.patient.values, "group_id": cov.group_id.values,
                                 "rank": k + 1, "ModelID": mid, "distance": ds[:, k].round(3),
                                 "weight": ws[:, k].round(4)}))
    top = pd.concat(top).sort_values(["group_id", "sample", "rank"])
    top = top.join(models[mcols], on="ModelID")
    top["tissue_match"] = top.OncotreeLineage.values == top.group_id.map(gtissue).values
    top.to_csv(RES / "patient_top_lines.tsv", sep="\t", index=False)

    # patient QC
    t1 = top[top["rank"] == 1].set_index("sample")
    qc = pd.DataFrame({
        "patient": cov.patient, "group_id": cov.group_id, "cohort": cov.cohort,
        "d1": ds[:, 0].round(3), "d10": ds[:, TOPK - 1].round(3), "n_eff_lines": n_eff.round(2),
        "nearest_line": t1.CellLineName.reindex(dist.index).values,
        "nearest_lineage": t1.OncotreeLineage.reindex(dist.index).values,
        "group_tissue": cov.group_id.map(gtissue).values,
        "frac_top10_tissue_match": top.groupby("sample").tissue_match.mean().reindex(dist.index).round(2).values,
        "estimate_score": cov.estimate_score, "absolute_purity": cov.absolute_purity, "blast_pct": cov.blast_pct,
    }, index=dist.index)
    qc["nearest_tissue_match"] = qc.nearest_lineage == qc.group_tissue
    # far from every line relative to the group, or low tumor content -> do not trust the match as-is
    q90 = qc.groupby("group_id").d1.transform(lambda x: x.quantile(0.9))
    qc["flag_far"] = qc.d1 > q90
    qc["flag_low_purity"] = (qc.absolute_purity < 0.3) | (qc.blast_pct < 60)
    qc.to_csv(RES / "patient_qc.tsv", sep="\t", index_label="sample")

    # purity diagnostic per group
    pur = []
    for gid, g in qc.groupby("group_id"):
        r_est = rho(g.d1, g.estimate_score); r_abs = rho(g.d1, g.absolute_purity); r_bl = rho(g.d1, g.blast_pct)
        pur.append(dict(group_id=gid, n=len(g), rho_d1_estimate=r_est[0], p_estimate=r_est[1],
                        rho_d1_absolute=r_abs[0], p_absolute=r_abs[1], n_absolute=r_abs[2],
                        rho_d1_blast=r_bl[0], n_blast=r_bl[2], n_low_purity=int(g.flag_low_purity.sum())))
    pur = pd.DataFrame(pur)
    pur["fdr_absolute"] = np.nan
    ok = pur.p_absolute.notna()
    if ok.any():
        pur.loc[ok, "fdr_absolute"] = bh(pur.loc[ok, "p_absolute"])
    pur.to_csv(RES / "purity_diagnostic.tsv", sep="\t", index=False, float_format="%.4g")

    # line coverage
    t10 = top
    cov_l = t10.groupby("ModelID").agg(n_patients_top10=("sample", "nunique"), n_groups=("group_id", "nunique"),
                                       n_rank1=("rank", lambda r: int((r == 1).sum())))
    tg = t10.groupby(["ModelID", "group_id"]).size().rename("n").reset_index().sort_values("n", ascending=False)
    cov_l["top_groups"] = tg.groupby("ModelID").apply(
        lambda x: ", ".join(f"{a} {b}" for a, b in zip(x.group_id.head(4), x.n.head(4))), include_groups=False)
    cov_l = models[mcols].join(cov_l, how="right").sort_values("n_patients_top10", ascending=False)
    allm = pd.DataFrame(index=dist.columns).join(models[mcols]).join(cov_l[["n_patients_top10", "n_groups", "n_rank1"]])
    allm = allm.fillna({"n_patients_top10": 0, "n_groups": 0, "n_rank1": 0})
    allm.sort_values("n_patients_top10", ascending=False).to_csv(RES / "line_coverage.tsv", sep="\t",
                                                                  index_label="ModelID")
    print(f"lines never in any patient's top {TOPK}: {(allm.n_patients_top10 == 0).sum()} of {len(allm)}")

    # group summary
    n_tissue_lines = models.loc[dist.columns].OncotreeLineage.value_counts()
    rows = []
    for gid, g in qc.groupby("group_id"):
        wm = weights.loc[g.index].mean().sort_values(ascending=False)
        t1g = g.nearest_line.value_counts()
        rows.append(dict(group_id=gid, n_patients=len(g), tissue=gtissue.get(gid),
                         n_lines_in_tissue=int(n_tissue_lines.get(gtissue.get(gid), 0)),
                         median_d1=round(g.d1.median(), 2), median_n_eff=round(g.n_eff_lines.median(), 1),
                         pct_nearest_tissue_match=round(100 * g.nearest_tissue_match.mean(), 1),
                         mean_frac_top10_tissue_match=round(g.frac_top10_tissue_match.mean(), 2),
                         n_distinct_nearest=len(t1g),
                         top_nearest=", ".join(f"{k} {v}" for k, v in t1g.head(5).items()),
                         top_by_mean_weight=", ".join(f"{models.CellLineName.get(k)} {v:.3f}"
                                                      for k, v in wm.head(5).items())))
    gs = pd.DataFrame(rows).merge(pur[["group_id", "rho_d1_estimate", "rho_d1_absolute", "fdr_absolute"]])
    gs.to_csv(RES / "group_summary.tsv", sep="\t", index=False)

    # lesion check (TCGA)
    les_t, covered = tumor_lesions(set(qc.loc[qc.cohort.str.startswith("TCGA"), "patient"]))
    les_l = line_lesions()
    les_l = les_l[les_l.ModelID.isin(dist.columns)]
    line_has = les_l.groupby("lesion").ModelID.apply(set).to_dict()
    pt2s = qc.reset_index().set_index("patient")["sample"]
    topsets = t10.groupby("sample").ModelID.apply(list)
    rows = []
    for (gid, lesion), g in les_t.assign(group_id=les_t.patient.map(qc.reset_index().set_index("patient").group_id)) \
            .dropna(subset=["group_id"]).groupby(["group_id", "lesion"]):
        gp = set(qc.loc[(qc.group_id == gid) & qc.patient.isin(covered), "patient"])
        car = set(g.patient) & gp
        if len(car) < 5:
            continue
        lines_with = line_has.get(lesion, set())
        frac = lambda pts: np.mean([np.mean([m in lines_with for m in topsets[pt2s[p]]]) for p in pts]) if pts else np.nan
        tis = gtissue.get(gid)
        tis_lines = models.index[(models.OncotreeLineage == tis) & models.index.isin(dist.columns)]
        rows.append(dict(group_id=gid, lesion=lesion, n_profiled=len(gp), n_carriers=len(car),
                         carriers_top10_frac_with=round(frac(car), 3),
                         noncarriers_top10_frac_with=round(frac(gp - car), 3),
                         tissue_lines_frac_with=round(np.mean([m in lines_with for m in tis_lines]), 3) if len(tis_lines) else np.nan,
                         all_lines_frac_with=round(np.mean([m in lines_with for m in dist.columns]), 3),
                         pct_carriers_any_top10_with=round(100 * np.mean([any(m in lines_with for m in topsets[pt2s[p]])
                                                                          for p in car]), 1)))
    lc = pd.DataFrame(rows).sort_values(["group_id", "n_carriers"], ascending=[True, False])
    lc.to_csv(RES / "lesion_check.tsv", sep="\t", index=False)

    print(f"tau = {tau:.3f}; median effective lines per patient {np.median(n_eff):.1f}")
    print(gs[["group_id", "n_patients", "tissue", "n_lines_in_tissue", "median_d1", "pct_nearest_tissue_match",
              "mean_frac_top10_tissue_match", "rho_d1_absolute", "top_nearest"]].to_string(index=False))
    print(lc.to_string(index=False))


if __name__ == "__main__":
    main()
