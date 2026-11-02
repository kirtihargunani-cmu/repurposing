#!/usr/bin/env python
"""Where healthy tissue lands among cell lines, per disease cohort (docs/notes/SPEC_celligner_normals.md, D62 on main).

GTEx samples: the v2 healthy reference of the disease side (main: results/step4_v2/gtex_reference.tsv, use == TRUE;
one sample per donor per tissue, colon / kidney sub-region filter), all donors, read with `git show main:...`.
TCGA adjacent normals: the existing adjacent arm (scripts/72, CONTROL=adjacent; all 737 projected in one run).

Projection (as scripts/72): the main run's Celligner model; all main-run tumors + a chunk of GTEx samples as ONE
batch, transform(compute_cPCs=False) against the main run's cPC-regressed lines. Celligner centers a batch on its own
mean, so normals must stay a small share: chunks of <= CHUNK samples (default 350, ~3% of the batch), balanced by
tissue. A fixed seeded set of N_ANCHOR GTEx samples is added to every chunk to measure chunk-to-chunk drift
(accept: same nearest line >= 90% and nearest-distance Spearman >= 0.98 across chunks).

Weights: per sample, softmax over negative distance with the main run's tau (median over tumors of (d10 - d1) / ln 10),
in hybrid mode as scripts/73 (the group's tissue-key lines plus any line closer than the sample's nearest key line);
cancer models only. A group's normal weight vector = mean over its normals; the tumor vector = the disease line-weight
matrix used for drug matching (results/celligner/disease_line_weights_long.tsv).

Per group (>= MIN_N normals): nearest-line distance tumors vs normals (AUC = P(normal farther)); shared weight
sum_l min(w_tumor, w_normal) (1 = same lines in the same proportions); lineage composition (weights summed by DepMap
OncotreeLineage) and its total variation distance; organ-driven lines (tumor weight >= 5% and normal weight >= half of
it) and cancer-specific lines (tumor weight >= 5%, normal weight < 1%).
Outputs: results/celligner/normal_projection_{gtex_v2_samples,group,lineage,drift}.tsv;
         data/processed/celligner/normal_projection_gtex_v2_{dist,weights}.tsv.gz
"""
import io
import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
PROC = ROOT / "data/processed/celligner"
RES = ROOT / "results/celligner"
CHUNK = int(os.environ.get("CHUNK", "350"))
N_ANCHOR, MIN_N, TOPK, SEED = 50, 8, 10, 0


def load_toil(samples, genes):
    tpm_f = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
    use = [s for s in dict.fromkeys(samples) if s in header]
    chunks = []
    for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + use, chunksize=5000, index_col=0):
        ch.index = ch.index.str.split(".").str[0]
        chunks.append(ch[ch.index.isin(genes)])
    x = pd.concat(chunks)
    return np.log2(np.power(2.0, x) - 0.001 + 1).clip(lower=0).T  # log2(TPM+0.001) -> log2(TPM+1)


def hybrid_weights(d, key_lines, tau):
    """d: samples x lines distances; key lines + any line closer than the sample's nearest key line."""
    near_key = d[key_lines].min(axis=1).values
    keep = d.columns.isin(key_lines)[None, :] | (d.values < near_key[:, None])
    dv = np.where(keep, d.values, np.inf)
    w = np.exp(-(dv - np.nanmin(dv, axis=1, keepdims=True)) / tau)
    w /= w.sum(axis=1, keepdims=True)
    return pd.DataFrame(w, index=pd.Index(d.index, name="sample"), columns=pd.Index(d.columns, name="ModelID"))


def main():
    from celligner import Celligner

    t0 = time.time()
    rng = np.random.default_rng(SEED)
    ref = pd.read_csv(io.StringIO(subprocess.run(["git", "show", "main:results/step4_v2/gtex_reference.tsv"], cwd=ROOT,
                                                 capture_output=True, text=True, check=True).stdout), sep="\t")
    ref = ref[ref.use]
    models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0)
    main_dist = pd.read_csv(PROC / "dist_tumor_model.tsv.gz", sep="\t", index_col=0)
    cancer = [m for m in main_dist.columns if m in set(models.index[models.OncotreePrimaryDisease != "Non-Cancerous"])]
    ds = np.sort(main_dist[cancer].values, axis=1)
    tau = float(np.median(ds[:, TOPK - 1] - ds[:, 0]) / np.log(10))

    c0 = Celligner().load(PROC / "model.pkl")
    genes, tumors, lines = list(c0.ref_input.columns), list(c0.target_input.index), list(c0.ref_input.index)
    qc = pd.read_csv(RES / "patient_qc.tsv", sep="\t", index_col=0)
    ref = ref.assign(groups=ref.groups.str.split(";")).explode("groups").rename(columns={"groups": "group_id"})
    ref = ref[ref.group_id.isin(set(qc.group_id))]
    gtex = sorted(set(ref["sample"]))
    x = load_toil(tumors + gtex, set(genes))[genes]
    gtex = [s for s in gtex if s in x.index]
    tis = ref.drop_duplicates("sample").set_index("sample").tissue.reindex(gtex)

    # tissue-balanced chunks + fixed anchors
    anchors = list(rng.choice(gtex, N_ANCHOR, replace=False))
    rest = [s for s in tis.sort_values(kind="stable").index if s not in set(anchors)]
    n_chunks = int(np.ceil(len(rest) / CHUNK))
    chunks = [rest[i::n_chunks] for i in range(n_chunks)]
    print(f"{len(gtex)} GTEx samples in {n_chunks} chunks of ~{len(chunks[0])} + {N_ANCHOR} anchors; "
          f"normals {100 * (len(chunks[0]) + N_ANCHOR) / (len(tumors) + len(chunks[0]) + N_ANCHOR):.1f}% of each batch; tau {tau:.3f}", flush=True)

    dist_parts, drift, tum_d1 = [], [], None
    for i, ch in enumerate(chunks):
        c = Celligner().load(PROC / "model.pkl")
        c.ref_input = c.combined_output.loc[lines, genes]
        samples = tumors + ch + anchors
        c.transform(x.loc[samples], compute_cPCs=False)
        c.computeMetricsForOutput(model_ids=lines, tumor_ids=samples)
        d = c.tumor_CL_dist[cancer]
        md = main_dist.loc[tumors, cancer]
        rep = dict(chunk=i + 1, n_normals=len(ch) + N_ANCHOR,
                   tumor_d1_spearman_vs_main=round(pd.Series(d.loc[tumors].min(axis=1)).rank().corr(md.min(axis=1).rank()), 4),
                   tumor_same_nearest_vs_main=round((d.loc[tumors].idxmin(axis=1) == md.idxmin(axis=1)).mean(), 4))
        if i == 0:
            tum_d1 = d.loc[tumors]
        a = d.loc[anchors]
        drift.append((rep, a))
        dist_parts.append(d.loc[ch] if i else d.loc[ch + anchors])
        print(rep, f"{(time.time() - t0) / 60:.1f} min", flush=True)

    a0 = drift[0][1]
    pd.concat([a.assign(chunk=r["chunk"]) for r, a in drift]).round(4).to_csv(
        PROC / "normal_projection_anchor_dist.tsv.gz", sep="\t", index_label="sample")
    soft = lambda a: (lambda e: e / e.sum(axis=1, keepdims=True))(np.exp(-(a.values - a.values.min(axis=1, keepdims=True)) / tau))
    w0 = soft(a0)
    drows = []
    for rep, a in drift:
        rep.update(anchor_same_nearest_vs_chunk1=round((a.idxmin(axis=1) == a0.idxmin(axis=1)).mean(), 4),
                   anchor_d1_spearman_vs_chunk1=round(a.min(axis=1).rank().corr(a0.min(axis=1).rank()), 4),
                   # weight-level stability: per anchor, shared softmax weight with chunk 1 (1 = identical)
                   anchor_median_shared_weight_vs_chunk1=round(float(np.median(np.minimum(soft(a), w0).sum(axis=1))), 4))
        drows.append(rep)
    drift_df = pd.DataFrame(drows)
    drift_df["accept"] = (drift_df.anchor_same_nearest_vs_chunk1 >= 0.9) & (drift_df.anchor_d1_spearman_vs_chunk1 >= 0.98)
    drift_df.to_csv(RES / "normal_projection_drift.tsv", sep="\t", index=False)
    dist = pd.concat(dist_parts)
    dist.round(4).to_csv(PROC / "normal_projection_gtex_v2_dist.tsv.gz", sep="\t", index_label="sample")

    # weights (hybrid, per group key)
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    cm = models.loc[cancer]
    def key_lines(gid):
        k = keys.loc[gid]
        m = cm[cm.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        return list(m.index)
    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    adj = cov[(cov.role == "adjacent_normal") & cov.primary_per_patient & cov.group_id.isin(set(qc.group_id))]
    adj_dist = pd.read_csv(PROC / "normal_control_adjacent_dist.tsv.gz", sep="\t", index_col=0)[cancer]
    adj = adj[adj["sample"].isin(adj_dist.index)]
    tw = pd.read_csv(RES / "disease_line_weights_long.tsv", sep="\t")
    lin = cm.OncotreeLineage

    smp_rows, grp_rows, lin_rows, wl = [], [], [], []
    for gid in sorted(set(ref.group_id)):
        kl = key_lines(gid)
        gs = [s for s in ref.loc[ref.group_id == gid, "sample"] if s in dist.index]
        if not gs:
            continue
        W = hybrid_weights(dist.loc[gs], kl, tau)
        wl.append(W.stack().rename("weight").reset_index()
                  .query("weight >= 1e-4").assign(group_id=gid))
        d1n = dist.loc[gs].min(axis=1)
        nl = dist.loc[gs].idxmin(axis=1)
        smp_rows.append(pd.DataFrame({"sample": gs, "group_id": gid, "tissue": tis.reindex(gs).values, "d1": d1n.round(3).values,
                                      "nearest_line": cm.CellLineName.reindex(nl).values, "nearest_lineage": lin.reindex(nl).values}))
        if len(gs) < MIN_N:
            continue
        gt = qc.index[qc.group_id == gid]
        d1t = tum_d1.loc[gt].min(axis=1)
        wt = tw[tw.group_id == gid].set_index("ModelID").weight.reindex(cancer).fillna(0)
        wn = W.mean(axis=0).reindex(cancer).fillna(0)
        arms = {"gtex": wn}
        ga = adj.loc[adj.group_id == gid, "sample"].tolist()
        if len(ga) >= MIN_N:
            arms["adjacent"] = hybrid_weights(adj_dist.loc[ga], kl, tau).mean(axis=0).reindex(cancer).fillna(0)
        r = dict(group_id=gid, n_gtex=len(gs), n_adjacent=len(ga), n_tumors=len(gt),
                 median_d1_tumor=round(d1t.median(), 2), median_d1_gtex=round(d1n.median(), 2),
                 auc_gtex_farther=round(mannwhitneyu(d1n, d1t).statistic / (len(d1n) * len(d1t)), 3))
        lt = wt.groupby(lin).sum()
        for arm, w in arms.items():
            lw = w.groupby(lin).sum()
            r[f"shared_weight_{arm}"] = round(np.minimum(wt, w).sum(), 3)
            r[f"lineage_tvd_{arm}"] = round(0.5 * (lt - lw.reindex(lt.index).fillna(0)).abs().sum(), 3)
            if arm == "gtex":
                big = wt[wt >= 0.05]
                r["organ_driven_lines"] = ", ".join(f"{cm.CellLineName[m]} ({wt[m]:.2f}/{w[m]:.2f})" for m in big.index if w[m] >= 0.5 * wt[m])
                r["cancer_specific_lines"] = ", ".join(f"{cm.CellLineName[m]} ({wt[m]:.2f})" for m in big.index if w[m] < 0.01)
                r["top_gtex_lines"] = ", ".join(f"{cm.CellLineName[m]} {v:.2f}" for m, v in w.sort_values(ascending=False).head(3).items())
            lin_rows.append(pd.DataFrame({"group_id": gid, "arm": arm, "lineage": lt.index.union(lw.index)})
                            .assign(tumor_weight=lambda q: q.lineage.map(lt).fillna(0).round(4),
                                    normal_weight=lambda q: q.lineage.map(lw).fillna(0).round(4)))
        r["top_tumor_lines"] = ", ".join(f"{cm.CellLineName[m]} {v:.2f}" for m, v in wt.sort_values(ascending=False).head(3).items())
        grp_rows.append(r)

    pd.concat(smp_rows).to_csv(RES / "normal_projection_gtex_v2_samples.tsv", sep="\t", index=False)
    pd.concat(wl).round(6).to_csv(PROC / "normal_projection_gtex_v2_weights.tsv.gz", sep="\t", index=False)
    grp = pd.DataFrame(grp_rows).sort_values("shared_weight_gtex")
    grp.to_csv(RES / "normal_projection_group.tsv", sep="\t", index=False)
    L = pd.concat(lin_rows)
    L = L[(L.tumor_weight > 0) | (L.normal_weight > 0)]
    L.to_csv(RES / "normal_projection_lineage.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(drift_df.to_string(index=False))
    print(grp[["group_id", "n_gtex", "n_adjacent", "auc_gtex_farther", "shared_weight_gtex", "lineage_tvd_gtex"]
              + [c for c in ["shared_weight_adjacent", "lineage_tvd_adjacent"] if c in grp]].to_string(index=False))
    print(f"NORMAL PROJECTION DONE in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
