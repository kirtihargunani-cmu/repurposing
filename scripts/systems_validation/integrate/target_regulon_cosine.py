#!/usr/bin/env python
"""Target-restricted nested cosine: does a drug oppose ITS OWN TARGET'S regulon, as that regulon is
wired in a given disease - and in WHICH disease does it do so most?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

MOTIVATION (user, 2026-11-19). D149 showed regulon rewiring across diseases is signal, not noise
(same-cohort noise floor Jaccard 0.552 vs 0.029 across tissues). The first nested-cosine attempt
averaged over ~1,000 regulons per cohort and simply re-derived the production score (Spearman 0.77
with rev) while the anti-estrogen control scattered from the 9th to the 78th percentile - a generic
proliferation signature from a thousand regulators drowned the target-specific one.

This version asks the narrower question the hypothesis actually makes: take the drug's ANNOTATED
TARGET, take THAT regulator's regulon AS WIRED IN THIS COHORT, and ask whether this drug opposes
those genes more than other drugs do.

  disease vector  v_d = tumour-vs-normal z over the target's regulon genes in this cohort
  drug vector     v_x = drug response z over THE SAME genes, in a cohort-matched line
  cos             = cosine(v_d, v_x);  NEGATIVE = opposes the disease pattern gene by gene

The headline is not the raw cosine but `opposing_pct`: this drug's rank among EVERY drug scored
against the SAME regulon in the SAME cohort. That controls for regulons which are easy to oppose -
a generic cytotoxic opposing everything scores no better than chance on this statistic.

The repurposing read-out is then: for one drug and one target, which cohorts rank it highest.

Output: results/systems_validation/integrated/target_regulon_cosine.tsv.gz + run json
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/integrated"
PROC = ROOT / "data/processed/lincs"
MIN_GENES = 15
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def net_for(cohort, ovr, available):
    """The full set has 29 PARENT networks only - no subtype splits. Fall back to the parent."""
    for cand in ([ovr[cohort]] if cohort in ovr else []) + \
                (["COADREAD"] if cohort in ("TCGA-COAD", "TCGA-READ") else []) + \
                [cohort.split("-", 1)[1], cohort.split("-", 1)[1].split("_")[0],
                 ovr.get(cohort, "").split("_")[0]]:
        if cand in available:
            return cand
    return None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    pid2name = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)
    name2pids = annot.groupby("cmap_name").pert_id.apply(list).to_dict()
    units = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    u1 = units[units.k == 1][["cmap_name", "genes", "orient"]].drop_duplicates()
    M = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype="float32", mode="r",
                  shape=(len(idx), len(genes)))
    gpos = {g: i for i, g in enumerate(genes.gene_id)}
    # FULL regulator set (TFs + kinases + signalling). The TF-only set makes only ~5% of annotated
    # drug targets regulators; the full set makes ~55% (D150), because oncology targets are kinases
    # and receptors, not transcription factors.
    E = pd.read_csv(PROC / "trace_aracne/aracne_full_regulon_edges.tsv.gz", sep="\t")
    REG = {n: {r: list(d.gene_id) for r, d in g.groupby("regulator")} for n, g in E.groupby("network")}
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohorts = sorted({p.name.split("__")[0] for p in (ROOT / "data/processed/step4_v2").glob("*__gtex.tsv.gz")})

    rows = []
    for cohort in cohorts:
        net = net_for(cohort, ovr, REG)
        regs = REG.get(net) if net else None
        if not regs:
            continue
        DE = pd.read_csv(ROOT / f"data/processed/step4_v2/{cohort}__gtex.tsv.gz", sep="\t")
        dz = DE.set_index("gene_id").z.dropna()
        lines = lw[(lw.group_id == cohort) & lw.in_key].ModelID.unique()
        ri = idx.index[idx.ModelID.isin(lines)].to_numpy()
        if len(ri) < 50:
            continue
        X = np.asarray(M[ri], dtype=np.float32)
        meta = idx.loc[ri].reset_index(drop=True)
        meta["cmap_name"] = meta.pert_id.map(pid2name)
        # which annotated targets are regulators in THIS cohort's network
        tg = u1[u1.genes.isin(regs.keys())]
        for target, sub in tg.groupby("genes"):
            gl = [g for g in regs[target] if g in gpos and g in dz.index]
            if len(gl) < MIN_GENES:
                continue
            cols = np.array([gpos[g] for g in gl])
            vd = dz[gl].values.astype(np.float64)
            nd = np.linalg.norm(vd)
            if nd == 0:
                continue
            S = X[:, cols].astype(np.float64)
            nx = np.linalg.norm(S, axis=1)
            ok = nx > 0
            cos = np.full(len(S), np.nan)
            cos[ok] = (S[ok] @ vd) / (nx[ok] * nd)
            # percentile among ALL drugs scored against this same regulon (lower cos = more opposing)
            valid = np.flatnonzero(ok)
            order = np.argsort(cos[valid])
            pct = np.empty(len(valid))
            pct[order] = 100 * (1 - np.arange(len(valid)) / max(len(valid) - 1, 1))
            pmap = dict(zip(valid, pct))
            for drug in sub.cmap_name.unique():
                pids = set(name2pids.get(drug, []))
                hit = np.flatnonzero(meta.pert_id.isin(pids).values)
                for h in hit:
                    if h not in pmap:
                        continue
                    rows.append(dict(cohort=cohort, network=net, target=target, cmap_name=drug,
                                     pert_id=meta.pert_id.iloc[h], ModelID=meta.ModelID.iloc[h],
                                     n_genes=len(gl), cos=round(float(cos[h]), 5),
                                     opposing_pct=round(float(pmap[h]), 2),
                                     n_drugs_scored=len(valid),
                                     orient=sub[sub.cmap_name == drug].orient.iloc[0]))
        print(f"{cohort:22s} {net:12s} rows so far {len(rows):,}  [{time.time()-t0:.0f}s]", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "target_regulon_cosine_full.tsv.gz", sep="\t", index=False)
    print(f"\n{len(df):,} (cohort, target, drug, line) rows; {df.cohort.nunique()} cohorts, "
          f"{df.target.nunique()} targets, {df.cmap_name.nunique()} drugs")
    (OUT / "target_regulon_cosine_full_run.json").write_text(json.dumps(dict(
        step="target-restricted nested regulon cosine",
        statistic="opposing_pct = rank among all drugs scored against the SAME regulon in the SAME cohort",
        min_genes=MIN_GENES, n_rows=int(len(df)), n_cohorts=int(df.cohort.nunique()),
        n_targets=int(df.target.nunique()), n_drugs=int(df.cmap_name.nunique()),
        runtime_s=round(time.time() - t0, 1)), indent=1))
    print("wrote", OUT / "target_regulon_cosine_full.tsv.gz")


if __name__ == "__main__":
    sys.exit(main())
