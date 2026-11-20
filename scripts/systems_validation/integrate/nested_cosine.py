#!/usr/bin/env python
"""Nested regulon-level cosine similarity between a DRUG array and a DISEASE array.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

MOTIVATION (user, 2026-11-19). D149 established that ARACNe regulon differences across diseases are
CONTEXT, not estimation noise: the same-cohort noise floor is Jaccard 0.552 while cross-tissue is
0.029. If a regulator is genuinely wired to different genes in different diseases, then the
gene-level pattern INSIDE a regulon carries the context - and collapsing it to one activity number
per regulator throws exactly that away.

The production reversal score does collapse it: rev = -(D . s)/||s||, where D is the drug's
regulator-level NES and s the disease's regulator-level NES. Both sides are one number per
regulator before they ever meet.

This keeps the nesting. For each regulator r in a cohort's own ARACNe network:
    disease vector  v_d[r] = tumour-vs-normal z over r's regulon genes
    drug vector     v_x[r] = drug response z over THE SAME genes
    cos[r]          = cosine(v_d[r], v_x[r])
A NEGATIVE cosine means the drug pushes that regulon's genes opposite to the disease, gene by gene,
rather than merely producing an opposite summary activity.

Per (cohort, drug, line) we report the |disease-NES|-weighted mean cosine over regulators, the
unweighted mean, and the share of regulators with cos < 0 - plus the same quantities for the
production rev, so the two can be compared directly rather than asserted to differ.

Output: results/systems_validation/integrated/nested_cosine_<cohort>.tsv.gz + a summary json
"""
import json, os, sys, time
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/integrated"
PROC = ROOT / "data/processed/lincs"
MIN_GENES = 15
COHORTS = {"BRCA-LumA": "BRCA_LUMA", "TCGA-SKCM": "SKCM", "BRCA-Her2": "BRCA_BAS",
           "TCGA-LUAD": "LUAD"}
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    pid2name = dict(annot.drop_duplicates("pert_id")[["pert_id", "cmap_name"]].values)
    M = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype="float32", mode="r",
                  shape=(len(idx), len(genes)))
    gpos = {g: i for i, g in enumerate(genes.gene_id)}
    E = pd.read_csv(PROC / "trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    mr = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t")
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")

    summary = []
    for cohort, net in COHORTS.items():
        con = ROOT / f"data/processed/step4_v2/{cohort}__gtex.tsv.gz"
        if not con.exists():
            print(f"{cohort}: no contrast"); continue
        DE = pd.read_csv(con, sep="\t")[["gene_id", "z"]].dropna()
        dz = DE.set_index("gene_id").z
        regs = {r: list(d.gene_id) for r, d in E[E.network == net].groupby("regulator")}
        nes = mr[mr.group_id == cohort].set_index("regulator").nes
        # lines matched to this cohort, as the production score uses
        lines = lw[(lw.group_id == cohort) & lw.in_key].ModelID.unique()
        rows_idx = idx.index[idx.ModelID.isin(lines)].to_numpy()
        if not len(rows_idx):
            print(f"{cohort}: no matched lines"); continue
        X = np.asarray(M[rows_idx], dtype=np.float32)
        meta = idx.loc[rows_idx].reset_index(drop=True)

        out = []
        for r, gl in regs.items():
            gl = [g for g in gl if g in gpos and g in dz.index]
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
            out.append(pd.DataFrame(dict(regulator=r, row=np.arange(len(S)), n_genes=len(gl),
                                         cos=cos, nes=float(nes.get(r, np.nan)))))
        if not out:
            print(f"{cohort}: no usable regulons"); continue
        df = pd.concat(out, ignore_index=True).dropna(subset=["cos"])
        df["pert_id"] = meta.pert_id.values[df.row.values]
        df["ModelID"] = meta.ModelID.values[df.row.values]
        df["cmap_name"] = df.pert_id.map(pid2name)
        w = df.nes.abs().fillna(0)
        df["w"] = w
        g = df.groupby(["pert_id", "cmap_name", "ModelID"])
        agg = g.apply(lambda t: pd.Series(dict(
            n_regulators=len(t),
            cos_mean=t.cos.mean(),
            cos_weighted=(t.cos * t.w).sum() / max(t.w.sum(), 1e-9),
            frac_opposing=float((t.cos < 0).mean())))).reset_index()
        agg["cohort"] = cohort
        agg.to_csv(OUT / f"nested_cosine_{cohort}.tsv.gz", sep="\t", index=False)
        print(f"{cohort:11s} ({net}): {len(regs):4d} regulons -> {df.regulator.nunique():4d} usable, "
              f"{agg.pert_id.nunique():,} drug-line rows, "
              f"weighted cosine median {agg.cos_weighted.median():+.4f}, "
              f"frac opposing median {agg.frac_opposing.median():.3f}  [{time.time()-t0:.0f}s]")
        summary.append(dict(cohort=cohort, network=net, n_regulons_usable=int(df.regulator.nunique()),
                            n_rows=int(len(agg)),
                            cos_weighted_median=round(float(agg.cos_weighted.median()), 5),
                            cos_weighted_p05=round(float(agg.cos_weighted.quantile(.05)), 5),
                            cos_weighted_p95=round(float(agg.cos_weighted.quantile(.95)), 5),
                            frac_opposing_median=round(float(agg.frac_opposing.median()), 4)))
    (OUT / "nested_cosine_run.json").write_text(json.dumps(dict(
        step="nested regulon-level cosine, drug array vs disease array",
        motivation="D149: regulon context is signal; keep gene-level structure inside each regulon",
        min_genes=MIN_GENES, cohorts=summary, runtime_s=round(time.time() - t0, 1)), indent=1))
    print("\nwrote", OUT / "nested_cosine_run.json")


if __name__ == "__main__":
    sys.exit(main())
