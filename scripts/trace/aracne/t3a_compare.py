#!/usr/bin/env python
"""TRACE-ARACNe stage 4: the CollecTRI vs ARACNe comparison for the two real pairs.

Replays the identical L5 computation on both frozen footprints, so any difference is the
regulon and nothing else. Also reports how much the two networks' footprints overlap.
Output: results/lincs/trace_aracne/network_comparison.tsv
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
PAIRS = [("fulvestrant", "ACH-000019", "MCF7", "fulvestrant_ESR1_BRCA_LumA"),
         ("vemurafenib", "ACH-000219", "A375", "vemurafenib_BRAF_TCGA_SKCM")]


def main():
    genes = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    M = np.memmap(ROOT / "data/processed/lincs/drug_line_profiles_bing_specific.f32",
                  dtype="float32", mode="r", shape=(len(idx), len(genes)))
    gpos = {g: i for i, g in enumerate(genes.gene_id)}
    out = []
    for name, model, line, key in PAIRS:
        fps = {}
        for src, d in (("collectri", "trace"), ("aracne", "trace_aracne")):
            fp = pd.read_csv(ROOT / f"results/lincs/{d}/predicted_footprint_{key}.tsv", sep="\t")
            fps[src] = fp[fp.expected_sign != 0].set_index("gene_id").expected_sign
        both = fps["collectri"].index.intersection(fps["aracne"].index)
        agree = (float((fps["collectri"][both].values == fps["aracne"][both].values).mean())
                 if len(both) else float("nan"))
        pid = set(annot[annot.cmap_name == name].pert_id)
        rows = idx[(idx.ModelID == model) & (idx.pert_id.isin(pid))].sort_values("n_sig", ascending=False)
        prow = int(rows.index[0])
        line_rows = idx.index[idx.ModelID == model].to_numpy()
        for src in ("collectri", "aracne"):
            s = fps[src]
            cols = np.array([gpos[g] for g in s.index if g in gpos])
            e = np.array([v for g, v in s.items() if g in gpos], dtype=float)
            z = np.asarray(M[prow][cols], dtype=float)
            mean = float((np.tanh(z / 2) * e).mean())
            Z = np.asarray(M[line_rows][:, cols], dtype=float)
            means = (np.tanh(Z / 2) * e).mean(axis=1)
            out.append(dict(pair=key, drug=name, line=line, network=src, n_signed=len(s),
                            n_covered=len(e),
                            concordance=round(float((np.sign(z) == e).mean()), 4),
                            mean_signed_fz=round(mean, 6),
                            generic_percentile=round(float((means < mean).mean() * 100), 1),
                            overlap_with_other=len(both), sign_agreement_on_overlap=round(agree, 4)))
    df = pd.DataFrame(out)
    df.to_csv(ROOT / "results/lincs/trace_aracne/network_comparison.tsv", sep="\t", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
