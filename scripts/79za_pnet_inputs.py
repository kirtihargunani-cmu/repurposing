#!/usr/bin/env python
"""D86 perturbation-network pilot: ARACNe inputs from LINCS drug profiles within one cell line.

Rule fixed 2026-11-02 before running (docs/DECISIONS.md D86). Lines MCF7, PC3, A375; arms raw and line-specific
drug x line consensus profiles (79c), one profile per compound (BRD id); genes = the 978 L1000 landmark genes;
regulators = landmark genes that are Lambert TFs (the v2 TF regulator list, data/processed/step5/aracne_inputs_v2/).
Folds: compounds grouped by mechanism (Repurposing Hub MOA, else LINCS MOA, else compound name), groups assigned to
5 folds greedily by size (largest first, ties in seeded random order) so drugs sharing a mechanism leave together;
the same group -> fold map is used in every line. Inputs per line x arm: "all" (every compound) and f1..f5 (every
compound outside that fold).
Outputs: data/processed/lincs/pnet/<LINE>_<arm>/<set>.exp.tsv (gene x compound, first column "gene"),
         data/processed/lincs/pnet/regulators.txt; results/lincs/pnet_folds.tsv (pert_id, ModelID, cmap_name, group, fold)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
OUT = PROC / "pnet"
LINES = {"MCF7": "ACH-000019", "PC3": "ACH-000090", "A375": "ACH-000219"}
ARMS = {"raw": "drug_line_profiles_bing.f32", "specific": "drug_line_profiles_bing_specific.f32"}
N_FOLDS, SEED = 5, 1


def main():
    cov = pd.read_csv(ROOT / "results/lincs/cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    ids = cov.drop_duplicates("cell_iname").set_index("cell_iname").ModelID
    for k, v in LINES.items():
        assert ids[k] == v, (k, ids[k], v)
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    lm = np.where(genes.feature_space.eq("landmark").values)[0]
    assert len(lm) == 978
    gid = genes.gene_id.values[lm]
    tf = set(pd.read_csv(ROOT / "data/processed/step5/aracne_inputs_v2/BRCA.regulators_tf.txt", header=None)[0])
    regs = sorted(set(gid) & tf)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "regulators.txt").write_text("\n".join(regs) + "\n")
    print(f"landmark genes {len(gid)}; regulators (landmark Lambert TFs) {len(regs)}")

    ann = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t", dtype=str)[["pert_id", "cmap_name", "hub_moa", "lincs_moa"]]
    ix = idx.reset_index().rename(columns={"index": "row"})
    ix = ix[ix.ModelID.isin(LINES.values())].merge(ann, on="pert_id", how="left")
    ix["group"] = ix.hub_moa.fillna(ix.lincs_moa).fillna(ix.cmap_name).fillna(ix.pert_id).str.lower()
    sizes = ix.groupby("group").size()
    rng = np.random.default_rng(SEED)
    order = sorted(sizes.index, key=lambda g: (-sizes[g], rng.random()))
    load, fold_of = np.zeros(N_FOLDS), {}
    for g in order:
        f = int(np.argmin(load)); fold_of[g] = f + 1; load[f] += sizes[g]
    ix["fold"] = ix.group.map(fold_of)
    ix[["pert_id", "ModelID", "cmap_name", "group", "fold"]].to_csv(ROOT / "results/lincs/pnet_folds.tsv", sep="\t", index=False)
    print("profiles per fold:", ix.fold.value_counts().sort_index().to_dict(), "| mechanism groups:", len(sizes))
    anti = ix[ix.cmap_name.isin(["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene"])]
    print("anti-estrogen folds:", anti.groupby("cmap_name").fold.unique().to_dict())

    for arm, f in ARMS.items():
        X = np.fromfile(PROC / f, dtype=np.float32).reshape(len(idx), len(genes))
        for line, mid in LINES.items():
            d = OUT / f"{line}_{arm}"; d.mkdir(exist_ok=True)
            sub = ix[ix.ModelID == mid]
            for s in ["all"] + [f"f{k}" for k in range(1, N_FOLDS + 1)]:
                keep = sub if s == "all" else sub[sub.fold != int(s[1:])]
                M = pd.DataFrame(X[keep.row.values][:, lm].T, index=gid, columns=keep.pert_id.values)
                M.index.name = "gene"
                M.to_csv(d / f"{s}.exp.tsv", sep="\t", float_format="%.5f")
            print(f"{line} {arm}: {len(sub)} compounds")


if __name__ == "__main__":
    main()
