#!/usr/bin/env python
"""Drug x cell-line consensus profiles on landmark + best-inferred genes, for regulator scoring (runs in `celligner`).

Input: data/processed/lincs/level5_bing_hiq.npy (+ _sigs, genes; scripts/76b). Active signatures as everywhere on the
drug side: is_hiq == 1 and TAS >= 0.4 (user, 2026-10-09). One profile per drug x line = mean of its active signatures
(MODZ z-scores) over doses / times, as in scripts/77.
Two versions:
  raw       the consensus profiles
  specific  per line, the line's generic drug response g_l (mean raw profile over its drugs, unit vector) is projected
            out of each profile, as in scripts/77. Needed for regulator scoring (user 2026-10-12): every active drug
            lowers proliferation genes, and in tumor-derived ARACNe networks lineage regulators carry proliferation
            genes as negative targets (e.g. FOXA1 / GATA3 / XBP1 in breast), so on raw profiles VIPER reads any
            cytostatic drug as raising them; anti-estrogens in MCF7 then score as anti-reversing in luminal breast
            cancer although they lower GREB1 / MYB / XBP1 at the gene level.
Genes are mapped to the GENCODE v23 Ensembl IDs of the regulons (results/step5/id_map.tsv, primary_for_hgnc): by
Entrez ID, else by Ensembl ID; genes that map to no or to an already used v23 gene are dropped.
Outputs: data/processed/lincs/drug_line_profiles_bing.f32 (raw) and drug_line_profiles_bing_specific.f32
         (float32, profile-major: profiles x genes, C order),
         drug_line_profiles_index.tsv (pert_id, ModelID, n_sig), drug_line_profiles_genes.tsv (v23 gene_id, symbol)
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
TAS_MIN = 0.4


def main():
    X = np.load(PROC / "level5_bing_hiq.npy", mmap_mode="r")
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    genes = pd.read_csv(PROC / "level5_bing_genes.tsv", sep="\t", dtype={"gene_id": str})

    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t", dtype={"entrez_id": str})
    idm = idm[idm.primary_for_hgnc == True]
    by_entrez = idm.dropna(subset=["entrez_id"]).drop_duplicates("entrez_id").set_index("entrez_id").gene_id
    v23 = genes.gene_id.map(by_entrez)
    v23 = v23.fillna(genes.ensembl_id.where(genes.ensembl_id.isin(idm.gene_id)))
    keep = v23.notna() & ~v23.duplicated()
    print(f"genes: {len(genes)}; mapped to v23 {int(v23.notna().sum())}; kept (unique) {int(keep.sum())}")
    gpos = np.where(keep.values)[0]

    act = sig[sig.tas >= TAS_MIN]
    grp = act.groupby(["pert_id", "ModelID"]).indices
    idx = pd.DataFrame(list(grp.keys()), columns=["pert_id", "ModelID"])
    idx["n_sig"] = [len(v) for v in grp.values()]
    rows = act.index.values
    P = np.empty((len(idx), len(gpos)), dtype=np.float32)
    for i, v in enumerate(grp.values()):
        P[i] = np.asarray(X[np.sort(rows[v])][:, gpos], dtype=np.float64).mean(axis=0)
    assert not np.isnan(P).any()
    P.tofile(PROC / "drug_line_profiles_bing.f32")
    S = P.astype(np.float64)
    for ix in idx.groupby("ModelID").indices.values():
        g = S[ix].mean(axis=0)
        g /= np.linalg.norm(g)
        S[ix] = S[ix] - np.outer(S[ix] @ g, g)
    S.astype(np.float32).tofile(PROC / "drug_line_profiles_bing_specific.f32")
    idx.to_csv(PROC / "drug_line_profiles_index.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_id": v23[keep].values, "symbol": genes.gene_symbol[keep].values,
                  "feature_space": genes.feature_space[keep].values}).to_csv(
        PROC / "drug_line_profiles_genes.tsv", sep="\t", index=False)
    print(f"profiles: {len(idx)} ({idx.pert_id.nunique()} drugs x {idx.ModelID.nunique()} lines) x {len(gpos)} genes")


if __name__ == "__main__":
    main()
