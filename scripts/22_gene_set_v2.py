#!/usr/bin/env python
"""Step 3 (v2 rebuild, D46 / D54): the v2 gene set read by every v2 step (Steps 4, 5, 6, 10).

Start: protein-coding genes of GENCODE v23 (results/step3/genes_protein_coding_v23.tsv).
Drop:  gene mappability < 0.8 (Saha & Battle 2018; GRCh38 / GENCODE v26 resource in data/raw/crossmap/,
       matched on unversioned Ensembl IDs), and readthrough loci marked drop == TRUE in
       data/interim/readthrough_genes_v23.tsv.
Keep:  genes without a mappability value (flagged `no_mappability_value`).
Outputs (results/step3/):
  genes_protein_coding_v23_v2.tsv   kept genes: gene_id, symbol, mappability, no_mappability_value
  genes_v2_dropped.tsv              dropped genes with the reason(s)
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/step3"
MIN_MAP = 0.8


def main():
    pc = pd.read_csv(OUT / "genes_protein_coding_v23.tsv", sep="\t")
    pc["gene_id"] = pc.gene_id.str.split(".").str[0]
    gm = pd.read_csv(ROOT / "data/raw/crossmap/hg38_gene_mappability.txt.gz", sep="\t", header=None, names=["id", "map"])
    gm["gene_id"] = gm.id.str.split(".").str[0]
    mapv = gm.drop_duplicates("gene_id").set_index("gene_id")["map"]
    rt = pd.read_csv(ROOT / "data/interim/readthrough_genes_v23.tsv", sep="\t")
    rt = set(rt.loc[rt["drop"].astype(str).str.upper() == "TRUE", "gene_id"].str.split(".").str[0])

    pc["mappability"] = pc.gene_id.map(mapv)
    pc["no_mappability_value"] = pc.mappability.isna()
    low = pc.mappability < MIN_MAP
    thr = pc.gene_id.isin(rt)
    pc["reason"] = [";".join(r for r, f in (("mappability<0.8", a), ("readthrough", b)) if f) for a, b in zip(low, thr)]
    keep, drop = pc[~(low | thr)], pc[low | thr]
    keep[["gene_id", "symbol", "mappability", "no_mappability_value"]].to_csv(
        OUT / "genes_protein_coding_v23_v2.tsv", sep="\t", index=False)
    drop[["gene_id", "symbol", "mappability", "reason"]].to_csv(OUT / "genes_v2_dropped.tsv", sep="\t", index=False)
    print(f"protein-coding {len(pc)}; dropped {len(drop)} (mappability<0.8 {int(low.sum())}, readthrough {int(thr.sum())}, "
          f"both {int((low & thr).sum())}); kept {len(keep)} ({int(keep.no_mappability_value.sum())} without a mappability value)")


if __name__ == "__main__":
    main()
