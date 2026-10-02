#!/usr/bin/env python
"""Step 5a: one ID mapping table for the whole project.

Step 4 signatures are keyed on unversioned GENCODE v23 Ensembl IDs (2015 symbols); CollecTRI /
DoRothEA use current HGNC symbols and aracne.networks uses Entrez. Every v23 gene in the extended
universe is mapped to the current HGNC record, trying in order:
  ensembl  - HGNC ensembl_gene_id == v23 Ensembl ID
  symbol   - v23 symbol is a current approved symbol
  prev     - v23 symbol is a previous symbol of exactly one approved gene
  alias    - v23 symbol is an alias of exactly one approved gene
Outputs:
  results/step5/id_map.tsv      gene_id, v23_symbol, biotype, hgnc_id, symbol, entrez_id, map_method
  results/step5/tf_universe.tsv Lambert 2018 'Is TF? = Yes' with mapped Ensembl / symbol / Entrez
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
OUT = ROOT / "results/step5"


def explode(hgnc, col):
    """symbol-like column (pipe-separated) -> {name: hgnc_id} keeping names that point to one gene."""
    x = hgnc[["hgnc_id", col]].dropna().assign(name=lambda d: d[col].str.split("|")).explode("name")
    counts = x.name.value_counts()
    return x[x.name.map(counts) == 1].set_index("name").hgnc_id


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(ROOT / "results/step3/genes_extended_v23.tsv", sep="\t")
    hgnc_f = sorted((RAW / "hgnc").glob("hgnc_complete_set_*.txt"))[-1]
    hgnc = pd.read_csv(hgnc_f, sep="\t", dtype=str, low_memory=False)
    hgnc = hgnc[hgnc.status == "Approved"]
    hgnc["entrez_id"] = hgnc.entrez_id.str.replace(r"\.0$", "", regex=True)

    by_ens = hgnc.dropna(subset=["ensembl_gene_id"]).drop_duplicates("ensembl_gene_id").set_index("ensembl_gene_id").hgnc_id
    by_sym = hgnc.set_index("symbol").hgnc_id
    by_prev, by_alias = explode(hgnc, "prev_symbol"), explode(hgnc, "alias_symbol")

    hid = pd.Series(pd.NA, index=genes.index, dtype="object")
    method = pd.Series("unmapped", index=genes.index)
    for name, lookup, key in [("ensembl", by_ens, genes.gene_id), ("symbol", by_sym, genes.symbol),
                              ("prev", by_prev, genes.symbol), ("alias", by_alias, genes.symbol)]:
        hit = hid.isna() & key.isin(lookup.index)
        hid[hit] = key[hit].map(lookup)
        method[hit] = name
    rec = hgnc.set_index("hgnc_id")
    out = pd.DataFrame({
        "gene_id": genes.gene_id, "v23_symbol": genes.symbol, "biotype": genes.biotype,
        "hgnc_id": hid, "symbol": hid.map(rec.symbol), "entrez_id": hid.map(rec.entrez_id),
        "map_method": method})
    # two v23 genes can land on one HGNC record (e.g. merged loci): `primary_for_hgnc` marks the one used
    # when mapping symbol/Entrez-keyed regulons back to Ensembl (Ensembl match first, then protein-coding)
    dup = out.hgnc_id.notna() & out.duplicated("hgnc_id", keep=False)
    out["hgnc_shared"] = dup
    rank = out.map_method.map({"ensembl": 0, "symbol": 1, "prev": 2, "alias": 3}).fillna(9) * 10 \
        + (out.biotype != "protein_coding")
    out["primary_for_hgnc"] = out.hgnc_id.notna() & ~out.assign(_r=rank).sort_values("_r").duplicated("hgnc_id").reindex(out.index)
    out.to_csv(OUT / "id_map.tsv", sep="\t", index=False)

    tf = pd.read_csv(RAW / "tf_lambert/DatabaseExtract_v_1.01.csv")
    tf = tf[tf["Is TF?"] == "Yes"].rename(columns={"Ensembl ID": "gene_id", "HGNC symbol": "lambert_symbol"})
    tf = tf[["gene_id", "lambert_symbol", "DBD"]].merge(out[["gene_id", "symbol", "entrez_id"]], on="gene_id", how="left")
    tf.to_csv(OUT / "tf_universe.tsv", sep="\t", index=False)

    print(f"HGNC release file: {hgnc_f.name}")
    print("map_method by biotype:\n", pd.crosstab(out.biotype, out.map_method).to_string())
    print(f"v23 symbol differs from current HGNC symbol: {(out.symbol.notna() & (out.symbol != out.v23_symbol)).sum()}")
    print(f"genes sharing an HGNC record: {dup.sum()}")
    print(f"Lambert TFs: {len(tf)}; in v23 universe: {tf.symbol.notna().sum()}")


if __name__ == "__main__":
    main()
