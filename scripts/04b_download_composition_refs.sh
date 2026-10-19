#!/usr/bin/env bash
# Tissue-composition references (SPEC_tissue_composition.md, D61). All open access.
#   HPA single-cell types   Human Protein Atlas v24, rna_single_cell_type.tsv (nTPM per gene per single-cell type)
#   MCP-counter             gene lists (genes.txt) from github.com/ebecht/MCPcounter at a pinned commit; the R package
#                           is installed from the same commit (its default arguments fetch these lists at run time,
#                           so the pinned local copy is passed explicitly)
# Output: data/raw/composition/; sha256 appended to data/raw/checksums.sha256
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/raw/composition
MCP_SHA=b6eac73e91c246fcff0bb1a5c68a816cd588fc48   # ebecht/MCPcounter HEAD, pinned 2026-10-18
HPA_URL=https://v24.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip
mkdir -p $R
[ -s $R/rna_single_cell_type.tsv.zip ] || curl -fsSL --retry 5 -o $R/rna_single_cell_type.tsv.zip $HPA_URL
[ -s $R/rna_single_cell_type.tsv ] || unzip -o -q $R/rna_single_cell_type.tsv.zip -d $R
[ -s $R/mcpcounter_genes.txt ] || curl -fsSL --retry 5 -o $R/mcpcounter_genes.txt \
  https://raw.githubusercontent.com/ebecht/MCPcounter/$MCP_SHA/Signatures/genes.txt
Rscript=~/miniforge3/envs/repurposing/bin/Rscript
$Rscript -e "if (!requireNamespace('MCPcounter', quietly = TRUE)) remotes::install_github('ebecht/MCPcounter@$MCP_SHA', subdir = 'Source', upgrade = 'never')"
for f in $R/rna_single_cell_type.tsv.zip $R/mcpcounter_genes.txt; do
  grep -q " $f\$" data/raw/checksums.sha256 2>/dev/null || sha256sum $f >> data/raw/checksums.sha256
done
echo "composition references: $(ls $R | tr '\n' ' ')"
