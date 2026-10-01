#!/usr/bin/env bash
# Create the conda env, then install ESTIMATE from R-Forge (not on CRAN/Bioc).
set -euo pipefail
cd "$(dirname "$0")/.."
MAMBA=~/miniforge3/bin/mamba
$MAMBA env create -y -f environment.yml
$MAMBA run -n repurposing Rscript -e '
  install.packages("estimate", repos = "https://R-Forge.R-project.org", type = "source")
  for (p in c("limma","edgeR","DESeq2","recount3","biomaRt","org.Hs.eg.db","estimate"))
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  cat("R packages OK\n")'
$MAMBA run -n repurposing python -c 'import pandas, openpyxl, yaml, pronto; print("Python packages OK")'
$MAMBA env export -n repurposing > environment.lock.yml
echo "ENV DONE"
