#!/usr/bin/env bash
# Steps 5-7 sources: Lambert TF list, HGNC complete set (ID mapping), TCGA MC3 mutations (PanCanAtlas),
# TCGA GISTIC2 thresholded gene-level copy number (Xena TCGA hub). CollecTRI / DoRothEA are pulled in R (decoupleR) by
# scripts/40_regulons.R. Appends sha256 to data/raw/checksums.sha256.
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/raw
get() { mkdir -p "$1"; curl -fsSL --retry 5 -C - -o "$1/${3:-$(basename "$2")}" "$2"; }
get $R/tf_lambert https://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv
get $R/hgnc https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt hgnc_complete_set_$(date +%Y%m%d).txt
get $R/tcga_pancan https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc mc3.v0.2.8.PUBLIC.maf.gz
get $R/tcga_pancan https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz
(cd $R && find tf_lambert hgnc tcga_pancan/mc3.v0.2.8.PUBLIC.maf.gz tcga_pancan/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz \
   -type f | sort | xargs sha256sum | sed 's#  #  ./#') > $R/.new.sha256
cat $R/checksums.sha256 $R/.new.sha256 | sort -u -k2 -o $R/checksums.sha256 && rm $R/.new.sha256
echo "REGULON SOURCES DONE"
