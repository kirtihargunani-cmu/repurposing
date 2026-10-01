#!/usr/bin/env bash
# Pull every source in config/datasets.yaml into data/raw/<source>/ and record sha256 + versions.
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/raw
get() {  # get <dir> <url> [outname]
  mkdir -p "$1"; local out="$1/${3:-$(basename "$2")}"
  curl -fL --retry 5 -C - -o "$out" "$2"
}

TOIL=https://toil.xenahubs.net/download
get $R/solid_toil $TOIL/TcgaTargetGtex_rsem_gene_tpm.gz &
get $R/solid_toil $TOIL/TcgaTargetGtex_gene_expected_count.gz &
get $R/solid_toil $TOIL/TcgaTargetGTEX_phenotype.txt.gz
get $R/solid_toil $TOIL/probeMap/gencode.v23.annotation.gene.probemap

# BeatAML: pin to the current commit
SHA=$(curl -fsS -H 'Accept: application/vnd.github.sha' https://api.github.com/repos/biodev/beataml2.0_data/commits/main)
echo "$SHA" > $R/aml_beataml2.COMMIT
for f in beataml_waves1to4_norm_exp_dbgap.txt beataml_probit_curve_fits_v4_dbgap.txt \
         beataml_wv1to4_clinical.xlsx beataml_drug_families.xlsx; do
  get $R/aml_beataml2 https://github.com/biodev/beataml2.0_data/raw/$SHA/$f
done

get $R/mondo https://github.com/monarch-initiative/mondo/releases/download/v2026-09-01/mondo.obo
get $R/uberon http://purl.obolibrary.org/obo/uberon/basic.obo uberon-basic.obo
get $R/gencode https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz

# TCGA Pan-Cancer Atlas clinical (TCGA-CDR, Liu 2018) + molecular subtypes, for subtype splits
get $R/tcga_pancan https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp
get $R/tcga_pancan https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGASubtype.20170308.tsv.gz

mkdir -p $R/oncotree
curl -fs https://oncotree.mskcc.org/api/versions > $R/oncotree/versions.json
get $R/oncotree "https://oncotree.mskcc.org/api/tumorTypes?version=oncotree_latest_stable" tumorTypes.json

wait
grep -m1 '^data-version' $R/uberon/uberon-basic.obo > $R/uberon/VERSION || true
(cd $R && find . -type f ! -name '*.sha256' ! -name 'download.log' | sort | xargs sha256sum) > $R/checksums.sha256
echo "DOWNLOAD DONE"
