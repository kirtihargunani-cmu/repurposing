#!/usr/bin/env bash
# Step 3 covariate sources: GTEx donor/sample annotations, TCGA plates + LAML blasts + ABSOLUTE purity,
# TARGET clinical spreadsheets. Appends sha256 to data/raw/checksums.sha256.
set -euo pipefail
cd "$(dirname "$0")/.."
R=data/raw
get() { mkdir -p "$1"; curl -fsSL --retry 5 -o "$1/${3:-$(basename "$2")}" "$2"; }

GT=https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files
get $R/gtex $GT/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
get $R/gtex $GT/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

get $R/tcga_pancan https://tcga.xenahubs.net/download/TCGA.LAML.sampleMap/LAML_clinicalMatrix
get $R/tcga_pancan https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5 TCGA_mastercalls.abs_tables_JSedit.fixed.txt

# TCGA RNA-seq aliquot barcodes (plate = 6th field) from the GDC API; snapshot, dated in the filename
mkdir -p $R/tcga_gdc
curl -fsS 'https://api.gdc.cancer.gov/files' -H 'Content-Type: application/json' -d '{
  "filters":{"op":"and","content":[
    {"op":"=","content":{"field":"cases.project.program.name","value":"TCGA"}},
    {"op":"=","content":{"field":"data_type","value":"Gene Expression Quantification"}},
    {"op":"=","content":{"field":"experimental_strategy","value":"RNA-Seq"}}]},
  "fields":"cases.samples.portions.analytes.aliquots.submitter_id",
  "format":"TSV","size":20000}' > $R/tcga_gdc/rnaseq_aliquots_$(date +%Y%m%d).tsv

# TARGET clinical (GDC open-access XLSX), by file UUID
while read -r uuid name; do
  get $R/target_clinical "https://api.gdc.cancer.gov/data/$uuid" "$name"
done <<'LIST'
3bf97830-865d-4e66-aaf1-d5bac81f119a TARGET_AML_ClinicalData_Discovery_20230720.xlsx
5f8b7137-c1e1-4191-aee5-a5ea55d32ca7 TARGET_AML_ClinicalData_Validation_20230720.xlsx
e58771ff-8501-4100-9c0a-c6dfa7dfb503 TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement_20230720.xlsx
19889559-1020-4fe3-b81e-b2eaa3f5b16c TARGET_ALL_ClinicalData_Phase_II_Discovery_20230727.xlsx
c6a63e40-f4e5-4002-9b87-3e6f3b3dfab4 TARGET_ALL_ClinicalData_Phase_II_Validation_20230727.xlsx
304aefc8-6a70-4896-9f54-2605d7e15630 TARGET_ALL_ClinicalData_Phase_I_20230727.xlsx
e20963a6-fc93-430d-b297-3da4f4da6326 TARGET_ALL_Blast_Data_20230727.xlsx
f5bfc8a5-b903-4418-b633-62234388a635 TARGET_NBL_ClinicalData_Discovery_20230523.xlsx
cdb0d02c-43d0-4d50-931d-ec2dbdf71a24 TARGET_NBL_ClinicalData_Validation_20230322.xlsx
d3002b03-71e8-4395-90ab-cd5c1020efe1 TARGET_WT_ClinicalData_Discovery_20230322.xlsx
5e3502fb-7211-4716-8355-8eaaa90c0c43 TARGET_WT_ClinicalData_Validation_20230322.xlsx
LIST

(cd $R && find gtex tcga_gdc target_clinical tcga_pancan/LAML_clinicalMatrix tcga_pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt \
   -type f | sort | xargs sha256sum | sed 's#  #  ./#') > $R/.new.sha256
cat $R/checksums.sha256 $R/.new.sha256 | sort -u -k2 -o $R/checksums.sha256 && rm $R/.new.sha256
echo "COVARIATE DOWNLOAD DONE"
