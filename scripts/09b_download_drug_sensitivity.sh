#!/usr/bin/env bash
# Drug sensitivity for validating the reversal scores (user 2026-10-12): independent of LINCS (viability, not
# expression). PRISM Repurposing 20Q2 (Figshare 20564034; Broad IDs = LINCS pert_id, DepMap IDs = ModelID):
# secondary screen = 8-dose curves (AUC, IC50; ~1,450 drugs x ~480 lines), primary screen = one dose (2.5 uM,
# replicate-collapsed log fold change; ~4,500 drugs x ~580 lines). md5 = Figshare supplied_md5.
# GDSC release 8.5 (Sanger, cog.sanger.ac.uk; fitted dose response 27Oct23, GDSC1 and GDSC2; lines by COSMIC /
# Sanger model ID -> DepMap Model.csv). md5 = S3 ETag.
set -euo pipefail
cd "$(dirname "$0")/.."
P=data/raw/prism_20q2
G=data/raw/gdsc_8.5
mkdir -p $P $G
while read -r id name md5; do
  [ -s "$P/$name" ] || curl -fsSL --retry 5 -o "$P/$name" "https://ndownloader.figshare.com/files/$id"
  echo "$md5  $P/$name" | md5sum -c --quiet -
done <<'LIST'
36794613 prism-repurposing-20q2-secondary-screen-readme 47419770d1edfcc36e0b11465f7d1fd0
36794589 prism-repurposing-20q2-secondary-screen-cell-line-info.csv 4c3a5fbee769e8384064a8fbe4b4ed7d
36794595 prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv eb0c7ee3ccb9e480148c71fcdc97312e
36794619 prism-repurposing-20q2-secondary-screen-replicate-collapsed-treatment-info.csv 162da467eb97a67abdbdbee7af101091
36794577 prism-repurposing-20q2-primary-screen-readme.txt 88f8769237b8c9f739afab3e938a6d60
36794565 prism-repurposing-20q2-primary-screen-cell-line-info.csv 580b9f0b5118d44f473ffa4efac7e0c2
36794580 prism-repurposing-20q2-primary-screen-replicate-collapsed-logfold-change.csv 8f3eb7c32f4e4e35d6c623e1605f1c1d
36794583 prism-repurposing-20q2-primary-screen-replicate-collapsed-treatment-info.csv 1e69e34208fcc5c5c4118f1598d2ad89
LIST
B=https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5
while read -r name md5; do
  [ -s "$G/$name" ] || curl -fsSL --retry 5 -o "$G/$name" "$B/$name"
  echo "$md5  $G/$name" | md5sum -c --quiet -
done <<'LIST'
GDSC1_fitted_dose_response_27Oct23.xlsx d3db01db37e4bd461641a7bc9231012b
GDSC2_fitted_dose_response_27Oct23.xlsx a57b46470f325ae5800dca24c72c19f5
screened_compounds_rel_8.5.csv 1f8a7c36e02f3298704063f104b5ff81
LIST
(cd data/raw && find prism_20q2 gdsc_8.5 -type f | sort | xargs sha256sum | sed 's#  #  ./#') > data/raw/.new.sha256
cat data/raw/checksums.sha256 data/raw/.new.sha256 | sort -u -k2 -o data/raw/checksums.sha256 && rm data/raw/.new.sha256
echo "DRUG SENSITIVITY DONE"
