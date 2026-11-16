#!/usr/bin/env bash
# CTRPv2 cell-line drug sensitivity, as a functional label that this project has NOT previously read.
#
# Why a third screen: D128-A0's prior-use audit found PRISM 20Q2 was already read by 79j / 79k / 79m / 09b /
# d94_decide.py and cited in D94, D97, D102c, D119, D120, D122 and A43, and GDSC 8.5 only appears alongside it.
# Neither can serve as an independent outcome label for D128 stage A5. CTRP v2 has never been read here (checked
# by grep over scripts/, docs/ and config/: the only hits are the D128 spec naming it as a candidate).
#
# Source: the original Broad / CTD^2 "CTRPv2.0_2015_ctd2_ExpandedDataset" archive. The NCI CTD^2 path
# (ctd2-data.nci.nih.gov/Public/Broad/...) now redirects to studycatalog.cancer.gov and is retired; the
# CloudFront object below is the mirror used by the BHK lab's CTRPv2 Snakemake pipeline
# (github.com/BHKLAB-DataProcessing/CTRPv2-Pharmacoset_Snakemake, config/pipeline.yaml) and serves the same zip.
# md5 verified on download 2026-11-15.
#
# Cite: Rees et al. Nat Chem Biol 2016; Seashore-Ludlow et al. Cancer Discov 2015; Basu et al. Cell 2013.
# Use follows the CTD^2 data release policy (ocg.cancer.gov/programs/ctd2/using-ctd2-data).
#
# Join keys (verified 2026-11-15):
#   compound  v20.meta.per_compound.txt  broad_cpd_id -> first 13 chars = LINCS pert_id
#   line      v20.meta.per_cell_line.txt ccl_name -> DepMap Model.csv StrippedCellLineName -> ModelID
#   value     v20.data.curves_post_qc.txt area_under_curve (also apparent_ec50_umol), joined via experiment_id
#             -> v20.meta.per_experiment.txt -> master_ccl_id
# Coverage against this project: 303 of our scored pert_ids, 70 of our 73 LINCS lines; 229,880 AUC values for
# those compounds across 887 lines (median 817 lines per compound).
set -euo pipefail
cd "$(dirname "$0")/.."
C=data/raw/ctrp_v2
mkdir -p $C

ZIP=CTRPv2.0_2015_ctd2_ExpandedDataset.zip
URL=https://d1hp0ufb0xxisr.cloudfront.net/CTD2/Broad/$ZIP
MD5=e764f6074a624d515d1b9d75e8b306e7

[ -s "$C/$ZIP" ] || curl -fsSL --retry 5 -o "$C/$ZIP" "$URL"
echo "$MD5  $C/$ZIP" | md5sum -c --quiet -

# Extract only what the analysis needs. The four per_cpd_* tables are 1.8 GB of per-well raw readings and are
# not unpacked; re-run unzip on the archive if they are ever wanted.
unzip -o -q -j "$C/$ZIP" \
  v20.data.curves_post_qc.txt \
  v20.meta.per_compound.txt \
  v20.meta.per_cell_line.txt \
  v20.meta.per_experiment.txt \
  v20._COLUMNS.txt \
  v20._README.txt \
  -d "$C"

echo "CTRP v2 ready in $C:"
ls -1 "$C"
