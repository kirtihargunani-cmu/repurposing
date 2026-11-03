#!/usr/bin/env bash
# D92 (rule fixed 2026-11-03, docs/DECISIONS.md): breast cell-type networks through the pipeline, all outputs tagged _d92
# (production untouched). Run after scripts/d92/d92_inputs.py and VERSION=v2 REGSET=tf scripts/43_run_aracne.sh BRCA_BASX.
# Threads: <= 24 in total.
set -euo pipefail
cd "$(dirname "$0")/../.."
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
EXTRA=data/processed/step5/regulons_aracne_project_v2_tf_d92.rds
MAP=results/step5/d92_group_net.tsv
OVR=results/step6/aracne_project_v2_tf_d92
G=BRCA-LumA,BRCA-LumB,BRCA-Basal,BRCA-Her2
NETS=$(cut -f2 $MAP | tail -n +2 | sort -u | paste -sd,)
echo "networks: $NETS"
# 1. regulons for the new network(s) (BRCA_LUM already in the production regulon set)
NEW=$(for n in ${NETS//,/ }; do [ -d data/processed/step5/aracne_v2_tf/$n ] && [ "$n" != BRCA_LUM ] && echo $n; done | paste -sd' ')
VERSION=v2 REGSET=tf OUT_TAG=_d92 $RS scripts/44_aracne_regulons.R $NEW
# 2. disease side: breast groups on their cell-type network, and on pooled BRCA with the same seeds (comparison)
REGULON=aracne_project VERSION=v2 REGSET=tf RES_TAG=_d92 EXTRA_RDS=$EXTRA NET_MAP=$MAP ONLY_GROUPS=$G GROUP_ONLY=1 SEED_GROUP=1 N_CORES=4 \
  $RS scripts/50_viper.R
REGULON=aracne_project VERSION=v2 REGSET=tf RES_TAG=_d92pooled ONLY_GROUPS=$G GROUP_ONLY=1 SEED_GROUP=1 N_CORES=4 $RS scripts/50_viper.R
# 3. proliferation-linked regulators (P1) of the new network(s)
VERSION=v2 REGSET=tf EXTRA_RDS=$EXTRA OUT_TAG=_d92extra $RS scripts/79h_prolif_regulators.R
# 4. drug activity on the cell-type networks; other networks linked from the _d85 run
VERSION=v2 REGSET=tf PROFILE=specific RUN_TAG=_d92 EXTRA_RDS=$EXTRA NET_MAP=$MAP MR_OVERRIDE=$OVR ONLY_NETS=$NETS N_CORES=2 \
  $RS scripts/79d_drug_viper.R
D=data/processed/lincs/drug_viper_specific_v2_tf_d92
for f in data/processed/lincs/drug_viper_specific_v2_tf_d85/*; do b=$(basename "$f"); [ -e "$D/$b" ] || ln -s "$(realpath "$f")" "$D/$b"; done
# 5. rankings
export OPENBLAS_NUM_THREADS=4
VERSION=v2 REGSET=tf RUN_TAG=_d92 NET_MAP=$MAP MR_OVERRIDE=$OVR $PY scripts/79e_drug_reversal.py
VERSION=v2 REGSET=tf RUN_TAG=_d92 NET_MAP=$MAP MR_OVERRIDE=$OVR P1_EXTRA=results/lincs/prolif_regulators_v2_tf_d92extra.tsv \
  $PY scripts/79i_split_reversal.py
echo D92 PIPELINE DONE
