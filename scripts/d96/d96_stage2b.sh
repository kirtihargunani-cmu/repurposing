#!/usr/bin/env bash
# Resume of d96_stage2.sh from the D96 Step 6 runs (stage 2 stopped there: the map was written after its first use; map now written first).
# consensus rankings (groups + subgroups) for both; downstream worktrees (meta, t85).
set -euo pipefail
cd "$(dirname "$0")/../.."
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
X92=data/processed/step5/regulons_aracne_project_v2_tf_d92b.rds
X96=data/processed/step5/regulons_aracne_project_v2_tf_d96all.rds
REGULON=aracne_project VERSION=v2 REGSET=tf RES_TAG=_d96new EXTRA_RDS=$X96 NET_MAP=results/step5/d96_group_net.tsv \
  ONLY_GROUPS=PCPG-PHEO,UCEC-EEC GROUP_ONLY=1 SEED_GROUP=1 N_CORES=2 $RS scripts/50_viper.R > logs/d96_viper.log 2>&1
REGULON=aracne_project VERSION=v2 REGSET=tf STRATA=results/step10_v2/stratum_patients.tsv RES_TAG=_d96 EXTRA_RDS=$X96 NET_MAP=results/step5/d96_group_net.tsv \
  ONLY_GROUPS=PHEO-kinase-signaling,PHEO-pseudohypoxia-Wnt GROUP_ONLY=1 SEED_GROUP=1 N_CORES=2 $RS scripts/50_viper.R > logs/d96_viper_strata.log 2>&1
$PY scripts/d96/d96_prep.py > logs/d96_prep.log 2>&1
VERSION=v2 REGSET=tf PROFILE=specific RUN_TAG=_d96 EXTRA_RDS=$X96 NET_MAP=results/step5/d96_group_net.tsv MR_OVERRIDE=results/step6/aracne_project_v2_tf_d96 \
  ONLY_NETS=PCPG_PHEO,UCEC_EEC N_CORES=2 $RS scripts/79d_drug_viper.R > logs/d96_79d.log 2>&1
P=data/processed/lincs
for f in $P/drug_viper_specific_v2_tf_d93/*; do b=$(basename "$f"); [ -e $P/drug_viper_specific_v2_tf_d96/$b ] || ln -s "$(realpath "$f")" $P/drug_viper_specific_v2_tf_d96/$b; done
mkdir -p $P/drug_viper_specific_collectri_d96; for f in $P/drug_viper_specific_collectri_d85/*; do b=$(basename "$f"); [ -e $P/drug_viper_specific_collectri_d96/$b ] || ln -s "$(realpath "$f")" $P/drug_viper_specific_collectri_d96/$b; done
# consensus rankings, groups then subgroups, split
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
arm() { local tag=$1; shift
  env "$@" VERSION=v2 REGULON=consensus RUN_TAG=$tag $PY scripts/79e_drug_reversal.py > logs/d9x_79e$tag.log 2>&1
  env "$@" VERSION=v2 REGULON=consensus RUN_TAG=$tag STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d9x_79e_strata$tag.log 2>&1
  env "$@" VERSION=v2 REGULON=consensus RUN_TAG=$tag STRATA=1 $PY scripts/79i_split_reversal.py > logs/d9x_79i$tag.log 2>&1; }
arm _d95 NET_MAP=results/step5/d95_group_net.tsv MR_OVERRIDE=results/step6/aracne_project_v2_tf_d95 P1_EXTRA=results/lincs/prolif_regulators_v2_tf_d95extra.tsv &
arm _d96 NET_MAP=results/step5/d96_group_net.tsv MR_OVERRIDE=results/step6/aracne_project_v2_tf_d96 MR_OVERRIDE_STRATA=results/step6/aracne_project_v2_tf_strata_d96x \
    P1_EXTRA=results/lincs/prolif_regulators_v2_tf_d96all.tsv &
wait
echo "rankings done"
bash scripts/d94/d94_downstream.sh meta _d95 > logs/d95_down.log 2>&1 &
bash scripts/d94/d94_downstream.sh t85 _d96 > logs/d96_down.log 2>&1 &
wait
echo "STAGE2 DONE"
