#!/usr/bin/env bash
# D94 upstream (rule fixed 2026-11-03): subgroup re-scoring and consensus rankings (groups + subgroups) for CUR / LUM / SPLIT.
# Tags: CUR _d85, LUM _d94lum, SPLIT _d94split. Production outputs untouched. <= 24 threads.
set -euo pipefail
cd "$(dirname "$0")/../.."
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
P=data/processed/lincs
# 1. subgroups whose parent moved network (same map in both arms: no breast subgroups)
SUBS=LGG-noncodel-LGr2,LGG-noncodel-LGr3,PHEO-kinase-signaling,PHEO-pseudohypoxia-Wnt,TGCT-SEM-C1,TGCT-SEM-C2
OPENBLAS_NUM_THREADS=1 REGULON=aracne_project VERSION=v2 REGSET=tf STRATA=results/step10_v2/stratum_patients.tsv RES_TAG=_d94 \
  EXTRA_RDS=data/processed/step5/regulons_aracne_project_v2_tf_d93.rds NET_MAP=results/step5/d93_group_net.tsv \
  ONLY_GROUPS=$SUBS GROUP_ONLY=1 SEED_GROUP=1 N_CORES=6 $RS scripts/50_viper.R
# 2. activity folders: ARACNe from _d93 (+ BRCA_LUM from _d92b for LUM), CollecTRI unchanged (_d85)
for arm in lum split; do
  for set in v2_tf collectri; do
    D=$P/drug_viper_specific_${set}_d94$arm; mkdir -p $D
    src=$([ $set = v2_tf ] && echo $P/drug_viper_specific_v2_tf_d93 || echo $P/drug_viper_specific_collectri_d85)
    for f in $src/*; do b=$(basename "$f"); [ -e "$D/$b" ] || ln -s "$(realpath "$f")" "$D/$b"; done
  done
done
for f in BRCA_LUM.f32 BRCA_LUM.regulators.txt; do ln -sf "$(realpath $P/drug_viper_specific_v2_tf_d92b/$f)" $P/drug_viper_specific_v2_tf_d94lum/$f; done
[ -e $P/drug_viper_specific_collectri_d85/ALL.f32 ] || { echo "missing collectri _d85"; exit 1; }
# 3. consensus rankings, groups then subgroups (79e), split scores (79i, STRATA=1 does both)
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
run_arm() {  # tag map override
  local tag=$1 map=$2 ovr=$3 env=()
  [ -n "$map" ] && env=(NET_MAP=$map MR_OVERRIDE=$ovr MR_OVERRIDE_STRATA=results/step6/aracne_project_v2_tf_strata_d94 P1_EXTRA=results/lincs/prolif_regulators_v2_tf_d93extra.tsv)
  env "${env[@]}" VERSION=v2 REGULON=consensus RUN_TAG=$tag $PY scripts/79e_drug_reversal.py > logs/d94_79e$tag.log 2>&1
  env "${env[@]}" VERSION=v2 REGULON=consensus RUN_TAG=$tag STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d94_79e_strata$tag.log 2>&1
  env "${env[@]}" VERSION=v2 REGULON=consensus RUN_TAG=$tag STRATA=1 $PY scripts/79i_split_reversal.py > logs/d94_79i$tag.log 2>&1
  echo "arm $tag done"
}
run_arm _d85 "" "" &
run_arm _d94lum results/step5/d94_lum_group_net.tsv results/step6/aracne_project_v2_tf_d94lum &
run_arm _d94split results/step5/d94_split_group_net.tsv results/step6/aracne_project_v2_tf_d94split &
wait
echo "D94 UPSTREAM DONE"
