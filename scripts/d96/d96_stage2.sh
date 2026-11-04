#!/usr/bin/env bash
# D95 / D96 stage 2 (after logs/d95_d96_stage1.sh): D95 drug activity + combination; D96 regulons, Step 6, activity;
# consensus rankings (groups + subgroups) for both; downstream worktrees (meta, t85).
set -euo pipefail
cd "$(dirname "$0")/../.."
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
X92=data/processed/step5/regulons_aracne_project_v2_tf_d92b.rds
# D95: activity on BRCA_LUM / BRCA_BAS with the four breast groups' filters, then combine
for t in lum bas; do N=$([ $t = lum ] && echo BRCA_LUM || echo BRCA_BAS)
  VERSION=v2 REGSET=tf PROFILE=specific RUN_TAG=_d95$t EXTRA_RDS=$X92 NET_MAP=results/step5/d95_all_$t.tsv MR_OVERRIDE=results/step6/aracne_project_v2_tf_d95$t \
    ONLY_NETS=$N N_CORES=1 $RS scripts/79d_drug_viper.R > logs/d95_79d_$t.log 2>&1 &
done; wait
OPENBLAS_NUM_THREADS=4 $PY scripts/d95/d95_combine.py > logs/d95_combine.log 2>&1
# D96: regulons of the two new networks, combined regulon file, Step 6 for PHEO / EEC and PHEO subgroups, P1, activity
VERSION=v2 REGSET=tf OUT_TAG=_d96 $RS scripts/44_aracne_regulons.R PCPG_PHEO UCEC_EEC > logs/d96_44.log 2>&1
$RS -e 'a <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf_d93.rds"); b <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf_d96.rds"); saveRDS(c(a, b[setdiff(names(b), names(a))]), "data/processed/step5/regulons_aracne_project_v2_tf_d96all.rds")'
X96=data/processed/step5/regulons_aracne_project_v2_tf_d96all.rds
VERSION=v2 REGSET=tf EXTRA_RDS=$X96 OUT_TAG=_d96extra0 $RS scripts/79h_prolif_regulators.R > logs/d96_79h.log 2>&1
$PY -c "import pandas as pd; x=pd.read_csv('results/lincs/prolif_regulators_v2_tf_d96extra0.tsv',sep='\t'); x[x.net.isin(['PCPG_PHEO','UCEC_EEC'])].to_csv('results/lincs/prolif_regulators_v2_tf_d96extra.tsv',sep='\t',index=False)"
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
