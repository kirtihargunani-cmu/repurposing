#!/usr/bin/env bash
# D94: gated arm (cohort gate: BRCA stays pooled, 7 cohorts on cell-type networks) upstream, then downstream for all arms.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
E=(NET_MAP=results/step5/d94_gated_group_net.tsv MR_OVERRIDE=results/step6/aracne_project_v2_tf_d93
   MR_OVERRIDE_STRATA=results/step6/aracne_project_v2_tf_strata_d94 P1_EXTRA=results/lincs/prolif_regulators_v2_tf_d93extra.tsv)
env "${E[@]}" VERSION=v2 REGULON=consensus RUN_TAG=_d94gated $PY scripts/79e_drug_reversal.py > logs/d94_79e_d94gated.log 2>&1
env "${E[@]}" VERSION=v2 REGULON=consensus RUN_TAG=_d94gated STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d94_79e_strata_d94gated.log 2>&1
env "${E[@]}" VERSION=v2 REGULON=consensus RUN_TAG=_d94gated STRATA=1 $PY scripts/79i_split_reversal.py > logs/d94_79i_d94gated.log 2>&1
echo "gated upstream done"
( bash scripts/d94/d94_downstream.sh cur _d85 > logs/d94_down_cur.log 2>&1; bash scripts/d94/d94_downstream.sh lum _d94lum > logs/d94_down_lum.log 2>&1 ) &
( bash scripts/d94/d94_downstream.sh split _d94split > logs/d94_down_split.log 2>&1; bash scripts/d94/d94_downstream.sh gated _d94gated > logs/d94_down_gated.log 2>&1 ) &
wait
echo "D94 REST DONE"
