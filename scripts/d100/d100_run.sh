#!/usr/bin/env bash
# D100: early and late rankings through the production pipeline (D97 networks by default), tags _d100early / _d100late.
set -euo pipefail
cd "$(dirname "$0")/../.."
RS=~/miniforge3/envs/repurposing/bin/Rscript; PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=4 $PY scripts/d100/d100_profiles.py > logs/d100_profiles.log 2>&1
D=data/processed/lincs/d100
for a in early late; do
  ( VERSION=v2 REGSET=tf PROFILE=specific PROFILE_FILE=$D/${a}_specific.f32 RUN_TAG=_d100$a N_CORES=10 $RS scripts/79d_drug_viper.R > logs/d100_79d_tf_$a.log 2>&1
    VERSION=v2 REGULON=collectri PROFILE=specific PROFILE_FILE=$D/${a}_specific.f32 RUN_TAG=_d100$a N_CORES=1 $RS scripts/79d_drug_viper.R > logs/d100_79d_col_$a.log 2>&1 ) &
done; wait
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
for a in early late; do
  ( VERSION=v2 REGULON=consensus RUN_TAG=_d100$a $PY scripts/79e_drug_reversal.py > logs/d100_79e_$a.log 2>&1
    VERSION=v2 REGULON=consensus RUN_TAG=_d100$a STRATA=1 $PY scripts/79e_drug_reversal.py > logs/d100_79e_strata_$a.log 2>&1
    VERSION=v2 REGULON=consensus RUN_TAG=_d100$a STRATA=1 PROFILE_PREFIX=$D/$a $PY scripts/79i_split_reversal.py > logs/d100_79i_$a.log 2>&1 ) &
done; wait
echo "D100 RUN DONE"
