#!/usr/bin/env bash
# D105b downstream without the literature search (user 11-04: review the new pairs before searching the literature).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
step() { echo "$(date +%T) $1"; }
step "79y_standin"; $PY scripts/79y_standin.py > logs/d105b_79y_standin.log 2>&1
step "79y2 counterfactual"; THREADS=8 SBX=$PWD/logs/sbx_borrowed_cf bash scripts/79y2_borrowed_counterfactual.sh > logs/d105b_79y2.log 2>&1
for s in 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page 79z_standin_page 79q_fda_page; do step "$s"; $PY scripts/$s.py > logs/d105b_$s.log 2>&1; done
step "D105B REST DONE"
