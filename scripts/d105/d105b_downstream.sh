#!/usr/bin/env bash
# D105b downstream refresh (79i2 default = D105a FDA main targets): 79i2 + 79i3 on production and the D100 early / late tables, D100 timing, then every output
# that reads spec_mech_pct (as scripts/d102/d102_downstream.sh). <= 24 threads.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
step() { echo "$(date +%T) $1"; }
for t in _v2_consensus _v2_consensus_d100early _v2_consensus_d100late; do
  step "79i2 + 79i3 $t"; SFX=$t $PY scripts/79i2_target_presence.py > logs/d105b_79i2$t.log 2>&1; SFX=$t $PY scripts/79i3_receptor_lines.py > logs/d105b_79i3$t.log 2>&1
done
step "d100 timing"; $PY scripts/d100/d100_decide.py > logs/d105b_d100_decide.log 2>&1
for s in 79p_fda_ranking 79o_literature 79y_standin; do step "$s"; $PY scripts/$s.py > logs/d105b_$s.log 2>&1; done
step "79y2 counterfactual"; THREADS=8 SBX=$PWD/logs/sbx_borrowed_cf bash scripts/79y2_borrowed_counterfactual.sh > logs/d105b_79y2.log 2>&1
for s in 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page 79z_standin_page 79q_fda_page; do step "$s"; $PY scripts/$s.py > logs/d105b_$s.log 2>&1; done
step "D105B DOWNSTREAM DONE"
