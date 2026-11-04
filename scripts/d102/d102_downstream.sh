#!/usr/bin/env bash
# D102c downstream refresh: gate summaries, D100 timing on gated early / late rankings, then every output that reads
# spec_mech_pct (same order as scripts/d97/d97_switch.sh step 5). <= 24 threads.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
step() { echo "$(date +%T) $1"; }
step "79i2 production (summaries)"; $PY scripts/79i2_target_presence.py > logs/d102_79i2.log 2>&1
mkdir -p results/lincs/d100/pre_d102; cp -n results/lincs/d100/*.tsv results/lincs/d100/pre_d102/ || true
for a in early late; do step "79i2 d100$a"; SFX=_v2_consensus_d100$a $PY scripts/79i2_target_presence.py > logs/d102_79i2_$a.log 2>&1; done
step "d100 timing"; $PY scripts/d100/d100_decide.py > logs/d102_d100_decide.log 2>&1
for s in 79p_fda_ranking 79o_literature 79y_standin; do step "$s"; $PY scripts/$s.py > logs/d102_$s.log 2>&1; done
step "79y2 counterfactual"; THREADS=8 SBX=$PWD/logs/sbx_borrowed_cf bash scripts/79y2_borrowed_counterfactual.sh > logs/d102_79y2.log 2>&1
for s in 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page 79z_standin_page 79q_fda_page; do step "$s"; $PY scripts/$s.py > logs/d102_$s.log 2>&1; done
step "D102 DOWNSTREAM DONE"
