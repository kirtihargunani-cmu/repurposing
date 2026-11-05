#!/usr/bin/env bash
# Resume d105b_rest.sh after the sandbox 79i2 step failed (missing fda_matches.tsv link; fixed in the 79y2 sh).
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python; SBX=$PWD/logs/sbx_borrowed_cf
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
step() { echo "$(date +%T) $1"; }
step "sandbox 79i2 + 79i3"; (cd "$SBX" && "$PY" scripts/79i2_target_presence.py > i2.log 2>&1 && "$PY" scripts/79i3_receptor_lines.py > i3.log 2>&1)
step "79y2"; SBX="$SBX" $PY scripts/79y2_standin_vs_borrowed.py > logs/d105b_79y2_cmp.log 2>&1
for s in 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page 79z_standin_page 79q_fda_page; do step "$s"; $PY scripts/$s.py > logs/d105b_$s.log 2>&1; done
step "PRESENT_ONLY 79x"; PRESENT_ONLY=1 $PY scripts/79x_rankings_page.py > logs/d105b_79x_present.log 2>&1
step "D105B RESUME DONE"
