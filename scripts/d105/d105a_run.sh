#!/usr/bin/env bash
# D105 evaluation run: copies of the production split tables (suffix _v2_consensus_d105), then 79i2 + 79i3 with
# TARGETS=d105a on the copies, then d105_eval.py. Production tables are not written.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4
for k in drug stratum; do cp results/lincs/${k}_reversal_split_v2_consensus.tsv results/lincs/${k}_reversal_split_v2_consensus_d105a.tsv; done
SFX=_v2_consensus_d105a TARGETS=d105a $PY scripts/79i2_target_presence.py > logs/d105a_79i2.log 2>&1
SFX=_v2_consensus_d105a TARGETS=d105a $PY scripts/79i3_receptor_lines.py > logs/d105a_79i3.log 2>&1
TAG=d105a $PY scripts/d105/d105_eval.py
