#!/usr/bin/env bash
# D132 step 4: re-rank with the corrected activities using the PRODUCTION 79e_drug_reversal.py, unmodified.
# RUN_TAG=_d132corr makes it read the corrected drug-viper sandbox and write *_v2_consensus_d132corr.tsv.
# Production tables are not touched.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=${PY:-$HOME/miniforge3/envs/celligner/bin/python}
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
echo "=== groups ==="
env VERSION=v2 REGULON=consensus RUN_TAG=_d132corr "$PY" scripts/79e_drug_reversal.py
echo "=== subgroups ==="
env VERSION=v2 REGULON=consensus STRATA=1 RUN_TAG=_d132corr "$PY" scripts/79e_drug_reversal.py
