#!/usr/bin/env bash
# D132 step 2: score the null profiles with the PRODUCTION 79d_drug_viper.R, unmodified.
# PROFILE_FILE / INDEX_FILE (D90, D125-A0) point it at the null matrix; RUN_TAG (D85) keeps the output in a sandbox.
set -euo pipefail
cd "$(dirname "$0")/../.."
R=${R:-$HOME/miniforge3/envs/repurposing/bin/Rscript}
D=data/processed/lincs/d132
export N_CORES=8
for NULL in replicate rotation; do
  for SET in tf collectri; do
    if [ "$SET" = tf ]; then ENVV="VERSION=v2 REGSET=tf"; else ENVV="REGULON=collectri"; fi
    echo "=== $NULL / $SET ==="
    env $ENVV PROFILE=specific \
        PROFILE_FILE="$D/${NULL}_profiles.f32" INDEX_FILE="$D/${NULL}_index.tsv" \
        RUN_TAG="_d132${NULL:0:3}" N_CORES=8 "$R" scripts/79d_drug_viper.R
  done
done
