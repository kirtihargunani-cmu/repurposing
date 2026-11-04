#!/usr/bin/env bash
# D94 downstream (rule fixed 2026-11-03): one isolated git worktree per arm (../repurposing_d94_<arm>), production file
# names filled with the arm's consensus outputs, then the downstream steps. data/ is linked to the main checkout (the
# downstream steps only add to append-only caches under data/raw); results/ is the worktree's own, so production results
# are never written. Usage: d94_downstream.sh <arm> <tag>   e.g. split _d94split
set -euo pipefail
ARM=$1; TAG=$2
MAIN=$(cd "$(dirname "$0")/../.." && pwd -P); WT=$(dirname "$MAIN")/repurposing_d94_$ARM
PY=~/miniforge3/envs/celligner/bin/python
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
cd "$MAIN"
[ -d "$WT" ] || git worktree add --detach "$WT" HEAD
cd "$WT"
[ -L data ] || { rm -rf data; ln -s "$MAIN/data" data; }
mkdir -p logs
rsync -a --ignore-existing "$MAIN/results/" results/        # untracked inputs (e.g. .gz) from the main checkout
L=results/lincs
for b in drug_reversal_v2_consensus drug_reversal_summary_v2_consensus drug_reversal_controls_v2_consensus \
         drug_reversal_split_v2_consensus split_controls_v2_consensus split_overview_v2_consensus prolif_regulators_p2_v2_consensus \
         stratum_reversal_v2_consensus stratum_reversal_summary_v2_consensus stratum_reversal_controls_v2_consensus stratum_reversal_split_v2_consensus; do
  src="$MAIN/$L/${b}${TAG}.tsv"; [ -f "$src" ] || { echo "missing $src"; exit 1; }
  cp "$src" "$L/$b.tsv"
done
for s in 79w_receptor_gate 79p_fda_ranking 79l_known_use 79n_dependency 79m_gate_sensitivity 79x_rankings_page; do
  echo "$(date +%T) $s"; $PY scripts/$s.py > logs/d94_$s.log 2>&1 || { echo "FAILED $s (see $WT/logs/d94_$s.log)"; exit 1; }
done
echo "D94 DOWNSTREAM DONE $ARM"
