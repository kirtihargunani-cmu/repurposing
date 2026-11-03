#!/usr/bin/env bash
# D86 perturbation-network pilot: ARACNe-AP on LINCS drug profiles within a line (inputs: scripts/pnet/pnet_inputs.py).
# Same settings as scripts/43_run_aracne.sh: MI threshold p = 1e-8 (seed 1); 100 bootstraps (seeds 1..100); DPI on;
# consolidation (Bonferroni-corrected Poisson p < 0.05). Regulators: data/processed/lincs/pnet/regulators.txt.
# One network per <LINE>_<arm>/<set> (set = all, f1..f5) in data/processed/lincs/pnet/<LINE>_<arm>/net_<set>/.
# Threads: THREADS per network x PARALLEL networks; default 8 x 3 = 24 (user: under half the machine's 64 threads).
# Resumable: finished bootstraps are kept, finished networks skipped.
set -euo pipefail
cd "$(dirname "$0")/../.."
D=${PNET_DIR:-data/processed/lincs/pnet}   # PNET_DIR: other network sets (e.g. pnet_depmap, 2026-11-03)
JAVA=${ARACNE_JAVA:-~/miniforge3/envs/aracne/bin/java}
JAR=${ARACNE_JAR:-~/tools/ARACNe-AP/dist/aracne.jar}
N_BOOT=100; THREADS=${THREADS:-8}; PARALLEL=${PARALLEL:-3}
run_one() {
  local exp=$1 d; d="$(dirname "$exp")/net_$(basename "$exp" .exp.tsv)"
  mkdir -p "$d"
  [ -s "$d/network.txt" ] && { echo "$(date +%F' '%T) $d already done"; return; }
  ls "$d"/miThreshold_* >/dev/null 2>&1 || \
    $JAVA -Xmx8G -jar $JAR -e "$exp" -o "$d" --tfs $D/regulators.txt --pvalue 1E-8 --seed 1 --calculateThreshold > "$d/threshold.log" 2>&1
  local have; have=$(ls "$d"/bootstrapNetwork_*.txt 2>/dev/null | wc -l)
  for ((i = have + 1; i <= N_BOOT; i++)); do
    $JAVA -Xmx8G -jar $JAR -e "$exp" -o "$d" --tfs $D/regulators.txt --pvalue 1E-8 --seed $i --threads $THREADS >> "$d/bootstraps.log" 2>&1
  done
  $JAVA -Xmx8G -jar $JAR -o "$d" --consolidate > "$d/consolidate.log" 2>&1
  echo "$(date +%F' '%T) $d done: $(($(wc -l < "$d/network.txt") - 1)) edges"
}
export -f run_one; export D JAVA JAR N_BOOT THREADS
ls $D/*/*.exp.tsv | xargs -P $PARALLEL -I{} bash -c 'run_one {}'
echo "$(date +%F' '%T) ALL DONE"
