#!/usr/bin/env bash
# Step 5d: ARACNe-AP networks per cohort (inputs from scripts/42_aracne_inputs.py).
#   MI threshold p = 1e-8 (seed 1); 100 bootstraps (seeds 1..100); DPI on; consolidation with
#   Bonferroni-corrected Poisson p < 0.05 (ARACNe-AP defaults). Resumable: finished bootstraps are kept.
# Usage: scripts/43_run_aracne.sh [NET ...]   (no args = queue below, 2 cohorts at a time)
#   VERSION=v2: inputs from aracne_inputs_v2/ (v2 gene set, rebuild spec D51), networks to aracne_v2/; the queue
#     adds BRCA_LUM (luminal-only breast network, D65).
#   REGSET=tf (with VERSION=v2): Lambert TFs only as regulators (<NET>.regulators_tf.txt), networks to aracne_v2_tf/ (D64).
# Tools: ARACNe-AP commit 30967ea built with OpenJDK 17 (conda env "aracne"); see docs/PROGRESS.md.
set -euo pipefail
cd "$(dirname "$0")/.."
SFX=$([ "${VERSION:-v1}" = v2 ] && echo _v2 || echo "")
IN=data/processed/step5/aracne_inputs$SFX
OUT=data/processed/step5/aracne$SFX
REGF=regulators
if [ "${REGSET:-full}" = tf ]; then
  [ "$SFX" = _v2 ] || { echo "REGSET=tf needs VERSION=v2" >&2; exit 1; }
  REGF=regulators_tf; OUT=${OUT}_tf
fi
JAVA=~/miniforge3/envs/aracne/bin/java
JAR=~/tools/ARACNe-AP/dist/aracne.jar
N_BOOT=100; THREADS=28; PARALLEL=2
# option (c): Steps 8-10 cohorts first, then small cohorts, then the rest (largest last)
QUEUE="BRCA COADREAD $([ "$SFX" = _v2 ] && echo BRCA_LUM) UCS KICH ACC WT TGCT NBL GBM PCPG PAAD ESCA UCEC SARC KIRP CESC LIHC OV BLCA \
STAD SKCM PRAD LUSC THCA LGG LUAD HNSC KIRC"

run_net() {
  local net=$1 d=$OUT/$1
  mkdir -p "$d"
  [ -s "$d/network.txt" ] && { echo "$(date +%F' '%T) $net already done"; return; }
  echo "$(date +%F' '%T) $net start"
  ls "$d"/miThreshold_* >/dev/null 2>&1 || \
    $JAVA -Xmx40G -jar $JAR -e $IN/$net.exp.tsv -o "$d" --tfs $IN/$net.$REGF.txt \
      --pvalue 1E-8 --seed 1 --calculateThreshold > "$d/threshold.log" 2>&1
  local have; have=$(ls "$d"/bootstrapNetwork_*.txt 2>/dev/null | wc -l)
  for ((i = have + 1; i <= N_BOOT; i++)); do
    $JAVA -Xmx60G -jar $JAR -e $IN/$net.exp.tsv -o "$d" --tfs $IN/$net.$REGF.txt \
      --pvalue 1E-8 --seed $i --threads $THREADS >> "$d/bootstraps.log" 2>&1
    echo "$(date +%F' '%T) $net bootstrap $i/$N_BOOT"
  done
  $JAVA -Xmx60G -jar $JAR -o "$d" --consolidate > "$d/consolidate.log" 2>&1
  echo "$(date +%F' '%T) $net done: $(($(wc -l < "$d/network.txt") - 1)) edges"
}
export -f run_net; export IN OUT REGF JAVA JAR N_BOOT THREADS

nets=${*:-$QUEUE}
printf '%s\n' $nets | xargs -P $PARALLEL -I{} bash -c 'run_net {}'
echo "$(date +%F' '%T) ALL DONE"
