#!/usr/bin/env bash
# D94 network uncertainty (rule fixed 2026-11-03): consolidate each network from two disjoint halves of its 100 bootstraps
# (bootstrap files ordered by creation time = seed order: 1-50 and 51-100), as 43_run_aracne.sh consolidates (ARACNe-AP
# defaults, Bonferroni Poisson p < 0.05). Outputs data/processed/step5/aracne_v2_tf_halves/<NET>/{h1,h2}/network.txt;
# Jaccard in scripts/d94/d94_reproducibility.py. PARALLEL consolidations at once (single-threaded each).
set -euo pipefail
cd "$(dirname "$0")/../.."
JAVA=${ARACNE_JAVA:-~/miniforge3/envs/aracne/bin/java}; JAR=${ARACNE_JAR:-~/tools/ARACNe-AP/dist/aracne.jar}
NETS="BRCA CESC ESCA LGG PCPG SARC TGCT UCEC BRCA_LUM $(awk -F'\t' 'NR>1{print $1}' results/step5/d93_networks.tsv | paste -sd' ')"
one() {
  local n=$1 src=data/processed/step5/aracne_v2_tf/$1 out=data/processed/step5/aracne_v2_tf_halves/$1
  mapfile -t f < <(ls -tr $src/bootstrapNetwork_*.txt)
  [ ${#f[@]} -eq 100 ] || { echo "$n: ${#f[@]} bootstraps"; return 1; }
  for h in 1 2; do
    d=$out/h$h; [ -s $d/network.txt ] && continue; mkdir -p $d
    for i in $(seq $(( (h-1)*50 )) $(( h*50-1 ))); do ln -sf "$(realpath ${f[$i]})" $d/; done
    $JAVA -Xmx8G -jar $JAR -o $d --consolidate > $d/consolidate.log 2>&1
  done
  echo "$(date +%T) $n done"
}
export -f one; export JAVA JAR
printf "%s\n" $NETS | xargs -P ${PARALLEL:-6} -I{} bash -c 'one {}'
echo "HALVES DONE"
