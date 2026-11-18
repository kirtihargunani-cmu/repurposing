#!/bin/bash
# Gate F2 — validate the no-good-cut enumerator against fixtures with a KNOWN number of
# equally-good explanations. The spec requires this: an ambiguity count nobody else checks is
# only trustworthy if it reproduces cases where the true answer is known by construction.
#
# INTERPRETERS: R and python both from /home/kirti/miniforge3/envs/sysval
set -u
cd /home/kirti/projects/repurposing
E=/home/kirti/miniforge3/envs/sysval
D=$(mktemp -d); OUT=results/systems_validation/synthetic
mkdir -p "$OUT"

$E/bin/Rscript -e '
suppressPackageStartupMessages(library(CARNIVAL))
D <- commandArgs(trailingOnly=TRUE)[1]
o <- defaultCbcSolveCarnivalOptions(); o$solverPath <- "/home/kirti/miniforge3/envs/sysval/bin/cbc"
o$outputFolder <- D; o$workdir <- D; o$keepLPFiles <- TRUE
mk <- function(name, net, meas, pert=NULL) {
  dir.create(file.path(D,name), showWarnings=FALSE)
  o2 <- o; o2$outputFolder <- file.path(D,name); o2$workdir <- file.path(D,name)
  if (is.null(pert)) generateLpFileCarnival(measurements=meas, priorKnowledgeNetwork=net, carnivalOptions=o2)
  else generateLpFileCarnival(perturbations=pert, measurements=meas, priorKnowledgeNetwork=net, carnivalOptions=o2)
}
n <- function(s,i,t) data.frame(source=s, interaction=i, target=t, stringsAsFactors=FALSE)
# 1 unique: a single chain, one explanation
mk("unique", n(c("X","A"), c(1,1), c("A","B")), c(B=1), c(X=1))
# 2 two interchangeable sources
mk("two", n(c("S1","S2","H","H"), c(1,1,1,1), c("H","H","T1","T2")), c(T1=1,T2=1))
# 3 three interchangeable sources
mk("three", n(c("S1","S2","S3","H","H"), c(1,1,1,1,1), c("H","H","H","T1","T2")), c(T1=1,T2=1))
# 4 two INDEPENDENT binary choices -> 2 x 2 = 4 explanations
mk("four", n(c("A1","A2","HA","B1","B2","HB"), c(1,1,1,1,1,1), c("HA","HA","TA","HB","HB","TB")),
   c(TA=1,TB=1))
' "$D" >/dev/null 2>&1

declare -A EXPECT=( [unique]=1 [two]=2 [three]=3 [four]=4 )
echo "fixture    expected  found  status      core  union"
pass=0; total=0
for f in unique two three four; do
  lp=$(ls "$D/$f"/*.lp 2>/dev/null | head -1)
  if [ -z "$lp" ]; then echo "$f: NO LP GENERATED"; continue; fi
  j=$($E/bin/python scripts/systems_validation/fault_isolation/enumerate_optima.py \
        --lp "$lp" --max-solutions 40 --out "$D/$f.json" >/dev/null 2>&1; cat "$D/$f.json")
  got=$($E/bin/python -c "import json,sys;d=json.loads(sys.argv[1]);print(d['n_alternative_optima'])" "$j")
  st=$($E/bin/python -c "import json,sys;d=json.loads(sys.argv[1]);print(d['status'])" "$j")
  core=$($E/bin/python -c "import json,sys;d=json.loads(sys.argv[1]);print(d['n_core'])" "$j")
  uni=$($E/bin/python -c "import json,sys;d=json.loads(sys.argv[1]);print(d['n_union'])" "$j")
  exp=${EXPECT[$f]}; total=$((total+1))
  if [ "$got" = "$exp" ]; then mark="PASS"; pass=$((pass+1)); else mark="FAIL"; fi
  printf "%-10s %8s  %5s  %-11s %4s  %5s   [%s]\n" "$f" "$exp" "$got" "$st" "$core" "$uni" "$mark"
  cp "$D/$f.json" "$OUT/gate_f2_enum_$f.json" 2>/dev/null
done
echo
echo "GATE F2: $pass/$total fixtures recovered the known number of equally-good explanations"
