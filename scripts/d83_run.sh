#!/usr/bin/env bash
# D83 comparison driver: after the three group-VIPER arms (main, 50_viper.R REGSET=tf_<arm> SEED_GROUP=1 GROUP_ONLY=1)
# finish: known-biology controls (51), sync to this worktree, drug activity (79d), proliferation regulators (79h),
# consensus ranking (79e) and split (79i) per arm, then the comparison (79s). Logs: logs/d83_*.log.
set -euo pipefail
cd "$(dirname "$0")/.."
MAIN=$HOME/projects/repurposing; R=$HOME/miniforge3/envs/repurposing/bin/Rscript; PY=$HOME/miniforge3/envs/celligner/bin/python
ARMS="base sig1 sig2"
for a in $ARMS; do until grep -q "VIPER DONE" $MAIN/logs/viper_tf_$a.log; do grep -q "Execution halted" $MAIN/logs/viper_tf_$a.log && { echo "VIPER $a failed"; exit 1; }; sleep 30; done; done
echo "$(date +%T) VIPER done"
for a in $ARMS; do (cd $MAIN && REGULON=aracne_project_v2_tf_$a $R scripts/51_compare_mr.R > logs/d83_controls_$a.log 2>&1); done
for a in $ARMS; do d=aracne_project_v2_tf_$a; mkdir -p results/step6/$d; cp $MAIN/results/step6/$d/{group_mr,regulator_filter,group_summary}.tsv results/step6/$d/; done
echo "$(date +%T) controls + sync done"
for a in $ARMS; do VERSION=v2 REGSET=tf_$a N_CORES=12 $R scripts/79d_drug_viper.R > logs/d83_drug_viper_$a.log 2>&1 & done; wait
echo "$(date +%T) drug viper done"
for a in $ARMS; do VERSION=v2 REGSET=tf_$a $R scripts/79h_prolif_regulators.R > logs/d83_prolif_$a.log 2>&1 & done; wait
echo "$(date +%T) prolif done"
for a in $ARMS; do VERSION=v2 REGULON=consensus ARACNE_SET=tf_$a $PY scripts/79e_drug_reversal.py > logs/d83_ranking_$a.log 2>&1 & done; wait
echo "$(date +%T) ranking done"
for a in $ARMS; do VERSION=v2 REGULON=consensus ARACNE_SET=tf_$a $PY scripts/79i_split_reversal.py > logs/d83_split_$a.log 2>&1 & done; wait
echo "$(date +%T) split done"
until [ -s results/lincs/knockdown_benchmark_compare_d83.tsv ]; do sleep 30; done
$PY scripts/79s_d83_compare.py > logs/d83_compare.log 2>&1
echo "$(date +%T) COMPARE DONE"
