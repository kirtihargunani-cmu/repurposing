#!/usr/bin/env bash
# Runbook for the v2 drug side (docs/PIPELINE_DRUG_SIDE_V2.md explains every step, its decisions and outputs).
# Each step is one command with its switches, in pipeline order; the step scripts hold the analysis.
# Usage:
#   scripts/run_v2_drug.sh list            list the steps
#   scripts/run_v2_drug.sh <step> [...]    run steps in the order given; stops at the first failure
#   scripts/run_v2_drug.sh show <step>     print the command without running it
# Logs: logs/v2_drug_<step>.log in this worktree. Env overrides: PY (celligner python), R (Rscript), MAIN (main checkout).
# The disease side (main branch, scripts/run_v2.sh) must have run first; sync_disease copies its outputs here.
# Celligner itself (70-74) and the LINCS profiles (76, 76b, 79a, 79c) do not change with v2 and are not repeated here.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-$HOME/miniforge3/envs/celligner/bin/python}
R=${R:-$HOME/miniforge3/envs/repurposing/bin/Rscript}
MAIN=${MAIN:-$HOME/projects/repurposing}

SYNC="set -e; for d in aracne_project_v2_tf collectri_v2 aracne_project_v2_tf_strata collectri_v2_strata; do mkdir -p results/step6/\$d; cp $MAIN/results/step6/\$d/{group_mr,regulator_filter,group_summary}.tsv results/step6/\$d/; done; mkdir -p results/step4_v2 results/step10_v2; cp $MAIN/results/step4_v2/signature_summary.tsv results/step4_v2/; cp $MAIN/results/step10_v2/{strata,stratum_patients}.tsv results/step10_v2/"

STEPS=(
  "sync_disease|Copy the disease-side inputs (master regulators, filters, subgroups) from the main checkout|$SYNC"
  "drug_viper_tf|Drug x cell-line regulator activity, ARACNe TF-only regulons|VERSION=v2 REGSET=tf $R scripts/79d_drug_viper.R"
  "drug_viper_collectri|Drug x cell-line regulator activity, CollecTRI (D75)|REGULON=collectri $R scripts/79d_drug_viper.R"
  "drug_viper_dorothea|Drug x cell-line regulator activity, DoRothEA A-C|REGULON=dorothea_ABC $R scripts/79d_drug_viper.R"
  "prolif_regulators_tf|Proliferation-linked regulators, ARACNe TF-only|VERSION=v2 REGSET=tf $R scripts/79h_prolif_regulators.R"
  "prolif_regulators_collectri|Proliferation-linked regulators, CollecTRI|REGULON=collectri $R scripts/79h_prolif_regulators.R"
  "knockdown_subset|LINCS TF knockdown signatures (benchmark input)|$PY scripts/79m_knockdown_subset.py"
  "knockdown_benchmark|Knockdown benchmark of regulon sets (reported, not used; D75)|$R scripts/79n_knockdown_benchmark.R"
  "subgroup_line_weights|Celligner line weights for the subgroups|VERSION=v2 $PY scripts/74b_stratum_line_weights.py"
  "subgroup_target_gate|Target check on each subgroup's own tumors|STRATA=1 VERSION=v2 $PY scripts/79b_target_gate.py"
  "transfer_confidence|Borrowed-drug transfer confidence and low_transfer flag (D77)|VERSION=v2 $PY scripts/79l_transfer_confidence.py"
  "ranking|Drug ranking per disease: consensus, NES-weighted, flags, evidence status (D74-D78)|VERSION=v2 REGULON=consensus $PY scripts/79e_drug_reversal.py"
  "ranking_subgroups|Drug ranking per subgroup, with sibling specificity|VERSION=v2 REGULON=consensus STRATA=1 $PY scripts/79e_drug_reversal.py"
  "split|Proliferation vs mechanism split, diseases and subgroups (spec_mech, spec_sib_mech)|VERSION=v2 REGULON=consensus STRATA=1 $PY scripts/79i_split_reversal.py"
  "data_gaps|Where good-quality lineage-matched data run out, per group and subgroup (C24)|$PY scripts/79r_data_gaps.py"
  "literature|PubMed / trial links for each disease's top 10 by specificity and by reversal|$PY scripts/79o_literature.py"
  "fda_ranking|FDA-approved drug rankings per disease and subgroup (Drugs@FDA)|$PY scripts/79p_fda_ranking.py"
  "literature_fda|PubMed / trial links for each cohort's top 3 FDA drugs, by FDA name|MODE=fda $PY scripts/79o_literature.py"
  "fda_page|Browsable FDA rankings page (results/lincs/fda_rankings_v2.html)|$PY scripts/79q_fda_page.py"
)

field() { local name=$1 k=$2; for s in "${STEPS[@]}"; do IFS='|' read -r n d c <<< "$s"; [[ $n == "$name" ]] && { [[ $k == d ]] && echo "$d" || echo "$c"; return 0; }; done; return 1; }

case "${1:-list}" in
  list) for s in "${STEPS[@]}"; do IFS='|' read -r n d c <<< "$s"; printf '  %-28s %s\n' "$n" "$d"; done ;;
  show) shift; for n in "$@"; do field "$n" c || { echo "unknown step: $n" >&2; exit 1; }; done ;;
  *) mkdir -p logs
     for n in "$@"; do
       cmd=$(field "$n" c) || { echo "unknown step: $n (scripts/run_v2_drug.sh list)" >&2; exit 1; }
       echo "$(date +%F' '%T) START $n" | tee -a logs/v2_drug_runbook.log
       if bash -c "$cmd" > "logs/v2_drug_$n.log" 2>&1; then echo "$(date +%F' '%T) OK    $n" | tee -a logs/v2_drug_runbook.log
       else echo "$(date +%F' '%T) FAIL  $n (see logs/v2_drug_$n.log)" | tee -a logs/v2_drug_runbook.log; exit 1; fi
     done ;;
esac
