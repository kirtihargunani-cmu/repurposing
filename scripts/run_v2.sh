#!/usr/bin/env bash
# Runbook for the v2 disease-side pipeline (docs/PIPELINE_V2.md explains every step, its decisions and outputs).
# Each step is one command with its switches, in pipeline order; the step scripts hold the analysis.
# Usage:
#   scripts/run_v2.sh list             list the steps (name, what it does)
#   scripts/run_v2.sh <step> [...]     run one or more steps in the order given; stops at the first failure
#   scripts/run_v2.sh show <step>      print the command without running it
# Logs: logs/v2_<step>.log. Env overrides: R (Rscript), PY (python), N_CORES.
# Long steps: aracne_build / aracne_build_tf take days (run on a large machine; resumable).
# Drug side: celligner branch, scripts/run_v2_drug.sh.
set -euo pipefail
cd "$(dirname "$0")/.."
R=${R:-$HOME/miniforge3/envs/repurposing/bin/Rscript}
PY=${PY:-$HOME/miniforge3/envs/repurposing/bin/python}
N_CORES=${N_CORES:-16}
ST=results/step10_v2/stratum_patients.tsv

STEPS=(
  "covariates|Step 3: covariates incl. GTEx collection variables|$PY scripts/20_covariates.py"
  "gene_set|Step 3: v2 gene set (no mappability < 0.8, no readthrough; D46)|$PY scripts/22_gene_set_v2.py"
  "gtex_reference|Step 4a: GTEx reference (sub-region filter D59, anchor and collection covariates D57)|$R scripts/32_gtex_reference.R"
  "signatures|Step 4b: tumor-vs-normal signatures, voomByGroup (D69)|VERSION=v2 N_CORES=$N_CORES $R scripts/30_signatures.R"
  "signature_checks|Step 4: checks C1 C3 C4 C6 C11 C12|N_CORES=$N_CORES $R scripts/37_step4_v2_checks.R"
  "aracne_inputs|Step 5: ARACNe inputs per cohort (+ BRCA_LUM, D65)|VERSION=v2 $PY scripts/42_aracne_inputs.py"
  "aracne_build|Step 5: ARACNe-AP networks, full regulator set (~2.5 days)|VERSION=v2 bash scripts/43_run_aracne.sh"
  "aracne_build_tf|Step 5: ARACNe-AP networks, TF-only set (D64; ~0.7 days)|VERSION=v2 REGSET=tf bash scripts/43_run_aracne.sh"
  "regulons|Step 5: regulons, cross-mappable edges removed (full set)|VERSION=v2 $R scripts/44_aracne_regulons.R"
  "regulons_tf|Step 5: regulons, cross-mappable edges removed (TF-only)|VERSION=v2 REGSET=tf $R scripts/44_aracne_regulons.R"
  "edge_support|Step 5: bootstrap support per edge (full set)|VERSION=v2 $R scripts/45_edge_support.R"
  "edge_support_tf|Step 5: bootstrap support per edge (TF-only)|VERSION=v2 REGSET=tf $R scripts/45_edge_support.R"
  "edge_confounding|Step 5: edge strength after purity / quality / plate (full set)|VERSION=v2 $R scripts/46_edge_confounding.R"
  "viper_tf|Step 6a: group + patient VIPER, ARACNe TF-only (primary)|REGULON=aracne_project VERSION=v2 REGSET=tf N_CORES=$N_CORES $R scripts/50_viper.R"
  "viper_full|Step 6a: group + patient VIPER, ARACNe full set|REGULON=aracne_project VERSION=v2 N_CORES=$N_CORES $R scripts/50_viper.R"
  "viper_collectri|Step 6a: group + patient VIPER, CollecTRI|REGULON=collectri VERSION=v2 N_CORES=$N_CORES $R scripts/50_viper.R"
  "viper_dorothea|Step 6a: group + patient VIPER, DoRothEA A-C|REGULON=dorothea_ABC VERSION=v2 N_CORES=$N_CORES $R scripts/50_viper.R"
  "null_test|Step 6a: covariate-adjusted null test (D73; result: not adopted)|$R scripts/64_null_adjust_test.R"
  "reference_swap|Step 6a: GTEx vs adjacent reference comparison|N_CORES=$N_CORES $R scripts/65_reference_swap.R"
  "purity_tf|Step 6b: NES-level purity correction, TF-only|REGULON=aracne_project_v2_tf $R scripts/56_purity_correct.R"
  "patient_tf|Step 6b: patient VIPER with purity + quality correction, TF-only|VERSION=v2 REGSET=tf N_CORES=$N_CORES $R scripts/61_quality_correct.R"
  "patient_full|Step 6b: patient VIPER with purity + quality correction, full set|VERSION=v2 N_CORES=$N_CORES $R scripts/61_quality_correct.R"
  "step6_checks|Step 6: checks C2 C5 C7 C9 C10, v1 vs v2 full|SET=full $R scripts/63_step6_v2_checks.R"
  "step6_checks_tf|Step 6: checks C2 C5 C7 C9 C10, v2 full vs TF-only|SET=tf $R scripts/63_step6_v2_checks.R"
  "compare_mr|Step 7: master-regulator summaries + known-biology controls, all v2 sets|for s in aracne_project_v2_tf aracne_project_v2 collectri_v2 dorothea_ABC_v2 aracne_project_v2_tf_strata collectri_v2_strata dorothea_ABC_v2_strata; do REGULON=\$s $R scripts/51_compare_mr.R || exit 1; done"
  "compare_sets|Step 7: three regulon sets compared; D18 confirmed regulators (TF arm, full, subgroups)|VERSION=v2 REGSET=tf $R scripts/59_compare_regulon_sets.R && VERSION=v2 $R scripts/59_compare_regulon_sets.R && VERSION=v2 REGSET=tf STRATA=1 $R scripts/59_compare_regulon_sets.R"
  "modules_tf|Step 7: regulator modules, TF arm|NES_DIR=data/processed/step6/aracne_project_v2_tf_patient $R scripts/58_regulator_modules.R"
  "edge_confounding_tf|Step 5: edge strength after purity / quality / plate (TF-only)|VERSION=v2 REGSET=tf $R scripts/46_edge_confounding.R"
  "recurring_module|Step 7: recurring-module diagnostics (C8)|NES_DIR=data/processed/step6/aracne_project_v2_patient $R scripts/62_recurring_module.R"
  "clustering|Step 8: consensus clustering, active build (D30)|N_CORES=$N_CORES $R scripts/80_consensus_clustering.R"
  "clustering_no_c8|Step 8: consensus clustering without the recurring module|EXCLUDE_MODULE=recurring N_CORES=$N_CORES $R scripts/80_consensus_clustering.R"
  "consistency|Step 8: subgroup keep rule metrics (D71)|N_CORES=$N_CORES $R scripts/84_subgroup_consistency.R"
  "consistency_no_c8|Step 8: keep rule metrics, no-module clustering|STEP8_SUFFIX=_no_c8 N_CORES=$N_CORES $R scripts/84_subgroup_consistency.R"
  "strata|Step 10: subgroup table (D79)|$R scripts/85_step10_strata_v2.R"
  "strata_signatures|Step 10: subgroup signatures vs the parent's normals|VERSION=v2 STRATA=$ST N_CORES=$N_CORES $R scripts/30_signatures.R"
  "strata_viper_tf|Step 10: subgroup VIPER, ARACNe TF-only|REGULON=aracne_project VERSION=v2 REGSET=tf STRATA=$ST N_CORES=14 $R scripts/50_viper.R"
  "strata_viper_collectri|Step 10: subgroup VIPER, CollecTRI|REGULON=collectri VERSION=v2 STRATA=$ST N_CORES=14 $R scripts/50_viper.R"
  "strata_viper_dorothea|Step 10: subgroup VIPER, DoRothEA|REGULON=dorothea_ABC VERSION=v2 STRATA=$ST N_CORES=14 $R scripts/50_viper.R"
  "strata_variability|Step 10: variability inside subgroups|$R scripts/86_subgroup_variability.R"
  "strata_characterize|Step 10: what defines each subgroup|$R scripts/87_subgroup_characterize.R"
  "status|Show which outputs exist for the active build|$R scripts/00_status.R"
)

field() { local name=$1 k=$2; for s in "${STEPS[@]}"; do IFS='|' read -r n d c <<< "$s"; [[ $n == "$name" ]] && { [[ $k == d ]] && echo "$d" || echo "$c"; return 0; }; done; return 1; }

case "${1:-list}" in
  list) for s in "${STEPS[@]}"; do IFS='|' read -r n d c <<< "$s"; printf '  %-24s %s\n' "$n" "$d"; done ;;
  show) shift; for n in "$@"; do field "$n" c || { echo "unknown step: $n" >&2; exit 1; }; done ;;
  *) mkdir -p logs
     for n in "$@"; do
       cmd=$(field "$n" c) || { echo "unknown step: $n (scripts/run_v2.sh list)" >&2; exit 1; }
       echo "$(date +%F' '%T) START $n: $cmd" | tee -a logs/v2_runbook.log
       if bash -c "$cmd" > "logs/v2_$n.log" 2>&1; then echo "$(date +%F' '%T) OK    $n" | tee -a logs/v2_runbook.log
       else echo "$(date +%F' '%T) FAIL  $n (see logs/v2_$n.log)" | tee -a logs/v2_runbook.log; exit 1; fi
     done ;;
esac
