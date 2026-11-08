#!/usr/bin/env bash
# Runbook for the v2 drug side (docs/PIPELINE_DRUG_SIDE_V2.md explains every step, its decisions and outputs).
# Each step is one command with its switches, in pipeline order; the step scripts hold the analysis.
# Usage:
#   scripts/run_v2_drug.sh list            list the steps
#   scripts/run_v2_drug.sh <step> [...]    run steps in the order given; stops at the first failure
#   scripts/run_v2_drug.sh show <step>     print the command without running it
# Runs from main or from the celligner worktree. Logs: logs/v2_drug_<step>.log in this checkout. Env overrides: PY (celligner python), R (Rscript), MAIN (main checkout).
# The disease side (main branch, scripts/run_v2.sh) must have run first; in the celligner worktree, sync_disease copies its outputs
# there (on main it does nothing).
# Celligner itself (70-74) and the LINCS profiles (76, 76b, 79a, 79c) do not change with v2 and are not repeated here.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-$HOME/miniforge3/envs/celligner/bin/python}
R=${R:-$HOME/miniforge3/envs/repurposing/bin/Rscript}
MAIN=${MAIN:-$HOME/projects/repurposing}

SYNC="set -e; if [ \"$(realpath $MAIN)\" = \"$(pwd -P)\" ]; then echo 'on main: disease-side outputs already here, nothing to copy'; exit 0; fi; for d in aracne_project_v2_tf collectri_v2 aracne_project_v2_tf_strata collectri_v2_strata; do mkdir -p results/step6/\$d; cp $MAIN/results/step6/\$d/{group_mr,regulator_filter,group_summary}.tsv results/step6/\$d/; done; mkdir -p results/step4_v2 results/step10_v2; cp $MAIN/results/step4_v2/signature_summary.tsv results/step4_v2/; cp $MAIN/results/step10_v2/{strata,stratum_patients}.tsv results/step10_v2/"

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
  "receptor_gate|Receptor evidence gate: ChEMBL binding, assembly, receptor in the cohort's matched lines; selectivity window|$PY scripts/79w_receptor_gate.py"
  # D117: gates and flags after the split, in this order (a partial rerun must resume from the first step it changes)
  "fda_matches|Drugs@FDA matching (fda_matches.tsv, read by the target-presence gate); 79p also writes provisional FDA lists|$PY scripts/79p_fda_ranking.py"
  "target_presence|Target presence in the typical tumor, FDA main targets (79i2: D102a/c, D105b, D110, D113), three timing variants|for t in _v2_consensus _v2_consensus_d100early _v2_consensus_d100late; do SFX=\$t $PY scripts/79i2_target_presence.py; done"
  "receptor_lines|Receptor present in the matched lines (79i3: D104b/c), three timing variants|for t in _v2_consensus _v2_consensus_d100early _v2_consensus_d100late; do SFX=\$t $PY scripts/79i3_receptor_lines.py; done"
  "timing_decision|Early / late timing scores vs aggregated (D100, reported)|$PY scripts/d100/d100_decide.py"
  "evidence_flags|Evidence flags: dose, inactive matched line per drug (D116), generic mechanism; openFDA label queue (79f3)|STAGE=flags $PY scripts/79f3_evidence_flags.py"
  "plate_confined|Plate-confined activity flag (79f4, D109a); must follow evidence_flags|$PY scripts/79f4_plate_confined.py"
  "fda_docs|Drugs@FDA label / review passages for drugs without a label Cmax (79f5; readings are a manual subagent pass, cmax_docs/)|$PY scripts/79f5_fda_docs.py"
  "exposure|Cmax from verified label / document readings -> low_exposure (79f3, D109b, D111); must follow plate_confined|STAGE=exposure $PY scripts/79f3_evidence_flags.py"
  "tiers|Final tiers and evidence status after all gates and flags (79i4, D117)|$PY scripts/79i4_tiers.py"
  "data_gaps|Where good-quality lineage-matched data run out, per group and subgroup (C24)|$PY scripts/79r_data_gaps.py"
  "literature|PubMed / trial links for each disease's top 10 by specificity and by reversal|$PY scripts/79o_literature.py"
  "fda_ranking|FDA-approved drug rankings per disease and subgroup (Drugs@FDA)|$PY scripts/79p_fda_ranking.py"
  "literature_fda|PubMed / trial links for each cohort's top 3 FDA drugs, by FDA name|MODE=fda $PY scripts/79o_literature.py"
  "literature_fda_all|PubMed / trial links for every FDA drug row of the rankings page, by FDA name (cached, resumable)|MODE=fda_all $PY scripts/79o_literature.py"
  "trials_fetch|Trials of each cohort's top 10 FDA drugs: stop reasons, completed-trial reading queue (79o2, D108)|STAGE=fetch $PY scripts/79o2_trial_outcomes.py"
  "trials_merge|Merge verified trial readings and stopped-trial roles (79o2; readings are a manual subagent pass, d108/readings, d108/roles)|STAGE=merge $PY scripts/79o2_trial_outcomes.py"
  "standin|Closest-compound stand-ins for borrowed / weak FDA drugs, with validation (79y)|$PY scripts/79y_standin.py"
  "standin_vs_borrowed|Stand-in vs the drug's own borrowed score: borrowed counterfactual in a sandbox (79e, 79i) + comparison (79y2)|PY=$PY scripts/79y2_borrowed_counterfactual.sh"
  "fda_page|Browsable FDA rankings page (results/lincs/fda_rankings_v2.html)|$PY scripts/79q_fda_page.py"
  "rankings_page|Browsable FDA drug x disease rankings with all scores, flags, receptors, literature (disease_drug_rankings_v2.html)|$PY scripts/79x_rankings_page.py"
  "standin_page|Browsable LINCS compounds tested in each cohort's lines with their closest FDA drugs (standin_compounds_v2.html)|$PY scripts/79z_standin_page.py"
  "known_use|Known-use recovery: do a cancer's approved drugs rank high for it? Hub + ChEMBL indications (79l, C28)|$PY scripts/79l_known_use.py"
  "dependency|CRISPR dependency of the drugs' targets in each cohort's lines: yardstick, gate test (79n, C29)|$PY scripts/79n_dependency.py"
  "target_role|Target-role score (dependency, alteration, over-activity of the main target), decoded; known-use analysis (79t, D120; reported only)|$PY scripts/79t_target_role.py"
  "gate_sensitivity|PRISM / GDSC selective killing vs the gated rankings: gate, matched vs borrowed, stand-ins (79m, C28)|$PY scripts/79m_gate_sensitivity.py"
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
