#!/usr/bin/env bash
# THE formal end-to-end pipeline, post ARACNe / VIPER (D186, 2026-11-22).
#
# User, 2026-11-22: "make sure that all the most recent changes are updated to produce the
# results. i do not want to be working with results that haven't gone through the correct
# filtering steps so let's formalize our pipeline with the decisions we have made so far and
# run it end to end post the ARACNe/VIPER builds."
#
# WHY THIS EXISTS. run_rebuild_chain.sh (D156) covered 79e -> 79i -> 79i2 -> 79t3 -> 79p plus the
# D179 nested arrays, and stopped there. Everything built since - the window arm (window_reversal
# -> 79j -> 79k -> 79e3 -> 79x2), the labelling arm (79aa) and the two registered tests (79ac
# D183, 79ad D184) - lived in NO runner and had only ever been run by hand. That is exactly how a
# result gets produced from a table that skipped a filtering step.
#
# WHAT IS UPSTREAM AND IS **NOT** RERUN HERE (the "post ARACNe/VIPER" boundary):
#   44_aracne_regulons.R     ARACNe networks
#   50_viper.R               disease VIPER (group + patient + strata)   -> group_mr.tsv
#   79d_drug_viper.R         drug x cell-line AND condition-level VIPER -> drug_viper_*_cond
#   celligner line weights, 79b target gate, 79w receptor gate, 79f3/79f4 flags, 79y stand-ins
# The last group is upstream-stable BY THE SAME ARGUMENT D156 USED: those stages depend on drug
# annotation, ChEMBL binding and Cmax caches, not on the regulons, so a regulon-arm rebuild does
# not invalidate them. Per the 2026-11-22 scope decision, the network / manual-reading stages
# (79o literature, 79o2 trials, 79f5 FDA docs, 79ab ChEMBL name fetch) are NOT run here either;
# they reuse their caches and stay manual.
#
# D186: _v2_tf (ARACNe-only, D160) IS PRODUCTION. 79i2 / 79i3 / 79p / 79t3 / 79i4 now default to
# it. D208 (2026-11-23, user: "remove the consensus arm completely"): the consensus arm is GONE -
# its outputs are deleted (recoverable from git history only) and every script refuses SFX/REGULON
# values naming it. The overactivity component no longer averages in CollecTRI NES.
#
#   scripts/run_pipeline.sh              # the whole thing
#   SKIP_GATE=1 scripts/run_pipeline.sh  # with a known, accepted failing pre-flight check
#   PAGES=1 scripts/run_pipeline.sh      # also rebuild the browsable HTML pages
#   FROM=windows scripts/run_pipeline.sh # resume at a stage group: disease|core|windows|labels|tests|pages
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
RES=results/lincs
PY=/home/kirti/miniforge3/envs/celligner/bin/python
R=/home/kirti/miniforge3/envs/repurposing/bin/Rscript

# D165: one thread per worker. Workers used 4-8 BLAS threads each and breached the <32 of 64
# budget; the cap belongs here, not in each script.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VERSION=v2 REGSET=tf REGULON=aracne_project N_CORES=${N_CORES:-6}

FROM=${FROM:-disease}
PAGES=${PAGES:-0}
declare -A ORDER=([disease]=1 [core]=2 [windows]=3 [labels]=4 [tests]=5 [pages]=6)
[[ -n "${ORDER[$FROM]:-}" ]] || { echo "bad FROM=$FROM" >&2; exit 2; }
want() { (( ${ORDER[$1]} >= ${ORDER[$FROM]} )); }

if [[ "${SKIP_GATE:-0}" != 1 ]]; then
  echo "== pre-flight: check_builds.py --strict"
  python scripts/check_builds.py --strict || {
    echo "ABORT: stale or crashed builds upstream. Fix, or re-run with SKIP_GATE=1 if accepted." >&2
    exit 1; }
fi

# D188 parameter gate. DELIBERATELY NOT under SKIP_GATE: the staleness gate is skipped routinely
# (d108's manual readings keep it red), and if a parameter check rode on the same switch it would
# be skipped routinely too - which is precisely how the TAS floor's state went unverified. A stale
# build is a bookkeeping problem; a build that used the wrong threshold is a wrong result.
# It needs pyyaml, so it runs under $PY rather than the bare interpreter.
if [[ "${SKIP_PARAM_GATE:-0}" != 1 ]]; then
  echo "== pre-flight: check_params.py (adopted thresholds)"
  $PY scripts/check_params.py || {
    echo "ABORT: a build used a parameter that disagrees with config/adopted_parameters.yaml." >&2
    echo "       Fix the build or the registry - do NOT run with a threshold in an unknown state." >&2
    echo "       SKIP_PARAM_GATE=1 overrides, and should be a deliberate, stated exception." >&2
    exit 1; }
fi

run() {  # run <label> <out...> -- <cmd...>
  local label="$1"; shift
  local outs=(); while [[ "$1" != "--" ]]; do outs+=(--out "$1"); shift; done; shift
  echo; echo "== $label  $(date +%H:%M:%S)"
  scripts/runlog.sh "${outs[@]}" --note "D186 formal end-to-end pipeline, ARACNe-only (D160/D186)" -- "$@" || {
    echo "ABORT: $label failed. See the ledger: python scripts/check_builds.py --runs" >&2; exit 1; }
}

# ---------------------------------------------------------------- disease nested arrays (D179)
# Reads step5/step6 only, so it is independent of everything on the drug side and goes first.
# EDGES=1 also writes the regulon edge table that 79e3 needs further down.
if want disease; then
run "50b contrib groups" "results/step6/aracne_project_v2_tf/regulon_contributions.tsv.gz" \
                  "data/processed/step5/regulons_aracne_project_v2_tf_edges.tsv.gz" \
                  -- env EDGES=1 STRATA=0 $R scripts/50b_regulon_contributions.R
run "50b contrib strata" "results/step6/aracne_project_v2_tf_strata/regulon_contributions.tsv.gz" \
                  -- env STRATA=1 $R scripts/50b_regulon_contributions.R
fi

# ---------------------------------------------------------------------------- drug-side core
# 79i2 and 79i3 edit the split tables IN PLACE and idempotently; 79t3:218 reads target_source, so
# skipping 79i2 makes 79t3 die with AttributeError (D164). Order here is load-bearing.
if want core; then
run "79e groups"  "$RES/drug_reversal_v2_tf.tsv"          -- env STRATA=0 $PY scripts/79e_drug_reversal.py
run "79e strata"  "$RES/stratum_reversal_v2_tf.tsv"       -- env STRATA=1 $PY scripts/79e_drug_reversal.py
run "79i split"   "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env STRATA=1 $PY scripts/79i_split_reversal.py
run "79i2 presence" "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env SFX=_v2_tf $PY scripts/79i2_target_presence.py
run "79i3 receptor lines" "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env SFX=_v2_tf $PY scripts/79i3_receptor_lines.py
run "79t3 lead"   "$RES/layered_score_v2_v2_tf.tsv" \
                  -- env ESSENTIAL=1 PROLIF_WEIGHT=1 SPLIT_SFX=_v2_tf $PY scripts/79t3_lead_score.py
# D210: evidence flags were run by hand outside this runner, and D208 found they had covered only the
# consensus rows. Placement is forced both ways: 79f3's label queue goes through fda_rows(), which sorts
# by layered_pct, so it must follow 79t3 (a fresh run crashed with KeyError 'layered_pct' when it sat
# before 79t3 - the hand-run order had only worked on stale columns); 79i4 applies the flags, so they
# must precede it. 79t3 reads no flag. STAGE=exposure follows STAGE=flags; 79f4 adds plate_confined.
run "79f3 flags"    "$RES/evidence_flags_v2.tsv" -- env STAGE=flags $PY scripts/79f3_evidence_flags.py
run "79f3 exposure" "$RES/evidence_flags_v2.tsv" -- env STAGE=exposure $PY scripts/79f3_evidence_flags.py
run "79f4 plate"    "$RES/evidence_flags_v2.tsv" -- $PY scripts/79f4_plate_confined.py
# 79i4 edits the split tables IN PLACE (tier_mech / tier_spec_mech / tier_lead) and refreshes the
# summary tables' evidence_status; it writes no table of its own.
# D210: it now runs BEFORE 79p. 79p copies evidence_status from the summaries, so running it first
# stamped the stale pre-79i4 status onto every FDA row (26 cohorts 'borrowed only', 15 'matched'
# where D122 gives 'candidates'), and a second pass gave a different file. 79i4 reads no 79p output.
run "79i4 tiers"  "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env SFXS=_v2_tf $PY scripts/79i4_tiers.py
run "79p fda"     "$RES/fda_ranking_v2_tf.tsv" "$RES/fda_ranking_strata_v2_tf.tsv" \
                  -- env SFX=_v2_tf $PY scripts/79p_fda_ranking.py
# D210: target role (layer 1 alone; D127 reports it) reads the split tables and the ARACNe-only
# overactivity (D208). It was run by hand outside this runner.
run "79t target role" "$RES/target_role" -- $PY scripts/79t_target_role.py
# D179 pooled nested array. Reads 79p's FDA rows, so it comes after 79p; it stops if its
# regulator terms no longer sum to the stored rev, which is the check that it describes THIS build.
run "79e2 contrib drug" "$RES/regulon_contributions_v2_tf" \
                  -- $PY scripts/79e2_regulon_contributions.py
fi

# ------------------------------------------------------------------------------- window arm
# window_reversal reads group_mr.tsv + the condition-level profiles directly, NOT 79e's output, so
# it does not depend on the core block - but 79j/79k annotate the split tables, so the whole group
# runs after core. D182: 79e3 mirrors window_reversal, NOT 79e; following 79e would silently
# produce numbers that do not match the window table.
if want windows; then
run "window_reversal" "results/systems_validation/integrated/window_reversal.tsv.gz" \
                  -- $PY scripts/systems_validation/integrate/window_reversal.py
# D211: the 14 subgroups get their own window run (own MRs and matched lines, parent network). 79j reads it.
run "window_reversal strata" "results/systems_validation/integrated/window_reversal_strata.tsv.gz" \
                  -- env STRATA=1 $PY scripts/systems_validation/integrate/window_reversal.py
run "79j window cols" "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env SPLIT_SFX=_v2_tf $PY scripts/79j_window_scores.py
run "79k window long" "$RES/drug_ranking_by_window_v2_tf.tsv" \
                  -- env SPLIT_SFX=_v2_tf $PY scripts/79k_window_long.py
# Two runs. The default one is the per-window nested array for the FDA rows; ALL=1 NO_TARGETS=1
# writes the regulator vectors over EVERY scored drug, which is what 79ad consumes. The target
# table at ALL=1 scale would be ~95M rows and no test reads it, hence NO_TARGETS.
run "79e3 window contrib" "$RES/regulon_contributions_v2_tf/regulator_contributions_by_window.tsv.gz" \
                  -- env SFX=_v2_tf $PY scripts/79e3_window_contributions.py
# D194: COSINE=1 adds, per (cohort, drug, window, regulator), the cosine between the gene-level
# DISEASE contributions (50b) and the gene-level DRUG contributions. Both sides already carry each
# edge's own tfmode x likelihood, so the per-gene weights are maintained end to end rather than
# rebuilt from z-scores. 18.5M rows; 79ah reduces it to what the page can embed.
run "79e3 window contrib ALL" "$RES/regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz" \
                  -- env SFX=_v2_tf ALL=1 COSINE=1 NO_TARGETS=1 $PY scripts/79e3_window_contributions.py
# D211: subgroup gene-level contributions per window, the subgroup input to 79ar. COSINE=1 is required:
# 79e3 only computes rev_gene_contrib inside the COSINE branch (the first D211 run crashed without it).
run "79e3 window contrib ALL strata" "$RES/regulon_contributions_v2_tf/regulator_contributions_by_window_all_strata.tsv.gz" \
                  -- env SFX=_v2_tf ALL=1 COSINE=1 NO_TARGETS=1 STRATA=1 $PY scripts/79e3_window_contributions.py
run "79ah cosine summary" "$RES/contribution_cosine_per_drug.tsv.gz" "$RES/contribution_cosine_summary.tsv" \
                  -- $PY scripts/79ah_cosine_summary.py
# D211: THE PRODUCTION RANKING. Filter on combined dependency / alteration / overactivity (layered_T >= 50),
# rank the eligible pairs by gene-level reversal (rev_gene), separately per timepoint. Needs 79t3, 79j's
# unknown reasons and both 79e3 ALL tables.
run "79ar gene ranking (D211)" "$RES/gene_ranking_v2_tf.tsv" "$RES/gene_ranking_chemo_v2_tf.tsv" \
                  -- env SFX=_v2_tf $PY scripts/79ar_gene_ranking.py
fi

# ---------------------------------------------------------------------------------- labels
# 79aa labels the WHOLE scored universe (41 groups + 14 subgroups), graded by clinical phase, and
# carries the window columns - so it must follow 79j. It only labels: no score is read, nothing
# is fitted, no row is filtered. It reuses the D181b ChEMBL id cache (79ab is not run here).
if want labels; then
run "79aa label pairs" "$RES/labels/pair_labels_v2_tf.tsv.gz" \
                  -- env SFX=_v2_tf $PY scripts/79aa_label_pairs.py
fi

# ----------------------------------------------------------------------- registered tests
# Both are PRE-REGISTERED (docs/SPEC_window_imputation.md, docs/SPEC_approved_profile_concordance.md)
# and are rerun so their numbers describe THIS build, not the previous one. D183 returned CLOSED
# and is expected to keep returning CLOSED; a change in that verdict is a finding, not a bug fix.
if want tests; then
run "79ac window imputation (D183)" "$RES/window_imputation" \
                  -- $PY scripts/79ac_window_imputation.py
run "79ad approved concordance (D184)" "$RES/approved_concordance" \
                  -- $PY scripts/79ad_approved_concordance.py
run "79ae moa-disjoint shortlist (D185)" "$RES/approved_concordance/moa_disjoint_shortlist.tsv" \
                  -- $PY scripts/79ae_moa_disjoint_shortlist.py
# D192: approval RECOVERY per (cohort, window), each window's AUC against the pooled score on the
# SAME rows. Needs 79aa's labels and 79j's window columns, so it sits here, after both.
run "79ag window accuracy (D192)" "$RES/window_accuracy.tsv" \
                  -- env SPLIT_SFX=_v2_tf $PY scripts/79ag_window_accuracy.py
# D210: D127 is the drug side's primary validation (D209: one vote per target unit is its headline)
# and was never in this runner. Needs 79t3 (layered_pct), 79f3 flags and 79t target role.
run "D127 biomarker (D127/D207/D209)" "$RES/d127" -- $PY scripts/d127/d127_biomarker.py
fi

# ------------------------------------------------------------------------------------ pages
# Off by default: they are presentation, not results, and 79x / 79q / 79z have not been audited
# for the _v2_tf arm. PAGES=1 to include them.
if want pages && [[ "$PAGES" == 1 ]]; then
run "79x2 window page" "$RES/drug_rankings_by_window_v2_tf.html" \
                  -- env SPLIT_SFX=_v2_tf $PY scripts/79x2_window_page.py
fi

echo; echo "== pipeline complete $(date +%H:%M:%S).  Audit:"
python scripts/check_builds.py

# D187: per-stage runtime, as measured, normalised per cohort and per drug so the numbers
# transfer to a different disease set rather than only describing this one.
echo; echo "== stage runtimes"
python scripts/stage_timings.py --write
