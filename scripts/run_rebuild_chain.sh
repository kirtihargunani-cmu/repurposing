#!/usr/bin/env bash
# D156 rebuild, ranking half: 79e -> 79i -> 79t3 -> 79p, ARACNe-only (D160).
#
# Every stage runs under scripts/runlog.sh, so a crash or a kill is recorded and each
# output carries a .build_manifest.json (D163). The chain gates on
# check_builds.py --strict first: D162 was a downstream stage consuming an arm that had
# never been rebuilt, and the whole point is not to start that way again.
#
# Lands beside production, not on top of it: 79t3:292 (D139) refuses to overwrite
# layered_score_v2.tsv when SPLIT_SFX is non-default, so this writes
# layered_score_v2_v2_tf.tsv and fda_ranking_v2_tf.tsv.
#
#   scripts/run_rebuild_chain.sh            # run it
#   SKIP_GATE=1 scripts/run_rebuild_chain.sh  # run with a known, accepted failing check
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
RES=results/lincs
PY=/home/kirti/miniforge3/envs/celligner/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VERSION=v2 REGSET=tf REGULON=aracne_project N_CORES=6

if [[ "${SKIP_GATE:-0}" != 1 ]]; then
  echo "== pre-flight: check_builds.py --strict"
  python scripts/check_builds.py --strict || {
    echo "ABORT: stale or crashed builds upstream. Fix, or re-run with SKIP_GATE=1 if accepted." >&2
    exit 1; }
fi

run() {  # run <label> <out...> -- <cmd...>
  local label="$1"; shift
  local outs=(); while [[ "$1" != "--" ]]; do outs+=(--out "$1"); shift; done; shift
  echo; echo "== $label  $(date +%H:%M:%S)"
  scripts/runlog.sh "${outs[@]}" --note "D156 rebuild chain, ARACNe-only (D160)" -- "$@" || {
    echo "ABORT: $label failed. See the ledger: python scripts/check_builds.py --runs" >&2; exit 1; }
}

run "79e groups"  "$RES/drug_reversal_v2_tf.tsv"          -- env STRATA=0 $PY scripts/79e_drug_reversal.py
run "79e strata"  "$RES/stratum_reversal_v2_tf.tsv"       -- env STRATA=1 $PY scripts/79e_drug_reversal.py
run "79i split"   "$RES/drug_reversal_split_v2_tf.tsv" "$RES/stratum_reversal_split_v2_tf.tsv" \
                  -- env STRATA=1 $PY scripts/79i_split_reversal.py
run "79t3 lead"   "$RES/layered_score_v2_v2_tf.tsv" \
                  -- env ESSENTIAL=1 PROLIF_WEIGHT=1 SPLIT_SFX=_v2_tf $PY scripts/79t3_lead_score.py
run "79p fda"     "$RES/fda_ranking_v2_tf.tsv" "$RES/fda_ranking_strata_v2_tf.tsv" \
                  -- env SFX=_v2_tf $PY scripts/79p_fda_ranking.py

echo; echo "== chain complete $(date +%H:%M:%S).  Audit:"
python scripts/check_builds.py
echo
echo "HOLD HERE. Nothing is promoted: this arm sits beside production as *_v2_tf."
echo "Next: D156 acceptance tests 3 (known-use AUC, 79l_known_use.py C28) and 4 (D127"
echo "biomarker pairs), against data/backup_pre_d156_rebuild_20261120. Attribution is"
echo "confounded - the TAS floor and the regulon change landed together - so a consensus"
echo "arm may be needed to separate them."
