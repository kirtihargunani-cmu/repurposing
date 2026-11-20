#!/usr/bin/env bash
# Run-and-record wrapper. Wrap ANY build step so that crashes are recorded and every
# output directory carries provenance saying which run wrote it, with what inputs.
#
# Why this exists (D162, 2026-11-20): drug_viper_raw_v2_tf held 28 files dated 11-02
# while its sibling arm was rebuilt on 11-20. Nothing recorded that, 79e would have
# half-crashed and half-silently-NaN'd, and the only visible symptom was a file count
# that was off by two for an unrelated reason. Dates per ARM, not per script.
#
# Usage:
#   scripts/runlog.sh --out DIR_OR_FILE [--out ...] [--note "..."] [--expect N] -- cmd args...
#
#   --out     an output directory or file this run is responsible for. Repeatable.
#             After a successful run each DIR gets .build_manifest.json; for a FILE the
#             manifest lands beside it as <name>.build_manifest.json.
#   --expect  how many files the run should leave in the (first) --out dir. Recorded,
#             and checked by scripts/check_builds.py. Omit if it varies.
#   --note    free text, goes in the ledger and the manifest.
#
# Writes:
#   logs/BUILD_LOG.tsv      append-only ledger, two rows per run (start, end)
#   logs/runs/<run_id>.log  the run's combined stdout+stderr
#   <out>/.build_manifest.json  provenance, only on success
#
# Append-only is deliberate: a SIGKILL cannot rewrite a line, so a start row with no
# end row IS the crash record. check_builds.py reports those as `died`.
set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LEDGER="$ROOT/logs/BUILD_LOG.tsv"
RUNDIR="$ROOT/logs/runs"
mkdir -p "$RUNDIR"

OUTS=(); NOTE=""; EXPECT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out)    OUTS+=("$2"); shift 2 ;;
    --note)   NOTE="$2";    shift 2 ;;
    --expect) EXPECT="$2";  shift 2 ;;
    --)       shift; break ;;
    *) echo "runlog: unexpected argument '$1' (did you forget -- before the command?)" >&2; exit 2 ;;
  esac
done
[[ $# -gt 0 ]] || { echo "runlog: no command given" >&2; exit 2; }

RUN_ID="$(date +%Y%m%dT%H%M%S)-$$"
LOG="$RUNDIR/$RUN_ID.log"
GIT_SHA="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo nogit)"
GIT_DIRTY="$(git -C "$ROOT" status --porcelain 2>/dev/null | head -1 | grep -q . && echo dirty || echo clean)"
# name the run after the script file in the command; fall back to the executable
SCRIPT="$(basename "$1")"; for a in "$@"; do case "$a" in *.py|*.R|*.sh|*.Rscript) SCRIPT="$(basename "$a")";; esac; done
CMD="$*"

# the env vars that actually change what a build produces in this project
ENVKEYS=(VERSION REGSET REGULON PROFILE STRATA SFX SPLIT_SFX RUN_TAG ONLY_NETS ONLY_GROUPS
         NET_MAP MR_OVERRIDE EXTRA_RDS PKN_FILE OUT_SFX TF_SFX PROLIF_MODE BETA_TAS_MIN
         ESSENTIAL PROLIF_WEIGHT CONS N_CORES THREADS)
ENVSTR=""
for k in "${ENVKEYS[@]}"; do [[ -n "${!k:-}" ]] && ENVSTR="$ENVSTR${ENVSTR:+ }$k=${!k}"; done

tsv() { local IFS=$'\t'; printf '%s\n' "$*" >> "$LEDGER"; }
[[ -s "$LEDGER" ]] || tsv run_id event ts script status exit_code dur_s git_sha git_state pid outputs env expect note log

OUTSTR="$(IFS=,; echo "${OUTS[*]:-}")"
START_EPOCH=$(date +%s)
tsv "$RUN_ID" start "$(date -Iseconds)" "$SCRIPT" running "" "" "$GIT_SHA" "$GIT_DIRTY" "$$" "$OUTSTR" "$ENVSTR" "$EXPECT" "$NOTE" "logs/runs/$RUN_ID.log"

echo "runlog: run_id=$RUN_ID  script=$SCRIPT  env='$ENVSTR'  log=logs/runs/$RUN_ID.log" >&2

# snapshot pre-existing outputs so the manifest can say what it superseded
PRE="$(mktemp)"; : > "$PRE"
for o in "${OUTS[@]:-}"; do
  [[ -e "$o" ]] && find "$o" -maxdepth 1 -type f -printf '%p\t%s\t%T@\n' >> "$PRE" 2>/dev/null
done

"$@" > "$LOG" 2>&1
RC=$?
DUR=$(( $(date +%s) - START_EPOCH ))

if [[ $RC -eq 0 ]]; then
  STATUS=ok
  for o in "${OUTS[@]:-}"; do
    [[ -e "$o" ]] || continue
    RUN_ID="$RUN_ID" CMD="$CMD" SCRIPT="$SCRIPT" GIT_SHA="$GIT_SHA" GIT_DIRTY="$GIT_DIRTY" \
    ENVSTR="$ENVSTR" NOTE="$NOTE" EXPECT="$EXPECT" PRE="$PRE" DUR="$DUR" \
      python -I "$ROOT/scripts/_runlog_manifest.py" "$o" || echo "runlog: manifest failed for $o" >&2
  done
else
  STATUS=crash
  echo "runlog: CRASH rc=$RC after ${DUR}s -- last lines:" >&2
  tail -15 "$LOG" >&2
fi
rm -f "$PRE"

tsv "$RUN_ID" end "$(date -Iseconds)" "$SCRIPT" "$STATUS" "$RC" "$DUR" "$GIT_SHA" "$GIT_DIRTY" "$$" "$OUTSTR" "$ENVSTR" "$EXPECT" "$NOTE" "logs/runs/$RUN_ID.log"
exit $RC
