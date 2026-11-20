# Build provenance: logging crashes, rebuilds and overwrites

**Why this exists.** D162 (2026-11-20). `data/processed/lincs/drug_viper_raw_v2_tf/` held 28 files
dated 11-02 while its sibling arm `drug_viper_specific_v2_tf/` had been rebuilt on 11-20.
`scripts/79e_drug_reversal.py:219-220` reads **both** arms per network. Nothing in the repo recorded
that one arm had not been rebuilt; the only visible symptom was a file count that was off by two,
for an unrelated reason (two networks orphaned by D97). Had it not been chased, `79e` would have
crashed on 26 networks and silently written `rev_raw = NaN` for the other 12 -- and `rev_raw` is the
diagnostic D154/D157 read to show how much reversal survives without the proliferation subtraction.

The rule that follows: **check dates per arm, not per script.** A stage reported "complete" can have
an arm that never ran. This page is the machinery that enforces it instead of relying on memory.

## Run every build step through the wrapper

```
scripts/runlog.sh --out <dir|file> [--out ...] [--expect N] [--note "..."] -- <command>
```

Example -- the rebuild chain:

```
scripts/runlog.sh --out data/processed/lincs/drug_viper_raw_v2_tf --expect 76 \
  --note "D162: raw arm was 18 days stale" -- \
  env VERSION=v2 REGSET=tf REGULON=aracne_project PROFILE=raw N_CORES=6 \
      Rscript scripts/79d_drug_viper.R
```

It writes three things:

| what | where |
|---|---|
| append-only run ledger | `logs/BUILD_LOG.tsv` |
| the run's stdout+stderr | `logs/runs/<run_id>.log` (gitignored) |
| provenance for each `--out` | `<out>/.build_manifest.json` |

The ledger gets **two rows per run**, `start` and `end`. Append-only is deliberate: a `SIGKILL`
cannot rewrite a line, so *a start row with no end row is itself the crash record*. It captures the
env vars that actually change what a build produces (`VERSION`, `REGSET`, `REGULON`, `PROFILE`,
`STRATA`, `SPLIT_SFX`, `PKN_FILE`, `PROLIF_MODE`, ... see `ENVKEYS` in the wrapper), the git SHA and
whether the tree was dirty.

## The manifest answers "may I trust this directory?"

`.build_manifest.json` records the run id, finished-at, duration, command, env, git SHA, per-file
sizes and mtimes, and two fields that exist specifically because of D162:

- **`files_not_rewritten_by_this_run`** -- files that were present before the run and that the run
  did not touch. A non-empty list means a **partial rebuild**: older files are sitting in the
  directory under current-looking names. This is caught at write time.
- **`supersedes`** -- the previous manifest's run id, finish time and env. What this build replaced.

## Audit at any time

```
python scripts/check_builds.py           # everything
python scripts/check_builds.py --runs    # just crashes / killed runs
python scripts/check_builds.py --strict  # exit 1 on any ERROR (for a pre-flight gate)
```

| finding | meaning |
|---|---|
| `MIXED_AGE` | files in one output span > 24h. **The D162 signature.** A consumer reading all of them mixes builds. |
| `NOT_REWRITTEN` | the manifest says the latest run left pre-existing files untouched. Partial rebuild. |
| `COUNT` | file count != the run's `--expect`. |
| `UNRECORDED` | files changed after the manifest was written -- something wrote outside a logged run, so provenance is void. |
| `NO_MANIFEST` | registered output with no provenance at all. |
| `DIED` | ledger start row with no end row: killed, OOM, or the box rebooted. |
| `CRASH` | ledger end row with a non-zero exit code. |

`config/build_outputs.tsv` registers outputs whose freshness matters, with the file count a complete
build should leave. Anything carrying a manifest is audited whether registered or not.
`results/lincs` is deliberately **not** registered: it is an accumulating results directory, not a
single atomic build, so a wide mtime span there is normal rather than a staleness signal.

## Honest limits

- **It records, it does not prevent.** A build run outside the wrapper leaves no ledger row; the
  next `check_builds.py` flags it as `UNRECORDED` or `NO_MANIFEST`, which is detection after the
  fact, not a block. Run `--strict` before trusting a downstream stage.
- **Manifests backfilled during the D162 audit** (`run_id: retro-20261120`, on
  `drug_viper_specific_v2_tf` and `drug_viper_specific_collectri`) are **reconstructed from file
  mtimes**. They carry a note saying so. Their inputs were never verified; they assert when the
  files were written, not what produced them.
- **Content is not hashed**, only size and mtime. A file rewritten with identical bytes still reads
  as a fresh build, and a corrupt file of the right size reads as fine.
- `MAX_SPAN_H = 24` is a heuristic. A legitimately slow build that straddles a day boundary will
  trip `MIXED_AGE`; read the oldest/newest filenames in the message before acting.

## Related

- `config/builds.yaml` -- named immutable *disease-side* builds. Its principle, "never rename or
  rebuild a directory in place: add a new build instead", is the stronger guarantee. The drug-side
  VIPER directories **are** rebuilt in place, which is why they need this weaker, after-the-fact
  machinery.
- D162 in `docs/DECISIONS.md`.
