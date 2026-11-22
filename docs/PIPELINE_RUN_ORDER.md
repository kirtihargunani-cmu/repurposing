# The formal run order, post ARACNe / VIPER (D186)

**This document and `scripts/run_pipeline.sh` are the same thing in two forms.** If they disagree, the
script is right and this file is stale — fix it.

User, 2026-11-22: *"i do not want to be working with results that haven't gone through the correct
filtering steps so let's formalize our pipeline with the decisions we have made so far and run it end to
end post the ARACNe/VIPER builds."*

---

## 1. What this fixes

Two defects, both found while writing it.

**(a) Orphan stages.** `run_rebuild_chain.sh` (D156) stopped after `79e → 79i → 79i2 → 79t3 → 79p` plus the
D179 nested arrays. Everything built after it — the window arm (`window_reversal` → `79j` → `79k` → `79e3`
→ `79x2`), the labelling arm (`79aa`), and the two registered tests (`79ac` D183, `79ad` D184) — was in **no
runner at all** and had only ever been run by hand. A result produced that way cannot be shown to have
consumed a current upstream table.

**(b) Two arms, silently diverged.** Every script written from D170 onward defaulted to `SFX=_v2_tf`, while
`79i2`, `79i3`, `79p`, `79t3` and `79i4` still defaulted to `_v2_consensus`. `run_v2_drug.sh` drove the
consensus arm, which **D160 had already retired from scoring**. So the windows, the labels and D179–D185
were all computed on an arm the repo never called production.

**D186 resolves (b): `_v2_tf` (ARACNe-only, per D160) IS production.** The five defaults are flipped. The
consensus arm stays reachable with `SFX=_v2_consensus` and keeps its legacy un-suffixed output names, so it
is a frozen reference rather than a deletion.

## 2. The boundary — what is upstream and is NOT rerun

"Post ARACNe/VIPER" means these are taken as given:

| stage | why it is upstream |
|---|---|
| `44_aracne_regulons.R` | the networks themselves |
| `50_viper.R` | disease VIPER → `group_mr.tsv`, patient, strata |
| `79d_drug_viper.R` | drug × cell-line **and** condition-level VIPER |
| celligner line weights, `79b` target gate, `79w` receptor gate, `79f3`/`79f4` flags, `79y` stand-ins | **by D156's own argument:** these depend on drug annotation, ChEMBL binding and Cmax caches, not on the regulons, so a regulon-arm rebuild does not invalidate them |

Also **not** run, by the 2026-11-22 scope decision: `79o` (PubMed), `79o2` (trials), `79f5` (Drugs@FDA
docs), `79ab` (ChEMBL name fetch). They hit external services or need manual reading passes. They reuse
their caches and stay manual, invoked separately.

## 3. The order, and what each edge depends on

Stages are grouped; `FROM=<group>` resumes at one.

### disease — the D179 nested arrays
| stage | reads | writes |
|---|---|---|
| `50b` groups, `EDGES=1` | step5 / step6 only | `regulon_contributions.tsv.gz`, **`regulons_aracne_project_v2_tf_edges.tsv.gz`** |
| `50b` strata | step5 / step6 only | strata `regulon_contributions.tsv.gz` |

Independent of the whole drug side, so it goes first. `EDGES=1` is **load-bearing**: the edge table it
writes is what `79e3` needs five stages later.

### core — the drug-side ranking
`79e` groups → `79e` strata → `79i` → **`79i2`** → `79i3` → `79t3` → `79p` → `79i4` → `79e2`

- `79i2` and `79i3` edit the split tables **in place and idempotently**. `79t3:218` reads `target_source`,
  so skipping `79i2` makes `79t3` die with `AttributeError` — this is D164, found the hard way.
- `79i4` also edits the split tables in place (tiers) and writes no table of its own.
- `79e2` reads `79p`'s FDA rows, so it is last. It **stops if its regulator terms no longer sum to the
  stored `rev`**, which is the check that it describes *this* build rather than the previous one.

### windows
`window_reversal` → `79j` → `79k` → `79e3` → `79e3 ALL=1`

- `window_reversal` reads `group_mr.tsv` and the condition profiles **directly, not `79e`'s output**, so it
  has no dependency on the core block. It is sequenced after it only because `79j`/`79k` annotate the split
  tables.
- **D182: `79e3` mirrors `window_reversal`, NOT `79e`.** Conditions rather than drug × line profiles, in-key
  lines only (no direct/borrowed split), `other` and unrecorded times dropped. Following `79e` would
  silently produce numbers that do not match the window table.
- `79e3` runs **twice**: the default run is the per-window nested array over the FDA rows; `ALL=1
  NO_TARGETS=1` writes regulator vectors over every scored drug, which is what `79ad` consumes. The target
  table at `ALL=1` scale would be ~95M rows and no test reads it.

### labels
`79aa` — labels the whole scored universe (41 groups + 14 subgroups) graded by clinical phase, and carries
the window columns, **so it must follow `79j`**. It only labels: no score is read, nothing is fitted, no row
is filtered. Reuses the D181b ChEMBL id cache.

### tests — both pre-registered
`79ac` (D183, `SPEC_window_imputation.md`) → `79ad` (D184, `SPEC_approved_profile_concordance.md`) → `79ae`
(D185 shortlist).

They are rerun so their numbers describe this build. **D183 returned CLOSED and is expected to keep
returning CLOSED; a change in that verdict is a finding to investigate, not a bug to fix.** `79ad` reads
`regulator_contributions_by_window_all.tsv.gz` (`79e3 ALL=1`), `pair_labels_v2_tf.tsv.gz` (`79aa`) and the
split tables, so it is genuinely last.

### pages — off by default
`PAGES=1` adds `79x2`. Pages are presentation, not results, and `79x`/`79q`/`79z` have **not** been audited
for the `_v2_tf` arm; they are deliberately not in the chain yet.

## 4. The pre-flight gate

The chain gates on `check_builds.py --strict`, because D162 was a downstream stage consuming an arm that
had never been rebuilt. As of 2026-11-22 the audit is **5 error, 1 warn, 11 checks**, all pre-existing and
all in manual-reading territory:

```
d108                   MIXED_AGE / NOT_REWRITTEN / UNRECORDED   (manual trial readings, 11-06..11-20)
backup_pre_d177_...    UNRECORDED
run 20261120T164652    CRASH: 79t3_lead_score.py exit=1          (historical, superseded)
```

None touches a table this chain consumes, so the D186 run was launched with `SKIP_GATE=1`. **That is an
accepted exception, not a clean gate**, and it should be cleaned up rather than left as the normal way to
run.

## 5. Thread budget

`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=1` is exported by the runner, not left to each script: D165 was
workers taking 4–8 BLAS threads each and breaching the "< half of 64" budget. `N_CORES` defaults to 6.
