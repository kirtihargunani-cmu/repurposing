# Stage 0 — audit (mechanistic_incremental_v1)

Protocol `5cc2227` (spec Amendment 1 + `config/mechanistic_incremental_v1.yaml`). Run `20261123T163252-3315568`,
25 s. **Status: pass, no blockers.** No outcome of any experiment was computed in this stage.

**1. What did we need to establish?** That the inputs are the 2026-11-23 full run, that every table joins without
duplicating rows, and that the production scores decompose exactly into the per-regulator terms the experiments
will use. Any failure would have stopped the dependent stages.

**2. What was checked.**

| Check | Result |
|---|---|
| Pinned inputs carry the reference run's ids | 10 of 10 match (`inputs_manifest.json`, with sha256 of every file) |
| Upstream inputs the run consumed unchanged | 13 sources, 89 files hashed |
| D208 guard (`_v2_tf` default, consensus refused) | present in `fda_rows.py` and `79ar_gene_ranking.py` |
| Key uniqueness | all unique: window tables (17,460 and 2,402 rows), layered score (14,300), gene ranking (8,595), regulator terms (18,512,024 and 2,478,356) |
| Sum of regulator terms = stored reversal | every one of 21,588 cancer and 2,808 subgroup cells, max error 5.7e-5 |
| Same sum = timepoint table's reversal (stored to 5 dp) | all cells, max error 1.6e-5, none unmatched |
| Sum of gene-level terms = production `rev_gene` | 2,172 and 278 ranked cells, max relative error 2.1e-6 |
| D209 target-unit count | reproduced: 20 units, 16 correct |

**3. What the data offer the experiments** (`coverage.tsv`):

| Level | Window | Drug rows with window data | Rows with >= 2 conditions (split-half possible) | Ranked in production |
|---|---|---|---|---|
| Cancers (29 with a window run) | early | 9,951 | 2,293 | 547 |
| Cancers | late | 11,637 | 4,717 | 777 |
| Subgroups (8 with a window run) | early | 1,352 | 260 | 68 |
| Subgroups | late | 1,456 | 572 | 105 |

Eligible pairs with no timepoint data: 2,046 profiled but never in the cancer's matched lines, 1,336 in a unit with no
window run. They stay `unknown` and count as coverage, not as results.

**4. Found while auditing.** `79e3`'s split-half mode seeds each row with Python `hash()` of strings, which changes
between processes, so its half assignment is not reproducible run to run. This study seeds with `crc32(key) ^ seed`
and does not use that output. Worth fixing in production separately.

**5. What it means.** The inputs are the most recent run, the joins are safe, and the decomposition the experiments
rely on is exact. The split-half coverage (2,293 early and 4,717 late cancer rows) is enough for the reproducibility
arm of Experiment A in 14 early and 17 late cancers with >= 10 such drugs.
