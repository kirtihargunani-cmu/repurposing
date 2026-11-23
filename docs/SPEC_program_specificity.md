# SPEC — Does a specificity-decomposed program score rank better? (D201)

**Status: PRE-REGISTRATION. Written and committed before any ranking metric is computed.**
User, 2026-11-22: *"let's test specificity-adjusted opposition"* and, earlier, *"we don't want to
remove all generic activity and proliferation as of yet"* — so this **decomposes** rather than
deletes.

---

## 1. Two premises measured BEFORE writing this, one of which killed the original design

The candidate was first specified as *opposition residualised on ‖drug activity‖*. That is now
**withdrawn**, because the premise is false:

| score | median \|ρ\| vs true ‖drug activity‖ | cells ≥ 0.5 |
|---|---|---|
| program-sum | **0.054** | 0 / 30 |
| `rev_pct` | 0.059 | 0 / 30 |

**The ranking is not magnitude-confounded**, so residualising on magnitude would have been a no-op
and the test would have measured nothing. Recorded because the original design would have passed
or failed for reasons unrelated to its hypothesis — the D184 first-run failure mode.

The surviving premise is different and is measured:

| | median ρ with the proliferation-program score | cells ≥ 0.5 |
|---|---|---|
| **total** program score | **+0.440** (max 0.932, TCGA-LIHC late) | 10 / 23 |
| **clean** programs only | **−0.069** | 3 / 23 |

Magnitude and proliferation are **different confounds**: a weakly active drug can still be
strongly proliferation-aligned. D196 found PC1 (proliferation) ranks approved cytotoxics well in
COAD/PRAD/LUSC; that is the axis with something to remove, and it is separable by partition.

## 2. The decomposition

Programs are partitioned using the **objective** annotations from D199 — no hand labelling, no
gene lists read, thresholds inherited from D37:

- **artefact**: `r2_plate ≥` that cohort's own plate chance `+ 0.10`, or `|r_length| ≥ 0.5`
- **proliferation**: `|r_proliferation| ≥ 0.5`
- **clean**: everything else

Per (cohort, drug, window), summing `opposition` over each part:

```
S_total  = sum over ALL programs
S_prolif = sum over proliferation programs
S_clean  = sum over clean programs          <- the specific score
S_art    = sum over artefact programs       <- reported, never scored on
specificity = S_clean / (|S_clean| + |S_prolif|)
```

Nothing is deleted: all four are carried, and `S_prolif` stays available as its own ranking.

## 3. What is compared, and against what

Approval recovery (AUC, `label_phase ≥ 3`) per (cohort, window), **always on identical rows**, for:

| id | score |
|---|---|
| **B0** | `rev_pct` — the production score, the bar that matters |
| B1 | window-resolved `rev_<window>_pct` |
| B2 | `S_total` |
| **T1** | **`S_clean`** — the hypothesis |
| T2 | `S_prolif` — reported; expected to do WELL where approved sets are cytotoxic |
| T3 | `specificity` |

## 4. Acceptance — fixed here

1. **`S_clean` beats `rev_pct` (B0) by ≥ 0.03 median AUC across cells**, and in **≥ 60%** of them.
2. The gain is **not** explained by proliferation: `S_prolif` must not beat `S_clean`.
3. Holds in both windows reported separately, never pooled.

Fail → the decomposition does not improve ranking, written up as prominently as a pass.

**Pre-committed:** no threshold moved after seeing a number; the partition rule is fixed in §2 and
not tuned; cells with < 10 rows or < 3 positives are dropped and listed.

## 5. The caveat that bounds any result here

Approval recovery is **weak validation** and D192 showed it is bimodal — 9 cohorts significantly
above chance, 9 significantly below (UCEC-EEC z = −7.7), because some approved sets are cytotoxic.
A score that *loses* AUC by shedding proliferation may be **more** useful under
[[discovery-is-the-goal]], not less. So:

- the D127 biomarker-pair check is reported alongside and is the primary biological check;
- a drop in AUC for `S_clean` accompanied by a rise for `S_prolif` is reported as
  **confound separation**, not as failure.

This is why criterion 1 is a ranking claim and criterion 2 is an attribution claim: passing 1
without 2 would mean we rediscovered proliferation.

## 6. Output

```
results/lincs/programs/specificity_scores.tsv.gz   per cohort, drug, window: the four sums + specificity
results/lincs/programs/specificity_auc.tsv         per cohort, window: AUC for B0-B2, T1-T3
results/lincs/programs/specificity_run.json        acceptance verdict as registered
```

Script: `scripts/79al_program_specificity.py`.

## 7. Registration

The git commit adding this file is the registration. Later changes to §2–§4 go in a new dated
section, never an edit in place.
