# SPEC — Can a drug's missing timepoint window be imputed? (D183)

**Status: PRE-REGISTRATION. Written and committed before any model is fitted or any accuracy number is
computed.** User, 2026-11-22: *"let's spec the pipeline of imputing the missing window thoroughly to test if
it can work."*

Everything in sections 3–9 is fixed by this document. If a choice here turns out to be wrong, the finding is
reported with the choice that was registered, and a second run under a changed rule is labelled as such.

---

## 1. The question, and the decision it informs

Most drugs are scored in only one timepoint window. 13,332 (cohort, drug) rows carry a single window — 5,823
`early_only`, 7,509 `late_only` — against 4,128 with both. Those single-window rows are currently reported as
what they are, and D180 labels the rest `unknown` rather than guessing.

**Question: for a drug measured in one window, can the other window's reversal score be predicted accurately
enough to be used in a ranking?**

Three outcomes, all acceptable in advance:

- **Yes** — imputation is adopted for the rows that meet the accepted conditions, with the error attached.
- **No** — imputation is closed with a number, and `unknown` stays the answer (D180). This is a real result.
- **Only under conditions** — e.g. it works for concordant drug classes but not discordant ones; those
  conditions must be stateable *without* reference to the held-out values.

## 2. Why this is pre-registered

Three recent tests in this project came back null or refuting after the fact: the RGES benchmark (D174), the
trial-outcome test (D175/D176, a clean null on outcome with a strong gradient on volume), and D169, whose own
hypothesis the check refuted. The temptation to tune is highest right after a null, and imputation is a method
with many knobs. The acceptance rule is therefore fixed here, before the first fit.

**Pre-committed:** a negative result is written up with the same prominence as a positive one, and no
threshold in section 8 is moved to recover a number.

## 3. Data

Unit of analysis: **(cohort, drug, window) reversal score** from
`results/systems_validation/integrated/window_reversal.tsv.gz` (D155/D171), the same table the ranking
consumes through 79j.

| | value |
|---|---|
| rows with BOTH windows (train + test pool) | **4,128** |
| …distinct drugs | **939** |
| …distinct cohorts | **14** |
| …distinct tissue keys | **8** |
| rows to impute if adopted (`early_only` / `late_only`) | 5,823 / 7,509 |
| median conditions behind one window of a both-row | **1** |
| both-rows with ≥ 2 conditions in BOTH windows | **1,224** |
| cohorts per drug in the both-set (median / max) | 4 / 12 |

**Two structural problems, named here because they bound every result that follows.**

**(a) The target is noisy.** The median both-row rests on **one condition per window**. A quantity measured
once cannot be predicted better than its own reproducibility. Reporting a model against a perfect-prediction
ceiling of 1.0 would therefore understate it, and reporting it against nothing would overstate it. Section 5
estimates the ceiling first.

**(b) The training set is not the application set.** The both-window pool spans **14 cohorts and 8 tissue
keys**; the rows that would be imputed span 29 cohorts. Any adopted model extrapolates to cohorts and tissues
it never saw. Section 7 tests that directly rather than assuming it.

## 4. What is held out, and how

The quantity predicted is the **target window's `rev`**, in its native units, for a row where the other
window is observed. Both directions are tested separately and must each pass on their own:

- **late → early** (predict `rev_early` from `rev_late` + features)
- **early → late** (predict `rev_late` from `rev_early` + features)

They are not pooled: the windows have different coverage and different condition counts, and a single pooled
number would hide a one-sided failure.

## 5. The noise ceiling — computed BEFORE any model

From the **1,224** both-rows with ≥ 2 conditions in both windows, split that row's conditions within a window
into two halves at random, score each half with the D155 formula, and correlate half A against half B across
rows. Spearman–Brown steps that up to full length. This gives `rel_early` and `rel_late`: how well the score
agrees with itself.

**The attainable ceiling for predicting window Y from window X is `sqrt(rel_Y)`** — no predictor can exceed
the square root of the target's reliability except by chance.

Every model and baseline result in section 8 is reported **both raw and as a fraction of that ceiling**. The
ceiling is computed once, from the data only, with no model involved, and is frozen before fitting.

**If `rel` is so low that the ceiling sits near the identity baseline, the test is inconclusive rather than
negative**, and that must be said plainly — it would mean the instrument cannot support the question, not
that imputation fails.

## 6. Baselines — fixed before fitting

A method must beat the **best** of these, not the weakest. Stated in advance so the comparison cannot drift:

| id | baseline |
|---|---|
| **B0** | the cohort's median `rev` in the target window |
| **B1** | **identity: predict the target window from the observed window unchanged.** This is the "early ≈ late" assumption the pooled score implicitly makes. D155 already measured it: Spearman(early, late) = **+0.388**, **37%** sign-discordant. |
| **B2** | ordinary least squares on the observed window alone (per direction, fitted on train) |
| **B3** | the drug's mean score in the target window across its OTHER cohorts, where it has any |

B1 is the bar that matters. The pooled score is only defensible if B1 is good, and a model is only worth
adopting if it is better.

## 7. Splits — grouped, never random

A drug appears in a median of 4 cohorts (max 12), and cohorts share cell lines and a disease NES vector. A
random row split would place the same drug on both sides and measure memorisation.

1. **Primary: grouped 5-fold by `pert_id`.** No drug appears in both train and test.
2. **Secondary: leave-one-tissue-key-out, 8 folds.** Tests the extrapolation in 3(b). A model that passes (1)
   and fails (2) is not adopted for cohorts outside the training tissues.
3. **Reported, not gated:** performance on the 1,224 high-condition rows vs the rest, since the target is
   less noisy there.

Hyperparameters are chosen inside the training folds only, by nested CV. Test folds are touched once.

## 8. Models, metrics, acceptance

**Models.** In increasing order of capacity; all are compared against the same baselines and splits.

- **M1** linear: target ~ observed window.
- **M2** linear + context: observed window, `n_conditions` in each window, `n_lines`, the cohort's `n_mrs`,
  and the drug's `line_match`.
- **M3** gradient-boosted trees on M2's features plus drug annotation (MOA class, number of annotated
  targets, anticancer flag) and the cohort's tissue key.
- **M4** the regulator-level predictor: use the per-regulator contributions from `79e3`
  (`regulator_contributions_by_window.tsv.gz`) — predict the target window's contribution vector and sum it,
  rather than predicting the scalar. Richer, and the one with a mechanistic story: if the early→late map is
  regulator-specific, this is where it would show.

**Metrics**, per direction and per split:

| metric | why |
|---|---|
| Spearman(pred, actual) | rank agreement; comparable to the D155 numbers and to the ceiling |
| **sign agreement** | **the decision-relevant one.** 37% of both-rows are sign-discordant; a predictor with good magnitude and wrong sign would move drugs the wrong way in a ranking |
| MAE / sd(actual) | error against the spread of the thing being predicted |
| calibration | if the model emits an uncertainty, does predicted error track actual error? Required for any conditional adoption |

**ACCEPTANCE — all four must hold, per direction, on the primary split:**

1. Spearman beats **B1 by ≥ 0.10** in absolute terms.
2. Spearman reaches **≥ 0.60 of the section-5 ceiling**.
3. **Sign agreement ≥ 0.75.**
4. The leave-one-tissue-out split (7.2) degrades Spearman by **< 0.10** against the primary split.

Pass on 1–3 but fail 4 → adoption restricted to the 8 training tissue keys, and said so.
Fail 1 or 3 → **imputation is closed.** The `unknown` label (D180) stands and the write-up says so.

## 9. What this does NOT test

- Whether a *continuous* time curve could be fitted. It could not: LINCS has effectively two usable
  timepoints (6 h → 255 dose-fittable cells, 24 h → 1,768), and only **196** drug × line pairs have ≥ 2
  dose-fittable timepoints. A two-point curve is a slope, which the window analysis already reports.
- The dose axis. That is `SPEC_dose_response_bmd.md` and is independent of this.
- Whether an imputed score improves any ranking. That is a separate, later test with its own registration;
  passing section 8 earns the right to ask it, nothing more.

## 10. Outputs

```
results/lincs/window_imputation/
  reliability.tsv      section 5: rel_early, rel_late, ceilings, n rows behind each
  baselines.tsv        B0-B3 per direction per split
  models.tsv           M1-M4 per direction per split, every metric
  acceptance.tsv       the four criteria, per direction, pass/fail as registered
  predictions.tsv.gz   held-out predictions, for inspection and for the error distribution
  run.json             seeds, fold assignment hash, package versions, input file mtimes
```

Script: `scripts/79ac_window_imputation.py`. Fixed seed, folds written out, so the run reproduces exactly.

## 11. Registration

Committed before the first fit. The git commit that adds this file is the registration; any later change to
sections 3–9 is a new section in this document with its own date and reason, never an edit in place.
