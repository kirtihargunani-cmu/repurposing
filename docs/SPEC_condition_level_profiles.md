# SPEC — condition-level drug profiles: keep the dose–time structure

**Status:** design only. Nothing implemented, nothing in production changed.
**Proposed ID:** D153 · **Drafted:** 2026-11-19
**Raised by:** the user — *"I thought the pipeline uses the full dose response curve"* and
*"how can we keep the dose response curve per drug per gene"*.
**Read first:** `docs/SPEC_dose_response_bmd.md` §7, `docs/DECISIONS.md` D152.

## 1. What the pipeline does today

`scripts/79c_drug_profiles.py`, in its own words:

> *"One profile per drug × line = **mean of its active signatures (MODZ z-scores) over doses / times**"*

Concretely (`79c:46-52`): signatures are filtered to `is_hiq & TAS ≥ 0.4`, grouped by
`(pert_id, ModelID)`, and **averaged**. Dose and time are not carried into the output at all —
`drug_line_profiles_index.tsv` has three columns: `pert_id, ModelID, n_sig`.

Every drug-side quantity in the project is computed downstream of that average: regulator activity
(`79d`), reversal (`79e`), the split (`79i`), the receptor gate (`79w`), the knockdown benchmark
(`79n`), and the nested-cosine work of D150.

## 2. Why this matters, measured

**D152.** In MCF7 the ESR1 programme under anti-estrogens responds **−0.056 at 6 h** (correct) and
**+0.094 at 24 h** (inverted), flipping for all five drugs separately. Against the BRCA-LumA disease
vector the cosine is **−0.132 at 6 h** (opposes the disease, would rank high) and **+0.093 at 24 h**
(agrees, ranks low).

Two losses stack before any score is computed:

1. **The TAS floor deletes the correct timepoint.** TAS measures response *size*; early responses are
   small. Of 107 anti-estrogen 6 h signatures in MCF7, **1** clears TAS ≥ 0.4, against **72 of 264**
   at 24 h. The 6 h median TAS is 0.202 — genuine early responses below a threshold calibrated on
   late ones.
2. **The average then collapses whatever survived.** A drug that engages its target at 1 µM and one
   that moves nothing below 10 µM produce the same row.

**D136 §7** recorded the same defect from the dose side: **8,778 (compound, line) cells have ≥ 3
distinct doses** and the production pipeline uses none of that structure.

## 3. The proposed artifact

Change the row key from `(pert_id, ModelID)` to **`(pert_id, ModelID, pert_dose, pert_time)`**. Same
matrix layout, same gene columns, same float32 profile-major order.

```
drug_condition_profiles_raw.f32       conditions × genes
drug_condition_profiles_specific.f32  conditions × genes
drug_condition_profiles_index.tsv     pert_id, ModelID, pert_dose, pert_time, n_sig, median_tas
```

Sizing — storage was never the reason to collapse:

| variant | rows | matrix |
|---|---|---|
| current, one per drug × line | 9,678 | 0.39 GB |
| condition-level, `is_hiq & TAS ≥ 0.4` | 29,464 | **1.20 GB** |
| condition-level, `is_hiq & TAS ≥ 0.2` | 68,535 | 2.79 GB |
| condition-level, `is_hiq` only | 105,779 | 4.30 GB |

**Written alongside the existing files, never replacing them.** Production reads what it reads today;
the two can be compared on identical drugs.

## 4. The one piece that is not a mechanical change

`79c:55-59` projects out a per-line generic response: `g_l` is the mean profile over that line's
drugs, normalised, and removed from each profile. This is **load-bearing** — the script's own header
records that without it, *"on raw profiles VIPER reads any cytostatic drug as raising
[FOXA1 / GATA3 / XBP1], and anti-estrogens in MCF7 then score as anti-reversing in luminal breast"*.

The generic response at 6 h is **not** the generic response at 24 h: early responses are smaller and
differently shaped. A single `g_l` pooled over timepoints would subtract a predominantly-24 h vector
from 6 h profiles and manufacture a new artefact in the arm we are trying to fix.

**Rule: `g` is computed per `(ModelID, pert_time)`**, over that line-and-timepoint's own drugs, and
requires **≥ 30 drugs** in the stratum; conditions in a stratum below that floor are written to the
raw file, flagged `generic_removed = False`, and excluded from specific-arm analyses rather than
being projected against a poorly estimated vector. The per-stratum drug count is recorded in the
index.

## 5. Where aggregation moves to

Today: average profiles → score once. Proposed: **score each condition → summarise the scores.**

The summary is deliberately not fixed here; three are reported side by side so the choice is visible:
`mean over conditions` (reproduces today's behaviour and is the control), `most potent condition`,
and the **benchmark-dose / coherence** pair already built for O-3, which is the only one that can say
whether a drug's score is stable across its own dose series.

This also removes the need to special-case anti-estrogens: a drug whose score flips sign between 6 h
and 24 h becomes an **incoherent curve**, which is a reportable property, rather than an averaged
number that is quietly wrong.

## 6. Pre-stated acceptance tests

1. **Reproduction.** Averaging the condition-level specific profiles back over `(pert_id, ModelID)`
   with the `TAS ≥ 0.4` filter reproduces the current `drug_line_profiles_bing_specific.f32` to within
   float32 tolerance, **except** for the per-timepoint `g`. Any row differing beyond that is explained
   before proceeding. This is the gate: if we cannot reproduce today's artifact, the new one is wrong.
2. **Anti-estrogen direction.** With conditions kept, the ESR1 programme in MCF7 under anti-estrogens
   is negative at 6 h. Pre-stated from the D152 measurement; it is a confirmation that the plumbing
   preserves what Level 5 already shows, **not** independent evidence.
3. **No silent coverage loss.** The number of (drug, line) pairs representable is unchanged; only the
   number of rows grows. Pairs with a single condition must round-trip identically.
4. **Knockdown benchmark (`79n`) does not fall** by more than 2 points under the `mean over
   conditions` summary — the arm that should behave like today.

## 7. Scope and what this does not claim

- It does **not** change any score, threshold or ranking. It changes what is available to score.
- It does **not** fix the TAS floor. Keeping conditions is necessary but not sufficient: at
  `TAS ≥ 0.4` the 6 h anti-estrogen signatures are still absent. The floor is a separate
  pre-registration (D152), and this spec is what makes that one testable.
- The median (drug, line) cell has **one** dose and **one** timepoint, so the benefit is concentrated
  in the 8,778 cells with ≥ 3 doses and 2,736 with ≥ 2 timepoints, not spread evenly.
- Downstream consumers (`79d`, `79e`, `79i`, `79w`, `79n`, `79r`, `79y`, `d126/b1`) continue to read
  the existing files until a separate decision moves them.

## 8. Build order

1. `79c2_condition_profiles.py` writing the three artifacts, with `g` per `(ModelID, pert_time)`.
2. Acceptance test 1 (reproduction) — **stop if it fails**.
3. `79d` run per condition for a bounded cohort set, to produce regulator activity per condition.
4. Acceptance tests 2–4.
5. Only then: a separate proposal for which summary, if any, production should adopt.
