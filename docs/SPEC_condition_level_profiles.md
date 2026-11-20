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

---

## 9. Time windows: how conditions are grouped once kept (added 2026-11-19)

**Raised by the user:** a single continuous number for drug response is biologically wrong, because
half-life and pharmacodynamics mean the response is not stationary. Their first proposal was fixed
0–12 h / 12–24 h windows; their second was to locate the boundary **dynamically**, at the turning
point where a drug's regulon response reaches a global maximum.

### 9.1 What the data supports

Timepoints per (drug, line), production filter (TAS ≥ 0.4), 9,678 cells:

| distinct timepoints | cells | share |
|---|---|---|
| exactly 1 | 8,182 | **84.5%** |
| ≥ 2 | 1,496 | 15.5% |
| **≥ 3** | **267** | **2.8%** |
| ≥ 4 | 7 | 0.1% |

**The dynamic version cannot be primary.** Locating a turning point needs ≥ 3 timepoints, which
exists for 2.8% of cells. For the rest it would be interpolation between two points, not estimation.

A distinction worth keeping separate: the user wrote both *inflection in the dose-response curve* and
*global max in regulon response*. Those are different axes. **Dose** structure is rich — 8,778 cells
with ≥ 3 doses. **Time** structure is thin. If the turning point is in time, which is what half-life
implies, the dose series cannot locate it.

### 9.2 But the early/late split is well supported, and is how LINCS was designed

Timepoints are not scattered — **96% of signatures are 6 h or 24 h**, and among multi-timepoint cells
the dominant set is literally `(6, 24)` (1,025 cells), then `(6, 24, 48)` (231). There is essentially
no 12 h data in LINCS (360 signatures in total), so a 12 h boundary is not merely arbitrary, it is
unmeasurable. **The boundary is 6 h vs 24 h because that is the experiment that was run.**

Single-timepoint cells are also **not** uniformly late:

| | 6 h only | 24 h only | other |
|---|---|---|---|
| production filter | **1,818 (22.2%)** | 6,250 (76.4%) | 114 |
| all `is_hiq` | 11,177 (36.8%) | 19,027 (62.6%) | 188 |

A 6 h-only cell is a clean early observation, not a deficient one. The early share rising from 22% to
37% when the TAS floor is dropped is the D152 effect again: the floor removes early signatures for
being small.

So roughly **9,300 of 9,678 cells can be assigned to a window** — 1,818 early-only, 6,250 late-only,
1,256 with both.

### 9.3 The rule

1. **Every condition carries `window ∈ {early, late, other}`**, from `pert_time ≤ 6` / `≥ 24` /
   otherwise. Recorded in the index; nothing is averaged across windows.
2. **Primary analyses are WITHIN-CELL**: restricted to the 1,256 cells holding both an early and a
   late condition, where the same drug in the same line is compared against itself. This is where the
   anti-estrogen result already lives, and it is the only comparison free of selection effects.
3. **Across-cell comparisons are reported with the window as a stated covariate, never as a free
   comparison.** Early-only and late-only cells are *not* a random partition — which timepoint a drug
   was profiled at reflects the screening campaign it came from, so comparing an early-only drug
   against a late-only drug confounds the window with whatever drove that choice. Any table mixing
   them must carry the window column.
4. **The dynamic turning point is an arm, not the rule.** On the 267 cells with ≥ 3 timepoints,
   locate the global maximum of regulon response and ask whether it falls between 6 and 24 h. This
   uses the rare dense cells to *earn* the simple boundary applied to the sparse ones, rather than
   applying a dynamic rule the data cannot support. If the maximum consistently falls before 6 h,
   then 6 h is already a late window and the framing changes — that would be a finding, not a
   failure.

### 9.4 Limit on the known-answer control

The anti-estrogen case has 6 h and 24 h and nothing between. It is therefore sufficient to show the
response **reverses between 6 and 24 h**, and insufficient to say when. No claim about 9 h or 18 h is
supported, and none will be made.

---

## 10. Toxicity and proliferation must be removed PER TIMEPOINT (user, 2026-11-19)

The project removes generic growth-arrest / toxicity signal in **three** places. Keeping conditions
is pointless if any of them is still estimated on a pooled, timeless basis, because the magnitude of
growth arrest is itself a function of time: at 6 h a cytostatic drug has barely arrested anything,
at 24 h it has. A single correction applied to both windows **over-corrects early conditions and
under-corrects late ones**, and it does so in the direction that destroys exactly the early-window
signal D152 depends on.

| # | where | what is removed | currently estimated over | required |
|---|---|---|---|---|
| 1 | `79c` generic-response projection | the line's mean drug response `g_l` | all of a line's profiles, pooled over time | **per `(ModelID, pert_time)`** — done in `79c2`, 121 strata, ≥ 30-drug floor |
| 2 | `79h` `prolif_linked` flag, criterion **P2** | regulators whose activity tracks the proliferation score | all drug × line profiles, pooled | **per window** — a regulator may track growth arrest at 24 h and not at 6 h |
| 3 | `79i` proliferation adjustment `A_adj = A − β·π` | the growth-arrest component of every regulator's activity | OLS over all profiles, one β per (network, regulator) | **β per (network, regulator, window)** |

Item 1 is implemented. **Items 2 and 3 are not, and they are the larger effect**: item 1 removes a
direction from the profile, while item 3 removes a fitted component from every regulator's activity
and is what the mechanism score is actually built on.

### 10.1 Why a pooled β is wrong, concretely

`π` is the profile's proliferation score (HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT, 303 genes
on the LINCS profiles). Pooled over conditions, `π` is dominated by 24 h observations — 76.4% of
single-timepoint cells and 96% of all signatures are 6 h or 24 h with 24 h in the majority. So the
fitted `β` describes the late relationship between activity and growth arrest, and is then applied
to 6 h conditions where the relationship has a different slope.

For a drug whose mechanism signal is early and whose arrest signal is late — which is the
pharmacology the user raised, and the anti-estrogen case specifically — this subtracts a late-fitted
arrest component from an early profile that does not yet contain one. The residual is not
mechanism; it is mechanism minus an artefact.

### 10.2 The rule

1. **`π` is computed per condition**, never per drug × line.
2. **P2 (`79h`) is evaluated per window.** A regulator is proliferation-linked in a window if its
   activity tracks `π` within that window. The flag becomes `(network, regulator, window)`, and a
   regulator linked at 24 h but not 6 h is reported as such rather than flagged globally.
3. **The adjustment regresses within window**: `β(network, regulator, window)`, with the window's own
   `π`. Strata below a stated minimum condition count fall back to the pooled `β` and are **flagged**
   `beta_pooled = True`, never silently.
4. **`prolif_effect` stays on raw profiles**, as today, but is reported per window.

### 10.3 Pre-stated check

Fit `β` separately at 6 h and 24 h over the 1,256 cells holding both windows — the within-cell set of
§9.3, so the same drugs and lines appear on both sides. **If the two β distributions are
indistinguishable, the pooled adjustment was harmless and items 2–3 can stay as they are.** If they
differ materially, the pooled correction has been distorting every mechanism score in the project,
and the size of that distortion is the finding.

This is a cheap test and it is the one that decides whether items 2 and 3 need rebuilding at all. It
should run before any of it is implemented.
