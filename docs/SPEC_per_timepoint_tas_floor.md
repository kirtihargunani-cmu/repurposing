# SPEC — per-timepoint activity floor

**Status:** design only. Nothing implemented, no threshold adopted, production unchanged.
**Proposed ID:** D156 · **Drafted:** 2026-11-20
**Raised by:** the user — *"did we derive this transcriptional activity score threshold of 0.4 from
somewhere? can we lower it for the 6 hr timepoint?"*
**Read first:** `docs/DECISIONS.md` D152, D155; `docs/SPEC_dose_response_bmd.md` §3 (the same
calibrate-the-instrument procedure, already used once).

## 1. Where 0.4 came from

Nowhere measurable. `docs/notes/drug_side_decisions_and_caveats.md`, 2026-10-09:

> *"A drug signature counts as active if it is high quality (`is_hiq`) and TAS >= 0.4"*

No derivation, no sweep, no sensitivity analysis, no citation. It is a round number fixed by
decision and applied uniformly ever since — the same situation as `MIN_TARGETS = 10` before D138.3
measured it.

## 2. Why a uniform floor is wrong here

**TAS (transcriptional activity score) measures response MAGNITUDE.** Early responses are genuinely
smaller — that is pharmacology, not poor quality. Reproducibility is handled separately by `is_hiq`,
which these signatures already pass.

Measured, `trt_cp`, `is_hiq` only:

| timepoint | signatures | median TAS | share ≥ 0.4 | quantile that 0.4 sits at |
|---|---|---|---|---|
| 3 h | 2,230 | 0.240 | 17.7% | 0.823 |
| **6 h** | **20,294** | **0.264** | **24.3%** | **0.757** |
| 24 h | 66,615 | 0.289 | 35.5% | 0.645 |
| 48 h | 1,537 | 0.372 | 47.7% | 0.523 |

TAS rises monotonically with time. A single absolute floor therefore demands that a 6 h response be
as large as a 24 h one to count as real, and deletes it when it is not.

**The cost, measured (D155 follow-up):** of 35,793 (drug, line) cells, **16,567 (46.3%) have a 6 h
signature**, but only **3,254 (9.1%)** have one surviving TAS ≥ 0.4. **13,313 cells have 6 h data
that exists and is discarded** — 17,392 of 22,836 six-hour signatures, 76%, with median TAS 0.219
and p90 0.350, i.e. clustered just under the line.

This is why 17 of 20 anti-estrogen rows in D155 are `late_only`, and why the one row that kept its
6 h data put fulvestrant at the 88.3rd percentile in BRCA-LumA against 60.2nd late.

## 3. The proposed rule

**Set the floor at the same QUANTILE of each timepoint's own TAS distribution**, rather than the
same absolute value. "Active" then means the same thing relative to what is achievable at that
timepoint.

Anchoring on the status quo — 0.4 sits at the 0.645 quantile of the 24 h distribution — gives:

| timepoint | floor | signatures kept |
|---|---|---|
| 3 h | 0.298 | 791 of 2,230 |
| **6 h** | **0.331** | **7,202 of 20,294** |
| 24 h | 0.400 (unchanged) | 23,635 of 66,615 |
| 48 h | 0.511 | 545 of 1,537 |

Note the rule *raises* the 48 h floor as well as lowering the early ones. It is not a device for
admitting more data; it is a device for making the criterion comparable across time.

## 4. Why this is calibration and not fitting

Identical in form to the dose-response `k` procedure (`SPEC_dose_response_bmd.md` §3), which this
project has already run once:

1. **The quantity calibrated on is label-blind.** TAS is computed from the signature alone — no
   target annotation, no regulon, no disease contrast, no expected direction. It describes the
   instrument, not the result.
2. **The PROCEDURE is pre-registered, not the value.** The anchor quantile is declared before any
   score is recomputed, and every headline is reported as a **profile across anchor quantiles**
   (0.55, 0.645, 0.75), with a conclusion counting only if it holds throughout.
3. **Choosing the anchor because it makes a control pass stays forbidden.** The anchor is fixed by
   matching the existing 24 h behaviour, which is a status-quo choice, not an outcome-driven one.

## 5. Pre-stated acceptance tests

Lowering a quality floor is precisely how noise gets admitted, so the guards are the point.

1. **Knockdown benchmark (`79n`) must not fall by more than 2 points.** This is the test that the
   admitted signatures carry signal rather than noise. **Failing this fails the proposal**, whatever
   the coverage gain.
2. **Anti-estrogen direction in MCF7 at 6 h is negative** — pre-stated from D152's direct Level-5
   measurement (−0.056 over the ESR1 regulon; cosine −0.132 against the LumA disease vector).
   Confirmation that the admitted data behaves as the raw matrix already shows, **not** independent
   evidence.
3. **Known-use recovery does not fall** by more than 0.01 median AUC on either label set (C28).
4. **The D127 biomarker pairs do not degrade** — the project's primary check per
   `approval-is-weak-validation`.
5. **Report what was admitted:** the TAS distribution of newly included signatures per timepoint,
   and how many drugs move from `late_only` to `both` in the D155 window table.

## 6. The honest risk

17,392 six-hour signatures sit below 0.4. **Some are genuine early responses; some are drugs that do
nothing.** The floor exists to separate those, and a lower one admits more of the second kind. Test 1
is what decides whether the trade is worth it, and it is stated as a hard fail rather than a
consideration.

A second risk worth naming: the newly admitted signatures are systematically *weaker*, so any score
built on them has a lower signal-to-noise ratio than one built on 24 h data. Window-level comparisons
(D155) must therefore keep reporting the condition counts, and a 6 h hit resting on one weak
signature is not equivalent to a 24 h hit resting on six.

## 7. Scope

- Changes **which signatures count as active**. Changes no scoring formula, no regulon, no gate.
- Touches `scripts/77`, `79b`, `79c`/`79c2` — everywhere `TAS_MIN` is applied.
- Implemented as a per-timepoint lookup with the current uniform 0.4 as the default, so production
  behaviour is unchanged until the anchor is adopted.
- **Does not** address cells with no 6 h signature at all (53.7%). That is a genuine absence; see the
  user's extrapolation question, deferred pending D134, since imputing early response from chemically
  similar drugs would bake in the similarity confound D134 was written to control.
