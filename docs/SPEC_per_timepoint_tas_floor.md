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

1. ~~**Knockdown benchmark (`79n`) must not fall by more than 2 points.**~~ **WITHDRAWN 2026-11-20 —
   this guard was invalid.** `79n` scores **genetic** perturbations (`trt_sh.cgs`, `trt_xpr`) and
   never reads a compound signature or the TAS floor, so it cannot respond to this change at all. It
   would have passed unchanged however much noise was admitted, and would have been reported as
   evidence the change was safe.

   **REPLACEMENT (hard fail): drug-specific signal in the admitted signatures.**
   `scripts/systems_validation/observability/tas_floor_signal_test.py`. For each 6 h signature of
   drug *d* in line *l*, correlate it with *d*'s own accepted 24 h profile in *l* (`r_self`) and with
   random other drugs in *l* (`r_other`). `delta = r_self − r_other` is the drug-specific signal. If
   admitted signatures are noise, delta → 0 for them. **Pass requires the admitted band to retain
   ≥ 50% of the accepted band's median delta and ≥ 90% of signatures with delta > 0.**
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

---

## 8. Measured: how far the 6 h floor can go (2026-11-20)

The replacement guard, run on 2,221 six-hour signatures that have a same-drug 24 h reference in the
same line. Accepted band (TAS ≥ 0.4) median `delta` = **0.364**; that is the 100% reference.

| TAS band | n | delta | % positive | vs accepted |
|---|---|---|---|---|
| 0.700+ | 63 | 0.445 | 90.5% | **122%** |
| 0.500–0.700 | 504 | 0.404 | 96.2% | 111% |
| 0.400–0.500 | 377 | 0.309 | 97.1% | 85% |
| **0.331–0.400** | 318 | 0.275 | 95.6% | **76%** |
| **0.284–0.331** | 249 | 0.231 | 91.6% | **63%** |
| **0.250–0.284** | 163 | 0.216 | 93.3% | **59%** |
| 0.200–0.250 | 250 | 0.168 | 88.0% | 46% |
| 0.150–0.200 | 187 | 0.133 | 87.7% | 37% |
| < 0.150 | 110 | 0.098 | 83.6% | 27% |

Cumulatively, a 6 h floor at X:

| floor | 6 h signatures kept | median delta | % positive |
|---|---|---|---|
| 0.400 (current) | 42.5% | 0.364 | 96.2% |
| 0.331 | 56.8% | 0.334 | 96.0% |
| 0.284 | 68.0% | 0.316 | 95.3% |
| **0.250** | **75.4%** | **0.303** | **95.1%** |
| 0.200 | 86.6% | 0.283 | 94.2% |
| 0.150 | 95.0% | 0.271 | 93.6% |

**Signal degrades gracefully and never collapses.** Under the replacement pass rule — admitted band
≥ 50% of accepted delta and ≥ 90% positive — bands down to **0.250** qualify (59%, 93.3%); 0.200–0.250
fails on both counts (46%, 88.0%).

**Revised proposal: a 6 h floor of 0.25**, which keeps 75.4% of 6 h signatures against 42.5% today.
That is lower than the quantile-anchored 0.331 of §3, and it is justified by measured signal rather
than by matching the 24 h quantile. The quantile rule remains the mechanism for the other timepoints,
where no equivalent signal measurement has been made.

**One thing the table also shows about the status quo:** TAS ≥ 0.7 carries delta 0.445, *122%* of the
accepted band's median, so the current 0.4 floor already mixes strong and marginal signatures. It was
never a clean separation.
