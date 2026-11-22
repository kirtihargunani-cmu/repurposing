# SPEC — Do our novel top-ranked drugs share the mechanism profile of the approved ones? (D184)

**Status: PRE-REGISTRATION. Written and committed before any concordance number is computed.**
User, 2026-11-22: *"can we run some sort of test on the already FDA approved drugs for diseases. For
example, their profile is favorable for these diseases and do our top ranked drugs (other than these drugs)
also exhibit this profile."*

---

## 1. The question, and what it is NOT

For a disease, the drugs already approved for it are the closest thing we have to a positive control. They
each have a **mechanism profile** — the per-regulator contribution vector from D179 (`79e2`), i.e. which of
the disease's regulators the drug moves and in which direction, summing exactly to `rev`.

**Question: do our top-ranked drugs that are NOT approved for this disease have profiles resembling the
approved ones, beyond what reversal sign alone explains?**

This is distinct from the two tests already in the project:

| test | asks |
|---|---|
| known-use AUC (79l, C28) | do approved drugs *rank high*? |
| D127 biomarker pairs | does the score rank a drug higher in the cohort carrying its biomarker? |
| **this (D184)** | **conditional on ranking high, does a novel drug look mechanistically like what works?** |

**It cannot validate the ranking on its own, and the write-up must not claim it does.** Under
[[discovery-is-the-goal]] a drug that looks exactly like fulvestrant is a me-too, not a discovery. High
concordance is therefore *not* unambiguously good and low concordance is *not* unambiguously bad.

**What makes it informative is section 6's MOA control:** a novel drug that resembles the approved profile
*without sharing their annotated targets* is mechanistic convergence by a different route. That is the
quantity worth having, and it is the headline of this test.

## 2. The confound, measured before designing (not a result)

A contribution is `-nes_i x activity_i / norm`, and `nes` is shared by every drug in a cohort. So drugs could
be similar by construction. Measured on BRCA-LumA (143 drugs x 854 regulators):

| quantity | value | consequence |
|---|---|---|
| pairwise cosine, mean / median | 0.027 / 0.020 | profiles are **not** all alike; the metric has range |
| mean cosine, same-sign `rev` pairs | **0.111** | **the confound is the SIGN of reversal** |
| mean cosine, opposite-sign pairs | **−0.066** | |
| Spearman(cosine, min abs `rev`) | **0.010** | the **magnitude** of reversal does NOT drive similarity |
| PC1 / PC2 of normalised profiles | 25.3% / 15.5% | one dominant shared axis, not swamping |

**Design consequences, fixed here:** match on the *sign* of `rev`, not its magnitude; and report every result
twice, with and without the leading shared component projected out.

## 3. Data and scope

- **Profiles:** `results/lincs/regulon_contributions_v2_tf/regulator_contributions.tsv.gz` (D179), one
  `rev_contrib` vector per (cohort, drug) over the cohort's tested regulators. Cohorts only; no subgroup is
  in scope (they are in the table, but the reference sets below are too small there).
- **Reference set A (per cohort):** drugs with `label_phase >= 3` for that cohort from
  `results/lincs/labels/pair_labels_v2_tf.tsv.gz` (D181), restricted to `rev > 0`.
  **40 of 41 cohorts have >= 3 such drugs; 37 have >= 5.** Cohorts with < 3 after the `rev > 0` restriction
  are dropped and listed.
  Secondary, stricter: `label_phase == 4` (22 cohorts have >= 3). Reported separately, never pooled.
- **Test set T (per cohort):** the top `TOP_PCT` = 10% by `layered_pct` among gate-passing rows, **excluding
  every drug with any label (`label_phase > 0`) for that cohort** and excluding A by construction.
- **Minimum sizes:** a cohort contributes only if `|A| >= 3` and `|T| >= 10` after exclusions.

## 4. The statistic

For drug `d` and cohort `c`, with profile vector `x_d` (L2-normalised over the cohort's regulators):

```
conc(d) = cosine( x_d , centroid(A \ {d}) )
```

Leave-one-out for members of A, so an approved drug is never compared against a centroid containing itself.

Reported per cohort and pooled: median `conc` for **T**, for **A** (the ceiling this could reach), and for
the null of section 5.

## 5. The null — matched on reversal sign, by construction

The comparison is **not** T against all other drugs. For each cohort, the null set **N** is drugs that are:

- **not** labelled for this cohort (`label_phase == 0`), and
- **not** in T (i.e. outside the top 10% by `layered_pct`), and
- **matched to T on the sign of `rev`** — N is restricted to `rev > 0`, as A and T are.

`conc` is computed identically for N. The effect is `median conc(T) - median conc(N)`, and significance comes
from **10,000 permutations that shuffle the T / N assignment within a cohort**, pooled across cohorts by the
median of per-cohort effects. Seed fixed, written to `run.json`.

This answers "do top-ranked novel drugs look like the approved ones *more than other reversing drugs do*",
which is the only version of the question that is not answered by reversal sign alone.

## 6. The MOA control — the part that makes the test worth running

Repeat the whole of sections 4–5 on the subset of T whose drugs share **no annotated target and no MOA class**
with any member of A (targets from `79i2` `evidence_targets`, MOA from `hub_moa`).

- A positive effect that **survives** this restriction is mechanistic convergence by a different route —
  the result worth reporting.
- A positive effect that **disappears** means the test rediscovered drug classes, which target annotation
  alone would have given without any of our scoring. That must be said plainly.

## 7. Shared-axis sensitivity

Every number in sections 4–6 is computed twice:

1. on the raw normalised profiles;
2. after projecting out the leading principal component of the cohort's profile matrix (25.3% of variance in
   the LumA probe), which is the generic proliferation-like axis `spec_mech` already removes elsewhere.

A conclusion counts only if it **holds in both**. If the effect exists only in (1), it is the shared axis and
must be reported as such, not as mechanism concordance.

## 8. Acceptance — fixed in advance

Pre-registered, per reference definition (`phase >= 3` primary, `phase == 4` secondary):

1. **Effect:** median per-cohort `conc(T) - conc(N)` > 0 with permutation **p < 0.01**.
2. **Breadth:** the effect is positive in **>= 60% of contributing cohorts** (not carried by two or three).
3. **Robustness:** holds under section 7's PC1 removal.
4. **The one that matters:** holds on section 6's MOA-disjoint subset, with its own p < 0.05.

Pass 1–3, fail 4 → reported as "the ranking recovers known drug classes", explicitly **not** as mechanistic
convergence.
Fail 1 or 2 → novel top picks do **not** share the approved profile; reported as a null, with the caveat from
section 1 that this does not by itself condemn the ranking.

**Pre-committed:** no threshold here is moved after seeing a number, `TOP_PCT` is not tuned, and a null is
written up as prominently as a pass.

## 9. Known limitations, stated now

- **`conc` is bounded by how coherent A itself is.** Section 4 reports median `conc` within A as the ceiling.
  If approved drugs for a disease do not resemble *each other*, the test is inconclusive for that cohort and
  it is dropped from the pooled statistic with the list reported.
- Profiles are the **pooled** (timepoint-`unknown`) ones from D179. A per-window version is possible from
  D182's `79e3` tables once the window arm's coverage is adequate; it is not in scope here.
- A is small (median 10 at `phase >= 3`, 4 at `phase == 4`), so per-cohort centroids are noisy; this is why
  the pooled permutation across cohorts, not any single cohort, carries the result.
- An orphaned earlier attempt exists at `results/lincs/regulon_profile/` (no script committed, shared /
  specific decomposition, its own negative control failing at the 92nd percentile). **Nothing here is built
  on it** and its numbers are not comparable.

## 10. Outputs

```
results/lincs/approved_concordance/
  per_drug.tsv       cohort, drug, set (A/T/N), conc, conc_pc1removed, moa_disjoint, layered_pct, rev
  per_cohort.tsv     |A|, |T|, |N|, median conc each, effect, ceiling (median conc within A)
  acceptance.tsv     the four criteria, per reference definition, pass/fail as registered
  run.json           seed, TOP_PCT, cohorts dropped and why, input mtimes
```

Script: `scripts/79ad_approved_concordance.py`.

## 11. Registration

Committed before the first concordance number. The git commit adding this file is the registration; later
changes to sections 3–8 go in a new dated section, never an edit in place.


---

## 12. AMENDMENT, 2026-11-22 — window-matched only; the pooled version is withdrawn

**User, after the D183 result:** *"there is no point in running D184 as pooled profiles because it dilutes all
of the activity and we have a lot of single timepoint values so the pooled versus single timepoint would be a
confounding factor in itself. unless we want to test the approved for that disease at only matching
timepoints."*

**This is correct and sections 3–8 are amended accordingly. Nothing has been run under the original
version.**

### Why the pooled design is withdrawn

A pooled profile (D179) is the average over whatever timepoints a drug happened to be measured at. That is
**not the same quantity for every drug**: for a 24 h-only drug the pooled profile *is* its 24 h profile,
while for a both-window drug it is an average of two states that D183 has now shown are largely independent
(~80% of a window's reliable variance is absent from the other). So "pooled" carries a hidden covariate —
the timepoint composition behind it — and that covariate is not balanced between the approved set A and the
test set T. Any concordance measured that way could be composition rather than mechanism.

### The amended design

Everything in sections 4–8 is unchanged **except that it runs inside a single window**. The unit becomes
**(cohort, window)**, with `window` in {early, late}:

- **Profiles:** `regulator_contributions_by_window.tsv.gz` (D182, `79e3`) instead of the pooled table.
- **A:** approved (`label_phase >= 3`) drugs for the cohort **with a profile in that window** and `rev > 0`.
- **T:** top `TOP_PCT` by `layered_pct`, unlabelled for the cohort, **with a profile in that window**.
- **N:** unlabelled, non-top-ranked, `rev > 0`, **in that window**.
- Early and late are reported **separately and never pooled**, for the same reason D183 kept its two
  directions apart.

### The coverage cost, stated before running

| | pooled (withdrawn) | window-matched |
|---|---|---|
| cohorts with >= 3 reference drugs | 40 | **10 early / 11 late** |
| (cohort, window) cells with >= 5 | — | **21 of 25** |
| null pool per cell (median) | — | 167 unlabelled reversing drugs |

So the amended test covers roughly a quarter of the diseases. That is the price of removing the confound and
it is paid knowingly. **Section 8's acceptance criterion (2) is restated for the smaller n: the effect must
be positive in >= 60% of contributing (cohort, window) cells, counted per window.** No other threshold
changes.

### Prerequisite

`79e3` currently writes profiles for FDA-ranking drugs only (2,416 rows). The null set N needs vectors for
unlabelled drugs too, so it must first be re-run with `ALL=1` over all window rows. Regulator-level output
only; the target-level table at that scale is not needed for this test.
