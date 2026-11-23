# SPEC — Reliability-weighted gene covariance as a correction to `rev` (D205)

**Status: PRE-REGISTRATION. Written and committed before the score is evaluated against any
criterion.** User, 2026-11-22: *"yes run a complete clean test with the covariance."*

---

## 1. Why there is a fourth attempt, and what was wrong with the third

Three scoring attempts have returned null: D194 (contribution cosine), D201 (program score), D203
(`rev_gene`). **Two of those were implementation errors of mine, not tests of the idea**, and that
is the reason for this document rather than another quick pass.

The D203 error, found by decomposition after the user objected that gene-level information should
not perform worse than the regulator totals:

```
Σⱼ dⱼCⱼ  =  n·cov(d,C)  +  (Σⱼ dⱼ)(Σⱼ Cⱼ)/n
```

With a median regulon of **88 genes**, `rev_gene` keeps `rev`'s own term divided by 88 and
multiplies the covariance by 88. Measured: Spearman(`rev_gene`, covariance term) = **+0.809**.
`rev_gene` had largely *replaced* the regulator signal rather than adding to it. The user's
objection was correct.

**The isolated covariance term is the first positive signal in four attempts**, on an anchor that
uses no approval data and shares no inputs with the score:

| | Spearman with target role `T` |
|---|---|
| `rev` | **−0.124** |
| `rev_gene` | −0.015 |
| **covariance term alone** | **+0.059** (Δ +0.140, better in 21/25 cells) |

## 2. The score

Per regulator, the covariance term in reversal sign (already produced by `79e3 COSINE=1`):

```
covᵢ = rev_gene_contribᵢ − rev_contribᵢ / nᵢ        ( = −nᵢ·cov(dᵢ,Cᵢ)/‖s‖ )
```

The score keeps `rev` at **full weight** and adds the covariance as a correction weighted by how
reliable that regulon's gene vector actually is:

```
rev_cov = rev + Σᵢ wᵢ · covᵢ
```

**`wᵢ` is measured, not tuned.** From the D202 split-half table, `wᵢ` = the median `cos_half` over
drugs for that (cohort, window, regulator), clipped to [0, 1]; where a regulator has no split-half
estimate, `wᵢ` = that (cohort, window) cell's median. A regulon whose gene vector does not
reproduce contributes nothing; one measured from 10+ conditions (reliability 0.79) contributes
nearly fully.

**There is no λ and no threshold.** The only free choice — how much covariance to admit — is fixed
by measured reliability. An unweighted variant (`wᵢ ≡ 1`) is reported alongside purely to show
whether the weighting is doing work; it is **not** an alternative to be selected after the fact.

## 3. What it is judged on

**Primary anchors — independent of approval.** The D120 target-role components: `D` (CRISPR
dependency), `A` (alteration), `T` (combined). These use no approval data, and unlike
`layered_pct` they share **no inputs** with `rev` or the covariance — `rev` carries no target
annotation at all. `O` (over-activity) is reported but not gated: it is the component closest to
what `rev` already measures, and D203 found it moves the other way.

**Secondary, reported and NOT gated: approval recovery.** D192 measured it bimodal (9 cohorts
significantly above chance, 9 below) and D196 showed the bare proliferation axis beats the
mechanism score where approved sets are cytotoxic. A score that sheds proliferation can lose AUC
while being better. It is reported for completeness and cannot pass or fail this test.

## 4. Acceptance — fixed here

Per (cohort, window), ≥ 50 rows and ≥ 3 positives, **identical rows** for every score.

1. **Combined target role.** Spearman(`rev_cov`, `T`) exceeds Spearman(`rev`, `T`) by a median of
   **≥ 0.05** across cells, and in **≥ 60%** of them, sign-test **p < 0.05**.
2. **Dependency specifically.** The same for `D`, at **≥ 0.03** and ≥ 60%.
3. **Both windows.** The direction of (1) holds in early and late reported separately — a gain
   carried by one window alone is reported as that, not as a result.
4. **Not a proliferation proxy.** Spearman(`rev_cov`, proliferation-program score) must not exceed
   Spearman(`rev`, proliferation-program score) by more than **0.10**. A gain bought by aligning
   harder with proliferation is the confound returning (D196), not a win.

Fail any → the covariance correction is not adopted, written up as prominently as a pass.

**Pre-committed:** no threshold moves after a number is seen; `wᵢ` is the measured reliability and
is not re-specified; the unweighted variant is diagnostic only and cannot be promoted to the
result; cells below the size floor are dropped and listed.

## 5. What this does NOT claim

- Adoption. Passing earns a promotion proposal on its own terms (the D186 pattern), nothing more.
- That `rev_gene` or the program scores should be revisited — D203 and D201 stand as nulls.
- Anything about cross-cohort transfer (D197) or the window structure (D183/D191).

## 6. Outputs

```
results/lincs/interpretation/
  d205_rev_cov.tsv.gz        per cohort, drug, window: rev, cov, w-weighted cov, rev_cov
  d205_weights.tsv           per (cohort, window, regulator): the measured reliability wᵢ
  d205_acceptance.tsv        the four criteria, pass/fail as registered
  d205_run.json              counts, seeds, verdict
```

Script: `scripts/79ap_covariance_correction.py`.

## 7. Registration

The git commit adding this file is the registration. Later changes to §2–§4 go in a new dated
section, never an edit in place.
