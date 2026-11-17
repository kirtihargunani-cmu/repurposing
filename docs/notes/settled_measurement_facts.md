# Settled measurement facts — stop re-deriving these

Things this project has already measured, which keep being rediscovered and
re-litigated in new stages. **Read this before writing a spec or a stage brief.**
Each entry names the decision row that established it. If you think one is wrong,
check the row before re-running the measurement.

## 1. Imputed genes are accurate. "Not a landmark gene" is NOT a limitation.

**Established D129 (2026-11-16).** L1000 measures 978 landmark genes and infers the
rest. Held-out R² of an inferred gene on the 978 landmarks: **median 0.875**
(p10 0.834, p90 0.917). Control — each landmark predicted from the other 977:
**median 0.630**. Inferred genes are reconstructed *well*, and only ~10.8% fall
below the landmarks' own 10th percentile of replicate reproducibility.

**So do not flag "the target is not a landmark gene" as a limitation.** The pipeline
interpolates and extrapolates many quantities; this one is among the better-supported.
A spec that lists it as a caveat is adding noise, and reviewers will reasonably ask
why a 0.875 R² is being treated as missing data.

### The ONE narrow case where it genuinely fails — and why

Using **gene X's own imputed value to verify a perturbation of X**.

The imputation predicts X from its normal covariation with the landmarks. A targeted
knockout or knockdown of X **breaks exactly that covariation by construction** — X's
transcript is ablated while its usual correlates are not. The model has no channel
through which to see it.

Measured at D135-S0 (CRISPR, target's own z after knockout):

| target gene is… | median target z | QC pass (z ≤ −1) | genes |
|---|---|---|---|
| a **landmark** (measured) | **−3.90** | 81% | 28 |
| **inferred** | **−0.39** | 15% | 171 |

This is not evidence that imputation is poor. It is the one operation imputation
cannot support, because it asks the model to report a deviation from the very
relationship it was fitted on. Everywhere else, use inferred genes freely.

**Practical rule:** perturbation self-QC is available for the 48 landmark targets
(28 with CRISPR) and nobody else. For the rest use perturbation signature strength
or DepMap gene effect — not the target's own imputed transcript, and **not** an
exclusion.

## 2. Effective dimensionality — the real ceiling

**Established D130 / D130.1 / D131 / D133 (2026-11-16 to 17).**

- Whole BING gene space (10,112 genes): **n_eff ≈ 41**. All 978 landmarks: ≈ 57.
  A 905-gene footprint: ≈ 44. A random 478 genes: ≈ 39.
- Regulator activities within a cohort (~1,052 MRs): **n_eff ≈ 10**.
- **D133 established this is NOT regulon target overlap.** Overlap explains 16–26%;
  removing it moves n_eff only 9.7 → 12.7. Literally disjoint regulons still give
  10.8–14.6. CollecTRI overlaps *more* yet correlates *less*. Regulon-inference
  circularity accounts for a factor of ~2.3 inside a 30–60× collapse.
- Phrase it as **"this sum's variance behaves as if it were n_eff independent
  terms"**, never "the data hold only n_eff genes' worth of information". The
  correlation is real biology; that is why corrected nulls must *preserve* it.

## 3. Permutation nulls over genes were too narrow

**Established D130.** Both the gene-label and sign-shuffle nulls permuted genes as if
independent. Corrected, the D128 ESR1 control fails (z 2.81/2.92/1.52 → 1.09/0.83/0.61)
while **M itself is unchanged**. Naive nulls fire at ~18% where 5% is nominal;
corrected ones at 2–8%. This is a **calibration fix**, not a blanket deflation.
Do not rescale a published z by a rule of thumb — both shortcuts are *harsher* than
the measured null.

## 4. The reversal normaliser is wrong in scale but not in rank

**Established D131.** `rev = −Σwa / √(Σw²)` assumes independent regulators; the true
SD is `√(wᵀΣw)`, so `rev` is ~**10×** its claimed scale (median 10.21 over 55 cohorts,
range 4.10–15.46). But `norm` is one scalar **per cohort**, so every within-cohort
rank is invariant — verified 55 of 55 for `rev_pct`, `rev_z` and the tier mask.
**Shipped rankings and tiers are unaffected.** Cross-disease comparisons and absolute
readings of `rev` are not.

## 5. Data-shape facts that have caused repeated errors

- **The drug matrix is 73 cell lines**, not 87. `79c` keeps only TAS ≥ 0.4, taking
  `level5_bing_hiq_sigs.tsv.gz` from 90,988 signatures / 87 lines to 29,812 / 74.
  Any line or signature count taken from the pre-filter file is wrong for paired work.
  (Caused two errors, D134-S0 and D135-S0.)
- **Coverage is not scope.** A target "covered" by a perturbation arm is not a testable
  cell. CRISPR covers 208 targets but yields 2,350 *paired* (drug, gene, line) cells
  over 7 lines.
- **InChIKey first blocks do not detect analogs.** They hash connectivity, so two
  statins or two "-tinib" compounds differ. 307 compounds → 302 skeletons, 4 collapsing
  groups. Use Morgan/Tanimoto (RDKit 2026.09.2 is installed). Measured at D134:
  same-target pairs median Tanimoto **0.228** (p90 0.500) vs different-target **0.095**
  (p90 0.147); a 0.4 cut flags 17.8% of same-target pairs against 0.1% of others.
  Morgan measures scaffold, not pharmacophore — atorvastatin vs simvastatin is 0.075,
  so "analog-excluded" never means "same-class-excluded".

## 6. Things that are genuinely open

Not settled, do not cite as if they were: whether a drug's profile carries its own
target's signature (D135, not run); whether target identity is recoverable across
compounds (D134 S2, not run); whether M adds information over `layered_pct`
(needs an outcome label); the drug-side NES normalisation (D132, in progress).
