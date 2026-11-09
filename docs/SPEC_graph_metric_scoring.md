# Spec: graph-metric reversal scoring (D125) and the mechanism-consistency filter (D126)

Status: **drafted 2026-11-09, not run.** Two linked proposals from the scoring review of 2026-11-09. Decisions D125
and D126 in [`DECISIONS.md`](DECISIONS.md) once the rules below are fixed. Background: the methods review in
[`LITERATURE_SCORING.md`](LITERATURE_SCORING.md), the drug side in
[`PIPELINE_DRUG_SIDE_V2.md`](PIPELINE_DRUG_SIDE_V2.md).

## Question

`rev` compares a disease regulator vector `d` and a drug regulator vector `g` as a normalized signed dot product
over every tested regulator (D74). That treats regulators as independent. They are not: D91 found identity-linked
regulators running 7-30% of a network's tested regulators (BRCA 18%), D85 and the proliferation split (D80) show
shared axes dominating, and a correlated block of twenty regulators is currently counted twenty times.

Two questions follow.

**D125.** Does replacing the identity metric with a graph-derived one improve the score? That is,

    rev_L = - (d^T K g) / sqrt(d^T K d),   K = (L + lambda I)^-1

with `L` the graph Laplacian over the regulator set, in place of `rev = - (d^T g) / sqrt(d^T d)`.

**D126.** Can a signed signaling prior flag drugs whose observed regulator shift cannot be produced by their own
annotated target, as a reported flag beside `target_absent` and `receptor_mismatch`?

Why a graph metric rather than an empirical covariance: the pnet work (D86) and the latent-factor attempts found
that removing the dominant axes of tumor expression removed the regulators' own programs with them (92% of ESR1's
links run through them in breast). An empirical `Sigma^-1` estimated from the drug profiles risks the same thing,
because the correlation it would whiten is partly the signal. `L` makes the redundancy claim from topology instead,
which is a separate source.

## Why this is a one-parameter extension, not a new score

`K -> I / lambda` as `lambda -> inf`, and the normalization is scale-free, so `rev_L -> rev` in the limit. The
current score is nested inside the proposed one. Any loss is from variance at finite `lambda`, not from a different
claim, which is why `lambda` must be fixed on data independent of the adoption criteria (below).

## Part A: graph-metric reversal (D125)

### A1. Inputs

- Disease side: group and subgroup regulator NES, current production tables (Step 6 / 7 v2, `REGSET=tf`; CollecTRI).
- Drug side: `data/processed/lincs/drug_viper_specific_v2_tf/<NET>.f32`,
  `drug_viper_specific_collectri/ALL.f32`, unchanged.
- Regulon target sets, per regulon set and per network, as used to compute the NES.
- Optional secondary arm only: a signed prior (OmniPath or SIGNOR, release pinned in `config/datasets.yaml`).

### A2. Graph construction (primary arm: regulon overlap)

Nodes = the regulators tested in that network. Edge weight between regulators `i` and `j` = the overlap of their
target sets, as the Szymkiewicz-Simpson coefficient `|T_i ∩ T_j| / min(|T_i|, |T_j|)` (chosen over Jaccard because
regulon sizes differ by an order of magnitude and Jaccard would under-weight a small regulon nested in a large one).
Sign is ignored for the graph; the metric encodes redundancy, not direction.

Fixed before use: edges below 0.1 dropped; self-loops excluded; the graph is built per network and per regulon set,
then the two `K` matrices are applied to their own set's vectors before the D75 consensus average, so consensus
remains an average of two comparable scores.

Laplacian: symmetric normalized, `L_sym = I - D^-1/2 A D^-1/2`, eigenvalues in [0, 2], so `lambda` has the same
meaning across networks of different size and density.

This arm introduces **no new curated dependency**: the overlap graph is derived from the regulons already in
production, so the trust tier of the score is unchanged. That is the reason it is primary.

**Secondary arm (reported only):** the same construction with a signed signaling prior as the adjacency. This moves
the load-bearing step from the measured tier onto a curated one and is reported for comparison, never adopted by
this spec.

### A3. lambda, fixed before the adoption criteria are read

`lambda` in {0.1, 0.25, 0.5, 1, 2, 4, 8}, chosen as the value maximizing correct-direction rate on the **D84
perturbation benchmark** (LINCS shRNA, CRISPR and overexpression signatures, the 57-TF main analysis plus the
knockdown subset), then frozen. Rationale: that benchmark asks whether the regulator reading is right, which is the
quantity the metric is supposed to improve, and it shares no labels with the known-use, control or dominance
criteria used for adoption. If no `lambda` beats the `lambda = inf` limit on that benchmark by >= 1.0 points, stop
here and record D125 as not met without running the pipeline.

### A4. Downstream

`spec_L`, `rev_mech_L`, `spec_mech_L`, `spec_sib_mech_L` computed exactly as `79e` / `79i` compute their
counterparts, with `rev_L` substituted at the one place `rev` enters. Gates, flags, tiers and the D105b target
logic untouched, so the comparison is like-for-like.

Cost: `p` is 300-700 regulators per network, so `K` is one Cholesky per (network, regulon set), cached. Negligible
next to the existing VIPER runs.

### A5. Pre-registered decision rule

Arms: CUR (production) vs D125 at the frozen `lambda`. Criteria fixed before any result is read.

1. **Controls (lead criterion).** Anti-estrogens in Basal (negative control) median `spec_mech` percentile falls by
   >= 10 points vs CUR, and breast positive controls (anti-estrogens in LumA / LumB) do not fall. Same form and
   tolerance as D92 / D93.
2. **Dominance.** Drugs scored in >= 33 of 41 groups with `rev > 0` in >= 80% of them: rises by <= 0.5 points vs
   CUR. Standing tolerance from D93 / D94.
3. **Known use, non-inferiority only.** D118 metric (pooled within-cohort AUC, anticancer, gate-passing eligible
   rows): Hub and ChEMBL each fall by <= 0.01. Treated as a floor, not as evidence of benefit, because D90 and C28
   record that this yardstick cannot resolve changes of this size.
4. **Layer-3 refit (the informative test).** Refit the D121 within-cohort conditional logit with `spec_mech_L_pregate`
   in place of `spec_mech_pct_pregate`. D125 is informative if the fitted `v` exceeds 0.05 **and** drug-held-out
   within-cohort AUC rises by >= 0.01 on the fitting label set, with the other label set and the breast controls
   holding (the test D121 itself failed).

**Adopt** as the reported reversal score if 1, 2 and 3 are met. Changing the **lead** score requires criterion 4 as
well, and is a separate decision (D122 put the lead on `layered_pct`; this spec does not propose to move it).

The D84 benchmark is the `lambda` selection set and therefore **is not** an adoption criterion. Report its value,
do not score the decision on it.

### A6. Outputs

`results/lincs/d125/` (lambda sweep, benchmark table, per-arm metrics, decision), `*_v2_consensus_d125.tsv`,
`scripts/d125/`.

## Part B: mechanism-consistency filter (D126)

### B1. What it asks

For a drug with an annotated main target (D105b) and an observed regulator shift in the lines its score came from:
is there a signed path in a signaling prior from that target, at the direction the drug acts, to the regulators the
drug most strongly moved, whose net sign matches the observed direction?

Output per (drug, group): `mechanism_consistent_frac`, and the flag `mechanism_inconsistent`.

### B2. Primary method: signed path consistency (no solver)

For the top 25 regulators by `|activity|` in the drug's scoring-line consensus profile, find the shortest signed
path from the main target to that regulator in the prior, propagate the sign (drug action x product of edge signs),
and compare with the observed sign. `mechanism_consistent_frac` = consistent / reachable. Unreachable regulators are
excluded from the denominator and counted separately as `mechanism_coverage`.

Chosen as primary over CARNIVAL because the ILP would run over ~1,100 drugs x lines x 41 groups and needs a
commercial solver for interactome-scale instances; HiGHS or CBC at that size is not realistic.

**Secondary arm (reported only, one cohort):** CARNIVAL proper on breast, target as the perturbation node of known
sign, measurements = the drug's regulator NES, to check whether the path heuristic agrees with the full causal
reconstruction. If they disagree materially, the heuristic is not a stand-in and B must be rethought.

### B3. Pre-registered rule

1. **Known-failure control.** The flag must fire on anti-estrogens in MCF7, where the ARACNe reading of ESR1 is
   known to be wrong (C18, C25, D85). If it does not fire there, the filter does not work and D126 is closed.
2. **Degree confound.** Flag rate vs the main target's degree in the prior, Spearman `|rho| < 0.3`. At or above
   that, the flag is reading how well-studied the target is and is reported with that caveat, not used.
3. **Coverage floor.** `mechanism_coverage >= 0.4` for >= 60% of gate-passing FDA drugs, else the flag is too sparse
   to carry and is reported per drug only.

**Reported only, never ranks, gates or tiers**, following the D120 / D121 convention. Promotion to a gate needs its
own decision and its own controls.

Under R3 (LLM coordination notes): where the flag fires, the report must say the drug's observed regulator shift is
not explained by its annotated target in the prior used, and must name the prior and its release. This is a
statement about the prior, not proof the drug does not work.

### B4. Outputs

`results/lincs/d126/`, columns `mechanism_consistent_frac`, `mechanism_coverage`, `mechanism_inconsistent` joined
to the split tables; `scripts/d126/`.

## Prerequisite: the noise ceiling

Before either part, estimate the attainable ceiling: split the LINCS replicate-level signatures in half, run the
scoring path on each half, and compute the known-use AUC one half reaches against the other. Lim & Pavlidis 2021
(17% cross-phase compound recovery) and our own issue 2 in `LITERATURE_SCORING.md` imply this ceiling is low. It
decides how to read criterion 3, and whether a metric change can be expected to show up at all.

If the ceiling is near where `spec_mech` already scores, Part A can only be justified on controls and dominance,
and Part B becomes the more valuable half of this spec.

## Risks

- **Whitening removes signal.** The D86 / latent-factor prior. Mitigated by the nested `lambda -> inf` limit and the
  sweep; not eliminated. The overlap graph is the narrow claim (shared targets = partly redundant evidence), which
  is defensible in a way that whitening the full empirical covariance is not.
- **Tier inversion.** Part B introduces a curated signaling prior as a new dependency. Contained by keeping it
  flag-only and by the degree-confound check; the secondary arm of Part A is reported for the same reason.
- **Publication bias.** Both the prior and CollecTRI are biased to well-studied nodes. Criterion B2 tests this
  directly for Part B; Part A's primary arm avoids it by not using a curated graph.
- **Criterion 1 is a single drug class in a single tissue.** Anti-estrogens in breast are the only control where we
  know the right answer. Passing it is necessary, not sufficient.

## Open questions for the user

- Does criterion A5.4 (`v > 0.05`) read as the right bar, or should the layer-3 refit be reported rather than
  decision-bearing on its first run?
- Which signed prior for Part B, and does adding it need its own entry in `config/datasets.yaml` review?
- Should the noise ceiling run as its own numbered decision before D125, or inside it?
