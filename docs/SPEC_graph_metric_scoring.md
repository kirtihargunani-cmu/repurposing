# Spec: graph-metric reversal scoring (D125) and the mechanism-consistency filter (D126)

Status: **drafted 2026-11-09, amended the same day (user: layer-3 refit reported not decision-bearing on the
first run; CARNIVAL on breast; the noise ceiling feeds the metrics rather than sitting beside them); not run.** Two linked proposals from the scoring review of 2026-11-09. Decisions D125
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

## A0. The noise ceiling, and reporting against it (user 2026-11-09)

Runs first, and its output is used in the metrics rather than quoted beside them.

Lim & Pavlidis 2021 (17% cross-phase compound recovery) and issue 2 of `LITERATURE_SCORING.md` both imply `r` is low.

**Two estimators, both run (user 2026-11-09).** The original single estimator (random split-half of the
replicate-level signatures of every drug x line pair) does not survive a depth check: `cell_coverage.tsv` gives 126k
high-quality signatures over 48k scored compound x line pairs, a mean of **2.6 signatures per pair**. A pair needs
>= 6 replicates before two halves are each stable enough to correlate, so for most pairs the split is undefined, and
where it is defined on 2-3 signatures the halves measure the split procedure rather than the assay. That bites
hardest in the thin cohorts, which are exactly the ones per-cohort correction was introduced to protect. So:

**A0a. Paired-signature reliability (primary, covers every cohort).** For each drug x line pair with >= 2
independent high-quality level-5 signatures at different dose or time, score each signature separately through the
path and correlate the resulting `spec_mech_pct`. `r_c` = the per-cohort Spearman over gate-passing drugs. Material
already exists: `79c` collapses across these signatures today, so A0a needs one scoring pass per signature stratum
and **no profile or VIPER rebuild**. Pre-register which strata count as independent (dose and time fields of the
level-5 index) before running.

**A0b. Split-half on deep pairs (bound, partial coverage).** Restricted to drug x line pairs with >= 6 replicate-level
high-quality signatures, which exist mainly in the deep lines (MCF7 3,540 compounds, PC3 3,760, A375 3,565, A549
2,900, HA1E 3,154). This is the clean estimator, and it does require `79c` and `79d` to rerun per half — see the cost
note below — but only over the restricted subset.

**A0.0 gate, measured first.** Count the pairs and lines meeting A0b's >= 6 threshold from the level-4 index. If
fewer than 500 pairs qualify, or they fall in fewer than 3 lines, A0b is confined to the lines that do qualify and
reported as a line-level bound, not a cohort-level one.

**Reading the two together.** `r_A0b - r_A0a`, where both exist, is the dose-and-time share of the apparent noise.
Direction matters and runs against us: A0a conflates dose and time variation with noise, so it **under**estimates
`r`, which lowers `AUC_max` and therefore **inflates** attainment — it flatters the pipeline. A0b is the stricter
test. Rule fixed in advance: where A0b covers a cohort's scoring lines, its `r` sets the ceiling; elsewhere A0a's
value is used and flagged `r_optimistic`. If the two differ by >= 0.1 in `r`, say so wherever a ceiling-relative
number is reported.

**Reliability is per cohort, not pooled (user 2026-11-09).** `r_c` reported for each of the 41 groups and the 14
subgroups, with the median across cohorts as a summary only. Reason: cohorts differ by an order of magnitude in how
many LINCS-tested matched lines they have (D78: 14 of 41 `matched`), so reliability almost certainly varies with
that, and a single pooled `r` would over-correct the data-rich cohorts and under-correct the thin ones — the cohorts
whose results are already least trustworthy. Any per-cohort metric is corrected with that cohort's own `r_c`; a
pooled metric uses the median and says so.

**Split count for A0b, adaptive (user 2026-11-09).** Batches of 10 until the standard error of the median `r_c`
stabilises, floor 20, cap 100; report the final count and the SE with every `r`. Target: SE small relative to 0.05
in `r`, since `r` enters the metrics through `sqrt(r)`. If the cap is reached without stabilising, report the SE and
treat the ceiling as indicative only. A0a has no resampling loop and so no cap.

**Cost (measured 2026-11-09, corrected).** An earlier draft claimed the split could share the VIPER step. It cannot:
`79d` consumes `drug_line_profiles_bing_specific.f32` from `79c`, so two half-consensuses mean `79c` and `79d` rerun
per half — 2n full passes, not n cheap re-scorings. Timed at the real dimensions (48,274 drug x line profiles, 38
networks, median 990 scored regulators per network) the aREA core is ~2 min per pass, so 200 passes is ~7 core-hours,
under an hour on 8 cores, with ~1 TB read and ~1.5 TB written if intermediates are kept rather than streamed into
scoring. Those figures are an optimised-BLAS **floor**: R's `viper` carries per-regulator overhead not captured, so
scale up before planning. The cap therefore survives on compute; it was the replicate depth, not the timing, that
forced the two-estimator design.

**Use.** Convert AUC to Somers' D (`D = 2*AUC - 1`), which is the rank-correlation form, and report

    ceiling     AUC_max = (1 + sqrt(r)) / 2
    attainment  D_observed / sqrt(r)

`sqrt(r)` is the classical attenuation bound on the correlation any score can reach with truth when the score is
measured at reliability `r`. Attainment is the fraction of the attainable signal a scoring arm actually gets, and it
is the quantity criterion 3 is stated in. Worked example: at `r = 0.3`, `AUC_max` is about 0.77, so `spec_mech` at
D118's Hub 0.462 is attainment ~ -0.14 — i.e. below chance against an attainable 0.77, which would say the problem
is the score and not the ceiling. At `r = 0.05`, `AUC_max` is about 0.61 and the same number reads very differently.

**Assumptions, stated because they are load-bearing.** It treats the indication labels as noiseless (they are not;
repoDB and Hub coverage per cancer is thin, which depresses the observed AUC for reasons the split-half does not
see), and the AUC-to-D conversion assumes the attenuation bound carries over to a rank statistic on a binary
outcome, which is an approximation. So `AUC_max` is an upper bound on what the drug-side noise permits, not a
prediction of achievable performance. Report both raw and attainment-corrected numbers; never report attainment
alone. Per-cohort `r_c` from few drugs is itself noisy, so a cohort with fewer than 50 gate-passing drugs gets its
`r_c` reported with its CI and is corrected with the median instead, flagged `r_pooled`; this applies to both
estimators independently, so a cohort can be `r_pooled` on A0b and not on A0a.

**Consequence for this spec.** If `AUC_max` sits near where `spec_mech` already scores, Part A can only be justified
on controls and dominance, and Part B becomes the more valuable half. The estimate also gives the first honest
answer to how much of the pipeline's measured performance is capped by L1000 reproducibility, which is worth having
independently of D125.

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
3. **Known use, non-inferiority, ceiling-relative.** D118 metric (pooled within-cohort AUC, anticancer,
   gate-passing eligible rows), reported as **attainment** against the noise ceiling of section A0: Hub and ChEMBL
   attainment each fall by <= 0.02. Stated as a floor, not as evidence of benefit, because D90 and C28 record that
   this yardstick cannot resolve changes of this size; the ceiling is what makes "cannot resolve" a number rather
   than a judgement.
4. **Layer-3 refit — reported, not decision-bearing on this run (user 2026-11-09).** Refit the D121 within-cohort
   conditional logit with `spec_mech_L_pregate` in place of `spec_mech_pct_pregate`, and report the fitted `v`, the
   drug-held-out within-cohort AUC on both label sets, and the breast controls. Reason it is not decision-bearing:
   the fit would run on the same Hub labels whose D121 version failed to transfer to ChEMBL and to the breast
   controls, so one arm of evidence should not carry a change on its first run. Making it decision-bearing needs its
   own rule, with the bar (`v` and the held-out AUC gain) fixed in advance.

**Adopt** as the reported reversal score if 1, 2 and 3 are met. Criterion 4 is reported alongside and informs
whether a later decision on the **lead** score is worth writing; it does not affect adoption here (D122 put the lead
on `layered_pct`, and this spec does not propose to move it).

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

### B2. The prior (user 2026-11-09)

**OmniPath, filtered to signed and directed interactions**, with SIGNOR alone as a sensitivity arm.

Filter fixed before use: `consensus_direction` and `consensus_sign` true, and >= 2 supporting resources or >= 2
references. Where per-source signs disagree, the consensus sign is taken and the disagreement recorded per edge
(`sign_disputed`), so a later arm can drop disputed edges without a re-pull. **Pinning: cache the edge table, do not pin a version string (user 2026-11-09).** OmniPath is a service over moving
upstream resources, so a release label is a weaker guarantee than it is for a downloaded file. The filtered
signed-and-directed edge table is written to `data/raw/omnipath/` with a checksum and registered in
`config/datasets.yaml` alongside the OmniPath package version, the exact query parameters, and the retrieval date.
The analysis then depends on a file we hold, not on a service that answers differently next quarter. The registered
record states the **signed-and-directed subset size**, not OmniPath's headline interaction count. This is a stronger
guarantee than the other sources get, which is appropriate: it is the one curated dependency Part B introduces into
an otherwise measured-tier pipeline. The SIGNOR sensitivity arm is cached the same way.

Why not SIGNOR as primary: it is the better resource edge-for-edge (sign is the primary curated annotation, not a
derived one, and each entry is mechanism-level and PubMed-backed), but its coverage is the binding constraint here.
B3's floor is on `mechanism_coverage`, and signed *path* consistency needs the whole chain present, so a sparse
prior shrinks the denominator until the flag rests on a handful of regulators per drug. SIGNOR is contained in
OmniPath, so this is a filter setting rather than a different resource, and the sensitivity arm is the same code
path with a stricter filter — it buys the precision-vs-coverage comparison for one extra run.

What this costs: provenance is heterogeneous and aggregation compounds publication bias rather than correcting it,
because hubs accumulate edges across resources. Expect B3 criterion 2 to be harder to pass here than on SIGNOR, not
easier. The contingency below exists for that.

### B2a. Primary method: signed path consistency (no solver)

For the top 25 regulators by `|activity|` in the drug's scoring-line consensus profile, find the shortest signed
path from the main target to that regulator in the prior, propagate the sign (drug action x product of edge signs),
and compare with the observed sign. `mechanism_consistent_frac` = consistent / reachable. Unreachable regulators are
excluded from the denominator and counted separately as `mechanism_coverage`.

Chosen as primary over CARNIVAL because the ILP would run over ~1,100 drugs x lines x 41 groups and needs a
commercial solver for interactome-scale instances; HiGHS or CBC at that size is not realistic.

**Breast arm: CARNIVAL proper (user 2026-11-09).** Run on the breast cohorts, where the controls are and where the
known failure lives: target as the perturbation node of known sign (inhibitor -1, agonist +1), measurements = the
drug's regulator NES in the scoring lines, signed prior as the network. Scope fixed to keep the ILP tractable:
gate-passing FDA drugs with a curated ChEMBL mechanism target (D105b) in the four breast groups, one instance per
(drug, group). Solver: CPLEX or Gurobi under an academic licence if available, else the instance is reduced by
restricting the prior to the k-hop neighbourhood of the target (k = 3) before falling back to HiGHS; if neither
runs, record the arm as not executable rather than substituting the heuristic silently.

Read-out: agreement between CARNIVAL's per-drug consistency and `mechanism_consistent_frac` (Spearman, and
concordance on the flag). This is what licenses the heuristic outside breast. If they disagree materially
(Spearman < 0.5, or the flag disagrees on > 25% of drugs), the heuristic is not a stand-in, and either B runs on
CARNIVAL only where it is affordable or the flag is confined to breast.

### B2b. Degree-correction contingency (fixed now, applied only if B3.2 fails)

If flag rate correlates with the main target's degree at `|rho| >= 0.3`, tightening the filter toward SIGNOR is
**not** the remedy — curated resources are degree-biased too, so it would trade coverage for no gain on the
confound. Instead, score consistency against a degree-matched null: for each (target, regulator) pair, count the
signed paths of the same length from 100 nodes sampled to match the target's degree within +/- 20%, and report
`mechanism_consistent_frac` as the target's consistency minus the null mean, in null SD units. The flag then fires
on the corrected quantity.

Fixed in advance so the choice is not made after seeing which direction it moves the results. Applied only on a
B3.2 failure, and its own re-check is B3.2 re-run on the corrected score.

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

## Risks

- **Whitening removes signal.** The D86 / latent-factor prior. Mitigated by the nested `lambda -> inf` limit and the
  sweep; not eliminated. The overlap graph is the narrow claim (shared targets = partly redundant evidence), which
  is defensible in a way that whitening the full empirical covariance is not.
- **Tier inversion.** Part B introduces a curated signaling prior as a new dependency. Contained by keeping it
  flag-only and by the degree-confound check; the secondary arm of Part A is reported for the same reason.
- **Publication bias.** Both the prior and CollecTRI are biased to well-studied nodes. Criterion B2 tests this
  directly for Part B; Part A's primary arm avoids it by not using a curated graph.
- **A0a flatters the pipeline.** It conflates dose and time variation with assay noise, so `r` comes out low,
  `AUC_max` low, and attainment high. Wherever A0b does not cover a cohort, its ceiling is optimistic and flagged
  `r_optimistic`. Do not read an attainment number from A0a alone as evidence the score is close to its ceiling.
- **A0b may not reach cohort level at all.** It depends on the A0.0 count of pairs with >= 6 replicates, which is
  unmeasured; if those sit in a handful of deep lines, A0b is a line-level bound and most cohorts fall back to A0a.
- **Criterion 1 is a single drug class in a single tissue.** Anti-estrogens in breast are the only control where we
  know the right answer. Passing it is necessary, not sufficient.

## Open questions for the user

- Whether an academic CPLEX or Gurobi licence is available, since it decides whether the breast arm runs on the full
  prior or on the 3-hop reduction (B2).
- Which level-5 index fields define an independent stratum for A0a (dose only, time only, or both), fixed before
  the run.

Settled 2026-11-09: layer-3 refit reported, not decision-bearing (A5.4); CARNIVAL runs on breast (B2); the noise
ceiling runs first and feeds the metrics (A0); OmniPath signed-and-directed as the prior with SIGNOR as a
sensitivity arm, and the degree-correction contingency fixed in advance (B2, B2b); reliability per cohort with an
adaptive split count (A0); the OmniPath edge table cached with a checksum rather than pinned by version (B2); both
reliability estimators run, A0a primary and A0b as the stricter bound, after the replicate-depth check killed
split-half as a sole estimator (A0). The 100-split cap survives on compute (timed 2026-11-09).
