# Spec: graph-metric reversal scoring (D125) and the mechanism-consistency filter (D126)

Status (2026-11-10): **Both parts are closed. Part A is not met (A1 early stop); Part B is closed by its own B5.1 control.** B0 built and registered the prior (D126-B0); B1 ran on all 55 cohorts and the flag fires on 0 of 22 anti-estrogen control rows, so D126 is closed (D126-B1 result); the CARNIVAL arm is recorded as not executable. Nothing from this spec entered production. Plain-language write-up of both parts, with the charts: [`REPORT_graph_metric_scoring.html`](REPORT_graph_metric_scoring.html) (2026-11-11). A0 ran 11-09 (ceiling r 0.873, AUC_max 0.967; D125-A0, D125-A0.1). **A1 ran 11-10 and fired the early stop below: best gain +0.17 points against the +1.0 bar, so D125 is recorded as not met and A2 never ran** (D125-A1, D125-A1.1, D125-A1 result). B0-B3 (D126) are untouched. Originally: **specified, not run.** Drafted and settled 2026-11-09 (scoring review; all decisions below taken by the
user in that session). Consolidated from the drafting thread into run order; this file is the implementation
contract. Decisions D125 and D126 in [`DECISIONS.md`](DECISIONS.md). Background: the methods review in
[`LITERATURE_SCORING.md`](LITERATURE_SCORING.md), the drug side in
[`PIPELINE_DRUG_SIDE_V2.md`](PIPELINE_DRUG_SIDE_V2.md), reporting rules R3 / R4 in
[`LLM_COORDINATION_NOTES.md`](LLM_COORDINATION_NOTES.md).

**Pre-registration discipline.** Every threshold, tolerance and selection rule below is fixed now and is not to be
tuned once results are visible. Where a decision depends on a measurement we do not yet have, the rule names the
measurement and the action for each outcome. Arms whose criteria fail are recorded as failed, not adjusted.

## Question

`rev` compares a disease regulator vector `d` and a drug regulator vector `g` as a normalized signed dot product
over every tested regulator (D74). That treats regulators as independent. They are not: D91 found identity-linked
regulators running 7-30% of a network's tested regulators (BRCA 18%), D85 and the proliferation split (D80) show
shared axes dominating, and a correlated block of twenty regulators is currently counted twenty times.

**D125.** Does replacing the identity metric with a graph-derived one improve the score?

    rev_L = - (d^T K g) / sqrt(d^T K d),   K = (L + lambda I)^-1

with `L` the symmetric normalized Laplacian over the regulator set, in place of `rev = - (d^T g) / sqrt(d^T d)`.

**D126.** Can a signed signaling prior flag drugs whose observed regulator shift cannot be produced by their own
annotated target, as a reported flag beside `target_absent` and `receptor_mismatch`?

Why a graph metric rather than an empirical covariance: the pnet work (D86) and the latent-factor attempts found
that removing the dominant axes of tumor expression removed the regulators' own programs with them (92% of ESR1's
links run through them in breast). An empirical `Sigma^-1` estimated from the drug profiles risks the same thing,
because the correlation it would whiten is partly the signal. `L` makes the redundancy claim from topology instead,
which is a separate source.

**Why D125 is a one-parameter extension, not a new score.** `K -> I / lambda` as `lambda -> inf`, and the
normalization is scale-free, so `rev_L -> rev` in the limit. The current score is nested inside the proposed one.
Any loss is from variance at finite `lambda`, not from a different claim, which is why `lambda` must be fixed on
data independent of the adoption criteria.

## Run order

    A0.0  eligible-pair count for A0b                     -> gates A0b's scope
    A0a   paired-signature reliability, all cohorts       -> r_c, every cohort
    A0b   split-half reliability, deep pairs              -> r_c, stricter, partial coverage
    A1    lambda sweep on the D84 perturbation benchmark  -> lambda, or early stop
    A2    full-pipeline D125 arm vs CUR                   -> criteria 1-3, criterion 4 reported
    B0    OmniPath pull, filter, cache, checksum          -> the prior
    B1    signed path consistency, all cohorts            -> mechanism_consistent_frac, flag
    B2    CARNIVAL, breast cohorts                        -> agreement read-out
    B3    criteria

A0 runs first because its output is used inside the metrics, and because it recontextualises the existing D118
numbers regardless of what happens to D125.

# A0. The noise ceiling, and reporting against it

Lim & Pavlidis 2021 (17% cross-phase compound recovery) and issue 2 of `LITERATURE_SCORING.md` both imply `r` is
low. A0's output is used in the metrics, not quoted beside them.

## A0.0 Eligible-pair gate (measured first)

Count, from the level-4 index, the drug x line pairs with >= 6 high-quality replicate-level signatures, and the
lines they fall in. If fewer than 500 pairs qualify, or they fall in fewer than 3 lines, A0b is confined to the
qualifying lines and reported as a line-level bound, not a cohort-level one.

Why the gate exists: `cell_coverage.tsv` gives 126k high-quality signatures over 48k scored compound x line pairs,
a mean of **2.6 signatures per pair**. A pair needs >= 6 replicates before two halves are each stable enough to
correlate, so a plain split-half is undefined for most pairs and, at 2-3 signatures, measures the split procedure
rather than the assay. That bites hardest in the thin cohorts, which per-cohort correction exists to protect. Hence
two estimators rather than one.

## A0a Paired-signature reliability (primary, every cohort)

For each drug x line pair with >= 2 independent high-quality level-5 signatures at different dose or time, score
each signature separately through the path and correlate the resulting `spec_mech_pct`. `r_c` = the per-cohort
Spearman over gate-passing drugs.

No profile or VIPER rebuild: `79c` collapses across these signatures today, so A0a needs one scoring pass per
signature stratum. **Open item, fix before running:** which level-5 index fields define an independent stratum
(dose only, time only, or both).

## A0b Split-half reliability (stricter bound, partial coverage)

Restricted to drug x line pairs passing A0.0. These sit mainly in the deep lines (MCF7 3,540 compounds, PC3 3,760,
A375 3,565, A549 2,900, HA1E 3,154). This estimator does require `79c` and `79d` to rerun per half, over the
restricted subset only.

Split count, adaptive: batches of 10 until the standard error of the median `r_c` stabilises, floor 20, cap 100;
report the final count and the SE with every `r`. Target: SE small relative to 0.05 in `r`, since `r` enters the
metrics through `sqrt(r)`. If the cap is reached without stabilising, report the SE and treat the ceiling as
indicative only. A0a has no resampling loop and so no cap.

## A0c Reading the two together

`r_A0b - r_A0a`, where both exist, is the dose-and-time share of the apparent noise.

Direction matters and runs against us: A0a conflates dose and time variation with noise, so it **under**estimates
`r`, which lowers `AUC_max` and therefore **inflates** attainment. A0a flatters the pipeline; A0b is the stricter
test.

Rule, fixed in advance: where A0b covers a cohort's scoring lines, its `r` sets the ceiling; elsewhere A0a's value
is used and flagged `r_optimistic`. If the two differ by >= 0.1 in `r`, say so wherever a ceiling-relative number
is reported.

## A0d Per-cohort reliability

`r_c` is reported for each of the 41 groups and the 14 subgroups, with the median across cohorts as a summary only.
Cohorts differ by an order of magnitude in how many LINCS-tested matched lines they have (D78: 14 of 41 `matched`),
so reliability almost certainly varies with that, and a single pooled `r` would over-correct the data-rich cohorts
and under-correct the thin ones — the cohorts whose results are already least trustworthy.

Any per-cohort metric is corrected with that cohort's own `r_c`; a pooled metric uses the median and says so. A
cohort with fewer than 50 gate-passing drugs gets its `r_c` reported with its CI and is corrected with the median
instead, flagged `r_pooled`. This applies to both estimators independently, so a cohort can be `r_pooled` on A0b
and not on A0a.

## A0e Use in the metrics

Convert AUC to Somers' D (`D = 2*AUC - 1`), the rank-correlation form, and report

    ceiling     AUC_max    = (1 + sqrt(r)) / 2
    attainment  D_observed / sqrt(r)

`sqrt(r)` is the classical attenuation bound on the correlation any score can reach with truth when the score is
measured at reliability `r`. Attainment is the fraction of attainable signal a scoring arm gets, and is the quantity
criterion 3 is stated in.

Worked example: at `r = 0.3`, `AUC_max` is about 0.77, so `spec_mech` at D118's Hub 0.462 is attainment about -0.14
— below chance against an attainable 0.77, which would say the problem is the score and not the ceiling. At
`r = 0.05`, `AUC_max` is about 0.61 and the same number reads very differently.

**Assumptions, stated because they are load-bearing.** It treats the indication labels as noiseless (they are not;
repoDB and Hub coverage per cancer is thin, which depresses the observed AUC for reasons the reliability estimate
does not see), and the AUC-to-D conversion assumes the attenuation bound carries over to a rank statistic on a
binary outcome, which is an approximation. `AUC_max` is an upper bound on what drug-side noise permits, not a
prediction of achievable performance. Report raw and attainment-corrected numbers together; never attainment alone.

**Consequence.** If `AUC_max` sits near where `spec_mech` already scores, Part A can only be justified on controls
and dominance, and Part B becomes the more valuable half.

## A0f Cost

**Correction (11-09, D125-A0), measured on disk:** production has 9,678 drug x line profiles (29,812 active level-5
signatures; `drug_line_profiles_index.tsv`) and 40 ARACNe v2 TF networks plus CollecTRI, not 48,274 profiles and 38
networks; the 126k / 48k figures below count all high-quality signatures before the TAS >= 0.4 filter that `79c`
applies. Only level-5 data are on disk. A0 rebuilds only the changed pairs with each line's production g_l and splices
their VIPER columns into the production matrices (D125-A0 (4)); the dry run of `scripts/d125/a0b_splithalf.py`
reports the measured cost. The paragraph below is kept as drafted.

`79d` consumes `drug_line_profiles_bing_specific.f32` from `79c`, so two half-consensuses mean both rerun per half:
2n full passes for A0b, not n cheap re-scorings. (An earlier draft wrongly claimed VIPER could be shared.)

Timed at the real dimensions (48,274 drug x line profiles, 38 networks, median 990 scored regulators per network)
the aREA core is ~2 min per pass, so 200 passes is ~7 core-hours, under an hour on 8 cores, with ~1 TB read and
~1.5 TB written if intermediates are kept rather than streamed into scoring. Those are an optimised-BLAS **floor**:
R's `viper` carries per-regulator overhead not captured, so scale up before planning. The cap survives on compute;
replicate depth, not timing, forced the two-estimator design.

## A0g Outputs

`results/lincs/d125/a0/`: `reliability_a0a.tsv`, `reliability_a0b.tsv` (per cohort: `r_c`, CI, n drugs, n splits,
SE, `r_pooled`, `r_optimistic`), `eligible_pairs.tsv` (A0.0), `ceiling.tsv` (`AUC_max`, attainment per reported
metric). `scripts/d125/`.

# Part A: graph-metric reversal (D125)

## A1 Inputs

- Disease side: group and subgroup regulator NES, current production tables (Step 6 / 7 v2, `REGSET=tf`; CollecTRI).
- Drug side: `data/processed/lincs/drug_viper_specific_v2_tf/<NET>.f32`,
  `drug_viper_specific_collectri/ALL.f32`, unchanged.
- Regulon target sets, per regulon set and per network, as used to compute the NES.

## A2 Graph construction (primary arm: regulon overlap)

Nodes = the regulators tested in that network. Edge weight between regulators `i` and `j` = the overlap of their
target sets, as the Szymkiewicz-Simpson coefficient `|T_i n T_j| / min(|T_i|, |T_j|)`, chosen over Jaccard because
regulon sizes differ by an order of magnitude and Jaccard would under-weight a small regulon nested in a large one.
Sign is ignored: the metric encodes redundancy, not direction.

Fixed before use: edges below 0.1 dropped; no self-loops; the graph is built per network and per regulon set, and
each set's `K` is applied to its own vectors before the D75 consensus average, so consensus remains an average of
two comparable scores.

Laplacian: symmetric normalized, `L_sym = I - D^-1/2 A D^-1/2`, eigenvalues in [0, 2], so `lambda` means the same
thing across networks of different size and density.

This arm introduces **no new curated dependency** — the overlap graph comes from regulons already in production, so
the trust tier of the score is unchanged. That is why it is primary.

**Secondary arm (reported only, never adopted by this spec):** the same construction with the Part B signed prior as
the adjacency. It moves the load-bearing step onto a curated tier; reported for comparison only.

## A3 lambda, fixed before the adoption criteria are read

`lambda` in {0.1, 0.25, 0.5, 1, 2, 4, 8}, chosen as the value maximizing correct-direction rate on the **D84
perturbation benchmark** (LINCS shRNA, CRISPR and overexpression signatures: the 57-TF main analysis plus the
knockdown subset), then frozen.

Rationale: that benchmark asks whether the regulator reading is right, which is what the metric is supposed to
improve, and it shares no labels with the known-use, control or dominance criteria used for adoption.

**Early stop:** if no `lambda` beats the `lambda = inf` limit on that benchmark by >= 1.0 points, record D125 as not
met and do not run the pipeline arm.

> **Result 2026-11-10 — the early stop fired.** Consensus read, D84 main analysis, 1,774 signatures over 329 TFs:
> `lambda = inf` 0.6077; lambda 0.1 0.5727, 0.25 0.5958, 0.5 0.6037, 1 0.6048, 2 0.6094, 4 0.6077, 8 0.6065. Best
> **+0.17 points at lambda = 2**, bar +1.0. Monotone in the wrong direction; overexpression loses 10.3 points at
> lambda = 0.1 and never recovers the limit. CollecTRI's optimum is the limit itself. **A2 does not run.**
> Implementation choices and the identity amendment: D125-A1 and D125-A1.1. Outputs `results/lincs/d125/a1/`.

The D84 benchmark is the selection set and therefore **is not** an adoption criterion. Report its value; do not
score the decision on it.

## A4 Downstream

`spec_L`, `rev_mech_L`, `spec_mech_L`, `spec_sib_mech_L` computed exactly as `79e` / `79i` compute their
counterparts, with `rev_L` substituted at the one place `rev` enters. Gates, flags, tiers and the D105b target logic
untouched, so the comparison is like-for-like.

Cost: `p` is 300-1,200 regulators per network, so `K` is one Cholesky per (network, regulon set), cached. Negligible
beside the existing VIPER runs.

## A5 Pre-registered decision rule

Arms: CUR (production) vs D125 at the frozen `lambda`.

**Amended 11-09 by D125.1:** criterion 3 is reported only and no longer gates adoption; criterion 1 is widened to the
D127 biomarker-expectation test and remains the lead criterion. Criteria 2 and 4 are unchanged. The text below is as
drafted, with the amendment noted in place.

1. **Controls (lead criterion).** *(D125.1: widened.)* The **D127** pass rule is met on its primary set (49 instances,
   7 biomarker pairs: share >= 70% in the expected direction, sign test p < 0.05, median gap >= 10 percentile points),
   **and** the breast control does not regress: anti-estrogens in Basal (negative control) median `spec_mech` percentile
   falls by >= 10 points vs CUR, and breast positive controls (anti-estrogens in LumA / LumB) do not fall. Same form and
   tolerance as D92 / D93.
2. **Dominance.** Drugs scored in >= 33 of 41 groups with `rev > 0` in >= 80% of them: rises by <= 0.5 points vs
   CUR. Standing tolerance from D93 / D94.
3. **Known use — reported only (D125.1); it does not gate adoption.** D118 metric (pooled within-cohort AUC, anticancer,
   gate-passing eligible rows), reported as attainment against A0's ceiling: Hub and ChEMBL attainment each fall by
   <= 0.02. A floor, not evidence of benefit, because D90 and C28 record that this yardstick cannot resolve changes
   of this size; the ceiling is what makes "cannot resolve" a number rather than a judgement.
4. **Layer-3 refit — reported, not decision-bearing on this run.** Refit the D121 within-cohort conditional logit
   with `spec_mech_L_pregate` in place of `spec_mech_pct_pregate`; report the fitted `v`, drug-held-out
   within-cohort AUC on both label sets, and the breast controls. Not decision-bearing because the fit would run on
   the same Hub labels whose D121 version failed to transfer to ChEMBL and to the breast controls, so one arm of
   evidence should not carry a change on its first run. Promotion needs its own rule with the bar fixed in advance.

**Adopt** as the reported reversal score if 1 and 2 are met (D125.1: criterion 3 is reported, not required). Criterion 4 informs whether a later decision on the
**lead** score is worth writing; it does not affect adoption here (D122 put the lead on `layered_pct`, and this spec
does not propose to move it).

## A6 Outputs

`results/lincs/d125/`: lambda sweep, benchmark table, per-arm metrics, decision; `*_v2_consensus_d125.tsv`;
`scripts/d125/`.

# Part B: mechanism-consistency filter (D126)

## B1 What it asks

For a drug with an annotated main target (D105b) and an observed regulator shift in the lines its score came from:
is there a signed path in a signaling prior from that target, at the direction the drug acts, to the regulators the
drug most strongly moved, whose net sign matches the observed direction?

Output per (drug, group): `mechanism_consistent_frac`, `mechanism_coverage`, and the flag `mechanism_inconsistent`.

## B2 The prior

**OmniPath, filtered to signed and directed interactions**, with SIGNOR alone as a sensitivity arm.

Filter fixed before use: `consensus_direction` and `consensus_sign` true, and >= 2 supporting resources or >= 2
references. *(Corrected 2026-11-10, D126-B0: the service has no `consensus_sign` field. It returns
`consensus_stimulation` and `consensus_inhibition`, so the condition is implemented as exactly one of those two
being true. The intended filter is unchanged.)* Where per-source signs disagree, the consensus sign is taken and the disagreement recorded per edge
(`sign_disputed`), so a later arm can drop disputed edges without a re-pull.

**Pinning: cache the edge table, do not pin a version string.** OmniPath is a service over moving upstream
resources, so a release label is a weaker guarantee than for a downloaded file. The filtered signed-and-directed
edge table is written to `data/raw/omnipath/` with a checksum and registered in `config/datasets.yaml` alongside the
OmniPath package version, the exact query parameters and the retrieval date. The analysis then depends on a file we
hold. The registered record states the **signed-and-directed subset size**, not OmniPath's headline interaction
count. This is a stronger guarantee than the other sources get, which is appropriate: it is the one curated
dependency Part B introduces into an otherwise measured-tier pipeline. The SIGNOR arm is cached the same way.

Why not SIGNOR as primary: it is better edge-for-edge (sign is the primary curated annotation, not a derived one,
and each entry is mechanism-level and PubMed-backed), but coverage is the binding constraint here. B5's floor is on
`mechanism_coverage`, and signed *path* consistency needs the whole chain present, so a sparse prior shrinks the
denominator until the flag rests on a handful of regulators per drug. SIGNOR is contained in OmniPath, so this is a
filter setting rather than a different resource, and the sensitivity arm is the same code path with a stricter
filter.

What this costs: provenance is heterogeneous, and aggregation compounds publication bias rather than correcting it,
because hubs accumulate edges across resources. Expect B5 criterion 2 to be harder to pass here than on SIGNOR, not
easier. B4 exists for that.

## B3 Primary method: signed path consistency (no solver)

For the top 25 regulators by `|activity|` in the drug's scoring-line consensus profile, find the shortest signed
path from the main target to that regulator, propagate the sign (drug action x product of edge signs), and compare
with the observed sign. `mechanism_consistent_frac` = consistent / reachable. Unreachable regulators are excluded
from the denominator and counted as `mechanism_coverage`.

Chosen as primary over CARNIVAL because the ILP would run over ~1,100 drugs x lines x 41 groups and needs a
commercial solver at interactome scale.

**Breast arm: CARNIVAL proper.** Run on the breast cohorts, where the controls are and where the known failure
lives: target as the perturbation node of known sign (inhibitor -1, agonist +1), measurements = the drug's regulator
NES in the scoring lines, signed prior as the network. Scope fixed to keep the ILP tractable: gate-passing FDA drugs
with a curated ChEMBL mechanism target (D105b) in the four breast groups, one instance per (drug, group).

Solver: CPLEX or Gurobi under an academic licence (CMU is expected to have Gurobi; confirm before scheduling). If
unavailable, first try HiGHS on the full prior with a time limit and a reported MIP gap — an instance left
unsolved-to-optimality with its gap stated is more honest than a reduced one. Only then fall back to restricting the
prior to the target's k-hop neighbourhood (k = 3), and note that this changes the problem rather than the runtime:
CARNIVAL can no longer route through a distant explanatory node, so consistency is systematically lowered for
long-range mechanisms and must be reported as neighbourhood-restricted everywhere. If neither runs, record the arm
as not executable rather than substituting the heuristic silently.

Read-out: agreement between CARNIVAL's per-drug consistency and `mechanism_consistent_frac` (Spearman, and
concordance on the flag). This is what licenses the heuristic outside breast. If they disagree materially (Spearman
< 0.5, or the flag disagrees on > 25% of drugs), the heuristic is not a stand-in, and either B runs on CARNIVAL only
where affordable or the flag is confined to breast.

## B4 Degree-correction contingency (fixed now, applied only if B5.2 fails)

If flag rate correlates with the main target's degree at `|rho| >= 0.3`, tightening the filter toward SIGNOR is
**not** the remedy — curated resources are degree-biased too, so it would trade coverage for no gain on the
confound. Instead score consistency against a degree-matched null: for each (target, regulator) pair, count the
signed paths of the same length from 100 nodes sampled to match the target's degree within +/- 20%, and report
`mechanism_consistent_frac` as the target's consistency minus the null mean, in null SD units. The flag then fires
on the corrected quantity.

Fixed in advance so the choice is not made after seeing which direction it moves results. Its own re-check is B5.2
re-run on the corrected score.

## B5 Pre-registered rule

1. **Known-failure control.** The flag must fire on anti-estrogens in MCF7, where the ARACNe reading of ESR1 is
   known to be wrong (C18, C25, D85). If it does not fire there, the filter does not work and D126 is closed.

   > **Result 2026-11-10 — failed; D126 closed.** 0 of 22 scorable anti-estrogen rows in the four MCF7 breast
   > cohorts are flagged (7 of 22 without the coverage precondition). ESR1 is reached in none of them: it is never
   > among the top 25 regulators by |activity|, which are FOXM1, MXD3, PA2G4, E2F2, MYC, TEAD4. The check reads
   > growth arrest, not mechanism, and so endorses the failure it was built to catch. B5.3 also fails (coverage
   > median 0.08); B5.2 passes (rho 0.165).
2. **Degree confound.** Flag rate vs the main target's degree in the prior, Spearman `|rho| < 0.3`. At or above
   that, B4 applies.
3. **Coverage floor.** `mechanism_coverage >= 0.4` for >= 60% of gate-passing FDA drugs, else the flag is too
   sparse to carry and is reported per drug only.

**Reported only, never ranks, gates or tiers**, following the D120 / D121 convention. Promotion to a gate needs its
own decision and its own controls.

Under R3: where the flag fires, the report must say the drug's observed regulator shift is not explained by its
annotated target in the prior used, and must name the prior and its release. This is a statement about the prior,
not proof the drug does not work.

## B6 Outputs

`results/lincs/d126/`: `mechanism_consistent_frac`, `mechanism_coverage`, `mechanism_inconsistent`,
`sign_disputed` counts, CARNIVAL agreement table; columns joined to the split tables. `scripts/d126/`.

# Risks

- **Whitening removes signal.** The D86 / latent-factor prior. Mitigated by the nested `lambda -> inf` limit and the
  sweep; not eliminated. The overlap graph is the narrow claim (shared targets = partly redundant evidence), which
  is defensible in a way that whitening the full empirical covariance is not.
- **A0a flatters the pipeline.** It conflates dose and time variation with assay noise, so `r` comes out low,
  `AUC_max` low, and attainment high. Where A0b does not cover a cohort, the ceiling is optimistic and flagged
  `r_optimistic`. Do not read an A0a attainment number alone as evidence the score is near its ceiling.
- **A0b may not reach cohort level at all.** It depends on the A0.0 count, which is unmeasured; if the deep pairs
  sit in a handful of lines, A0b is a line-level bound and most cohorts fall back to A0a.
- **Tier inversion.** Part B introduces a curated signaling prior as a new dependency. Contained by keeping it
  flag-only, by the cache-and-checksum rule, and by the degree-confound check; Part A's secondary arm is reported
  for the same reason.
- **Publication bias.** Both the prior and CollecTRI are biased to well-studied nodes. B5.2 tests this directly;
  Part A's primary arm avoids it by not using a curated graph.
- **Criterion A5.1 is a single drug class in a single tissue.** Anti-estrogens in breast are the only control where
  we know the right answer. Passing it is necessary, not sufficient.
- **Citations unverified.** PANACEA, CARNIVAL and the OmniPath filter field names in this spec and in
  `LITERATURE_SCORING.md` were recorded without fetching the sources. Check them against the papers and the pinned
  release before building on them.

# Open items

- ~~Which level-5 index fields define an independent stratum for A0a.~~ Settled 11-09 (D125-A0): the (dose, time) cell.
  The same decision fixes A0.0 / A0b on active level-5 signatures (level-4 profiles are not on disk), the drugs entering
  r_c, and the spliced-profile implementation.
- ~~A0b stopping rule.~~ Settled 11-09 (D125-A0.1): SE(median r_c) <= 0.01 after >= 20 splits. The same amendment
  restricts r_c, for both estimators, to drugs whose every scoring line changes between the passes.
- ~~Whether the OmniPath release needs a `config/datasets.yaml` review entry of its own, given it is a service rather
  than a fixed file.~~ Settled 11-10 (D126-B0): yes, with a `review:` key recording that the resource carries no
  upstream release label, that the cached file is the dependency, and that a re-pull is a new decision.
- Gurobi academic licence availability (expected yes at CMU). Decides whether the breast CARNIVAL arm runs on the
  full prior, on HiGHS with a reported gap, or on the 3-hop reduction.
