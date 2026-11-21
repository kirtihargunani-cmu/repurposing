# SPEC — RGES benchmark (D174, pre-registration)

**Status:** pre-registration. Written and committed **before the comparison was run**, because a
benchmark we implement ourselves and then win is weak evidence otherwise — every implementation
choice is ours, and without fixing them in advance we would tune our way to a win without meaning to.
**Date:** 2026-11-20 · **Raised by:** the user — *"can we benchmark it against the other paper that
tried to use mechanistic reversal to return computationally generated hits for various cancer types?"*

## 1. The question

Our pipeline reverses the disease signature at the level of **master regulators** (ARACNe + VIPER).
The canonical published alternative — Chen et al., Reverse Gene Expression Score (**RGES**) — reverses
it at the level of **genes**, using the Lamb/CMap enrichment statistic. Both consume the same two
inputs we already have: a tumour-vs-normal signature and a drug's perturbational profile.

> **Does the regulon layer earn its complexity?** If gene-level RGES matches or beats our
> regulon-level reversal on our own evaluation, the ARACNe/VIPER machinery is not paying for itself
> and we should say so.

This is the baseline question. "3.84× chance" answers *better than random*; a reviewer asks *better
than something simple*.

## 2. What is held identical (so the only difference is gene vs regulon)

| | fixed to |
|---|---|
| disease signature | the contrast production uses per cohort (`results/step4_v2/signature_summary.tsv`, `status == ok`), gene-level `z` from `data/processed/step4_v2/<group>__<contrast>.tsv.gz` |
| drug profile | `drug_line_profiles_bing_specific.f32` — **the same matrix `79e` scores**; `raw` reported as a sensitivity arm only |
| cell lines | the same Celligner in-key lines and weights `79e` uses (`disease_line_weights_long.tsv`) |
| aggregation | weighted mean over lines, as `79e` |
| gene universe | the 10,112 genes of the LINCS profile matrix |

## 3. RGES, fixed now

1. Disease gene sets: among genes with `adj.P.Val < 0.05` **and** present in the LINCS gene
   universe, take the **top 150 by `z`** (up set) and the **bottom 150 by `z`** (down set). A cohort
   with fewer than 50 genes on either side is **excluded and reported**, not padded.
2. For each drug × line profile, rank all genes by their value, descending.
3. `ES_up`, `ES_down` = the Lamb/CMap Kolmogorov–Smirnov enrichment statistic of each set against
   that ranking.
4. `RGES = ES_up − ES_down`, **except** that `RGES = 0` when `ES_up` and `ES_down` share a sign
   (the published rule: a drug that moves both sets the same way is not reversing anything).
5. More **negative** RGES = stronger reversal. For comparability with our `rev` (higher = better)
   we report **`rges_rev = −RGES`** throughout. No other transformation.

## 4. The comparison, fixed now

Scored on the **D127 biomarker pairs**, unchanged: the same 49 pre-registered primary instances,
7 pairs, 4 classes, the same `high > low` direction test and sign test. D127 is the right venue
because each instance is the same drug in two cohorts, so every drug-level confound cancels — and
because our own reversal score already has a published number there to beat.

Reported side by side on **identical instances**:

| score | what it is | current value |
|---|---|---|
| `rges_rev` | gene-level reversal (this benchmark) | — |
| `spec_mech_pct_pregate` | our regulon reversal, mechanism-specific | **0.6939** (p = 0.0047) |
| `rev_pct` | our regulon reversal, raw | — |
| `layered_pct` | our lead score (target evidence, **not** reversal) | **0.7755** (p = 7e-05) |

**No pass/fail threshold is set, because this is a comparison and not a gate.** The result is
whichever ordering the numbers give. Three outcomes are named in advance so none can be presented
as a surprise:

- **RGES ≥ our reversal.** The regulon layer is not earning its complexity on this test. Report it.
- **RGES clearly below.** Evidence the regulon layer adds something, on one 49-instance test.
- **Both near 0.5.** Neither reversal score is carrying the ranking — consistent with D165, where
  the lead score was found never to read the reversal at all.

## 5. Declared limits, before results

- **49 instances, 4 drug classes.** This cannot separate close results. A gap of one or two
  instances is noise, as D166 found when 38/49 vs 36/49 turned on two 0.9-point margins.
- **Our implementation of someone else's method.** Chen et al. tuned RGES on their own data and
  their own sensitivity endpoint; a loss here is not proof their method is worse, only that it is
  worse *on our inputs under our fixed choices*. That asymmetry is unavoidable and must be stated
  wherever the result is quoted.
- **Target recovery (D158) is not part of this benchmark.** RGES produces no regulator ranking, so
  there is nothing to compare — this is a difference in what the methods output, not a result.
- **OncoTreat is deliberately excluded.** It is the same method family as ours (ARACNe + VIPER
  master-regulator reversal) and reimplementing a clinical pipeline from a paper would mostly
  measure how well we guessed its parameters. It is cited as provenance, not benchmarked.

## 6. Related

D127 (the test), D158/Route 1 (target recovery, not applicable here), D165 (the lead score never
reads the reversal), `docs/LITERATURE_SCORING.md` (RGES and OncoTreat prior art).
