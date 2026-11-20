# ARACNe regulon variation is mostly context, not noise — and a correction

**Date:** 2026-11-19 · **Decision row:** D149 · **Prompted by:** the user, who challenged my reading
of D138.1 · **Code:** `scripts/systems_validation/observability/regulon_reproducibility.py`
**Data:** `results/systems_validation/observability/regulon_reproducibility.{tsv,json}`

## The correction

D138.1 measured that the same target's footprint shares about **1%** of its genes between two ARACNe
(Algorithm for the Reconstruction of Accurate Cellular Networks) networks, and I reported that as
ARACNe failing to reproduce itself — evidence the regulons were unreliable.

**That reading was wrong.** The user's objection was that context-dependent targets are expected
biology: a regulator's accessible sites and available co-factors differ by cell type, so the same
drug *should* engage different downstream genes in different diseases. Low cross-tissue overlap is
then signal, not failure.

That is testable, and the control already existed.

## The control

`regulons_aracne_project_v2_tf_brca_controls.rds` holds **BRCA_SUB1/2/3** — three networks each built
from **769 randomly drawn BRCA tumours** (D81). Same tissue, same disease mix, same protocol. The
biology is held fixed, so any regulon disagreement between them is **estimation noise alone**. That
is the floor against which every other comparison has to be read.

## Result

Median pairwise Jaccard (shared genes ÷ union) of regulon membership:

| regulator | same cohort, biology fixed (**noise floor**) | nested breast cohorts | across all tissues | drop attributable to context |
|---|---|---|---|---|
| ESR1 | **0.490** | 0.056 | 0.010 | 0.480 |
| MYC | **0.553** | 0.170 | 0.026 | 0.528 |
| FOXA1 | **0.581** | 0.079 | 0.015 | 0.566 |
| GATA3 | **0.528** | 0.096 | 0.010 | 0.518 |
| SOX10 | **0.648** | 0.197 | 0.032 | 0.616 |
| JUN | 0.373 | 0.104 | 0.069 | 0.305 |
| ETS1 | **0.564** | 0.228 | 0.067 | 0.498 |
| E2F1 | **0.551** | 0.209 | 0.132 | 0.419 |
| **median** | **0.552** | 0.170 | **0.029** | — |

**Estimation noise alone produces ~0.55 Jaccard, not 0.03.** ARACNe recovers roughly half of a
regulon from an independent sample of the same tumours. The further fall to 0.17 between nested
breast cohorts and 0.03 across tissues is therefore **not the algorithm being unreliable** — it is
the regulon genuinely differing by context, exactly as the user argued.

## What this changes

**It supports the per-cohort network design.** Using the cohort-matched network is the right call,
not a weak compromise. And D139 — swapping the entire regulon set moved the reversal score by only
Spearman 0.99 — now reads as robustness of the activity readout rather than insensitivity to its
input.

**It does not rescue the causal map.** D147's failure rests on different evidence entirely: ARACNe
agrees with curated signs on **0.0% / 0.0% / 2.5%** of edges the curated network calls inhibiting,
contains the *reversed* edge more often than the forward one (170 vs 157), and CARNIVAL returned a
literally **empty** solution on it. None of that is about reproducibility. Co-expression cannot encode
inhibition however reproducible it is.

**It sharpens one real limit.** A ~0.55 ceiling means roughly **half of any single regulon is sampling
artefact**. For activity inference, which aggregates over the whole set, that averages out — which is
why the readout survives. For a claim of the form "these specific genes are ESR1's targets in this
tumour", it does not: about half that list would not survive a resample of the same patients.

So the defensible statements are:

- *"This drug acts through ESR1, and in this disease context ESR1's programme includes these genes"* —
  supported, with the resampling caveat on any individual gene.
- *"These genes are ESR1's targets"* as a context-free claim — not supported.
- *"Drug inhibits target, which inhibits regulator, which drives these genes"* — not supported, for
  the separate sign/direction reasons in D147.

## A follow-up worth taking seriously

If context-specific rewiring is **signal**, then *how* a regulon differs between diseases is itself
information about what a drug would do differently in each — which is closer to the repurposing
question than anything measured so far. The noise floor established here is what makes that
measurable: a rewiring is only interpretable if it exceeds ~0.55 Jaccard of expected resampling
churn. Recorded as revisit item **R-5**.
