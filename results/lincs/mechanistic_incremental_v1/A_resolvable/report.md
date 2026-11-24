# Experiment A — how much signal each processing step resolves

Protocol `5cc2227` (spec Amendment 1, A1.5). Run `20261123T163736-3317604`, 77 s, inputs = the 2026-11-23 full run
(Stage 0 pass). 7,842 drug-cancer-window rows split in half; 1,000 cluster-bootstrap resamples; 1,000 permutations.
Resolved = the 95% interval excludes the matched null (A1.7).

## What resolved

**Every processing step carries drug-specific signal, in both windows, for cancers and subgroups.** The two halves
of one drug's conditions agree far more than two different drugs in the same cancer, condition-count stratum and
window, and the signal sharpens as the profile is projected onto the cancer's ARACNe network:

| Step | Cancers, early | Cancers, late | Subgroups, early | Subgroups, late |
|---|---|---|---|---|
| S1 gene profile (same minus different-drug cosine) | **0.35** (0.33-0.38) | **0.56** (0.54-0.58) | **0.30** (0.24-0.36) | **0.62** (0.58-0.66) |
| S2 single-condition regulator activity | **0.37** (0.34-0.40) | **0.50** (0.48-0.52) | **0.29** (0.20-0.37) | **0.58** (0.54-0.62) |
| S3 window regulator activity | **0.43** (0.40-0.46) | **0.62** (0.59-0.65) | **0.35** (0.26-0.46) | **0.69** (0.62-0.75) |

Different-drug cosines sit at -0.03 to +0.07 throughout, so nearly all of the same-drug agreement is drug-specific.
S1 -> S3 adds +0.08 (early) and +0.06 (late) in cancers: aggregating over matched lines on the cancer's network
resolves drug identity better than the raw gene profile.

**The scores themselves are reproducible across drugs** (split-half Spearman-Brown reliability per cancer; the
square root is the best correlation any test of that score can reach):

| Score | Cancers, early (14 cancers) | Cancers, late (17) | Subgroups, late (4) |
|---|---|---|---|
| S4 regulator-level reversal `rev` (reference) | **0.65** (0.62-0.69) | **0.80** (0.75-0.82) | **0.86** (0.81-0.87) |
| S5 gene-level reversal `rev_gene` (production score) | **0.61** (0.58-0.69) | **0.74** (0.71-0.84) | **0.73** (0.17-0.85) |
| S6 `rev_gene` within the eligible set (production ranking) | **0.43** (0.31-0.59), 7 cancers | **0.62** (0.58-0.76), 13 cancers | 0.80, 2 subgroups |

Reliability ceilings for the production score: 0.78 early and 0.86 late in cancers. S6 is lower than S5 partly
because the eligible set is smaller and narrower (restriction of range), not only because of noise.

**The signal is mechanistic, by two labels that share no input with the filter or the score:**

| Label | Step | Cancers, early | Cancers, late | Subgroups, late |
|---|---|---|---|---|
| (b) Drug moves its on-target PROGENy pathway the expected way (share; null) | S1 | **0.76** (0.70-0.81); null 0.65 | **0.77** (0.73-0.82); null 0.64 | **0.79** (0.70-0.87); null 0.65 |
| (c) Drugs sharing a target unit vs other pairs (cosine gain) | S1 | **+0.28** (0.15-0.38) | **+0.38** (0.22-0.45) | **+0.50** (0.24-0.66) |
| | S3 | **+0.36** (0.24-0.49) | **+0.44** (0.28-0.53) | **+0.53** (0.29-0.78) |
| | S5 | **+0.30** (0.21-0.44) | **+0.36** (0.23-0.51) | **+0.50** (0.28-0.68) |

PROGENy: 1,266 early and 1,607 late labelled rows (172 and 216 ingredients), permutation p = 0.001 in both windows.
Shared target: 49-59 distinct target units in cancers. The production score's per-regulator vector (S5) keeps the
shared-target signal at about the strength of the regulator activity it comes from.

## Not resolved at current data depth

- **Primary-target rank in the drug's own regulator response** (label a): mean percentile 0.48-0.56 in cancers, intervals
  spanning 0.5. Only **14 distinct targets** are both a drug's primary target and a scored regulator (most drug targets
  are kinases or enzymes, not transcriptional regulators), so the interval is wide; the top-5% rate is at or slightly
  above chance (0.04-0.09 vs 0.05). This needs more regulator-type targets, or the downstream-pathway route of
  Experiment C, to resolve.
- **Subgroup early reliability (S4-S6)** and **subgroup S6 late**: 2 subgroups each with >= 10 split drugs, too few for
  an interval. Subgroup early PROGENy sits exactly at its null bound (0.71) despite p = 0.003.

## What this supports and what it does not

It establishes that the perturbation data carry reproducible, drug-specific, mechanism-consistent signal at every
step of the production scorer, in both windows, and that the production score is reliable enough to be studied
further (late reliability 0.74, ceiling 0.86). It does not establish efficacy, and it does not show that the score
ranks drugs better than another score; Amendment 1 does not ask that question. Not used: PRISM/GDSC (A1.4).

**Files:** `A_curve.tsv` (steps), `A_reliability_cells.tsv` (per cancer), `A_labels.tsv` (labels),
`A_halves.tsv.gz`, `A_progeny.tsv.gz`, `A_target_rank.tsv.gz`, `A_shared_target_pairs.tsv.gz`, `run.json`.
**Next:** Experiment B uses the split-half residuals from the 400 drug x line pairs with repeat conditions in both
windows as its noise model (A1.6).
