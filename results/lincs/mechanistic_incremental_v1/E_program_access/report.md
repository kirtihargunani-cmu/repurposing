# Experiment E — which disease programs can drugs reproducibly perturb?

Protocol `5cc2227` + E settings `6583b74` (frozen before any E computation). Programs = the D198 modules rebuilt in
this study's folder from the 2026-11-23 regulator contributions (`E0_programs.py`, production code with redirected
paths; production files untouched); patient-level modules (`58`, 2026-11-22) used as an upstream input. 29 cancers,
461 programs, 313,350 drug-program responses, 62 s. Explanatory only: no program score is proposed for ranking (D201).

## What resolved

**1. Every disease program responds to drugs reproducibly.** Program opposition split by conditions agrees across
drugs at a median Spearman-Brown reliability of **0.61** early (199 program-windows, 14 cancers) and **0.76** late (275,
17 cancers), the same level as the gene-level score in Experiment A (0.61 / 0.74). A drug's effect on a disease program
is a stable, measurable property at both timepoints.

**2. Responses are driven by the drug, not by the cell line.** After removing each drug's mean, cell line explains a
median **2.7%** (early, 165 programs) and **2.5%** (late, 218) of the remaining variance in program opposition. Pooling
a cancer's matched lines, as the production scorer does, loses little.

## Not resolved at current data depth

- **Disease programs are not more perturbable than random regulator sets of the same size.** Size- and
  |NES|-tertile-matched pseudo-programs are as reproducible (median 0.64 early, 0.77 late); 3 of 474 program-windows sit
  above the pseudo-program 95th percentile (1 clean early, 1 plate early, 1 plate late), fewer than the 5% expected by
  chance. Reproducibility here is the drug-identity signal of Experiment A carried into any regulator set, not a
  property of the disease-defined grouping.
- **Drugs whose target links into a program do not perturb it more.** Linked minus unlinked |opposition| percentile:
  -0.017 early (86 programs) and -0.014 late (130), intervals spanning 0; residualised on the drug's overall response
  magnitude, +0.016 and +0.012, also spanning 0; the same for clean programs only. Target-to-program links (primary
  target or its one-hop OmniPath targets among members) do not yet predict which programs a drug moves.
- **Plate-flagged programs are not stronger evidence:** their reliability (0.62 / 0.77) matches clean programs (0.58 /
  0.74), as the protocol required.

## What this supports

Program-level responses are reproducible and drug-driven, so programs remain a valid lens for interpreting a drug's
effect (D198), but at current depth they add no accessibility signal beyond the regulator-level response itself and
should not be used to rank (consistent with D196 and D201). Resolving program-specific accessibility would need more
repeated conditions per drug (median one signature per condition) or programs defined to match drug-response
covariation.

**Files:** `E_summary.tsv`, `E_repeatability.tsv`, `E_target_link.tsv`, `E_line_selectivity.tsv`,
`E_program_response.tsv.gz`, `programs/` (rebuilt programs and flags), `run.json`.
