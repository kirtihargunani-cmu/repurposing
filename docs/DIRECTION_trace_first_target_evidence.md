# Trace-first target mechanism evidence — design direction

**Status:** design note, not a new experimental result. **Repository state reviewed:** `main` at `9b2304bb` (November 17, 2026).

This document accompanies [the interactive testing history](testing_progression.html), [the mechanism trace procedure](PROCEDURE_mechanism_trace.md), and [settled measurement facts](notes/settled_measurement_facts.md). It records the decision to preserve evidence-level quantitative measurements but stop outputting a global numeric drug-suitability score.

## Why the output changes

D127's layered score correctly ordered 36/49 known primary development biomarker pairs. D128's target-anchored mechanism M ordered 16/49. D128 also found BRAF/MYC hub dominance; D130 corrected a falsely narrow gene-independence null; D131/D133 found strong regulator dependence. A single number hides mechanisms that break because a target is absent, a cell line is mismatched, a footprint is dominated by cytotoxicity, or exposure information is missing. It remains valid to retain historical numerical metrics as diagnostics and benchmarks; it is **not** valid to relabel them causal drug suitability.

## New evidence object 1: disease–target state (L2)

**Grain:** `(cohort, molecular_subgroup, target_gene_or_unit, defined_patient_population, version)`.

Retain separately: alteration prevalence **with actual numerator/denominator** (activating hotspot, loss-of-function, amplification, deletion, fusion, unknown function); individual-patient target/regulon expression or VIPER activity distributions; effect sizes and matched reference; functional dependency across *genotype/lineage-relevant* DepMap lines (explicitly not tumor-patient prevalence); independent driver and perturbation evidence; resistance/bypass, tumor purity and normal-tissue alternatives. A highly expressed target is not automatically a driver; a recurrent missense event is not automatically activating. Target-wide and biomarker-subgroup statements require different denominators.

## New evidence object 2: drug–target pharmacology (L1, exposure)

**Grain:** `(ingredient, molecule, target_unit, measured_assay, molecular_context, dose/exposure_source)`.

Retain ChEMBL direct mechanism/target source, agonist/antagonist/inhibitor/degrader action, family or complex membership and binding type. **Do not collapse Kd/Ki, biochemical/cellular IC50, and EC50** into one undocumented `binding efficacy` scalar. Keep pChEMBL's parent standard measurement type, assay, molecule and confidence. Missing polypharmacology evidence does not establish selectivity.

If a valid reversible occupancy model applies, use free concentration and a supported Kd and disclose assumptions; an IC50 is *not* automatically a Kd. D123 free fraction and Cmax are incomplete and plasma measurements do not prove tumor penetration. Occupancy and pharmacological efficacy are distinct. A missing Kd, fu or tumor exposure is `not_evaluable`, not a negative score.

The reviewed repository contains 372/497 numeric pChEMBL target-unit rows in `results/lincs/target_role/d121/units.tsv`, and an independently scoped 86/231 free-occupancy entries in `results/lincs/d123/occupancy.tsv`. These denominators are different populations, not two fractions to be added.

## Seven-link execution and immutable prediction

The [procedure](PROCEDURE_mechanism_trace.md) defines L1 drug→target, L2 target→disease, L3 target→regulators, L4 regulators→genes, L5 predicted→observed, L6 concentration/confound audit, L7 disease closure. Every link has `expected`, `observed`, `verdict`, `provenance`, `independence`, `concentration`, `alternatives`.

**Forward-only order:** commit the L1–L4 predicted target footprint and gene directions before reading drug-response or outcome data, then evaluate L5–L7 on those frozen predictions. No backfilling a mechanism from a favorable observation. A broken link makes the trace broken; no averaging, no final drug score. An unevaluable link remains unevaluable, not a zero.

T1 already froze fulvestrant→ESR1→BRCA-LumA in MCF7 and vemurafenib→BRAF→SKCM in A375 at `results/lincs/trace/t1_run.json`. L1–L4 freeze is *not* a test of response. For T2, use the frozen expected signs and show all disagreements, not a post-hoc set of supportive genes.

## Independence and confounds

- Signed target→TF paths and TF→gene sets share graph/regulon evidence; they are not two independent confirmations.
- VIPER and its target genes share an expression source; do not count them twice.
- ARACNe co-expression networks and disease signatures may share source tissue information.
- DepMap dependency is evidence from cancer cell lines, not patient prevalence.
- Report unique-gene attribution and hub share (D128 BRAF case ~65.1% structural MYC share); proliferation/cytotoxicity, polypharmacology, target-absent and line-context alternatives must remain explicit.
- D129 establishes inferred-gene quality is generally good (median held-out R² 0.875); **only** target-self-QC on the target's imputed mRNA is not supported.
- D130 correlation-preserving nulls are mandatory; do not independently shuffle correlated genes.

## Falsification, scope, next steps

Negative controls before any positive narrative: target absent in source line; RAF inhibitor in BRAF-WT/RAS-active context; random drug–cohort pair; cytotoxic response explained by generic death program; target with no anchored route. Report false-closure rate and first broken link **including failed traces**.

D134 target recoverability and D135 drug-vs-genetic perturbation have completed *S0 feasibility*, not measured target-recovery/concordance. D136 dose-response is design only. These may become orthogonal links, not excuses to treat a broken trace as successful.

**Next action:** implement L5–L7 for frozen T1 pairs; add typed disease-state and drug-binding adapters; run the five negative control classes; render each chain from evidence fields. Production ranking code and target scores are left unchanged by this documentation.