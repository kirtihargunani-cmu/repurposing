---
title: "Target-Mechanism Scoring: Full Implementation Plan"
subtitle: "D128 candidate | Repurposing project | Engineering and experimental validation"
author: "Implementation plan based on the D128 target-mechanism scoring specification"
date: "15 November 2026 (revised 15 November 2026: Section 0 added; workflow aligned with the project's main-branch practice)"
lang: en-US
geometry: "margin=0.85in"
fontsize: 10pt
linestretch: 1.12
colorlinks: true
linkcolor: blue
header-includes:
  - \usepackage{longtable}
  - \usepackage{booktabs}
  - \usepackage{array}
  - \usepackage{fancyhdr}
  - \pagestyle{fancy}
  - \fancyhf{}
  - \fancyfoot[C]{\thepage}
  - \renewcommand{\headrulewidth}{0pt}
  - \setlength{\emergencystretch}{3em}
---

# Executive summary

**Goal.** Extend the existing drug-repurposing pipeline with interpretable **target-anchored mechanism concordance (M)** and **target-local regulatory rescue (R)** scores. The implementation must retain the signed, weighted TF-to-gene network and individual gene-level LINCS responses, distinguish evidence quality from mechanistic disagreement, and establish whether these new measurements add value beyond the project's existing target-role, VIPER-reversal and layered rankings.

**Scope.** This is an implementation and experimental plan, **not** a report that D128 has been implemented or run. Work begins with read-only adapters for the project's already computed disease signatures, LINCS-specific profiles, ChEMBL target annotations, VIPER activities, and regulatory networks. All new outputs are sandboxed. The production lead ranking remains unchanged unless a later validation decision approves a change.

**Scientific distinction.** Four propositions should never be conflated:

- **T: Target relevance** - is the annotated target important in the disease? Already addressed by D120/D121/D122 and subsequent production refinements.
- **M: Mechanism concordance** - does the observed drug-response footprint match an independently annotated route from its target through downstream regulators to genes?
- **R: Regulatory rescue** - do the drug-responsive genes in that target-linked footprint move opposite the disease-associated direction?
- **V: Functional benefit** - does independent drug-sensitivity or phenotype evidence support therapeutic activity?

The first release reports **T, M, R, and evidence coverage separately**. No new all-in-one scoring formula is allowed until the component-level tests pass.

**Key previous results.** The D126 shortest-path check failed to flag 22/22 anti-estrogen/MCF7 control rows and had low mechanism coverage, so the new method anchors at the annotated target rather than the top 25 changed regulators. D127's current layered score passed 36/49 preselected biomarker-pair instances (73.47%), but the thyroid BRAF-like versus RAS-like pair was only 3/9. These are **seen development benchmarks**, not independent holdouts for D128.

**Precursor that bears directly on R (2026-11-12, commit `55c80d7`, `docs/PROGRESS.md`).** Gene-level reversal restricted to ESR1's regulon, with `tfmode`/`likelihood` kept, *inverted* the anti-estrogens in the ARACNe arm (ESR1 ranked them last in BRCA-Her2, LumB and LumA, and high only in ER-negative Basal); the CollecTRI arm showed no inversion but no lift. The diagnosed cause is the reference, not resolution: GTEx normal breast is ER-active, so ESR1 targets do not read as up in tumor. R (Section 3.3) is essentially that construction, so its behavior in breast is **predicted before any D128 run**. Section 0 (D0.1) and the A4 gate are written so that this prediction is recorded, not rediscovered and then tuned away.

## Deliverables at a glance

| Artifact | Required contents | Acceptance signal |
|:--|:--|:--|
| Coverage audit | Candidate drugs, cohorts, source lines, target/genotype status, network and L1000 coverage | All candidate inputs accounted for; noncomputable cases named |
| Target footprint graph | Target-unit to TF to gene paths, signs, confidence and ambiguity | Synthetic graphs have exactly reproducible paths/signs |
| Scoring module | M and R per drug-line-unit/cohort, with gene- and regulator-level attributions | Hand-calculated fixtures match; no double counting |
| Cohort aggregation | Reused Celligner/source-line logic; direct/borrowed flags | No silent denominator or weighting changes |
| Controlled benchmark | Matched positive/negative mechanisms, D127 replay and molecular holdouts | Results compared fairly against current production score |
| Robustness + decision | Assay, network, dose/time, degree-null and functional-screen analyses | Clear promote / retain as diagnostic / reject verdict |

# 0. Primary arm, frozen defaults and open decisions (revision 2026-11-15)

The first draft listed many sensitivity arms but named no primary arm, and left several weights undefined. That leaves room to report whichever arm passes. This section names **one primary arm** for every free choice. Everything else in the plan is a sensitivity arm, reported beside the primary and never substituted for it. The issues behind each entry are in `docs/notes/d128_plan_review.md`.

**How these become binding.** As with D125-A0, the choices are written into `docs/DECISIONS.md` as a **D128-A0** entry ("implementation choices, fixed before running"), together with any value A0 measures label-blind, such as coverage floors. That entry is committed and pushed **before** any M or R value is computed. Rows marked *user decision* stay open until the user decides them in that entry. Until then, A0 and A1 code may be written and run, because they compute no scores, but A2 does not run.

| ID | Choice | Primary arm | Sensitivity arms (reported only) | Basis | Status |
|:--|:--|:--|:--|:--|:--|
| P1 | Disease gene effect | Production Step 4 limma + voomByGroup (D69) columns of `data/processed/step4_v2/<group>__<contrast>.tsv.gz`: sign from `logFC`, magnitude from `z`, `adj.P.Val` for the gene set | `*_voom` columns (check C12) | Step 6 VIPER reads the `z` matrix (`50_viper.R:80`); `79t_target_role.py:115` reads `logFC` | Settled |
| P2 | Contrast per cohort | 79t's rule: `gtex` if its `signature_summary.tsv` status is `ok`, else `adjacent` | Breast: see D0.1 | Same rule as the production target-role score | Settled |
| P3 | Regulon network | **CollecTRI** (context-free, curated signs) for both M and R | ARACNe v2 TF cohort networks (`aracne_project_v2_tf`, subtype networks where production uses them) | ARACNe `tfmode` is a co-expression sign learned in tumors, not in the drug's line; the precursor's ARACNe ESR1 arm inverted | Proposed |
| P4 | Gene feature space | All 10,112 BING genes in the `79c` profiles (978 landmark + 9,134 best inferred) | Landmark-only; must agree in sign for promotion (Section 7.3) | Production 79c/79d use BING; landmark-only coverage is too thin for most footprints, and A0 measures by how much | Proposed |
| P5 | Drug profile | `79c` line-specific profiles (`drug_line_profiles_bing_specific.f32`) | Raw profiles (`drug_line_profiles_bing.f32`), required for the ESR1 control | Production uses specific profiles. If specific fails and raw passes, the verdict is "generic removal removes the target signal", **not** a pass | Proposed |
| P6 | Signaling hops | 0-2 hops, OmniPath primary (D126 cache) | 0, 1, 3 hops; SIGNOR-only | Cap from the draft; frozen here before any outcome is read | Proposed |
| P7 | M gene weight `wM` | **Unweighted:** each unique sign-resolved gene weighs 1 within a target unit | Edge-confidence weighted (ARACNe `likelihood`; CollecTRI weights where not ±1) | The draft's composite of confidence, assay reliability and annotation support had no definition; an unweighted primary has nothing to tune | Proposed |
| P8 | R gene set and `wD` | `G_ct` = target-linked genes with `adj.P.Val < 0.05` in the cohort's P2 contrast, weight 1 each (`sign(logFC)` therefore carries only significant genes) | `|z|` weight, capped at 10; all target-linked genes with `|z|` weight | Removes the full ±1 vote of null genes; label-free | Proposed |
| P9 | Bounded transform | `f(z) = tanh(z/2)` | `tanh(z/1)`, `tanh(z/3)` | A0 audits the distribution of the specific-profile values (they are residuals, not raw MODZ). The divisor may change **only** on label-blind distributional grounds, in the D128-A0 entry | Proposed |
| P10 | Minimum coverage | Kinase/receptor unit: >= 10 unique sign-resolved genes and >= 2 anchored TFs; direct TF unit: >= 10 regulon genes. Otherwise `NA` + `abstain_reason` | Floors x0.5 and x2 | Spec's suggested guards; A0 reports how many rows each floor admits before it is frozen | Proposed |
| P11 | Drug action sign | Reuse `79t_target_role.py`'s `LOWER`/`RAISE` sets and `orient()`. Rows `orient()` flags `direction_assumed` (e.g. `MODULATOR`) or `mixed` are **ambiguous for M** (E07), except under D0.2 | Production behavior (`direction_assumed` -> lower) | One action vocabulary across the project; E07 forbids a default inhibitor sign in M | Proposed; ESR1 SERMs: D0.2 |
| P12 | Compound key | `pert_id`, as in production (79e/79t rows). A0 lists compounds with several `pert_id`s (e.g. fulvestrant: 3 in MCF7); controls report each `pert_id` and their mean | | Never silently merge `pert_id`s | Proposed |
| P13 | Multi-unit drugs | Mean over scorable target units, unweighted | Potency-weighted (D121 units) | Draft Section 3.4; no best-unit selection | Proposed |
| P14 | Line aggregation | 79e's lines and Celligner weights, unchanged (direct/borrowed) | Direct-only; leave-one-line-out | Draft Section 3.4 | Settled |

## 0.1 Decisions that needed the user — SETTLED 2026-11-15 (row `D128-D0.1 / D0.2`)

**Both were decided as option (a).** D0.1: the production contrast (P2) is kept for R everywhere, including breast; the expected result is recorded in advance (uninformative on the CollecTRI primary, inverted on the ARACNe arm) and the mandatory anti-estrogen control is gated on **M alone**. D0.2: fulvestrant is the primary positive set, and the SERMs are a second, separately reported set scored as ESR1 antagonists in breast-lineage lines only. The options below are retained as the record of what was weighed; (b) and (c) under D0.1 are computed as reported-only arms.

**D0.1 - Reference for R in breast.** The P2 contrast for BRCA-LumA/LumB/Her2 is tumor vs GTEx normal breast, and the precursor shows ESR1 targets inverted against it. Options:

- **(a) Recommended: keep P2 for R everywhere and record a prediction.** The D128-A0 entry states in advance that R is expected to rank anti-estrogens low in the ER-positive cohorts against GTEx. The mandatory ESR1 control is then gated on **M only**, and R in breast is reported, with that prediction beside it. Rationale: R stays on the same reference as every production score, and the choice is not made after seeing which reference works.
- (b) Use adjacent normals for breast R (LumA 62, LumB 22; Her2 has only 9 and is skipped). Adjacent normal breast is also ER-active, so it may not fix the problem.
- (c) Use a within-tissue subtype contrast (LumA/LumB vs Basal) as the breast disease direction. This changes R's meaning in breast from "toward normal" to "away from the other subtype".

Under every option, (b) and (c) are computed as reported-only arms for breast.

**D0.2 - SERMs in the ESR1 control.** ChEMBL annotates fulvestrant on ESR1 as `ANTAGONIST|DEGRADER`. Tamoxifen, raloxifene, clomifene, toremifene and bazedoxifene are `MODULATOR` on ESR1 (production: lower, `direction_assumed`). Tamoxifen's one clean `ANTAGONIST` row is on ESRRG, not ESR1. MCF7 coverage is 18 fulvestrant signatures (3 `pert_id`s) and 3-27 per SERM (raloxifene 27, tamoxifen 10 across 3 `pert_id`s, bazedoxifene 5, clomifene 4, toremifene 3). Options:

- **(a) Recommended:** the **primary positive set is fulvestrant**. A second, separately reported positive set is the SERMs scored as ESR1 antagonists **in breast-lineage lines only**, under a written tissue-context rule (well established in breast). Elsewhere `MODULATOR` stays ambiguous under P11.
- (b) Fulvestrant only; SERMs abstain everywhere.
- (c) SERMs treated as antagonists everywhere, as production does.

There is no agonist sign-flip control: estradiol is not profiled in MCF7. A0 records this as a gap.

# 1. Current repository interfaces and reuse strategy

**Repository:** [kirtihargunani-cmu/repurposing](https://github.com/kirtihargunani-cmu/repurposing)  
**Parent specification:** `docs/SPEC_target_mechanism_scoring_benchmark.md` (draft PR #2, merged into main 2026-11-15). Where the spec and this plan number stages differently (in the spec `a1_score` is scoring; here A1 is the graph), **this plan's numbering is authoritative**.  
**Implementation location:** `scripts/d128/` (code and configuration; configuration is a Python module, `scripts/d128/config.py`, as in `scripts/d125/a0_config.py`), `tests/d128/`, `results/lincs/d128/`.  
**Review notes:** `docs/notes/d128_plan_review.md` (issues found before A0; Section 0 resolves them).

The existing data transformations are the source of truth. Rebuilding them increases computational cost and risks changing the definition of an experiment while claiming a like-for-like comparison.

| Current source | Reused information | Important contract |
|:--|:--|:--|
| `scripts/79c_drug_profiles.py` | Drug-by-line LINCS gene signatures, generic-response-removed | `drug_line_profiles_bing_specific.f32`, profile-major float32; index has `pert_id`, `ModelID`, `n_sig` |
| `scripts/79d_drug_viper.R` | Actual drug VIPER activities | Per-network `.f32` with `.regulators.txt`; verify layout and value alignment |
| `data/processed/step4_v2/<group>__<contrast>.tsv.gz`, `results/step4_v2/signature_summary.tsv` | Disease gene differential expression | Columns `gene_id, symbol, logFC, t, adj.P.Val, z` (limma + voomByGroup, primary); `*_voom` columns are check C12 only. Contrast chosen per P2 |
| `scripts/precursor/esr1_regulon_reversal.py`, `results/lincs/diagnostics/precursor_esr1_regulon/` | 11-12 ESR1 regulon gene-level reversal | Prior result that predicts R's breast behavior (D0.1); its consensus-over-lines code is a cross-check for A3 |
| `results/step6/` | Disease master-regulator NES | Match network, regulator and cohort exactly |
| `data/processed/step5/` and `results/step5/` | ARACNe, CollecTRI/DoRothEA regulons | Keep TF-to-gene sign (`tfmode`) and confidence (`likelihood`) separately |
| `scripts/79i2_target_presence.py` and ChEMBL cache | Main targets, compound mapping, mechanism targets | Preserve target units, complexes, annotations and unknown action |
| `scripts/79t_target_role.py`, `79t2_layers.py`, `79t3_lead_score.py` | Target-role components and current lead score | Current production includes D124b proliferation/essential-target behavior by default |
| `scripts/d126/b0_omnipath.py` | Signed, directed OmniPath/SIGNOR prior | Reuse cached files and checksums; do not refresh silently |
| `scripts/d127/d127_biomarker.py` | Baseline paired control instances | Replay identical pairs; mark all as previously seen |
| `scripts/77_lincs_extrapolation_test.py` | Original cross-cell-line validation | Optional transfer-fidelity diagnostic; not proof of clinical response |

**Non-goals:** do not re-download LINCS; do not modify `79c`, `79d`, `79e` or `79t3` in place; do not call a linear TF projection "VIPER"; do not treat high drug-induced transcriptional change as target engagement automatically; do not interpret regulatory rescue as cell killing.

# 2. Canonical representation and data contracts

Use **sparse normalized tables and matrix-backed arrays**, not duplicated nested JSON for every sample. A nested regulator array is a useful presentation layer; the computational representation should have separate immutable edges and context-dependent values.

Define the following objects:

$$
 B_{rg} = (s_{rg}, c_{rg}), \quad D_{cg}, \quad A_{cr}, \quad P_{dlg}, \quad Q_{dlr}, \quad U_{dt}.
$$

- `B`: TF (`r`) to gene (`g`) relation with sign `s` and confidence `c`.
- `D`: disease (`c`) gene effect, direction and reliability.
- `A`: disease regulator NES.
- `P`: drug (`d`) / source line (`l`) / gene response measured by LINCS Level-5 MODZ.
- `Q`: corresponding real drug VIPER regulator responses.
- `U`: annotated drug/target-unit relationships, action types and supporting evidence.
- `Path(t,r)`: signed signaling route from a drug target to a regulator, with hop count and ambiguity status.

## 2.1 Strict join keys

1. **Cohort:** `group_id`, optionally `stratum`, tissue key and contrast (`gtex`/`adjacent`). A project-defined broad subtype is not equivalent to a mutation-defined patient subgroup.
2. **Drug:** `pert_id` as LINCS profile identity, with explicit `cmap_name`, ChEMBL molecule/parent IDs and ingredient mapping. Ambiguous compounds remain marked, not automatically merged.
3. **Cell line:** DepMap `ModelID`, joined to LINCS cell name via existing coverage mapping. The profile index row represents `(pert_id, ModelID)`; duplicates are prohibited at the aggregated profile level.
4. **Gene:** project GENCODE-v23 stable gene ID for numerical joins; display symbols are secondary. Detect one-to-many symbol mappings and duplicated IDs.
5. **Regulator:** exact regulator identifier **plus network ID**. ARACNe per-disease network and the context-free CollecTRI network must never be merged by a bare gene symbol alone.
6. **Mechanism:** `target_unit_id`, component genes, action category and sign. A complex/family is a single target unit, not automatically several independent targets.

## 2.2 Proposed intermediate tables

| Table | Minimum fields / shape | Primary check |
|:--|:--|:--|
| `disease_genes.parquet` | `group_id,gene_id,logFC,z,t,adj_p,contrast` (P1, P2) | Contrast unique per unit; signed nonnull effect; values equal the step4_v2 files exactly |
| `disease_regulators.parquet` | `group_id,network,regulator,nes` | Aligned with drug network |
| `regulon_edges.parquet` | `network,regulator,gene_id,sign,confidence,source` | Sign ±1; confidence finite, positive |
| `drug_profiles_index.parquet` | `profile_idx,pert_id,ModelID,n_sig` | Number of rows equals binary matrix rows |
| `drug_genes.parquet` | `col_idx,gene_id,symbol,feature_space` | Column count and order match binary matrix |
| `drug_target_units.parquet` | `pert_id,unit_id,target_gene,action_sign,action_type,source` | Conflicted/unknown actions retained |
| `target_tf_paths.parquet` | `target_gene,network,tf,hops,path_sign,path_status` | Deterministic sign/path resolution |
| `line_weights.parquet` | `group_id,ModelID,weight,match_status,transfer_conf` | Existing direct/borrowed rules reproduced |

For every binary input, assert `(file bytes)/(4 bytes per float)` equals the expected product of dimensions. Read using `np.memmap` and sample a known row/column against the upstream producer. Add data hashes, row counts, software versions and Git commit SHA to `run_manifest.json`.

## 2.3 Crucial scale restriction

Tumor differential-expression **log fold change** and LINCS Level-5 MODZ **z-scores** are not the same physical quantity. **Never add them** as if they were directly calibrated post-treatment expression. The first implementation uses signed and bounded drug effects and a directional rescue score. An absolute normal-state distance or ``tumor + drug`` counterfactual requires a separate calibration study.

The LINCS BING gene set includes inferred features; report results both for the **978 directly measured landmarks** and for landmark-plus-inferred features, with reliability stratification. An effect from inferred targets is not equivalent to an independent measurement.

# 3. Scoring specification to implement

## 3.1 Build the target-anchored footprint

For a drug acting on a target unit `t`, identify routes:

`drug action -> target t -> signed signaling path -> TF r -> signed TF-gene edge -> gene g`.

For a direct TF target such as ESR1, allow `target = TF` at **zero signaling hops** even if ESR1 does not appear among the most changed drug VIPER regulators. For receptors and kinases, resolve signed paths in the frozen OmniPath graph, primarily at **0-2 signaling hops**; report 0/1/3-hop and SIGNOR-only sensitivity arms. The cap is an *exploratory proposal* and must be locked before new outcomes are read.

For each candidate gene:

$$
 e_{dt rg} = a_{dt}\;p_{tr}\;s_{rg},
$$

where `a` is the signed drug action, `p` the net path sign and `s` the TF-to-gene sign. If shortest valid paths disagree, mark the regulator/gene **ambiguous**, not whichever direction best matches the drug. Keep the conflicting paths for explanation.

**Coverage** is computed against the prespecified target-anchored gene universe *before* looking at whether drug responses are favorable. Report `n_target_units`, `n_anchored_tfs`, `n_genes_total`, `n_genes_direct`, `n_genes_inferred`, `n_sign_resolved`, `ambiguous_fraction`, and `hub_contribution`.

## 3.2 Mechanism concordance M

Let `z(d,l,g)` be the LINCS-specific gene response. Introduce an initially proposed bounded transform:

$$
 f(z) = \tanh(z / 2), \quad -1 \le f(z) \le 1.
$$

The divisor `2` is **not validated**; check its distribution in a label-blind A0 audit and freeze it before evaluating controls (P9). For each drug-line-target unit with sign-resolved genes:

$$
 M(d,l,t) = \frac{\sum_{g\in G_{dt}^{resolved}} w^M_{gt}\; e_{dtg}\;f(z(d,l,g))}{\sum_{g\in G_{dt}^{resolved}}w^M_{gt}}.
$$

In the primary arm `wM = 1` for every unique sign-resolved gene (P7); edge-confidence weighting is a sensitivity arm. Assay reliability enters through the landmark-only arm (P4), not through a weight. A gene linked through multiple TFs is counted **once per target unit**. Conflicting expectations are excluded from the signed numerator but counted in the ambiguity and coverage denominators. Empty or inadequately covered footprints produce `NA` with an explicit `abstain_reason`, not zero.

M is a measure of **transcriptional consistency with an annotated route**, **not** direct binding, pharmacological occupancy, or proof of causal engagement.

## 3.3 Target-local rescue R

Let `s(c,g) = sign(logFC(c,g))`, with disease-only weights `wD(c,g)` frozen independently from drug labels. For target-linked relevant genes:

$$
 R(d,c,l,t) = -\frac{\sum_{g\in G_{ct}} w^R_{cgt}\;s(c,g)\;f(z(d,l,g))}{\sum_{g\in G_{ct}}w^R_{cgt}}.
$$

Positive R means the target-linked gene response tends to oppose the disease direction. The gene set is chosen from the **precomputed footprint and disease data**, not after reading positive drug responses. Primary arm (P8): `G_ct` = target-linked genes with `adj.P.Val < 0.05` in the cohort's P2 contrast, `wR = 1`; `|z|`-weighted arms are sensitivity only. The reference for breast follows D0.1. Sign ambiguity can prevent M while a gene remains usable for directional R; record distinct coverage denominators for each score.

## 3.4 Aggregate without double counting

0. The drug key is `pert_id` (P12). Compounds with several `pert_id`s are reported per `pert_id`, and their mean is shown in controls; they are never merged silently.
1. Compute gene contributions once per `(drug, line, cohort, target unit, gene)`; a gene present in five TF regulons is not five independent observations.
2. Preserve per-regulator rows for explanation, but derive an *independently defined* unique-gene aggregate for each target unit.
3. For multi-target drugs, compute and retain each target-unit score. **Primary proposal:** mean over scorable units with unit weights fixed without outcome labels; use binding-potency weights only as a sensitivity arm. Do not select the best unit after observing rescue.
4. Aggregate source lines using the project's **existing** Celligner/direct/borrowed logic. Report `n_lines_used`, best line match, borrowed fraction, transfer confidence, source-line genotype and between-line heterogeneity. Missing line weights and missing M are tracked independently; no silent fallback to a confident result.
5. Preserve existing drug-level T and VIPER reversal in the final output alongside M and R. Do not mix them into a new lead score in the pilot.

## 3.5 Controls and null distributions

Use seeded nulls that preserve important structures: target-degree-matched alternative targets; shuffled TF-to-gene signs; gene-label permutations matched by regulon size/assay availability; unrelated drugs matched for transcriptional strength and cytotoxic/proliferation proxies. Score all nulls through the identical pipeline and report null median, spread, empirical percentile and hub sensitivity. No positive therapeutic label should enter these operations.

# 4. Incremental implementation work packages

All phases are **gated**. A later stage must not conceal a failing earlier stage. Tasks A0-A3 can be implemented modularly and developed in parallel after the shared schema is fixed, but biological results should be interpreted in order.

**Per-stage routine (the project's standing workflow; applies to every stage below):**

1. Before running a stage, add a `D128-<stage>` row to `docs/DECISIONS.md` recording the implementation choices, fixed before running. Commit and push.
2. Run the stage in the sandbox, with thread counts set explicitly and < 32 in total.
3. Add a `D128-<stage> result (date)` row to `DECISIONS.md` and a dated `docs/PROGRESS.md` entry, including failures and abstentions.
4. Update the plain-language drug-side doc (`docs/DRUG_SIDE_PLAIN_LANGUAGE.md`) and, where a reporting requirement arises, `docs/LLM_COORDINATION_NOTES.md`.
5. Commit and push to main, then merge main into `celligner` and push.
6. When a decision changes the drug-side design, bring the Claude Docs drug-side spec into step with the repo.

## A0 - Environment inventory and label-blind feasibility audit

**Files:** `config.py` (all Section 0 values, in one place), `io.py`, `a0_audit.py`, and optionally `a0_genotypes.py` (all under `scripts/d128/`).

**Actions**

1. Check required files and environments on the compute server. Record checksums, matrix shapes, effective memory demand, gene mapping rates, network names, target action coverage and source-line counts.
2. Implement read-only loaders with strict orientation and index verification. Use the existing ordered gene metadata; never sort the gene list.
3. Build all cohort x drug-class x compound x line x target-unit candidate rows, without consulting new functional outcomes.
4. Enumerate target/line molecular-state evidence (BRAF V600, EGFR activating mutations, ERBB2 copy number, RB1 loss); distinguish **known positive**, **known negative**, **unknown** and **mismatched**.
5. Quantify target-to-TF path and gene assay coverage on landmark-only and BING arms.
6. Report the label-blind quantities that Section 0 depends on:
   - the distribution of specific-profile values (P9);
   - the number of rows each coverage floor admits (P10);
   - landmark vs inferred share of each footprint (P4);
   - `orient()` flags per control compound (P11, D0.2);
   - `pert_id` multiplicity (P12);
   - SKBR3 ERBB2 copy number;
   - an RB1-null LINCS line, if one exists;
   - the A4 in-context status (Section 5).
7. Audit prior use of the external screens. List every script and decision that has read PRISM 20Q2 or GDSC 8.5 (known: `79j`, `79k`, `79m`; check D94) and whether that use touched cohort, compound or parameter choice. The A5 validation source is fixed from this list.
8. Freeze minimum scoreable coverage, the executable cohort contrasts and all Section 0 rows (including the user's D0.1/D0.2 calls) **only after this label-blind audit**, in the D128-A0 entry of `DECISIONS.md` plus a matching `scripts/d128/config.py` commit, tagged `d128-a0-freeze` and pushed.

**Outputs:** `coverage.tsv`, `eligibility.tsv`, `unmatched_ids.tsv`, `source_line_context.tsv`, `cohort_feasibility.tsv`, `profile_value_audit.tsv`, `external_screen_prior_use.tsv`, `run_manifest.json`.

**Pass gate:** 100% of candidate rows are labeled eligible, excluded or unscorable with a reason; no silent indexing ambiguities; proposed positive/negative controls have verified data availability; the D128-A0 entry is pushed with every Section 0 row marked settled. **If a drug or cohort is absent, do not invent a substitute.** No M or R value is computed before this gate.

## A1 - Regulatory graph construction and signed route indexing

**Files:** `scripts/d128/export_edges.R`, `scripts/d128/graph.py`, `scripts/d128/a1_build_paths.py`.

**Actions**

1. Export existing ARACNe-v2-TF and CollecTRI RDS networks with `tfmode` and `likelihood`, preserving network provenance and gene IDs.
2. Load cached D126 OmniPath primary and SIGNOR sensitivity graphs. Ensure signs represent regulatory effects, **not** differential gene expression.
3. Build a sparse target-to-TF path index, with zero-hop TF anchors, cycles handled, deterministic shortest-path signs, and explicit ambiguities.
4. Compose signs with TF-to-gene links; deduplicate shared downstream genes before scoring.
5. Cache footprints using target, network and source checksum keys so repeated cohort/drug computations reuse paths.

**Outputs:** `regulon_edges.parquet`, `target_tf_paths.parquet`, `target_gene_footprints.parquet`, `path_ambiguity.tsv`.

**Pass gate:** all directed/sign/path fixtures match hand calculations; ambiguous evidence remains ambiguous; no dependence on drug outcome direction.

## A2 - Mechanism and rescue score engine

**Files:** `scripts/d128/score.py`, `scripts/d128/a2_score.py`, `scripts/d128/nulls.py`.

**Actions**

1. Implement `score_mechanism()` and `score_rescue()` as pure functions with typed input arrays and explicit scoreability metadata.
2. Use bounded gene effects and predeclared weights. Return per-gene signed contributions, per-regulator explanation and per-target-unit scores.
3. Add unique-gene aggregation, ambiguous-sign masks, missing-data abstentions and measured-versus-inferred evidence splits.
4. Build null distributions that preserve degree, regulon size and transcriptomic response strength, and compute interpretable null statistics.
5. Profile chunked, vectorized execution across ~10,000 production drug-line profiles; use memory-mapped inputs and sparse operations.

**Outputs:** `scores_by_line.parquet`, `gene_contributions.parquet`, `regulator_contributions.parquet`, `target_units.parquet`, `nulls.parquet`, `score_run.json`.

**Pass gate:** numeric and invariance fixtures pass; all scores bounded to [-1,1] when defined; score of zero is distinguishable from `NA`; no target-only pathway with no genes gets a confident score.

## A3 - Cohort aggregation, uncertainty and baseline joins

**Files:** `scripts/d128/aggregate.py`, `scripts/d128/a3_aggregate.py`, `scripts/d128/baselines.py`.

**Actions**

1. Join source-line results back to the same direct/borrowed populations and Celligner weights used in `79e`; implement a diagnostic comparison against its line selection.
2. Handle multi-unit drugs with fixed aggregation, preserving each unit's score and the fraction of annotated units that remain unscorable.
3. Report uncertainty across biological source lines and dose/time signatures when available. Distinguish uncertainty from noise and from systematic genotype mismatch.
4. Join production `layered_pct` (including D124b behavior), `target_role_pct`, `spec_mech_pct_pregate`, and VIPER-based reversal. Never recompute the lead from an older implementation.
5. Produce rankings **only within a fixed, declared comparison population**. Show percentiles next to raw M and R so cross-cohort distributions are not confused.

**Outputs:** `scores_by_cohort.tsv`, `scores_by_cohort.parquet`, `line_heterogeneity.tsv`, `baseline_join_audit.tsv`, `abstentions.tsv`.

**Pass gate:** no unexpected loss/duplication in joins; existing production baseline columns reproduce source values; denominator and eligibility conventions are recorded.

## A4 - Controls and already-seen benchmark replay

**Files:** `scripts/d128/a4_controls.py`, `scripts/d128/a4_d127_compat.py`.

**Actions**

1. Run the mandatory **MCF7 anti-estrogen** positive mechanism test. D128 must explicitly include ESR1 target anchoring even where its drug VIPER activity is not in the top 25. Compare raw and generic-removed profiles and matched cytostatic controls. The positive set follows D0.2 (primary: fulvestrant, each of its 3 `pert_id`s). **The ESR1 gate is on M.** R in breast is reported against the prediction recorded under D0.1, not gated, unless the user's D0.1 call says otherwise. If M passes only on raw profiles (P5), record "generic removal removes the target signal", not a pass. Run this control first, before any other A4 item: it is the cheapest test of whether D128 can work at all.
2. Include eligible HER2-, BRAF- and EGFR-targeted compounds only after their LINCS line/genotype eligibility is confirmed.
3. Run mutation-/line-context negatives; especially RAF inhibitor paradox and inappropriate EGFR-wild-type/KRAS-active comparisons.
4. Replay exact D127 development pairs using their exact compound and cohort definitions. Report candidate M, R and coverage beside existing production lead percentile and target-role/reversal baselines.
5. Analyze failures by **mechanism** rather than just rank: missing network path, absent target expression, sign contradiction, inferred-only genes, cytotoxic signature, line mismatch, unknown genotype, polypharmacology.

**Outputs:** `control_results.tsv`, `d127_seen_pairs.tsv`, `d127_seen_summary.tsv`, `mechanism_case_reports/`.

**Pass gate:** high-information MCF7 ESR1 controls produce interpretable target-linked gene footprints (or a documented reason they cannot); negative controls do not become confident positives merely from common growth arrest. D127 is descriptive only; **do not optimize on its 49 previously seen instances**.

## A5 - Molecular holdouts and independent functional validation

**Files:** `scripts/d128/a5_genotype_pairs.py`, `scripts/d128/a5_external.py`, `scripts/d128/statistics.py`.

**Actions**

1. Before reading holdout scores, lock eligible mutation-defined patient strata and compounds. Priority: **BRAF-V600-mutant vs BRAF-WT/NRAS-driven melanoma** and **activating EGFR-mutant vs EGFR-WT/KRAS-mutant LUAD**. HER2 amplified vs nonamplified and RB1 intact vs deficient are secondary if sufficient profiles exist.
2. Rebuild only the **required molecular-stratum disease contrast and regulator scores** using existing Steps 4/6 functions as needed, carefully distinguishing patient strata from previously scored broad cohorts. No silent all-SKCM/all-LUAD substitution.
3. Compare the **same drug** across positive and negative strata, with identical scoring populations. Report raw M/R, percentile gaps, valid score coverage and evidence types.
4. Validate on independent drug sensitivity where matched compound/line/genotype data exist (PRISM, GDSC or CTRP). Use compound-by-line response and drug/class-aware resampling, not tumor-normal expression as a surrogate for functional sensitivity.
5. **Label the independence honestly.** A genotype holdout is new on the patient side only: the drug-side evidence (e.g. RAF inhibitors in A375) is the same set of profiles D127 and the earlier development work already read. Report these as "new patient strata, seen drug profiles", never as fully independent validation. The external screen follows the A0 prior-use audit: a screen whose earlier use touched cohort, compound or parameter choice is reported as development evidence, not as an independent label.
6. Keep functional data completely outside parameter fitting and cohort selection. Report cases where M or R cannot predict sensitivity (e.g. PARP/HRD synthetic lethality) as model-domain limitations, not discarded outliers.

**Outputs:** `genotype_holdout_instances.tsv`, `genotype_holdout_summary.tsv`, `external_viability.tsv`, `external_joins.tsv`, `holdout_manifest.json`.

**Pass gate:** truly unseen outcome labels, molecularly justified positive/negative classes and sufficient sample sizes; all borrowed/mismatched lines flagged. If these conditions are not met, the study remains exploratory.

## A6 - Robustness, interpretability and disposition

**Files:** `scripts/d128/a6_robustness.py`, `scripts/d128/report.py`, `docs/REPORT_target_mechanism_scoring.md`.

**Actions**

1. Compare ARACNe versus CollecTRI; OmniPath versus SIGNOR; hop limits 0/1/2/3; directly measured landmark versus inferred BING genes; raw versus generic-removed profiles.
2. Remove common proliferation/cell-death genes; assess whether M or R is dominated by a few hubs, target families or generic stresses.
3. Repeat per-line, leave-one-line-out and dose/time diagnostics; record mechanisms whose sign reverses across context.
4. Run signed-edge, degree-matched-target and gene-label nulls. Quantify whether extra complexity beats a simple global gene-reversal baseline and the actual production lead.
5. Publish a drug mechanism explanation table: drug/target unit, cohort, pathway, signed expected gene, observed response, disease direction, weighted contribution, assay origin and why a case abstained.
6. Issue an explicit **promote**, **keep as side diagnostic**, **needs better data**, or **reject** decision. Production modifications require a separate approval/PR and a frozen validation threshold.

**Outputs:** `sensitivity_arms.tsv`, `nulls.tsv`, `ablation.tsv`, `summary.json`, decision report, plots and manifest.

**Pass gate:** results survive relevant perturbations and provide incremental independent evidence. Failed controls or low coverage block adoption, regardless of a good pooled headline statistic.

# 5. Cohorts and test order

These are **candidate cohorts**, not a claim that all required drug signatures, genotypes and matching source lines are available. The A0 audit decides actual executability without consulting results.

| Wave | Cohort comparison | Mechanism | Positive evidence and hard negative |
|:--|:--|:--|:--|
| A1 | BRCA-LumA/LumB vs BRCA-Basal; MCF7 context | Direct ESR1-target transcriptional effects | SERMs/SERDs where scored; basal and unrelated cytostatic controls |
| A2 | BRCA-Her2 vs BRCA-Basal (amplification status checked) | ERBB2 receptor -> signaling -> TF | Eligible ERBB2 inhibitors versus nonamplified cases |
| A3 | BRAF V600 SKCM vs WT/NRAS SKCM | Driver-specific RAF/MAPK response | Eligible RAF inhibitors; paradoxical wild-type/RAS-active controls |
| A4 | EGFR-mutant LUAD vs WT/KRAS LUAD | EGFR driver -> downstream TF | Eligible EGFR TKIs; resistant or target-negative contexts |
| B1 | THCA BRAF-like-classical vs RAS-like-follicular | Subtle BRAF-vs-MEK distinction | **D127 failure case:** pooled BRAF/MEK 3/9; split drug classes |
| B2 | RB1-intact vs RB1-deficient breast | CDK4/6 pathway bypass | CDK4/6 drugs and resistance controls |
| B3 | BRCA/HRD-high vs HR-proficient ovarian | PARP synthetic lethality | PARP drugs; low transcriptomic R may coexist with efficacy |
| C | AR-driven prostate; ALK-fusion lung | Lineage and fusion requirements | Only if valid target, genotype and LINCS coverage |

**Data check made before A0 (2026-11-15, DepMap 24Q4 hotspot calls and the 79c profile index; A0 confirms and supersedes this):**

| Wave | In-context LINCS lines | Expected A0 status |
|:--|:--|:--|
| A1 | MCF7 (PIK3CA): fulvestrant, tamoxifen, raloxifene, bazedoxifene, toremifene, clomifene | Executable. No estradiol in MCF7 (no agonist sign-flip control) |
| A2 | SKBR3: lapatinib, neratinib, afatinib, canertinib, erlotinib | Executable if A0 confirms ERBB2 amplification in SKBR3 |
| A3 | A375, SKMEL5 (BRAF): vemurafenib, dabrafenib, PLX-4720, SB-590885; KRAS-mutant A549, HCC515 for the paradox negative | Executable on the drug side; the patient strata need MC3-based Step 4/6 rebuilds |
| A4 | No EGFR-sensitizing-mutant line (no HCC827, PC9 or HCC4006). H1975 (L858R + T790M, resistant to first-generation TKIs) has 3 profiles (afatinib, neratinib) | **Expected not computable in context.** EGFR positives are borrowed or mismatched (E16), so they are reported only and cannot pass a gate |
| B2 | MCF7 (RB1-intact): palbociclib | Needs an RB1-null LINCS line, which A0 looks for |
| B3 | MCF7: olaparib, rucaparib, talazoparib, veliparib | Descriptive (E11) |

TCGA genotype inputs on disk: `data/raw/tcga_pancan/` (MC3 MAF, GISTIC2 thresholded copy number, TCGA subtypes).

**First biological execution:** start with A1 (MCF7/ESR1) to establish basic target anchoring, then A2 for one-hop receptor signaling. Proceed to A3/A4 only when genotype-specific patient and source-line evidence has passed A0. B1 is a **deliberate failure/stress case**, not the primary calibration cohort. Do not collapse selective RAF and MEK inhibitors into one assumed identical mechanism.

# 6. Edge-case test inventory

These cases must be represented by deterministic fixtures *and*, where data exist, an experimental audit row. A fixture verifies software correctness; it does not validate biological truth.

| ID | Failure mode | Required implementation behavior |
|:--|:--|:--|
| E01 | TF target is absent from top changed VIPER regulators | Zero-hop anchor still expands correct TF targets |
| E02 | Target mRNA rises after kinase inhibitor | Do not interpret this alone as failure of kinase inhibition |
| E03 | Activating versus repressing TF edges | Correct sign composition through target -> TF -> gene |
| E04 | Two shortest paths produce opposite signs | Mark ambiguous, exclude from resolved M, retain coverage/attributions |
| E05 | Shared gene appears in many regulons | Count it once in aggregate; preserve regulator-level views |
| E06 | Multiple targets/complex components | Score units separately; no arbitrary best-unit cherry-pick |
| E07 | Unknown, mixed or contradictory MoA | Explicit unknown/ambiguous state; no default inhibitor sign |
| E08 | No route, too few targets, no observed genes | `NA` plus specific abstain reason, not numerical zero |
| E09 | BRAF-wild-type, RAS-active RAF paradox | Flag context; avoid false confident predicted inhibition |
| E10 | RB1-loss CDK4/6 bypass | Do not infer phenotype solely from CDK4/6 expression |
| E11 | PARP synthetic lethality | Distinguish mechanistic activity/viability from normal-state rescue |
| E12 | Generic cytostatic or cytotoxic drug | Compare stress-matched negative drugs and proliferation ablation |
| E13 | BING inferred genes dominate footprint | Report landmark-only and feature-space-specific results |
| E14 | Misjoined gene IDs, TF network or binary layout | Hard fail with diagnostics; no silently shuffled columns |
| E15 | Dose/time profiles disagree | Report stratum variability; avoid hiding sign switch in mean |
| E16 | Drug measured only in genotype-mismatched line | Confidence flag or abstain; never label as direct biological evidence |
| E17 | Very high-degree target or regulator | Degree-preserving null and maximum-hub-contribution report |
| E18 | Target alteration present in only a minority of tumors | Report carrier fraction and stratify; do not claim cohort-wide mechanism |
| E19 | Small/unequal cohort comparator populations | Explicit denominator, matched compound keys, paired uncertainty |
| E20 | No independent outcome or circular gene-level VIPER evidence | Separate exploratory consistency from orthogonal validation |

## 6.1 Test fixtures and assertions

Recommended test module layout:

- `tests/d128/test_io.py`: byte size/orientation, gene ordering, duplicate/unknown IDs, badly mapped TF symbols, missing network.
- `tests/d128/test_graph.py`: direct TF, inhibitor-to-activator, repression, cycles, equal-length opposite paths, shortest-only rule, multi-target complex.
- `tests/d128/test_scores.py`: exact hand-calculated M/R, boundedness, sign flip, neutral vs missing, shared-gene uniqueness, confidence scaling invariance.
- `tests/d128/test_abstain.py`: missing action, no paths, insufficient unique genes, inferred-only genes, unknown genotype, empty matrices.
- `tests/d128/test_aggregate.py`: line-weight identity, borrowed-vs-direct, one-unit and multiple-unit aggregation, missing-unit coverage.
- `tests/d128/test_nulls.py`: permutation reproducibility, degree preservation, null scoring does not alter real target set.
- `tests/d128/test_d127_compat.py`: exact D127 `(compound, high, low)` identity match and original production baseline columns.

pytest is not yet in the `repurposing` environment. Installing it (conda-forge, pinned and added to `environment.lock.yml`) is part of A0. At minimum each test includes a synthetic expected outcome, a fixed random seed where stochastic, and a failure message that names the dataset key being misinterpreted. Pipeline smoke tests should run on a tiny on-disk fixture without R, LINCS downloads, GitHub credentials or a GPU.

# 7. Evaluation design and decision thresholds

## 7.1 Baselines and required ablations

Compare **like for like** on the same compound/cohort set:

1. Current production **`layered_pct`** (including current D124b essential/proliferation behavior).
2. Existing `target_role_pct` (D120) where defined.
3. Existing VIPER-based `spec_mech_pct_pregate` and original reversal.
4. Global gene-signature reversal using the same specific LINCS profiles but without target anchoring.
5. Target-only signed mechanism concordance **M**.
6. Target-local directional rescue **R**.
7. M and R together **as a diagnostic panel**; only later, if justified, one frozen integration with the production score.
8. Degree-matched/size-matched target swaps, sign permutations, and stress-matched compounds.

All variants should report **coverage**; comparing a strict new method's 20 scorable drugs with the old score's 200 without accounting for missingness is not a valid ranking comparison.

## 7.2 Statistics

- Primary descriptive development outcome: per-pair correct high > low, median percentile gap, ties, eligible count; D127's 49 preselected instances are **already seen**.
- New holdouts: prespecify paired same-drug molecular contrast, drug/class-aware confidence intervals, and one primary test *before* score outputs are opened.
- Avoid pseudo-replication: compounds within the same mechanism class, overlapping cohorts and shared genes are not independent observations. Prefer class-balanced macro-averages and bootstrap by compound plus mechanism/cohort where feasible.
- Mechanistic control outcomes: resolved path coverage, sign concordance, specificity above degree- and stress-matched nulls, and direction consistency in the **relevant source-line genotype**.
- External functional screens: drug sensitivity effect size, AUC/AUPRC where labels support it, line/genotype-stratified analyses and confidence intervals. Chronos dependency used as a target-role input is **not** an independent validation label for the same target claim.
- Transfer quality: optionally reference the 77/77b held-out pair benchmark; always distinguish cross-line correlation from drug efficacy.

## 7.3 Proposed promotion gate

Do **not** hard-code a numerical improvement threshold as though it was already validated. Stage A0 may set coverage floors based on input completeness (without opening outcomes); the final independent holdout design must then freeze an effect-size bar, uncertainty criterion, and cohort list. The following qualitative conditions are mandatory:

- Known high-information controls have interpretable target-linked footprints; failed ESR1 detection **by M on the primary arm** (Section 0) despite adequate data blocks promotion.
- The landmark-only arm agrees in sign with the primary (all-gene) arm on the controls (P4).
- Context-specific negatives (RAF paradox, target absence and resistance where measurable) do not rise merely from broad proliferation suppression.
- There is demonstrable value **beyond** target annotations and the actual production `layered_pct`, on genuinely independent molecular/functional data.
- Results are not driven solely by inferred genes, the largest signaling hub, one source line, one dose, one drug class or one permissive prior.
- Every unscorable cohort/drug is reported; production scoring is untouched unless a separate approval is granted.

Possible verdicts: **promote as side metric**, **study further**, **insufficient evidence**, or **reject**. Adoption into `79t3` is a distinct later decision.

# 8. Implementation sequence, dependencies, and review boundaries

| Milestone | Work packages | Can proceed in parallel? | Deliverable / review |
|:--|:--|:--|:--|
| M0 | A0 environment, manifest, schema and coverage | Adapters can be divided by input source | Data inventory + locked feasibility table |
| M1 | A1 signed graph and regulon exporter | Alongside A2 scoring unit fixtures after schema freeze | Validated graph index; path explanations |
| M2 | A2 M/R engine and nulls | Null-generation and score engine independently | Synthetic tests + small smoke run |
| M3 | A3 baseline joins and source-line aggregation | Starts as soon as M/R interface is stable | Cohort scores + integrity checks |
| M4 | A4 mechanism controls and D127 replay | Controls can run separately from D127 compatibility | Diagnosed failures and seen-data comparisons |
| M5 | A5 genotype holdouts and functional screen joins | Genotype definition and screen ID mapping can develop in parallel | Frozen independent outcome comparison |
| M6 | A6 robustness, report, disposition | After M4 and sufficient M5 evidence | Decision report; optional later score-integration PR |

**Critical path:** input access and integrity -> footprint graph/coverage -> score correctness -> baseline-equivalent cohort aggregation -> mandatory controls -> independent holdout -> decision. GPU access is not necessarily the bottleneck; RAM, source files, correct mappings and experimental annotations matter more for sparse scoring. Use a GPU only if profiling demonstrates a benefit or for a later learned cross-line predictor.

**Commit sequence (main branch, the project's standing workflow; replaces the earlier PR decomposition):**

Work is committed directly to `main` and pushed after each step, then `main` is merged into `celligner` and pushed. Each block below is one or more commits, each preceded by its `D128-<stage>` entry in `DECISIONS.md` (Section 4, per-stage routine):

- A0: `scripts/d128/config.py`, `io.py`, `a0_audit.py`, manifests and input fixtures (no scoring); then the D128-A0 freeze (tag `d128-a0-freeze`).
- A1: network exporter, signed graph, `a1_build_paths.py` and graph tests.
- A2: M/R scoring, nulls, per-gene explanations, score tests.
- A3-A4: aggregation, baseline comparison and mandatory controls (ESR1 first).
- A5-A6: molecular and functional holdouts, robustness and final report (tag `d128-eval-freeze` before any holdout outcome is read).

Do not fold outcome-specific tweaks back into earlier methodological commits. A change after a freeze is a new, dated `DECISIONS.md` amendment (as D125-A0.1 did), never an edit of the frozen entry. Adoption into production (`79t3`) is a separate decision with its own entry and is outside this plan.

# 9. Command-line interface and output tree (proposed)

The following are **planned interfaces**; they are not executable commands in the repository yet.

All values come from `scripts/d128/config.py`; scripts take no configuration flags except the thread count. Run from the repo root with the `repurposing` environment, keeping total threads < 32.

```bash
PY=~/miniforge3/envs/repurposing/bin/python

# Initial audit, no new efficacy labels
$PY scripts/d128/a0_audit.py --threads 8

# Build/cached target-to-regulator-to-gene graph
~/miniforge3/envs/repurposing/bin/Rscript scripts/d128/export_edges.R
$PY scripts/d128/a1_build_paths.py --threads 8

# Mechanistic scores and source-line aggregation (only after the d128-a0-freeze tag)
$PY scripts/d128/a2_score.py --threads 16
$PY scripts/d128/a3_aggregate.py --threads 8

# Already-seen development controls (NOT independent validation); ESR1 first
$PY scripts/d128/a4_controls.py --threads 8
$PY scripts/d128/a4_d127_compat.py

# Only after the d128-eval-freeze tag
$PY scripts/d128/a5_genotype_pairs.py --threads 8
$PY scripts/d128/a5_external.py
$PY scripts/d128/a6_robustness.py --threads 16

# Offline software tests
$PY -m pytest -q tests/d128
```

Expected output structure:

```text
results/lincs/d128/
  run_manifest.json
  a0/coverage.tsv
  a0/eligibility.tsv
  a0/cohort_feasibility.tsv
  a0/unmatched_ids.tsv
  a1/target_tf_paths.parquet
  a1/target_gene_footprints.parquet
  a1/path_ambiguity.tsv
  a2/scores_by_line.parquet
  a2/gene_contributions.parquet
  a2/regulator_contributions.parquet
  a2/target_units.parquet
  a2/nulls.parquet
  a3/scores_by_cohort.tsv
  a3/line_heterogeneity.tsv
  a3/abstentions.tsv
  a4/control_results.tsv
  a4/d127_seen_pairs.tsv
  a5/genotype_holdout_instances.tsv
  a5/external_viability.tsv
  a6/sensitivity_arms.tsv
  a6/ablation.tsv
  a6/summary.json
  a6/decision_report.md
```

## 9.1 Minimum score-table schema

One row per eligible `(group_id, pert_id)` (optionally `stratum`), with the compound and cohort evidence metadata:

`group_id, stratum, pert_id, cmap_name, target_units_scored, target_units_total, M, R, M_status, R_status, abstain_reason, n_linked_genes, n_sign_resolved, n_landmarks, n_inferred, n_anchored_tfs, sign_ambiguity, hub_fraction, n_source_lines, direct_or_borrowed, line_match, transfer_confidence, target_genotype_status, layered_pct, target_role_pct, spec_mech_pct_pregate, manifest_hash`.

Separate wide explanations are insufficient: export the per-gene and per-path evidence with explicit signed expectations, measured z-score, disease direction, confidence and source provenance.

# 10. Operational reproducibility and safeguards

- **No writes to production:** new code is read-only with respect to upstream files; everything goes under the D128 output root. Fail if the output root points into production score filenames.
- **Determinism:** explicit seed for permutations, sorted IDs, pinned source checksums and a manifest describing cohort selection and frozen thresholds.
- **Privacy/compliance:** do not commit any restricted patient-level data, server credentials, SSH keys or large raw matrices to GitHub. Commit only code, aggregate summaries permitted by source terms and small nonidentifying fixtures.
- **Logging:** every stage logs elapsed time, matrix shapes, coverage/abstention counts and its upstream manifest fingerprint.
- **Idempotency:** rerunning on identical manifest yields identical files/summary; changed source hashes invalidate dependent caches.
- **Diagnostics before speed:** fail loudly on row-order errors; run a small source-line smoke test, then scale by chunking. Avoid a silent full RAM load of multigigabyte arrays.
- **Compute budget:** fewer than 32 threads in total across concurrent jobs (half the server's 64), set explicitly (`--threads`, `N_CORES`, `OMP_NUM_THREADS`).
- **What is committed:** code, `DECISIONS`/`PROGRESS`/plain-language docs, manifests and small summary tables (`.tsv`, `.json`). Large per-gene or per-path tables (`*.parquet`, `*.tsv.gz`) stay on the server: `*.gz` is already gitignored, and A0 adds `results/lincs/d128/**/*.parquet` to `.gitignore`.
- **Review cadence:** inspect the coverage audit before thresholds; inspect synthetic tests before biology; inspect mandatory controls before interpreting unseen holdouts; only then consider a production integration proposal.

# 11. Immediate next actions

1. ~~Pull the D128 spec branch~~ Done 2026-11-15: the spec is merged into main, and the input files named in Section 1 were confirmed on the server.
2. **User:** decide D0.1 (reference for breast R) and D0.2 (SERMs in the ESR1 control), and confirm or amend the proposed rows of Section 0.
3. Implement and run **A0 only** (the code may be written while step 2 is open), exporting the complete cohort-by-drug-by-line target feasibility matrix, plus explicit data gaps and the label-blind quantities in A0 action 6.
4. Review A0's measured coverage and source-line molecular matches. Write and push the D128-A0 entry, freezing every Section 0 row, then tag `d128-a0-freeze`, *without opening held-out outcome labels*.
5. Build the **ESR1 zero-hop footprint** and its synthetic test suite first; use that to validate graph semantics before adding receptor/kinase pathways.
6. Only after the A1/A2 software gates pass, execute the MCF7 control (M-gated; Section 0) and then the staged cohort benchmarks.

**Success at the end of this project is not a visually compelling ranked drug list; it is a calibrated, interpretable and falsifiable answer to whether target-anchored regulatory footprints improve mechanism fidelity and biomarker-specific drug prioritization over the existing layered score.**

# References and project links

- [D128 scientific specification (PR #2)](https://github.com/kirtihargunani-cmu/repurposing/pull/2)
- [Original LINCS extrapolation and activity-transfer extension (PR #1)](https://github.com/kirtihargunani-cmu/repurposing/pull/1)
- [Subramanian et al., Cell 2017: LINCS L1000](https://doi.org/10.1016/j.cell.2017.10.049)
- [Hatzivassiliou et al., Nature 2010: paradoxical RAF activation](https://doi.org/10.1038/nature08833)
- [Chapman et al., NEJM 2011: vemurafenib in BRAF V600E melanoma](https://doi.org/10.1056/NEJMoa1103782)
- [Mok et al., NEJM 2009: EGFR-mutant lung cancer/gefitinib](https://doi.org/10.1056/NEJMoa0810699)
- [CARNIVAL causal network inference](https://pmc.ncbi.nlm.nih.gov/articles/PMC6848167/)

