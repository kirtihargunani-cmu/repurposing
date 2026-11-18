# Research Specification: Causal Fault Isolation, Network Controllability, and Disease-State Transitions

**Status:** Draft for experimental review; no new experiments performed.
**Repository:** kirtihargunani-cmu/repurposing
**Created:** 2026-11-17
**Scope:** Three independently falsifiable investigations using existing tumor transcriptomics, gene-regulatory networks, LINCS drug and genetic perturbations, and drug-target annotations.
**Primary deliverable:** Validated, context-specific mechanistic hypotheses and explicit limits of inference—not a numerical drug-suitability ranking.

## 1. Scientific motivation

The pipeline maps disease and drug perturbations into a shared molecular feature space:

**Disease:** tumor-versus-reference differential expression -> regulatory/pathway representation -> candidate molecular causes and dependencies.

**Drug:** target annotation and binding evidence -> observed LINCS perturbation -> regulatory/pathway consequences -> comparison with disease-associated abnormalities.

The limitation is that transcriptional agreement does not establish a drug's causal mechanism or therapeutic suitability. A drug can broadly oppose a disease signature without acting on the mechanisms that maintain the pathology. Conversely, an effective molecular intervention might produce only a small, transient or non-transcriptional response.

Three complementary research questions follow:

| ID | Field | Research question |
|---|---|---|
| **F** | Industrial causal fault isolation | Which upstream molecular causes best explain the disease-associated transcriptional abnormalities? |
| **O** | Observability and controllability | Which molecular processes are measurably identifiable, and which are responsive to available interventions? |
| **S** | Ecological restoration and nonlinear dynamics | Can an intervention alter a process that maintains the pathological state, rather than merely make a snapshot look normal? |

The experiments are independent. A success in one is not evidence of success in another. No composite drug score is proposed.

## 2. Existing repository results and constraints

### 2.1 The present mechanism trace is not validated

The seven-link trace tested two intended positives and five negative controls using cohort-matched ARACNe. **Three of five negative controls closed (60% in this small set), and neither positive trace closed.** The network-specific gene footprints were also unstable: for the ESR1 case, the median pairwise Jaccard overlap across 41 ARACNe networks was **0.011**. These results bind this proposal: a plausible signed network route is **not** sufficient evidence of target-specific mechanism.

Sources: docs/PROCEDURE_mechanism_trace.md (D137–D138); results/lincs/trace_aracne/false_positive_rate.json, network_stability.tsv and network_comparison.tsv.

### 2.2 Measurement correlation and null calibration

D130 showed that independent-gene permutation nulls were too narrow; the revised nulls must preserve dependence. D131 and D133 demonstrated substantial covariance among regulator activities and limited effective variance dimensionality. Correlated genes or TFs cannot be counted as independent replications.

Sources: docs/notes/settled_measurement_facts.md; results/lincs/d130/, d131/, d133/.

### 2.3 Perturbational evidence exists, but usable scope is narrower than annotation scope

D135's feasibility audit identified **2,350 CRISPR-matched (drug, target, line) cells** over **7 lines**, but only 46–121 cells survived stringent illustrative exclusions, depending on single-target versus two-target definitions and perturbation strength. D134 found **71 target units** with eligible same-line, between-drug target-recovery comparisons. Neither feasibility audit proves target concordance or retrieval.

Sources: results/lincs/d135/d135_s0_run.json; results/lincs/d134/d134_s0_run.json.

### 2.4 Dose and time are available, but are not equivalent

D136 documented **8,778 compound–line cells with >=3 distinct doses** and **2,736 with >=2 timepoints**. These support dose-response and limited temporal questions, not automatic identification of stable attractors or post-treatment persistence. D136 remains a proposed analysis, not a demonstrated biological result.

Source: docs/SPEC_dose_response_bmd.md.

### 2.5 Other project-specific restrictions

- Do not reintroduce a claim that inferred BING genes are categorically unreliable. D129 measured strong held-out reconstruction. The narrow exception is using a gene's *own inferred transcript* to QC a perturbation of that gene; that self-QC is invalid.
- The production, TAS-filtered drug matrix contains **73 cell lines**. Do not substitute a pre-filter count.
- ARACNe-only is the current regulon decision; independently specified alternative network representations may be used *as sensitivity arms*, not cherry-picked to obtain a favorable result.
- Avoid post hoc removal of genes uniquely or multiply attributed to TFs, or recomputing a drug score after deleting a major downstream regulator. The current trace procedure treats transcriptional response holistically; regulatory concentration is described, not subtracted.
- Freeze test procedures before inspecting their holdout outcomes; label all existing D127/D128 controls as **seen/development**.

## 3. Shared experimental representation

### 3.1 Molecular measurements

Let x in R^G denote a gene-expression observation. For cohort c, a disease-associated observational contrast is

    Δx_c = x_disease,c - x_reference,c.

This contrast can reflect drivers, consequences, cell identity, immune/stromal composition, and reference differences. It is **not** automatically a vector of causal therapeutic targets.

### 3.2 Program representation

Let z = W x represent K regulatory or pathway program activities. W may be derived from the existing regulator activity pipeline, a fixed curated program set, or an independently validated low-dimensional representation. Alternative representations remain separate experimental arms; no assumption is made that a cohort-specific ARACNe graph is ground-truth causal wiring.

Use the existing gene features and VIPER measurements with source network, gene IDs, reference, uncertainty, and cohort fully recorded. Avoid interpreting numerical activity scaling without accounting for D131's covariance normalization issue.

### 3.3 Drug intervention

Define a drug intervention record as u_d = (T_d, A_d, P_d), including annotated molecular targets, type/sign of pharmacological action, and potency/exposure evidence. Retain compound, line, concentration and timepoint identifiers for the response Δx_(d,l,q,τ).

Separate **binding affinity (Kd/Ki)**, **functional IC50/EC50**, **transcriptional BMD**, **unbound exposure**, and **therapeutic effect**. A ChEMBL pChEMBL value alone does not identify which of these was assayed; do not substitute an IC50 for a Kd in an occupancy equation without justification.

### 3.4 Primary data and reuse

| Data | Existing project interface |
|---|---|
| Disease differential expression | Step 4 tumor–reference contrast files |
| Disease master regulators | Step 6 VIPER outputs |
| Regulatory networks | ARACNe; signed signaling priors where applicable |
| LINCS compound response | Level-5 profiles; 79c aggregated data only when appropriate |
| Genetic perturbations | D135 CRISPR and shRNA feasibility/paired data |
| Target activity/availability | D120/D121 and DepMap |
| Genomic disease context | TCGA MC3/GISTIC and molecular strata |
| Pharmacology | D105/D121 ChEMBL target units, D123 exposure |
| Line matching | Celligner weights and source-line genotype/context |
| Prior evaluation | D127–D139 and trace results |

All joins must be audited at the correct grain. Unavailable or biologically mismatched data must remain explicitly **not evaluable**.

## 4. Experiment F — Causal fault isolation

### 4.1 Question and hypothesis

**Question:** Can we distinguish upstream disease-maintaining mechanisms from downstream molecular abnormalities?

**Hypothesis F-H1:** A constrained causal-network method predicts the effects of independent perturbations of candidate disease drivers better than differential-expression magnitude, regulatory centrality, or the D128 signed-path baseline.

The method must not pass merely by reproducing known ChEMBL or TCGA annotations supplied as inference constraints.

### 4.2 F1 — Candidate causes

For predefined cohorts, extract candidates with supporting genomic or functional evidence:

- Activating driver variants, with variant-specific evidence (not generic missense status).
- Functionally consequential amplifications/deletions, with subgroup prevalence.
- Target dependencies in molecularly relevant cell lines.
- Regulatory programs with independently observed activation.

Retain patient denominator, event frequency, source, effect direction, and missingness. Expression alone identifies abnormality, not cause.

### 4.3 F2 — Signed causal-prior graph

Construct a directed, signed graph with node identity, positive/negative/unknown action, provenance, evidence level, tissue context and ambiguity. Use the existing OmniPath/SIGNOR and ARACNe inputs as **hypothesis priors**; distinguish curated direction from co-expression.

Explicitly flag high-degree hubs, contradictory paths, latent/unmeasured mediators and context-dependent edges.

### 4.4 F3 — Candidate inference

Evaluate, on the same eligible cases:

1. Differential-expression/target-abundance baseline.
2. Regulator-activity and centrality baseline.
3. Existing D128 path consistency and the currently failed trace as controls.
4. A constrained causal explanatory model, e.g. CARNIVAL-style signed optimization.

Use a complexity penalty or other *label-blind* structural constraint; do not choose the most favorable route after observing drug results.

### 4.5 F4 — Freeze forward intervention predictions

For every proposed upstream cause, record before reading its held-out intervention response:

- Expected affected regulators and direction.
- Expected whole-profile or prespecified program-level response.
- Alternative mechanistic explanations.
- Unmeasurable consequences and reasons.
- Perturbation identity and matching conditions.

Commit a prediction artifact, checksum and evaluation plan before evaluating its observations.

### 4.6 F5 — Independent validation and controls

Primary validation: held-out genetic or pharmacological perturbations in genotype-appropriate lines. The perturbation that validates the prediction must not have been used to infer the candidate route.

Required controls: strong differential expression without dependency, incompatible genotype, degree-matched alternative targets, randomized **drug–target–line** tuples, matched cytotoxicity/proliferation, and alternative supported network representations.

The D137 random-cohort control was confounded by preserving a pharmacologically compatible line. Do not repeat that mistake.

### 4.7 F endpoints

- Accuracy/coverage of intervention direction predictions on held-out observations.
- Mechanism recovery conditional on what was supplied versus withheld.
- False-positive/false-negative rates for biologically meaningful negative controls.
- Number and diversity of equally plausible causal explanations.
- Network-choice sensitivity and rates of abstention.
- Context-specific estimates with cluster-aware confidence intervals.

**Failure:** The method merely tracks degree, class frequency or original driver annotations; produces plausible explanations for negative controls; changes arbitrarily with the network; or does not improve out-of-sample intervention predictions.

**Deliverable:** One disease fault-isolation record per cohort and candidate driver, with predicted consequences, observed test, competing causes and verdict.

## 5. Experiment O — Observability and controllability

### 5.1 Question and hypotheses

**Question:** Which internal molecular processes are distinguishable in the assay, and which can be modified by available drugs in the tested context?

Use the local model

    dz/dt = A z + B u,
    y = C z + ε,

where z represents program state, A its dynamics, B intervention influence, and C the observation mapping. The repository does not currently identify these matrices. **Graph reachability is not proof of controllability.**

Hypotheses:

- **O-H1, empirical observability:** Target-associated interventions produce distinguishable transcriptional profiles.
- **O-H2, intervention responsiveness:** Independently measured target perturbations predict an appropriate component of the drug response.
- **O-H3, response specificity:** The apparent response is not explained solely by generic stress, chemical analogs or proliferation.

### 5.2 O1 — Target observability

Reuse D134 target recoverability. Ask whether a compound's target is more recoverable from its profile than chance **within the same cell line**. Preserve target-frequency-aware nulls, held-out compound identity, Morgan/Tanimoto chemical similarity stratification and expression-strength controls. Singleton targets remain outside the between-drug retrieval test, not biologically unobservable.

This is **empirical target identifiability**, not the formal observability of the complete cellular state.

### 5.3 O2 — Drug versus genetic perturbation

Reuse D135's same-line drug–CRISPR paired units as the primary experiment; consensus shRNA is a separate secondary arm. Drug inhibition need not be equivalent to target knockout when scaffolding roles, partial inhibition, feedback, complex membership or pharmacological selectivity differ.

The assay's best-inferred gene estimates may be used for global profiles. The target's own inferred expression must **not** be used to confirm its knockout/knockdown.

Freeze the drug-action expectations and test against strength-matched alternative target knockouts, not independent-gene permutations.

### 5.4 O3 — Concentration-response evidence

Reuse D136's full Level-5 dose series rather than 79c's collapsed drug-line means.

Ask whether target-associated programs respond reproducibly and directionally as dose increases, with gene-level replicate noise, censoring and non-monotone curves handled explicitly. Keep the gene-level, regulon-level and target-level evidence separated.

Timepoints are not interchangeable. A response measured at 6 h may differ mechanistically from 24 h.

### 5.5 O4 — Disease-program responsiveness

For an intervention with independent target evidence, test effects on a **prespecified disease-relevant program**, with the program defined without access to the intervention outcome.

Compare with matched off-targets and controls. Do not remove major downstream effectors merely because they dominate the response.

### 5.6 O5 — Formal control analysis, only where identifiable

Where dense, independent time-series interventions exist, fit small local dynamical models with explicit observation and input operators. Only then consider controllability Gramians or finite-horizon intervention analyses. Report uncertainty in A, B and C. No biological controllability claim may rest on topology alone.

### 5.7 O outputs and gates

Per target × drug × line × context, report:

- Recoverability status, chance-adjusted retrieval and chemical similarity sensitivity.
- Drug–genetic concordance with perturbation quality.
- Concentration-dependent effect with uncertainty and censoring.
- Relevant disease-program response.
- Exposures/affinities and limitations.
- Alternative targets and generic-response explanations.
- Abstention reason, if needed.

**Failure:** Target specificity vanishes against matched nulls; response is inconsistent across relevant conditions; predictions do not generalize; or formal control claims depend on arbitrary/unidentified dynamics.

**Important:** Poor transcriptomic observability does not imply that a target is therapeutically unimportant.

## 6. Experiment S — Nonlinear disease-state transitions

### 6.1 Question and constraints

**Question:** Does intervention change a mechanism maintaining a pathological regulatory state, or only transiently displace its expression measurements?

Cross-sectional tumor–normal expression and a single post-treatment LINCS profile do **not** establish alternate stable attractors. Dose variation is not temporal variation. Persistence requires observation after drug withdrawal or other independent functional evidence.

Run two arms: (i) known-answer simulations, and (ii) biological analysis only if an explicit longitudinal feasibility gate is passed.

### 6.2 Hypotheses

- **S-H1:** Some pathological programs exhibit regulatory feedback or compensatory responses.
- **S-H2:** Interventions at causal maintenance nodes produce temporal responses distinguishable from those at downstream indicator nodes.
- **S-H3:** Selected interventions produce persistent changes after treatment ends; this requires washout/recovery data.

### 6.3 S1 — Independent low-dimensional program definition

Define z(t) in R^K using prespecified, biologically interpretable programs. Preserve distinct proliferation, lineage, microenvironment and target-specific components. Do not define states using the future intervention outcome.

### 6.4 S2 — Synthetic benchmarks

Construct known dynamical systems with:

1. Stable linear response.
2. Negative feedback.
3. Bistable positive feedback.
4. Redundant or compensatory pathways.
5. A downstream feature with greater apparent differential expression than its causal driver.
6. Unobserved latent regulators and non-identifiable parameters.

Test both correctly specified models and models forced to abstain when the data are insufficient.

### 6.5 S3 — Model comparison

Prespecify progressively more complex models:

- Static response baseline.
- Linear state-space dynamics.
- Small Boolean network or switch model.
- Sparse nonlinear dynamics when data permit.

Prefer the simplest model that predicts held-out temporal observations. A complex model must outperform simpler baselines on independent timepoints/replicates, not only fit its training observations.

### 6.6 S4 — Temporal feasibility audit

From Level-5 metadata, inventory cell-line × compound × matched-dose time series:

- Distinct timepoints at comparable doses.
- Biological and technical replication.
- Untreated time-matched controls.
- Viability and cell-composition information.
- Washout or recovery samples.
- Sampling intervals relative to measured response times.

Do not proceed with biological attractor inference from only two arbitrary observations.

### 6.7 S5 — Held-out trajectories

Fit using a predefined subset of timepoints, predict later observations, and report errors by independent compound–line trajectory. Compare with no-dynamics and linear-response baselines. Preserve observational uncertainty and response heterogeneity.

### 6.8 S6 — Persistence and functional interpretation

If post-washout measurements exist, distinguish return to baseline, incomplete return, sustained displacement and evidence of an alternate stable state. In the absence of washout measurements, the verdict is **persistence not evaluable**.

A sustained expression change is not inherently beneficial. Evaluate independent viability, proliferation, differentiation, apoptosis or functional outcomes suitable to the disease.

### 6.9 S failure and deliverable

**Failure/abstention:** Multiple models fit but disagree on intervention effects; poor held-out temporal prediction; inadequate timepoints; viability/composition confounding; persistence inferred without recovery observations; or nonlinear complexity yields no predictive benefit.

**Deliverable:** A context-specific state-trajectory report with model uncertainty, reversibility evidence, compensatory pathways, functional interpretation and explicitly unanswerable questions.

## 7. Integrated interpretation: evidence, not a scalar drug score

| Component | F: Cause | O: Response/control | S: Temporal state |
|---|---|---|---|
| Inference object | Candidate disease-maintaining mechanism | Intervention–program relationship | Dynamical state trajectory |
| Validation | Held-out perturbation direction | Target retrieval, genetics, dose response | Held-out timepoints, washout |
| Key alternative | Downstream association | Generic stress/off-target effects | Transient response or composition |
| Permissible conclusion | Causal *hypothesis* independently supported | Intervention effect supported in specified context | Persistence/transition only if observed |

Illustrative interpretation, NOT a new biological finding:

- F identifies a BRAF V600-associated MAPK driver mechanism in a genetically defined melanoma subgroup.
- O finds genotype-appropriate, concentration-dependent drug effects consistent with BRAF pathway inhibition.
- S observes a treatment trajectory but lacks recovery data.

The combined report should state **mechanism supported; intervention response supported; persistence and clinical benefit not established**. It must not average these observations into a numerical suitability score.

## 8. Incremental implementation

### Stage 0 — Feasibility and frozen test contract

- Audit disease molecular strata, intervention pairs, network choice and genomic context.
- Audit genetic perturbation, matched dose/time, recovery, and independent functional endpoints.
- Freeze what counts as an independent unit, null and endpoint for each experiment.
- Inventory measurement noise and structural missingness without outcome-driven tuning.

**Exit:** Every planned claim has an independent validation data source or a declared abstention.

### Stage 1 — Synthetic-known-answer benchmarks

- Implement small ground-truth causal/dynamical systems.
- Include confounded, latent, non-identifiable and wrong-model examples.
- Test null calibration and graceful failure.

**Exit:** Methods detect true mechanisms where identifiable and do not hallucinate explanations when they are not.

### Stage 2 — F: Disease root-cause inference

- Pilot ESR1/BRCA-LumA, BRAF-V600/SKCM, EGFR-driver/LUAD, ERBB2-amplified/BRCA as development cases.
- Write and commit forward predictions *before* examining target perturbation outcomes.
- Use network-choice and genotype-incompatible controls.

**Exit:** Held-out intervention predictions improve on declared baselines without a rise in mechanistic false positives.

### Stage 3 — O: Empirical observability and intervention response

- Use the D134 and D135 test contracts rather than duplicate their pipelines.
- Incorporate D136 dose series where calibrated and available.
- Keep empirical observability, target concordance, target potency, disease-program response and binding/exposure distinct.

**Exit:** At least a bounded, clearly defined target subset exhibits reproducible specificity beyond matched nulls.

### Stage 4 — S: Biological temporal feasibility

- Audit matched-dose longitudinal observations and washout/recovery data.
- Run known-answer simulations regardless of biological coverage.

**Exit:** If biological trajectories are not identifiable, stop that arm without fitting an unconstrained attractor model.

### Stage 5 — Independent generalization

- Predefine genuinely held-out drugs, mechanisms, cell lines or disease strata.
- Report direct versus borrowed and genotype-matched versus mismatched separately.
- Include non-cancer disease extension only when reference tissues and validation outcomes are appropriate.

**Exit:** Improvement generalizes across independent units and withstands relevant negative controls.

### Stage 6 — Evidence integration

- Populate structured mechanism records with provenance, competing explanations, confidence bounds and known limitations.
- Render scientific prose strictly from the structured evidence.
- Make no production ranking or tier changes without a separate validated proposal.

## 9. Proposed repository structure

    docs/
      SPEC_systems_theoretic_mechanism_validation.md
    scripts/
      systems_validation/
        config.py
        io.py
        nulls.py
        fault_isolation/
          build_candidates.py
          infer_causes.py
          predict_perturbations.py
          evaluate.py
        observability/
          target_recovery.py
          genetic_concordance.py
          dose_response.py
          evaluate.py
        state_transitions/
          synthetic_systems.py
          temporal_audit.py
          fit_dynamics.py
          evaluate_trajectories.py
        integrate/
          evidence_records.py
          render_report.py
    tests/
      systems_validation/
        test_fault_isolation.py
        test_observability.py
        test_controllability.py
        test_state_transitions.py
        test_nulls.py
    results/
      systems_validation/
        feasibility/
        synthetic/
        fault_isolation/
        observability/
        state_transitions/
        integrated/

**Reuse** existing D134–D136 implementations and file contracts instead of forking duplicate code.

## 10. Statistical safeguards

1. **Development versus validation:** D127/D128 cases and trace controls are *seen*; do not rebrand reruns as confirmatory.
2. **Dependence-preserving nulls:** Do not independently shuffle correlated genes. Use appropriate compound/target/line permutations, degree/strength matching and covariance-aware resampling.
3. **No annotation leakage:** A supplied driver label cannot simultaneously serve as independent proof of driver recovery. Hold out the variable corresponding to the primary prediction.
4. **Correct sampling unit:** Genes, correlated TFs and technical replicates are not independent biological trials. Report clustered uncertainty by compound, target class, line, cohort and trajectory as appropriate.
5. **Freeze procedures:** Prespecify splits, baseline definitions, nulls, endpoints and abort conditions before examining held-out results.
6. **Measurement-range calibration:** If cutoffs must be calibrated, select them using predeclared label-blind noise/availability procedures, not by optimizing to successful mechanisms.
7. **Negative controls first:** Declare the full negative-control set and calculate false-positive behavior before interpreting positive mechanistic narratives.
8. **No false cross-platform arithmetic:** Tumor logFC/t statistics and LINCS MODZ z scores are different measurements; do not directly add them to infer post-treatment expression.
9. **No silent data borrowing:** Source-line genotype, tissue and exposure mismatch must appear in each record; unsupported transfer is unscorable.

## 11. Edge-case register

| Case | Interpretation required |
|---|---|
| Large DE but no functional driver evidence | Abnormality is not automatically causal |
| Driver activation without raised mRNA | Preserve genomic and functional evidence |
| Binding affinity with no relevant cellular response | Engagement is unestablished, not automatically absent |
| High-degree pathway hub | Graph reachability is not controllability |
| Contradictory network footprints | Report model disagreement |
| Polypharmacology | Preserve target-specific explanations |
| Cytotoxic whole-profile response | Avoid generic-stress mechanism attribution |
| Weak or unmatched LINCS line | Limit inference to actual experimental context |
| CRISPR knockout differs from inhibitor | Distinct intervention types; no assumed equivalence |
| Perturbed gene's own inferred transcript | Invalid direct knockout self-QC, not grounds to reject all inferred genes |
| Dose-dependent sign reversal | Preserve concentration context |
| One timepoint | Dynamics not identifiable |
| No washout/recovery | Persistence not testable |
| Near-normal transcription with poor viability | Not automatically therapeutic |
| Similar signatures across diseases | Shared causality not established |
| Unobserved driver/latent cause | Competing hypotheses or abstention |

## 12. Expected outcomes

**Best case:** A bounded set of mechanisms passes independent upstream-cause inference, observed target-specific perturbation, dose/time response, and relevant functional validation. The system generates falsifiable repurposing hypotheses rather than universal rankings.

**Realistic case:** Different targets pass different evidence links. Some genomically established drivers are transcriptionally difficult to observe; others yield measurable response without driver evidence. The available longitudinal data may not establish attractor changes. A valuable deliverable is the **domain of applicability**.

**Negative case:** Causal reconstruction adds nothing beyond supplied annotations, drug-specific observability fails, or temporal data cannot identify dynamics. Record this rather than adding an unvalidated score.

## 13. Immediate recommendation

**Start with the shared Stage 0 audit and Experiment F's known-answer synthetic benchmark.** In parallel, continue D134, D135 and D136 as the near-term evidence routes for Experiment O. Delay biological Experiment S until its longitudinal feasibility gate is satisfied.

**First milestone:** For BRAF-V600 melanoma and ESR1-associated breast cancer, demonstrate whether the method can (i) identify plausible upstream mechanisms, (ii) freeze predictions of their consequences, and (iii) distinguish true from convincing but incorrect explanations using independent perturbations.

## 14. Literature foundations

- Yang et al. (2012), *Progress in Root Cause and Fault Propagation Analysis of Large-Scale Industrial Processes*, Journal of Control Science and Engineering.
- Liu et al. (2019), *From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL*, npj Systems Biology and Applications, DOI: 10.1038/s41540-019-0118-z.
- Liu, Slotine and Barabási (2011), *Controllability of complex networks*, Nature, DOI: 10.1038/nature10011.
- Liu, Slotine and Barabási (2013), *Observability of complex systems*, PNAS, DOI: 10.1073/pnas.1215508110.
- Gates and Rocha (2016), *Control of complex networks requires both structure and dynamics*, Scientific Reports, DOI: 10.1038/srep24456.
- Angulo, Moog and Liu (2019), *A theoretical framework for controlling complex microbial communities*, Nature Communications, DOI: 10.1038/s41467-019-08890-y.
- *Network inference from perturbation time course data* (2022), npj Systems Biology and Applications, DOI: 10.1038/s41540-022-00253-6.

## 15. Final objective

Replace the narrow question

> Which compound produces the largest numerical reversal of a disease-associated signature?

with

> Which mechanisms maintain the disease-associated state, which interventions can modify those mechanisms, and what independent evidence supports the predicted consequences?

The intended output is a set of **falsifiable, context-dependent mechanism hypotheses** with explicit measurements, contradictory evidence, reproducibility requirements and limits of inference.

**Draft only.** No new experimental claim, positive mechanism verdict, production score change or drug reranking follows from this document.
