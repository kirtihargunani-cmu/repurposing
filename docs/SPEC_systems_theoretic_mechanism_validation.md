# Research Specification: Causal Fault Isolation, Network Controllability, and Disease-State Transitions

**Status:** Draft for experimental review; mathematical implementation and drug-repurposing-output amendment 2026-11-17; no new experiments performed.
**Audit status (agent, 2026-11-18):** independently audited in `docs/AUDIT_systems_theoretic_spec.md`.
All nine numeric claims in §2 reproduce exactly from their artifacts. Gates **F0 (4/4)**, **O0 (10/10)**
and **S0 (8/8)** now **pass** in the isolated `sysval` environment. One substantive defect found —
causal-ambiguity enumeration requires CPLEX; see §4.8a. Agent-inserted amendments below are marked
**[AGENT AUDIT 2026-11-18]** for review.
**Repository:** kirtihargunani-cmu/repurposing
**Created:** 2026-11-17
**Scope:** Three independently falsifiable investigations using existing tumor transcriptomics, gene-regulatory networks, LINCS drug and genetic perturbations, and drug-target annotations.
**Primary deliverable:** **Falsifiable, evidence-grounded proposals to repurpose existing drugs for a different disease or molecular subgroup, explicitly explaining why the disease target is relevant and how the drug could modify it.** Numerical scores are optional rather than required.

## 0. Primary output: drug repurposing proposals, with explicit why and how

**Research objective.** Starting from the mechanisms implicated in a **target disease**, identify an existing or previously developed drug **with a different established/development indication** whose molecular action could change one of those mechanisms. The required output is a *candidate nomination with a testable mechanism*, not a percentile or universal ranking.

Each proposal must answer:

1. **For what new disease and which patients?** Specify a target disease, molecular subgroup and cohort denominator. Record the fraction of patients with the candidate driver state; do not equate pathway-wide changes with a driver in every patient.
2. **Why this drug's target?** Identify upstream disease-maintaining evidence (functional variant, dependency, causal hypothesis), distinguishing driver effects from downstream, lineage and proliferation signatures. State whether the inferred driver was supplied to or discovered by the model.
3. **How would the drug act?** Name annotated drug target(s), mechanism (inhibition, activation, degradation, etc.), expected signed signaling/regulatory changes and the disease-associated program those changes could alter.
4. **Does the proposed mechanism appear in a relevant experiment?** Compare frozen predictions with LINCS gene-level profiles, **drug-side ARACNe/VIPER regulator activity outputs**, genetic perturbations, dose/time response and target pharmacology where eligible. Verify source-line genotype, tissue and drug exposure, rather than silently borrowing mismatched lines.
5. **Why might this interpretation be wrong?** Report viable alternative targets, off-target/cytotoxic effects, uncertain pathway directions, network-choice instability, context transfer, and contradictory observations.
6. **What would test it next?** State a falsifying target-specific experiment, the expected molecular readout, independent functional outcome, appropriate negative controls and the result that would reject the candidate.

### 0.1 Required structured candidate dossier

The atomic record is **new disease/subgroup × drug × target unit × observed source line**, preserving relevant doses and times. Each dossier requires:

| Field | Content |
|---|---|
| Drug and existing indication | Verified compound identifiers, approved/previously developed indications, and status; avoid claiming a known use is novel |
| Target disease and subgroup | Exact proposed *new* indication, relevant molecular subtype, patient numerator/denominator and reference tissue |
| Disease-driver rationale | Patient alteration prevalence, target activity, expression, dependency and limitations of causality |
| Drug-target hypothesis | Target(s), annotated pharmacological action, assay-specific affinity/potency, selectivity and exposure information |
| Frozen causal mechanism | Upstream target → predicted regulator/program response → proposed disease consequence |
| Drug-side measurement | Original LINCS whole-response arm and tumor-derived **ARACNe/VIPER** drug-response arm; source network/version, line, dose and time |
| Independent tests | Genetic-perturbation concordance, target recovery, concentration response and relevant biological phenotype when available |
| Contradictions and alternatives | All unresolved negative controls, polypharmacology, generic response, tissue mismatch and uncertainty |
| Falsification experiment | Explicit perturbation, controls, predicted direction and failure condition |
| Disposition | Mechanistically supported for experimental follow-up / preliminary hypothesis / contradicted / not evaluable |
| Optional prioritization | A score or ordered shortlist ONLY if calibrated, with stated reference population; not mandatory |

**Rendering order:** proposed drug and target disease → **why the target matters in this disease** → **how the drug might affect it** → what we observed → what conflicts with that interpretation → what to test next. All prose must be derived from structured, source-linked evidence. An impressive transcriptomic correlation cannot override a broken or untestable link.

**No-result is valid.** If the disease-driving target is not established, drug action or target exposure cannot be defended, relevant line evidence is lacking, or the mechanism fails prespecified negative controls, the output is **no supported repurposing proposal**. An unresolved mechanistic idea may be retained as explicitly *preliminary*, but not promoted.

### 0.2 Where each mathematical experiment contributes

| Experiment | Why/how contribution | Required output for a candidate | Primary falsification |
|---|---|---|---|
| **F — Causal fault isolation** | **Why:** distinguish upstream candidate drivers from descriptive expression abnormalities | Candidate disease mechanisms and actionable target(s), drug mapping, frozen downstream predictions | Alternative upstream causes and independent intervention results contradict the predicted mechanism |
| **O — Observability / controllability** | **How:** evaluate drug engagement, observable target-specific effects and intervention responsiveness | Drug → target → observed regulatory/program response, contextual relevance, potency/exposure limitations | Drug fails independent target-specific perturbation/response tests or is indistinguishable from generic stress |
| **S — Nonlinear state transitions** | **What happens after intervention:** test feedback, adaptation, reversal, persistence when temporal data exist | Predicted versus held-out trajectories, compensation and washout behavior if observed | Wrong trajectory predictions, feedback alternatives or absence of data for claimed persistence |

**S is conditional on longitudinal data and is not a hard requirement for every candidate.** No persistence conclusion can be made from a single post-treatment snapshot. These three experiments must be compared to existing gene-signature matching and the failed D137 mechanism trace in their **ability to produce correct, independently testable mechanism explanations**, not simply how many candidates they return.

### 0.3 Drug-side ARACNe implementation arm

The existing drug side uses tumor-derived, cohort-specific **ARACNe regulons applied to LINCS profiles** through VIPER (e.g. scripts/79d_drug_viper.R and data/processed/lincs/drug_viper_specific_v2_tf/). These are estimates of *drug-associated regulator activity*, not independently measured network causal effects and not networks trained directly from drugs.

- For **F**, freeze the causal mechanism and use drug-side regulator activities to test predicted TF/program changes, in parallel with the whole-gene LINCS arm and independent genetics.
- For **O**, use those regulator activity responses as a *primary empirical representation* of how an intervention affects disease-relevant programs, benchmarking against the original gene-level responses and D134/D135 controls.
- For **S**, when real matched time/dose profiles exist, compute VIPER separately at each sampled condition rather than treating aggregated 79d outputs as a dynamical time series.
- Preserve the exact cohort network and version. Given D137–D139 and the ESR1 anti-estrogen discrepancy, no ARACNe-derived mechanistic interpretation is accepted merely because it gives a favorable signed path. Report sensitivity to network choice and the independently tested whole-response prediction.

### 0.4 Validation criteria for a useful repurposing proposal

**Primary objective:** Improve the precision and generalizability of **mechanistically justified, new-indication drug proposals**, assessed against independent molecular perturbation and disease-functional response data with disease/subgroup, drug and target-class leakage controls. A proposal is supported only if the mechanism survives the required negative controls and has genuine independent intervention evidence; without that it is preliminary or not evaluable.

Report (i) nomination coverage; (ii) fraction of proposed pairs with adequately supported target, correct disease genotype and source-line context; (iii) held-out mechanistic prediction accuracy; (iv) false-positive rate on independent negative-control tuples; (v) disease-functional outcome when available; and (vi) number of candidates without enough data. Avoid using the total number of proposed candidates as a success metric.

A numerical prioritization score **may** be studied later to triage laboratory follow-ups, but must be separately calibrated and must not be presented as probability of clinical benefit. The required deliverable is **the drug–new-disease rationale and falsifiable mechanism**.

**Proposed output files:** results/systems_validation/repurposing/candidate_dossiers.jsonl, candidate_summary.tsv, frozen_predictions.jsonl, negative_controls.tsv, validation_summary.tsv, coverage.tsv and run_manifest.json. No production drug rankings are changed under this draft.

> **[AGENT AUDIT 2026-11-18] §0.5 — measured entry-criteria funnel.**
> The new-indication requirement is **feasible**, with this denominator. Computed from
> `results/lincs/known_use/labels.tsv`, where `lab_hub` / `lab_chembl` mark existing indications:
>
> | filter | pairs |
> |---|---|
> | cohort × FDA-ingredient rows | 7,570 (41 cohorts, 391 ingredients) |
> | minus known uses (Hub **or** ChEMBL) | 7,099 (−471) |
> | and gate-eligible | 5,408 |
> | **and source line matched, not borrowed** (§0.1) | **1,341** (266 drugs, 29 cohorts) |
> | and the drug is not already an anticancer | 618 (164 drugs, 20 cohorts) |
>
> Two cautions: coverage is skewed (median **5** candidates per cohort, range 1–132), and
> **12 of 41 cohorts have zero** eligible candidates — SARC-DDLPS/-LMS/-UPS, TARGET-WT, TCGA-ACC,
> -HNSC, -KICH, -KIRC, -KIRP, -THCA, TGCT-NSGCT, TGCT-SEM. For those the correct §0 output is
> "no candidate meets the entry criteria", stated up front. `labels.tsv` is the mechanical
> enforcement resource for "avoid claiming a known use is novel" (471 pairs, 97 drugs).


---

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

The three experiments independently contribute to **nominating, explaining, testing or rejecting drug–new-disease repurposing hypotheses**. A successful mathematical result in one experiment does not compensate for a failed intervention control in another. The final research artifact is a mechanistic candidate dossier, not a mandatory drug score. Quantitative prioritization remains an optional, separately validated research question.

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

### 2.6 Implementation audit findings (agent, 2026-11-18)

> **[AGENT AUDIT 2026-11-18]** Full detail in `docs/AUDIT_systems_theoretic_spec.md`.
>
> **(a) Mathematical gates pass.** F0 (CARNIVAL's own standard + inverse toys, both open solvers;
> vanilla matches the bundled `toy_result_expected_ex1` fixture at node-activity sign agreement
> **1.000**), O0 (**10/10** hand-computable fixtures), S0 (**8/8**; the Angulo invariant is
> symbolically **0** and numerically conserved to **<3.2e-10**). Manifests in
> `results/systems_validation/synthetic/`.
>
> **(b) Dependencies are installed in an isolated environment.** `sysval`
> (`environment.sysval.lock.yml`): CARNIVAL **2.20.0**, coin-or-cbc 2.10.13, r-lpSolve 5.6.23,
> sympy 1.14.0, scipy 1.17.1, networkx 3.6.1, R 4.5.3, Python 3.11.18. The pinned production envs
> `repurposing` and `celligner` were deliberately **not** modified.
>
> **(c) The production lead ranking is insensitive to the regulon.** D139: switching the drug side
> to ARACNe-only moves `rev` by Spearman 0.9905 and `spec_mech_pct` by 0.9866, but **`layered_pct`
> is bit-identical (Spearman 1.0000, all 41 groups)** — it is built from DepMap dependency, TCGA
> alteration and ChEMBL binding, and the regulon enters only via the `O` component where
> `o_source == NES`, **4.9%** of `gene_components` rows. §0's dossier deliverable makes this
> non-blocking, but the dossier↔production-ranking relationship is still unstated: does a supported
> dossier supersede, annotate, or ignore a `layered_pct` rank for the same pair?
>
> **(d) Temporal feasibility, measured.** From `data/raw/lincs2020/siginfo_beta.txt` (720,216
> `trt_cp` signatures, 201,494 compound × line cells): **2,038 cells (1.0%) have ≥3 timepoints**,
> **317** at `is_hiq` & TAS ≥ 0.4, 28 with ≥4. 97% of signatures are 6 h or 24 h. **`siginfo_beta`
> has no washout, recovery or drug-removal field at all.** Persistence (S-H3) is therefore not
> evaluable *by construction*, not by sample size — consistent with §0.2 and §6.9a making S
> conditional. S-math is unaffected and passes.
>
> **(e) Three-way environment split is a live hazard.** `sysval` has pandas 3.1.1 / numpy 2.4.6 /
> Python 3.11; `celligner` has pandas 1.4.4 / Python 3.10 and is the interpreter the drug-side
> pipeline runs under; `repurposing` holds R and has no SciPy. A `79i` run failed on
> `ModuleNotFoundError: scipy` under the wrong interpreter on 2026-11-17. Every
> `scripts/systems_validation/` file must declare its interpreter.

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


### 4.8 F-MATH — Direct implementation of CARNIVAL's integer-programming inference

**Source and implementation boundary.** Implement the published **CARNIVAL** software, rather than relabeling the D128 shortest-path procedure as CARNIVAL. CARNIVAL takes a signed/directed prior-knowledge network (PKN) and measured or inferred transcription-factor (TF) activities and solves an integer-linear-programming (ILP) causal-network-contextualization problem. Its two distinct modes are (a) **inverse CARNIVAL**, with no supplied perturbation target, and (b) **standard CARNIVAL**, conditioned on a specified perturbation target and its activating/inhibiting action. Their scientific claims differ:

- **F-inverse:** recover a plausible upstream cause *without giving the solver that cause*. Validate driver recovery only against withheld, independently annotated driver labels.
- **F-target-constrained:** fix the candidate target and intervention sign; evaluate whether the inferred downstream mechanism predicts an independent perturbation signature. Supplying BRAF as a constraint is **not** evidence that CARNIVAL rediscovered BRAF.
- Disease-side observational VIPER states and treatment-induced VIPER states are **different experimental contexts**. First reproduce the published perturbation-response use case; using tumor-versus-normal activities as CARNIVAL inputs is an explicitly labeled **adaptation**, evaluated separately.

**Published mathematical structure.** CARNIVAL represents candidate node activities and whether edges/nodes are selected using integer decision variables, with sign-consistency and connectivity constraints linking upstream nodes to TF observations. The optimizer trades TF/pathway-activity agreement against solution complexity. For orientation, the following is a **schematic experimental objective, not a transcription of the CARNIVAL source code**:

    min_{binary node/edge decisions, signed node activities}
        L(TF_activity_observed, TF_activity_explained)
        + λ * network_complexity
        + optional prior/activity penalties
    subject to:
        candidate edges ∈ frozen directed/signed PKN
        selected causal interactions obey action/sign constraints
        activated/inhibited nodes must satisfy the package's ILP logic
        imposed source target/sign (standard mode only)

**Do not manually substitute this schematic expression for the library objective.** Record the actual CARNIVAL package version, solver, configuration, TF measurement handling, regularization and termination status. The original implementation discretizes TF signs for ILP constraints and can use continuous TF/pathway activities in its objective. Inspect the installed package API and document the exact transformation; never treat a continuous VIPER NES as a binary sign without a reproducible rule.

> **[AGENT AUDIT 2026-11-18] Installed-API corrections.** Verified against CARNIVAL **2.20.0** in
> `sysval` (the channel version listing suggested 2.8.0 — confirm, as this section already instructs):
> 1. The entry points named below are the **legacy** interface. `runCARNIVAL(inputObj=, measObj=,
>    netObj=)` still exists, but the current API is **`runVanillaCarnival(perturbations, measurements,
>    priorKnowledgeNetwork, weights, carnivalOptions)`** and **`runInverseCarnival(measurements,
>    priorKnowledgeNetwork, weights, carnivalOptions)`**, configured by
>    `defaultLpSolveCarnivalOptions()` / `defaultCbcSolveCarnivalOptions()` /
>    `defaultCplexCarnivalOptions()`.
> 2. **CBC requires an absolute `solverPath`.** `Sys.which("cbc")` returns `""` when `Rscript` is
>    invoked by absolute path (the env `bin` is not on PATH), and CARNIVAL then aborts with
>    *"Invalid path to solver provided"*.
> 3. **The package's own inverse fixture is stale.** `toy_network_ex2.RData` ships columns
>    `Source/Effect/Target`; the 2.20.0 validator requires `source/interaction/target` (which
>    `toy_network_ex1.RData` already uses). The inverse toy fails out of the box until renamed —
>    a format-only change that must be recorded, not silently applied.
> 4. **Pre-intersect measurements with PKN nodes before calling.** If *any* measured node is absent
>    from the PKN, CARNIVAL aborts the whole call with *"None of your measurements are in prior
>    knowledge networks (PKN)"* — which is **factually wrong** when some do match (verified: the
>    matching node alone runs fine). Intersect first and log the dropped-TF count, or runs fail with
>    a misleading message and lost coverage goes unrecorded.

**Implementation choices:**

1. Add an R adapter \`scripts/systems_validation/fault_isolation/run_carnival.R\` to read a frozen \`tf_activity.tsv\` and \`signed_pkn.tsv\`, validate IDs/signs, invoke the **published** CARNIVAL R package, and serialize all competing feasible/optimal subnetwork solutions permitted by the solver.
2. Use the package's own toy input, measurement and PKN fixtures as the **first replication gate**. The documented R API includes \`runCARNIVAL(inputObj=..., measObj=..., netObj=...)\` for a target-constrained toy and \`runCARNIVAL(measObj=..., netObj=...)\` for an inverse toy. Confirm the installed API/version rather than pinning unverified function signatures. Prefer open-source **CBC** for substantive runs; \`lpSolve\` is acceptable for the package toy examples.
3. Export separate, frozen \`signed_pkn.tsv\` representations from the repo's curated signaling sources. Preserve source and target gene symbols/IDs, direction, sign (including unknown/conflicted), interaction provenance, degree, and tissue restrictions. ARACNe coexpression links must not silently become directed biochemical causal interactions.
4. Compute TF activities using the project's established VIPER outputs, and add an independent/cross-network arm to quantify dependence on the program representation. Keep disease TF activities separate from LINCS treatment TF activities; compare them in a **new** evaluation stage, never by pooling as one training response.
5. Predeclare the candidate upstream node universe and whether each candidate driver was supplied, hidden, or inferred. In inverse mode, prohibit the held-out label and phenotype evidence from influencing target selection, path choice, tuning or score interpretation.
6. Freeze ILP complexity weights, allowed TF-sign handling, PKN version, time limit, relative MIP gap and seed *before* held-out perturbations. Treat infeasible, timed-out, gap-not-closed and multiple-solution cases explicitly; no solution may be silently discarded for not matching the expected driver.
7. After committing the predicted regulator and gene consequences, evaluate them with independent, context-appropriate perturbations. All comparison genes/programs must be fixed before observing those responses.

**Known-answer mathematical fixtures** (hand-check each expected sign and graph, then compare against CARNIVAL output):

| Fixture | Ground truth | Required behavior |
|---|---|---|
| Four-node signed chain \`X -> A (+), A -> B (-), B -> C (+)\` | Activating X predicts \`A +, B -, C -\` | A target-constrained solution is sign-consistent; reversing the X perturbation reverses predicted effects |
| Two sources with indistinguishable measured TF effects | Driver cannot be uniquely recovered | Inverse mode reports alternative explanations / non-identifiability, not an invented unique driver |
| ↳ **[AGENT AUDIT 2026-11-18] TESTED — DOES NOT HOLD on open solvers** | `S1->H`, `S2->H`, `H->T1`, `H->T2`, `T1=T2=+1` | CARNIVAL inverse returned **`S1=100, S2=0`**: it picked one of two provably indistinguishable sources and reported **no** ambiguity. See §4.8a |
| Contradictory TF observations relative to available paths | No perfectly coherent explanation | Report residual disagreement or infeasibility, never invent a path |
| Disconnected measured TF | No connection in the frozen PKN | Do not fabricate reachability |
| Large high-degree decoy | Fewer observations uniquely linked to decoy | Recoverability is evaluated against degree-preserving and label-matched decoys |
| Hidden true cause | No driver label supplied | Recovery performance evaluated on held-out truth without any target constraint |

**Baselines and primary endpoints:** versus DE magnitude, VIPER activity/centrality, D128 signed-path tracing, and randomly degree-matched candidate causes, compare (i) held-out perturbation **signed-program accuracy**, (ii) ability to discriminate true versus false intervention hypotheses, (iii) false-positive rate on complete negative-control triples, (iv) causal ambiguity (number and diversity of comparably feasible subnetworks), (v) network-choice sensitivity, and (vi) solver feasibility/coverage. Driver-label accuracy is reported **only in the inverse arm**. Use distinct target classes and independent perturbation data as sampling units; preserve gene covariance.

**Gate F0 (published math):** the installed package reproduces its standard and inverse toy examples with auditable solutions. **Gate F1 (synthetic):** signs, ambiguity and impossible cases behave correctly. **Gate F2 (biological):** on genuinely held-out perturbation outcomes the method predicts intervention effects better than the frozen simple baselines, without a misleading increase in negative-control closures. An attractive inferred network by itself is **not a pass**.

**References:** [CARNIVAL paper](https://doi.org/10.1038/s41540-019-0118-z); [published package and vignettes](https://saezlab.github.io/CARNIVAL/); [R tutorial including both toy modes](https://rdrr.io/bioc/CARNIVAL/f/inst/doc/CARNIVAL.Rmd).

### 4.8a [AGENT AUDIT 2026-11-18] Causal-ambiguity enumeration requires CPLEX — amendment needed

> This section asserts, in several places, that alternative explanations will be enumerated:
> F-MATH §1 ("serialize **all competing feasible/optimal subnetwork solutions**"), the fixture table
> ("reports alternative explanations / non-identifiability"), §4.7 endpoint (iv) ("**causal ambiguity**
> — number and diversity of comparably feasible subnetworks") and §4.7's failure mode ("multiple
> indistinguishable explanations").
>
> **Solution-pool enumeration is a CPLEX-only capability in CARNIVAL 2.20.0.** Inspected the three
> option constructors:
>
> | solver | pooling fields exposed |
> |---|---|
> | `defaultLpSolveCarnivalOptions()` | **none** |
> | `defaultCbcSolveCarnivalOptions()` | **`poolrelGap` only** |
> | `defaultCplexCarnivalOptions()` | `poolrelGap, limitPop, poolCap, poolIntensity, poolReplace` |
>
> CPLEX is commercial and is not installed. **Confirmed empirically** by the fixture row added to the
> table above: a provably non-identifiable two-source case returned a single driver with no ambiguity
> signal.
>
> This matters because it is **precisely the failure this spec was written to prevent** — D137's
> lesson was that the procedure produced plausible explanations for negative controls, and measuring
> ambiguity was the countermeasure.
>
> **Three options, one of which must be chosen before F proceeds to biology:**
> 1. Obtain a **CPLEX academic licence** (free for academics; the only route to this spec as written).
> 2. **Enumerate externally** — re-solve with no-good cuts excluding each solution found. This is new
>    code and needs its own validation against a known-ambiguity fixture.
> 3. **Amend the spec**: state that ambiguity is not measurable on open solvers, drop endpoint (iv),
>    and label every F explanation "one of possibly many".
>
> **Not acceptable:** reporting CARNIVAL's single returned solution as *the* explanation.
>
> **Related, from the same F1 run:** contradictory observations are **silently zeroed**, not reported.
> Fixture `X->A(+), A-|B` with `A=+1, B=+1` returned `A=1, B=0` — no path invented (the spec's "never
> invent a path" holds), but no disagreement surfaced either. An F pipeline must detect dropped or
> zeroed measurements itself. The signed-chain fixture passed exactly: `X=+1 → A+100, B−100, C−100`,
> and reversing X reversed all three.

### 4.7a Candidate nomination from Experiment F

Experiment F is not complete when it outputs a plausible causal graph. It must nominate actionable **disease mechanisms and corresponding existing drug-target pairs**, with verified source versus proposed new indication, patient-subgroup driver evidence, frozen target-to-program intervention predictions and independent observations. A driver supplied to CARNIVAL cannot be reported as newly discovered. Deliver F_candidate_mechanisms.tsv and F_frozen_drug_predictions.jsonl, with alternative causes and withheld test IDs.

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


### 5.8 O-MATH — Structural matching, Kalman tests and finite-horizon intervention effects

**Source and implementation boundary.** The Liu–Slotine–Barabási network papers motivate **structural** control/sensing from directed topology. Classical Kalman tests apply to a **specified numerical linear dynamical model**, not directly to a transcriptomic adjacency matrix. Gates and Rocha supply the essential warning that dynamics can change a graph-only conclusion. Implement the mathematical procedures as *distinct* arms and compare them with actual perturbation data.

**State-space model and dimensions:**

    z_dot(t) = A z(t) + B u(t),        y(t) = C z(t) + ε(t)

    z ∈ R^K  (predefined molecular programs)
    A ∈ R^(K×K)  (program-to-program dynamics; convention A[i,j] = j -> i)
    B ∈ R^(K×M)  (M independently actuated perturbation inputs)
    C ∈ R^(P×K)  (observed assay outputs)
    u ∈ R^M, y ∈ R^P

The gene-to-program transformation \`W\` in §3 is **not automatically** a validated observation matrix \`C\`; specify the measurement model and its uncertainty. \`B\` is *not* defined by a ChEMBL target annotation alone. Signs and target engagement must be calibrated to the actual perturbation context where empirical control claims are intended.

**O-MATH-1: Structural controllability (topology-only baseline).**

1. Convert the frozen directed graph to a bipartite representation with \`K\` left (sources) and \`K\` right (destinations) copies. Add a bipartite edge \`j_left -> i_right\` for every permitted dynamical \`A[i,j]\` entry.
2. Find a maximum-cardinality matching (e.g., NetworkX bipartite matching / Hopcroft–Karp). The right-side unmatched nodes are the candidate *minimum-driver-node* locations under the **generic structural controllability** assumptions. If every right-side node is matched, at least one external input is still required for a nonempty system.
3. Test accessibility from proposed actuators (all relevant state nodes must be reachable through the directed input/state graph). Distinguish **input placement at arbitrarily selected nodes** from the real available drug target set. A hypothetical minimum number of unconstrained actuators is not a count of available repurposable drugs.
4. Compute the corresponding structural-observability sensor baseline using the correctly transposed/dual graph formulation; then check it with numerical observability tests.
5. Vary network edge definitions, unknown signs, self-loops, and source-network choice **as explicit model arms**. Structural controllability means generic feasibility for almost all choices of allowed nonzero weights under linear assumptions, not guaranteed control for the biological parameters.

**O-MATH-2: Numerical Kalman rank tests for an identified linear system.**

    Kc(A,B) = [ B, A B, A² B, …, A^(K-1) B ]
    Ko(A,C) = [ C; C A; C A²; …; C A^(K-1) ]

    fully controllable  <=> rank(Kc) = K
    fully observable    <=> rank(Ko) = K

Implement using small predeclared K, fixed numerical tolerances, SVD singular-value diagnostics, and a condition-number report. A full numerical rank under noise or an ill-conditioned \`Kc\` is not automatically *practical* controllability.

**O-MATH-3: Finite-horizon reachability and control energy.** When A and B are identified, calculate

    Wc(T) = ∫_0^T exp(A t) B Bᵀ exp(Aᵀ t) dt

and on the reachable subspace, for the unconstrained LTI minimum-energy problem,

    δ = z_target - exp(A T) z_initial
    E_min = δᵀ Wc(T)^(-1) δ

where the ordinary inverse is permitted only if \`Wc(T)\` is nonsingular. Otherwise report unreachable directions / the reachable-subspace condition and do **not** present a pseudoinverse energy as evidence that an unreachable target is attainable. This expression assumes unconstrained inputs and an exact LTI model; bounded drug doses, PK, irreversible pharmacology and toxicity require an additional constrained optimization.

For measured output objectives \`H z\` (e.g., a disease-specific program subset), explicitly distinguish **output controllability** from controlling all K internal states. Do not interpret a large graph-theoretic driver-node set as proof that a particular therapeutic program can be steered.

**O-MATH-4: Identify A, B and C only with sufficient perturbation/time information.** Implement a data-availability gate and compare a low-dimensional linear state-space fit with simple static and null models. Possible fitting methods include subspace identification or regularized state-space regression, but choose and freeze one only after an identifiability audit. A one-time LINCS endpoint or tumor differential expression alone cannot identify a continuous-time \`A\`. The current D134 target retrieval and D135 genetic concordance are *empirical evidence arms*, not estimates of Kalman matrices.

> **[AGENT AUDIT 2026-11-18] All 10 verified — 10/10 behave exactly as stated** (`scripts/systems_validation/observability/gate_o0_fixtures.py`, numpy 2.4.6 / scipy 1.17.1 / networkx 3.6.1). Kalman ranks 2/2, 1 for the wrong actuator and 1 for the wrong sensor; the directed path's unmatched driver node is `[1]`; the 3-cycle has a perfect matching and no unmatched node; a non-upstream actuator leaves rank 1 of 3; weak coupling gives full rank with a 1e-8 singular value; `Wc` is singular on unreachable directions. The bipartite construction and every stated answer are correct.

**Required hand-computable regression tests:**

| Fixture | A, B, C or graph | Correct answer |
|---|---|---|
| Two-state controllable/observable chain | \`A=[[0,1],[0,0]]\`, \`B=[[0],[1]]\`, \`C=[[1,0]]\` | \`rank(Kc)=2\`, \`rank(Ko)=2\` |
| Wrong actuator on the same chain | \`B=[[1],[0]]\` | \`rank(Kc)=1\` |
| Wrong sensor on the same chain | \`C=[[0,1]]\` | \`rank(Ko)=1\` |
| Directed path of three nodes | \`1->2->3\` | Unmatched input node is 1 in the structural maximum matching |
| Fully matched directed cycle | \`1->2->3->1\` | No unmatched node, but **at least one** input remains necessary |
| Reachability does not imply rank | Chain with one actuator only at a non-upstream node | At least one state unreachable; full control fails |
| Exact rank but poor conditioning | Weak inter-node coupling | Rank may pass while small singular values/control energy show impracticality |
| Missing dynamics | Graph only; no time-series identification | Structural output permitted; numerical biological controllability **not evaluable** |

**Implementation contract:** Provide \`structural_matching.py\`, \`kalman.py\`, \`gramians.py\`, \`identify_linear.py\` and \`evaluate.py\` under \`scripts/systems_validation/observability/\`. Use NumPy/SciPy/NetworkX with version-pinned environment metadata. The structural and numerical findings must have **different output fields**: \`structural_candidate_inputs\`, \`numeric_rank\`, \`min_singular_value\`, \`finite_horizon_energy\`, \`identification_status\`, and \`observed_target_response\`.

**Gate O0:** all rank/matching/toy energy checks exactly pass. **Gate O1:** structural results are stable enough under justified graph uncertainty to be meaningful as *hypotheses*. **Gate O2:** empirical D134/D135/D136 tests provide independent target specificity; a graph-only result never suffices. **Gate O3:** only identified and held-out-validated dynamic models receive numerical controllability/observability claims. If data cannot support \`A\`, \`B\`, \`C\`, record an untested dynamical claim rather than estimating arbitrary matrices.

**References:** [Liu et al., Nature 2011](https://doi.org/10.1038/nature10011); [Liu et al., PNAS 2013](https://doi.org/10.1073/pnas.1215508110); [Gates & Rocha 2016](https://doi.org/10.1038/srep24456).

### 5.7a Mechanistic drug-response evidence from Experiment O

Experiment O must take F's proposed **drug × new-disease × target** candidates and assess how the drug could modify the nominated program in the relevant context. Compare gene-level LINCS responses with the drug-side ARACNe/VIPER response as two prespecified, non-independent representations of the *same perturbation*. Supplement with D135 genetic evidence, D134 target recovery, D136 dose-response and annotated potency/exposure when eligible. Deliver O_candidate_engagement.tsv and O_program_prediction_checks.tsv with both supportive and contradictory observations. Do not equate a measurable response with therapeutic efficacy.

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

> **[AGENT AUDIT 2026-11-18] This audit is already done; here is the answer.** From
> `data/raw/lincs2020/siginfo_beta.txt` (720,216 `trt_cp` signatures, 201,494 compound × line cells):
>
> | distinct timepoints per (compound, line) | all cells | at `is_hiq` & TAS ≥ 0.4 |
> |---|---|---|
> | ≥ 2 | 47,216 | 1,811 |
> | **≥ 3** | **2,038 (1.0%)** | **317** |
> | ≥ 4 | 226 | 28 |
> | ≥ 5 | 24 | — |
>
> 97% of signatures sit at 6 h or 24 h. **`siginfo_beta` contains no washout, recovery or
> drug-removal field** — the only time/dose columns are `pert_time`, `pert_dose`, `pert_type`. LINCS
> measures continuous exposure to a fixed endpoint.
>
> **Therefore S-H3 (persistence) is not evaluable by construction, not by sample size**, and §6.8's
> rule forces `persistence_status = not_evaluable` for every compound. §0.2 and §6.9a already make S
> conditional, so this is documentation rather than a blocker. S-math is unaffected and passes (8/8).
> If temporal work is wanted, the better-powered question is the 6 h vs 24 h divergence in
> `SPEC_dose_response_bmd.md` §8 — 1,811 usable cells rather than 317.

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


### 6.10 S-MATH — Reproduce Angulo et al.'s autonomous-element example and test nonlinear accessibility

**Source and implementation boundary.** First reproduce the **published three-species generalized Lotka–Volterra (GLV) counterexample** from Angulo, Moog and Liu (2019). Its state variables are microbial population abundances, **not** gene-expression z-scores. The GLV model is a **mathematical benchmark**; any later molecular-regulatory model must have independently justified state variables and interactions. Distinguish accessibility (locally available motion directions) from global controllability, attractor stability, and persistence.

**S-MATH-1: Exact three-state benchmark** (Fig. 2 of Angulo et al.).

For positive species abundances \`x1,x2,x3\`, use

    x1_dot = x1 * (-1 + x3)
    x2_dot = x2 * ( 1 - x3)
    x3_dot = x3 * (-0.5 + 1.5*x3) + u3(t)

where only x3 is actuated, with the paper also considering impulsive interventions on x3. Define the **autonomous element**

    ξ(x) = x1 * x2.

For every admissible control \`u3(t)\` and every trajectory where the equations apply:

    ξ_dot = x2*x1_dot + x1*x2_dot
          = x1*x2*(-1+x3) + x1*x2*(1-x3) = 0.

Therefore \`x1(t)*x2(t) = x1(0)*x2(0)\` and targets with a different product are unreachable by interventions acting **only** on x3. This is an **exact algebraic property**, not an estimated statistical association.

Add a second control \`u1(t)\` directly to x1:

    x1_dot = x1 * (-1 + x3) + u1(t)
    ξ_dot = x2 * u1(t),

so the invariant is no longer preserved; for an impulse, \`Δξ = x2 * Δx1\` at the instant of the jump. Removing this invariant is **necessary evidence of expanded accessibility**, not by itself a guarantee that every biologically relevant state is globally reachable under positivity, amplitude, dose or time constraints.

**Simulation tests:** use \`scipy.integrate.solve_ivp\` and independently hand-check the invariant. For \`x(0)=(2,3,1)\`, x1*x2 remains 6 under multiple fixed, pulse and time-varying x3 controls, within tolerances established by solver convergence and time-step refinement. The nominal desired state \`x_target=(1,1,1)\` has product 1, so a controller acting only on x3 must **never** claim it is reached. Demonstrate a changed product after actuating x1. Test both continuous and impulsive control arms and positivity guards separately.

**S-MATH-2: Local nonlinear accessibility via vector fields and Lie brackets.**

For a control-affine system

    z_dot = f(z) + Σ_j g_j(z) u_j(t)

compute the span of available control directions and their successive Lie brackets; for example, under one convention

    [f,g](z) = Dg(z) f(z) - Df(z) g(z).

A full-rank accessibility distribution at a regular point is an appropriate **local accessibility** criterion under the relevant smoothness and control assumptions; it is not a certificate of global steering, stable-state transitions or bounded-input feasibility. State the bracket sign convention (rank does not depend on its overall sign).

- Implement symbolic differentiation for the published 3-state toy using SymPy and cross-check with finite-difference Jacobians.
- In the one-actuator toy, demonstrate the invariant's differential \`dξ\` annihilates the accessible directions; the distribution cannot have full rank on the positive domain.
- With independent actuation of x1 and x3, show the additional field \`g1=(1,0,0)^T\` allows full local rank in the positive domain for the published toy: \`g3=(0,0,1)^T\`, \`g1\`, and \`[f,g3]\` are linearly independent when x2>0.
- Reject conclusions based on singular points or numerical rank thresholds chosen after observing a favorable result.

**S-MATH-3: Disease-regulatory adaptation, explicitly not the paper's GLV model.** Construct small causal circuits with biological meanings and known dynamics (e.g., negative feedback, positive-feedback bistability, and a compensating parallel pathway). Start with dimension 2–5. Candidate model forms include signed Boolean networks and sparse nonlinear ODEs of the form

    z_dot = f(z; θ) + B(z; θ) u(t),     y(t) = C(z(t)) + ε(t).

Use biologically motivated saturating regulation where appropriate. Do **not** fit \`x_i*(r_i + Σ_j a_ij*x_j)\` directly to gene expression merely because the ecology paper used GLV for abundances.

Freeze synthetic ground truth (A/B dynamics, feedback, initial conditions, measurement map, perturbation schedule and noise) **before** testing a model. Compare a static baseline, a linear state-space model, a small Boolean switch and a sparse nonlinear model, with matching training/holdout protocols.

**S-MATH-4: Distinguish model properties that are often conflated.**

- **Equilibrium:** solve \`f(z*)=0\` under no drug (or a specified constant input), and assess local stability from Jacobian eigenvalues; a real-part-negative Jacobian spectrum supports local asymptotic stability for an isolated hyperbolic equilibrium.
- **Bistability:** demonstrate at least two stable equilibria with distinct attraction basins in the *identified model*; multiple transcriptomic clusters do not prove this.
- **Accessibility:** a property of allowed control fields and their Lie brackets, with stated domain.
- **Reachability:** existence of an admissible input taking one initial state to a desired set in finite time. Show an explicit feasible trajectory or report unknown; do not infer global reachability from a local rank test.
- **Persistence:** measured alteration after drug withdrawal. This cannot be established from continued treatment or a simulated attractor alone.
- **Clinical benefit:** independent disease-relevant functional outcome, not proximity to a reference expression vector.

**S-MATH-5: Biological identification gate.** First audit actual LINCS/other time-series granularity (same drug, source line, concentration, actual elapsed time, replicates and matched controls). Hold out whole trajectories or late timepoints; compare predictive performance against linear and no-dynamics baselines. A dataset with insufficient time resolution, no independent perturbations, or unmeasured recovery must **abstain** from claims of attractor changes or persistence. Keep dose-response BMD (D136) separate from dynamical parameter identification.

> **[AGENT AUDIT 2026-11-18] The Angulo rows verified — 8/8 exact** (`scripts/systems_validation/state_transitions/gate_s0_angulo.py`, sympy 1.14.0). `d(x1·x2)/dt` simplifies **symbolically to 0** under x3-only control and to **exactly `x2·u1`** with x1 actuation; numerically `|x1·x2 − 6| < 3.2e-10` from `x(0)=(2,3,1)` under constant, pulse and sinusoidal controls; the target `(1,1,1)` has product 1 and is correctly unreachable; accessibility rank is **2** with one actuator and **3** with two at `x=(2,3,1)`; `dξ` annihilates the one-actuator distribution. The stated Lie-bracket convention and ranks hold. The remaining rows (synthetic circuits, washout classification) are not yet implemented.

**Required synthetic unit and property tests:**

| Fixture | Ground truth | Required result |
|---|---|---|
| Angulo 3-state system, x3 control only | \`d(x1*x2)/dt=0\` exactly | Symbolic zero and numerically conserved invariant under multiple controls |
| Desired product differs from initial product | Unreachable under x3 control | Reachability optimizer may not claim success |
| Add x1 actuation | \`d(x1*x2)/dt=x2*u1\` | Invariant broken; confirm changed product |
| Positive-domain Lie bracket test | One actuator rank deficient; two actuators locally full rank | Correct rank and convention |
| Linear stable system | One known equilibrium and basin | Correct equilibrium and Jacobian stability result |
| Double-well / positive-feedback switch | Two attractors with distinct basins | Detect both; avoid spurious third attractor |
| Compensatory parallel pathway | Inhibit one branch, response rebounds through other | Reproduce known time course and missed effect by static baseline |
| Sparse 2-timepoint observations | Dynamics not identifiable | Model returns **not evaluable**, not a confident attractor |
| Time course during treatment only | Persistence unobserved | \`persistence_status=not_evaluable\` |
| Washout trajectory | Reversible versus persistent known by design | Correctly classify response without reading training trajectory twice |

**Implementation contract:** \`scripts/systems_validation/state_transitions/{angulo_glv.py,accessibility.py,synthetic_circuits.py,fit_dynamics.py,temporal_audit.py,evaluate_trajectories.py}\`; tests for exact invariants, solver convergence, local rank, equilibria and recovery. NumPy, SciPy and SymPy are proposed dependencies; pin their versions when implementing. Record numerical tolerances and solver methods in the run manifest.

**Gate S0 (published math):** reproduce the exact Angulo Fig. 2 invariant with both one- and two-input arms. **Gate S1 (synthetic):** distinguish reversible, persistent and non-identifiable cases under held-out validation. **Gate S2 (biological):** run only if matched longitudinal data support state identification. **Gate S3 (persistence):** no claim of irreversible/persistent disease-state change without actual washout/recovery observations and appropriate functional context.

**References:** [Angulo et al. 2019](https://doi.org/10.1038/s41467-019-08890-y); [open-access full text, Fig. 2 and methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6401173/).

### 6.9a Temporal qualifications for repurposing nominations

Experiment S adds temporal predictions to **eligible** candidate dossiers: rebound, compensation, persistence, or return toward a reference program. Use per-dose/time gene data and condition-specific VIPER inference, not collapsed drug-line averages. S is marked **not evaluable** when longitudinal and washout data are insufficient, without automatically rejecting otherwise supported disease–drug mechanisms. Deliver S_candidate_trajectories.tsv and a proposed next time-course/washout test.

## 7. Integrated interpretation: drug repurposing proposals and mechanism traces

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

The combined report must **name the proposed drug and new disease/subgroup, explain why the target is disease-relevant and how the drug could change its pathological pathway, then state what independent experiment supports or falsifies that proposal**. An appropriate finding may be **mechanism supported; intervention response supported; persistence and clinical benefit not established**. Do not require or infer a numerical suitability score.

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

### Stage 6 — New-indication candidate nomination and mechanism-evidence integration

- Generate the §0 candidate dossier per *drug × proposed new disease/subgroup*, documenting the drug's existing indication, disease driver, annotated action, predicted route, LINCS/ARACNe and genetic observations, pharmacology, contradictions and next falsification experiment.
- Test whether proposals survive complete negative controls and predict independently held-out interventions or functional outcomes. Report preliminary, contradicted and not-evaluable candidates without forcing positive nominations.
- Render **why/how mechanistic explanations** strictly from structured evidence and preserve uncertainty, without suggesting clinical efficacy or treatment.
- Optionally consider validated prioritization scores as a later separate workstream; no production ranking or tier changes are authorized by this draft.


### Mathematical replication checkpoints (required before the corresponding biological stage)

The **source-method replication** and the **biological adaptation** are separate deliverables. Passing a paper's toy benchmark does not establish that its premises hold for tumors.

| Method | Source-method checkpoint | Biological applicability gate | Failure mode explicitly detected |
|---|---|---|---|
| CARNIVAL ILP | Run both built-in standard and inverse toy cases; record optimal network, feasibility, solver gap | Program activities and signed priors match the intended context; independent perturbations available | Inferred driver merely supplied as input; multiple indistinguishable explanations |
| Maximum matching | Reproduce directed path and directed cycle driver-node cases | Reliable directed causal edge semantics and plausible drug-actuation sites | Topology-only control treated as observed effect |
| Kalman matrices/Gramian | Reproduce 2-state rank and finite-horizon reachability cases | Enough time-resolved interventions to identify A/B/C | Arbitrary dynamics, ill conditioning or unreachable target |
| Angulo GLV + accessibility | Prove & numerically test x1*x2 invariant; demonstrate breakdown with x1 actuation | Small, independently specified regulatory dynamics with adequate temporal data | Global state transition inferred from local accessibility |
| State persistence | Washout synthetic true/false classification | Biological recovery or washout observations | Treatment-only response interpreted as persistent state change |

**Test order:** source-method toy replication -> synthetic adversarial cases -> source-specific feasibility -> predictions frozen -> held-out perturbation/trajectory evaluation -> model-comparison decision. Freeze source-package versions, Git SHA, random seeds, network versions, numerical tolerances and solver options before each test. A method that fails its mathematical replication is not eligible for biological inference; a method that passes mathematics but fails biological controls remains a **negative or restricted-domain** scientific result.

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

> **[AGENT AUDIT 2026-11-18] Files that now exist, and two gaps in the tree above.**
>
>     scripts/systems_validation/
>       fault_isolation/gate_f0_carnival_toys.R    # Gate F0 - PASS (4/4)
>       fault_isolation/gate_f1_adversarial.R      # Gate F1 - see 4.8a
>       observability/gate_o0_fixtures.py          # Gate O0 - PASS (10/10)
>       state_transitions/gate_s0_angulo.py        # Gate S0 - PASS (8/8)
>     results/systems_validation/synthetic/
>       gate_f0_carnival.json  gate_f1_adversarial.json
>       gate_o0_fixtures.json  gate_s0_angulo.json
>     environment.sysval.lock.yml                  # pinned env for all of the above
>
> **Gap 1 — no `tests/` entries for the math fixtures** that Gates F1/O0/S0 depend on; the tree lists
> `tests/systems_validation/test_*.py` but not the gate fixtures themselves.
> **Gap 2 — no interpreter declaration.** Three environments are in play with incompatible versions
> (`sysval` pandas 3.1.1 / Python 3.11; `celligner` pandas 1.4.4 / Python 3.10, which is what the
> drug-side pipeline runs under; `repurposing` for R, no SciPy). Every file in this tree must state
> which interpreter it runs under. A `79i` run failed on `ModuleNotFoundError: scipy` under the wrong
> Python on 2026-11-17 — this hazard is live, not theoretical.
>
> Also missing from the tree: the `run_carnival.R` adapter that F-MATH §1 item 1 requires, and the
> **per-condition VIPER computation** that §0.1/§0.3/§6.9a imply. `79c` keeps one consensus profile
> per (drug, line) at TAS ≥ 0.4 and discards the dose/time structure, so the dossier pipeline needs
> its own activity computation from Level 5 rather than a join onto existing `79d` outputs.

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

**Best case:** The system nominates existing drugs for *specific new diseases or molecular subgroups* with independently validated upstream drivers, interpretable drug-target mechanisms, experimentally observed downstream effects, plausible exposure and explicit falsification experiments. It generates scientifically defensible **why/how drug repurposing proposals**, not merely ranked compounds.

**Realistic case:** A limited number of disease-specific repurposing proposals merit further experiments, while many candidates are preliminary, contradicted or lack suitable data. Some drivers are difficult to observe transcriptionally and longitudinal evidence is sparse. The valuable output is a **mechanistic domain-of-applicability map** and proposals with explicit evidence gaps.

**Negative case:** Causal reconstruction adds nothing beyond supplied annotations, the drug response cannot be attributed to the predicted target, or context-matched evidence is absent. Report **no supported repurposing proposals** rather than produce convincing narratives or unvalidated scores.

## 13. Immediate recommendation

> **[AGENT AUDIT 2026-11-18] Partly superseded.** Gates F0, O0 and S0 **have been run and pass**;
> dependencies are installed and locked. The remaining order is: (1) decide the §4.8a CPLEX question,
> which blocks F from biology; (2) answer the §2.6(c) dossier↔ranking relationship, which determines
> whether any of this can change an output; (3) Gate F1's remaining adversarial fixtures (high-degree
> decoy, hidden true cause) and the unimplemented S-MATH rows; (4) export the frozen `signed_pkn.tsv`
> from `data/raw/omnipath/omnipath_signed_directed.tsv.gz`. Nothing biological should start before (2).

**Start with the shared Stage 0 audit and Experiment F's known-answer synthetic benchmark.** In parallel, continue D134, D135 and D136 as the near-term evidence routes for Experiment O. Delay biological Experiment S until its longitudinal feasibility gate is satisfied.

**First milestone:** Using BRAF-V600 melanoma and ESR1-associated breast cancer as *seen development mechanisms*, determine whether the approach can (i) identify disease-relevant actionable upstream targets, (ii) map these to existing drugs with potentially different indications, (iii) freeze predictions, and (iv) issue a **drug → new disease → why → how → falsifying experiment** dossier or a justified no-proposal conclusion. These seen cases are not independent proof of repurposing efficacy.

## 14. Literature foundations

- Yang et al. (2012), *Progress in Root Cause and Fault Propagation Analysis of Large-Scale Industrial Processes*, Journal of Control Science and Engineering.
- Liu et al. (2019), *From expression footprints to causal pathways: contextualizing large signaling networks with CARNIVAL*, npj Systems Biology and Applications, DOI: 10.1038/s41540-019-0118-z.
- Liu, Slotine and Barabási (2011), *Controllability of complex networks*, Nature, DOI: 10.1038/nature10011.
- Liu, Slotine and Barabási (2013), *Observability of complex systems*, PNAS, DOI: 10.1073/pnas.1215508110.
- Gates and Rocha (2016), *Control of complex networks requires both structure and dynamics*, Scientific Reports, DOI: 10.1038/srep24456.
- Angulo, Moog and Liu (2019), *A theoretical framework for controlling complex microbial communities*, Nature Communications, DOI: 10.1038/s41467-019-08890-y.
- CARNIVAL software/documentation: https://saezlab.github.io/CARNIVAL/ ; package vignette: https://rdrr.io/bioc/CARNIVAL/f/inst/doc/CARNIVAL.Rmd .
- *Network inference from perturbation time course data* (2022), npj Systems Biology and Applications, DOI: 10.1038/s41540-022-00253-6.

## 15. Final objective

Replace the narrow question

> Which compound produces the largest numerical reversal of a disease-associated signature?

with

> Which mechanisms maintain the disease-associated state, which interventions can modify those mechanisms, and what independent evidence supports the predicted consequences?

**The intended output is a set of explicit, falsifiable drug-repurposing proposals for *new* diseases or molecular subgroups.** Each must explain **why** an existing drug's target matters in the new pathology and **how** its pharmacological action is expected to affect the disease-maintaining network, supported by independent measurements, competing explanations and a next experimental test. **Drug scores are not mandatory.**

**Draft only.** No new experimental claim, positive drug-repurposing verdict, clinical treatment recommendation, production score change or drug reranking follows from this document.
