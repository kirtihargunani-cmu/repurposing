# Draft spec: Target-mechanism scoring and incremental validation (D128 candidate)

**Status:** exploratory draft, 2026-11-14. **Not an adopted score or a pre-registered confirmatory test.** Freeze cohort eligibility, score variants, controls, and acceptance thresholds in a dated commit before accessing new holdout outcomes. This document specifies what to test; it reports no new benchmark run.

## 0. Why this experiment, and what is already established

The existing drug-side pipeline already has (a) tumor-versus-normal differential expression and VIPER master-regulator (MR) activities, (b) weighted/signed ARACNe, CollecTRI and DoRothEA regulons, (c) LINCS Level-5 cancer-line compound signatures, generic-response removal, cell-line weighting and drug VIPER, (d) target role / layered scores and (e) biomarker-pair benchmarks. Do **not** rebuild preprocessing. Preserve the live score and exports; run candidate scoring in a separate sandbox.

Relevant precedents:
- D120 measures whether the *annotated target* is dependent, altered or overactive in the cancer. It is a target-relevance score, not proof the drug engaged the target.
- D121/D122 create the current **layered_pct** lead ranking with target units, direction and relative potency. This is the primary production comparator; do not silently replace it.
- D125 regulator-overlap Laplacian did not clear the specified improvement bar (+0.17 vs +1.0 percentage point), so graph complexity alone is not validated.
- D126 signed shortest-path consistency **failed** its controls: flagged **0/22** anti-estrogen/MCF7 breast rows; median mechanism coverage ~0.08; no qualifying drug cleared its broad coverage floor. It selected the top 25 *changed* regulators, missing ESR1, and should not be rerun unchanged.
- D127 primary paired-biomarker test of **layered_pct** passed in its pre-declared set: **36/49**, share **0.7347**, sign-test p **0.0007**, median gap **15.1 percentile points**. But individual classes differ: thyroid BRAF-like vs RAS-like 3/9; SKCM vs LUSC 8/9; LUAD vs LUSC EGFR 5/8; HNSC vs LUSC EGFR 8/8; HER2 vs basal 5/6. All D127 pairs and results are *seen* development benchmarks for any new score. Do not label a rerun on them independent validation.

Working scientific hypothesis: target-anchored, signed and confidence-weighted **regulon member responses** contain useful, drug-specific information beyond aggregate regulator activity and target annotation, improving mechanistic fidelity and biomarker-specific ranking without mistaking common growth arrest for on-target action.

### Distinguish four propositions

1. **Target relevance (T):** the disease depends on a target or has an actionable molecular alteration. Already measured by D120/122; mutation/fusion and target status are separate from mRNA expression.
2. **Target-mechanism engagement (M):** observed drug perturbation has a footprint plausibly connected to its annotated target and action, evaluated in the actual source-line context.
3. **Disease-local rescue (R):** those target-linked genes/regulons move opposite to the disease-associated direction.
4. **Therapeutic phenotype (V):** independent viability/response evidence supports benefit. R and M do not imply V, particularly for synthetic lethality, cytotoxicity and resistance.

Report T, M, R, and V separately until an integration rule clears an external validation criterion. Missing M is **not** zero M; output unscorable/abstain with a reason.

## 1. Core inputs and frozen data contract

| Layer | Existing project source | Contract |
|---|---|---|
| Tumor gene signatures | results/step4_v2; TCGA/TARGET vs GTEx or specified adjacent normal | Gene ID, cohort/group, logFC, moderated statistic or confidence, contrast/provenance; no imputed treatment response |
| Disease regulator activity | results/step6/aracne_project_v2_tf and collectri_v2; strata equivalents | Same regulator and network identifiers used on drug side |
| Regulator-target edges | data/processed/step5/regulons_aracne_project_v2_tf.rds; results/step5/regulons_core.rds | For each edge: target, TF, sign (tfmode), likelihood/confidence, network, target inclusion |
| Drug gene responses | data/processed/lincs/drug_line_profiles_bing_specific.f32; drug_line_profiles_index.tsv; drug_line_profiles_genes.tsv, built by 79c | MODZ differential-expression *z-score* profile, profile-major float32; explicit gene ID and line/drug index; raw arm only as sensitivity |
| Drug regulator responses | data/processed/lincs/drug_viper_specific_v2_tf/<network>.f32; drug_viper_specific_collectri/ALL.f32, built by 79d | Real VIPER, not a renamed linear projection; verify row layout and network/TF match |
| Source-line context | results/celligner/disease_line_weights_long.tsv, stratum weights, config/celligner_tissue_keys.tsv | Reuse direct/borrowed split, weights and transfer-confidence flags |
| Target and action | D105b evidence_targets, drug_annotations.tsv, D120 ChEMBL direct-interaction mechanism/target units | Unit/complex membership, target source, known action, selectivity/potency confidence |
| Signed signaling prior | cached OmniPath and SIGNOR D126 tables | Signed, directed edges, provenance, degree; no inference from undirected links |
| Existing labels | results/lincs/d127 and D122 baseline tables | Development-only controls, not independent evidence |
| Independent functional labels | PRISM/GDSC/CTRP if available and compound/line joins verified | Not used to select model variants, calibrate gene signs, or define disease signatures |

Data audit **before scoring**: check SHA/version, gene IDs (GENCODE v23/Entrez/symbol), duplicates, row/column shapes, network naming, number of targets in the directly measured 978 L1000 landmarks versus BING inferred genes, missing signs, target annotation confidence, drug-line counts, dose/time coverage, null/constant profiles, and target-genotype status in the source line. Record every exclusion, not just the remaining sample size.

**Scale constraint:** disease tumor-vs-normal logFC (or t-statistic) and LINCS MODZ z-scores are not on a common additive scale. Never calculate tumor logFC + LINCS z-score as a post-treatment counterfactual. The primary new readout is rank/sign based; a squared-distance rescue score is permitted only in a separately validated cross-platform calibration experiment.

**Assay reliability:** primary mechanism-readout comparison should report a landmark-only arm and a BING arm. BING improves regulon coverage but contains inferred targets; stratify directly measured vs inferred and avoid treating them as equally independent measurements.

## 2. Representation: drug-target anchored hierarchical regulons

Store sparse edges separately, then assemble regulator objects for explanation. A gene can be targeted by many TFs and a drug can have many targets.

- **B[r,g]:** signed TF to gene regulatory relationship, with an independent confidence weight; immutable between disease and drug contexts unless an explicitly labeled context-specific network is supplied.
- **D[c,g]:** disease-associated signed gene effect and confidence/reliability in cohort c.
- **A[c,r]:** disease regulator NES.
- **P[d,l,g]:** drug-specific LINCS gene-response z-score in line l, with dose/time/replicate metadata where available.
- **Q[d,l,r]:** the *actual VIPER* drug-induced regulator response in line l.
- **U[d,t]:** target unit t, known sign of drug action (+ for agonism, - for inhibition/degradation), binding/relevance provenance; unknown or mixed action kept explicitly.
- **Path[t,r]:** signed, directed mechanistic connection from drug target to downstream regulator, including number of hops, ambiguous signs, prior/source/degree.

The scoring unit is **(drug, cohort, source-line/profile, target unit)**; cohort aggregation reuses existing weights and records how many independent drug-line observations contributed. Output every target/regulator/gene contribution and its provenance.

### 2.1 Primary simple signed target-local gene rescue (R)

Let disease direction be s(c,g) = sign(logFC[c,g]), and wD(c,g) a frozen, bounded confidence/importance weight derived only from the disease contrast (e.g., moderated-statistic magnitude and repeatability). Define bounded drug response b(d,l,g) = tanh(P[d,l,g] / 2). The divisor 2 is an *initial proposed fixed scale*; check its empirical distribution in a label-blind input audit before freezing the confirmatory run. Do not tune it on drug efficacy labels.

For the set G(d,c) of unique disease-relevant genes supported by at least one target-anchored regulon, with gene weight wG(g) derived from edge confidence, anchor reliability, wD and direct-measurement confidence:

    R(d,c,l) = - sum_{g in G} wG(g) * s(c,g) * b(d,l,g) / sum_{g in G} wG(g)

R is in [-1,1], positive when genes move opposite the disease-associated direction. Each **unique gene is counted once** regardless of how many regulons include it. Report a companion weighted **per-regulon** score for explanation, but do not sum these naïvely into the global score.

If G is empty or lacks a pre-specified minimum of independently measured/reliable genes, **abstain**; do not silently return a neutral zero. Overly aggressive drug responses can yield large R even when cytotoxic; require M and viability checks for interpretation.

### 2.2 Primary target-mechanism concordance (M)

Anchor the graph at the *drug's annotated target*, **not** the top 25 affected regulators as in D126. If the target is itself a transcription factor, permit the 0-hop TF anchor. For kinase/receptor targets, use signed paths to TFs through a cached signaling prior, then signed TF→gene regulons. Main arm: signed paths at most two signaling hops, with 0/1/3-hop and SIGNOR-only sensitivity arms. Ambiguous shortest-path signs, unknown action types, missing edges and mixed unit signs are **not known contradictions**; treat them as ambiguous or unscorable and count them.

For unambiguous target→TF→gene predictions, the expected gene response sign is

    expected(d,t,r,g) = action(d,t) × signed_path(t,r) × tfmode(r,g)

Compare expected direction with the observed drug gene response using its bounded magnitude and confidence. Aggregate over **unique genes** (conflicting expectations become ambiguous, not selected opportunistically). This creates M(d,c,l) in [-1,1] for reachable, resolvable, measured genes; also report path coverage, sign-ambiguity fraction, median hop distance, number of independent anchored TFs, and largest single hub's contribution.

A high M indicates a drug-induced footprint consistent with the supplied *annotated route*; it is **not** binding/target engagement proof. A low M with very low coverage must not be interpreted as negative mechanism evidence. Use a target-degree matched null to detect hub-driven false positives.

### 2.3 Additional observations and ablations

- **Regulator-level rescue:** existing VIPER Q against disease A, both signed/weighted and existing spec_mech; this is a distinct baseline, not a reconstructed gene-level result.
- **Target-specificity:** compare M for the annotated target against (i) degree-matched shuffled targets, (ii) drugs with unrelated targets but matched TAS and viability-proxy response, and (iii) shuffled edge signs and gene labels preserving regulon sizes.
- **Coverage:** fraction of relevant, independently observed gene mass with known direction; fraction of expected target units represented; loss from landmark-only versus BING.
- **Agreement:** M versus target-proximal regulator Q and M versus R. Agreement can be mechanistically informative but circular if regulator Q derives from the same member genes; provide held-out-target or orthogonal-pathway sensitivity.
- **No score mixing yet:** keep (T, M, R, existing reversal) as a panel. Only after the pilot define one frozen integration arm, e.g. stratified/rank-transformed R added to the existing layered score with a *fixed* coefficient selected on development cohorts and evaluated on held-out cohorts. Always report layered_pct and R-only/M-only ablations.

## 3. Cohort design: start where mechanism, negatives and data are testable

This is **not** a simple cancer-type ranking benchmark. Some target-drug pairs require mutation, amplification, fusion or functional pathway status, not merely tissue membership. Do not assert target genotype from expression alone.

| Wave | Cohort / proposed contrast | Positive drugs/expectation | Negative or hard control | Scientific question / principal hazard |
|---|---|---|---|---|
| **A1** | **BRCA-LumA / BRCA-LumB vs BRCA-Basal**, with MCF7 source-line evidence annotated | ESR1 antagonists / SERDs/SERMs: fulvestrant, tamoxifen-class where measured; ER target regulon suppression | Same anti-estrogens in Basal; unrelated cytostatic drug in MCF7 | Direct TF-target footprint and **known D126 ESR1 miss**; media/ligand and SERM context matter |
| **A2** | **BRCA-Her2 vs BRCA-Basal**, with ERBB2 amplification/protein status checked; LumA as secondary comparator | Small-molecule ERBB2 drugs, e.g. lapatinib if it has an eligible LINCS signature | ERBB2 nonamplified tumors; dual EGFR/HER2 and proliferation-only confounding | Receptor→MAPK/PI3K→TF→gene propagation; PAM50 HER2-enriched is not equivalent to confirmed amplification |
| **A3** | **TCGA-SKCM BRAF V600E/K vs BRAF WT / NRAS-driven SKCM** (derive molecular strata); existing SKCM vs LUSC is development-only | Selective BRAF inhibitors (vemurafenib/dabrafenib where available) in mutant setting | RAS-driven/BRAF-wild-type setting, where some RAF inhibitors paradoxically activate MAPK | Mutation-conditioned sign reversal, target expression != mutant target action; distinguish BRAF from MEK inhibitors |
| **A4** | **TCGA-LUAD EGFR exon19del/L858R vs EGFR WT/KRAS mutant** (derive molecular strata); existing LUAD vs LUSC is development-only | EGFR TKIs, e.g. gefitinib/erlotinib/osimertinib if eligible and the relevant mutation is present | KRAS-mutant/EGFR-wild-type, resistant EGFR variants where annotated | Receptor-driver specificity, resistance, context transfer and sparse target-proximal expression footprints |
| **B1 challenge** | **THCA-BRAF-like-classical vs THCA-RAS-like-follicular** | BRAF-selective vs MEK inhibitor analyses **separately** | RAS-like disease is **not** expected to be universally MEK-insensitive | Explicit D127 failure (3/9 when pooled BRAF/MEK); distinguish mechanistic RAF paradox from shared downstream MAPK dependence |
| **B2 challenge** | **BRCA ER-positive, RB1-intact vs RB1-deficient/high-p16** (molecular subset, not LumA vs Basal alone) | Palbociclib/ribociclib/abemaciclib where profiled | RB1 loss, bypass via cyclin E/CDK2 | Functional pathway bypass despite drug hitting CDK4/6; resistance not predictable from CDK expression alone |
| **B3 challenge** | **TCGA-OV / BRCA-mutant or HRD-high vs HR-proficient**, with clinical context/repair status | PARP inhibitors, e.g. olaparib/niraparib when transcriptional data exist | HR-proficient or known acquired BRCA reversion if measured | Synthetic lethality and DNA-damage phenotype need not resemble return-to-healthy expression; do not interpret weak R as drug failure |
| **C extensions** | **TCGA-PRAD** androgen receptor-driven vs AR-independent/neuroendocrine; **ALK-fusion LUAD** vs fusion-negative; other tissues once coverage audited | AR antagonists; ALK TKIs only with fusion evidence | Lineage switches, receptor absence, target fusions not captured by mutation-only annotations | Test generalization beyond breast/MAPK when cohorts and compounds are genuinely scorable |

**First implementation set = A1, A2, A3, A4.** A1/A2 reuse existing project cohorts; A3/A4 need new patient genotype strata and may be **blocked by data availability**. A3/A4 must not be replaced silently by all-SKCM/all-LUAD. B1 is the **first stress test**, not a validation success requirement.

Important practical limits: antibody therapeutics are not interchangeable with small-molecule LINCS compound signatures; no data → not testable. For an EGFR or BRAF mutant patient cohort, a wild-type source LINCS cell line does not establish the corresponding on-target response. Mark line-target genotype as known, likely, unknown, or mismatched; perform in-context and borrowed analyses separately.

## 4. Incremental test ladder and explicit exit conditions

### Stage 0 — availability, joins, and **blind** coverage audit

Before reading any new functional outcomes: build one row per candidate cohort × drug class × annotated compound × source line, with available network, target direction, patient genotype, line genotype, gene coverage, dose/time, TAS, expected action, and matched/borrowed status.

Report: available positive and negative drugs, number of lines, 978 vs BING coverage, % signs resolved, target-degree distributions and fraction of source lines genotype-matched. **If a class has no scored compounds, no target-anchored path, or no measured downstream targets, label not-computable.** Do not fabricate zeros or quietly substitute another drug/class. Freeze a minimum coverage criterion after this label-blind audit and before outcomes. Suggested pilot guardrail: >= 10 usable target-linked genes and >= 2 independently anchored TFs for a kinase pathway; direct TF drugs may qualify with >= 10 regulon targets. These are proposed guards, not previously validated cutoffs.

### Stage 1 — deterministic unit / adversarial fixtures

Synthetic tests with expected exact outcomes:
1. Complete opposite-direction rescue; concordant and discordant on-target response.
2. All effect signs reversed; invariant to uniformly scaling edge weights.
3. Shared gene in multiple TF regulons counted only once; TF-level explanation preserved.
4. Missing/unknown target, contradictory drug action, multiple targets with opposing mechanisms → ambiguous/abstain.
5. No reachable target-linked TF; isolated target; zero/nonfinite profile and missing gene → unscorable, no crash.
6. Direct ESR1 transcription-factor anchor is included even when ESR1 absent from top 25 changed VIPER TFs.
7. BRAF WT + RAS-driven paradoxical MAPK activation: expected negative control not falsely marked as a clean BRAF inhibition footprint.
8. Landmark-only vs BING missingness, duplicate symbols/Ensembl mappings, and inferred-gene evidence flags.
9. Dose/time heterogeneity: do not average strong and weak or opposing doses without recording it.
10. Gene label shuffle, edge-sign shuffle, target swap, and degree-preserving null give near-null signals; test with seeded randomization.
11. Drug nonselectivity, multiple ChEMBL molecules / family targets and a protein complex do not multiply one molecular interaction.
12. Null results retain missingness, do not spuriously pass a score threshold.

### Stage 2 — direct LINCS mechanism sanity checks

Choose **eligible** source-line controls after audit: MCF7 anti-estrogen as mandatory, melanoma A375/BRAF only if the line's target genotype and compounds are verified, and at least one eligible HER2/EGFR example if possible. Measure per drug × line:
- M, R and coverage for target versus degree-matched alternative targets.
- Relevant downstream gene contributions and regulator Q; known-direction consistency on directly measured genes.
- Concordance across independent dose × time strata; separately test raw and generic-removed profiles.
- Difference between target-local and genome-wide/proliferation response; test against TAS-/stress-/viability-matched non-target drugs.
- Significance by drug- and line-aware resampling, not an indiscriminate per-gene t-test.

**Stop mechanistic promotion** if the MCF7 anchor still cannot detect the ESR1 footprint *despite adequate observed coverage*. Diagnose whether biology is absent from the assay versus network/score failure before changing cutoffs.

### Stage 3 — re-score original D127 pairs (development/compatibility only)

Read the same exact drug/cohort instances as scripts/d127/d127_biomarker.py; report **layered_pct**, D120 target_role, current spec_mech, candidate R, candidate M and coverage with paired high-minus-low scores. Use the original D127 primary/seens labels explicitly; no new independence claim. Preserve its inclusion convention: do not exclude the expected-low cohort because target_absent/local_route/low_exposure flagged it (report end-to-end separately).

Report per-class and per-pair results plus pooled descriptive metrics. **Do not tune the new score on the 49 pairs then call its D127 performance confirmatory.** Treat thyroid failure and the MCF7 control as diagnostic tests; rank changes alone are insufficient.

### Stage 4 — new within-tissue molecular holdouts (first substantive validation)

Derive genotype/amplification-stratified patient subsets and confirm each signature and network have adequate support. Primary candidate cases: A3 BRAF mutant vs WT/NRAS SKCM and A4 activating EGFR-mutant vs EGFR WT/KRAS-mutant LUAD. Secondary: ERBB2-amplified vs nonamplified BRCA and RB1-intact vs deficient breast.

For each compound measured/scored in both cohorts, compare:
- Direction and size of **within-cohort percentile gap** (high minus low).
- Same-drug paired score difference, with equal scoring population definitions.
- Mechanism-vs-role attribution: does M/R explain an improvement **beyond** D120 alterations and D122 role?
- Same tissue and, when possible, same source-line family to minimize generic tissue effects.
- Gene overlap, mutation/fusion status and line genotype; stratify direct/borrowed instead of pooling.

Designate new holdouts *before reading their score outcomes*. If prior cohort/compound exploration already used a pair, mark it seen and choose a genuinely new independent validation cohort.

### Stage 5 — external phenotype checks and negative controls

Intersect eligible drugs/lines with PRISM, GDSC or CTRP, with independently verified compound and cell-line IDs. Use dose–response AUC/IC50 or viability response in **biomarker-defined** sensitive vs resistant lines; adjust for lineage and broad drug toxicity. Do not use DepMap Chronos dependency both as a model input and as proof of independent target relevance.

Negative controls include: target absent or genotype incompatible, pharmacologically inactive context, unrelated drug with matched stress intensity, incorrect annotated target, reversed mechanism sign, highly promiscuous compound, pure cytotoxic/proliferation signature and a drug with high R but no mechanistic concordance. Include at least one **failure by design** (BRAF RAF paradox, ESR1 top-25 failure, RB1 loss, and PARP transcription/phenotype mismatch).

### Stage 6 — robustness and prospective generalization

Sensitivities: ARACNe vs CollecTRI/DoRothEA; OmniPath vs SIGNOR; 0/1/2/3-hop; landmark vs BING; leave-one-source-line-out; direct vs borrowed; leave-drug-out/leave-class-out; dose/time matching; exclude proliferation/viability-linked genes; random graph degree controls; decompose target units and binding confidence. Compare unweighted vs confidence-weighted and target-local vs all-regulon rescue. Do not promote any arm that changes sign or passes only with inferred genes or a high-degree hub unless explained and independently corroborated.

## 5. Primary readouts, baselines, inference and failure criteria

### Baselines

B0 D122 **layered_pct** lead score; B1 D120 **target_role_pct** where defined; B2 existing **spec_mech_pct_pregate** / VIPER reversal; B3 global gene-level reversal using the same LINCS specific profiles but **without target anchoring**; B4 target-only enrichment/annotated-signaling footprint without disease rescue; B5 degree/size-matched random target+edge nulls. Candidate arms: **R-only**, **M-only**, and **target-local R with M diagnostic**. Do not compare a new combined score only to a weaker historical baseline.

### Report for every cohort and mechanism class

- Eligibility counts, unavailable fraction, target-sign coverage, measured/inferred fraction, and exclusions with reason.
- Within-pair correct direction and median high-minus-low *percentile* gap; paired bootstrap CIs at compound and mechanism-class level.
- Regulator and target-gene sign concordance; quantitative transfer correlation **relative to per-representation reproducibility** if using PR #1's held-out pair evaluation. Cross-line fidelity is not proof of tumor efficacy.
- Functional validation on external viability data (where available): paired sensitivity difference, class-aware AUC/AUPRC, effect size, CIs.
- Top changed drug rankings and decomposition by target, regulator and member gene, **including which expected drug negatives move up**.
- Stratified performance for line match, transfer confidence, target genotype, L1000 feature space and cytotoxicity proxy.
- Degree/null z-scores and results with proliferation genes removed. Flag potential target-annotation and pathway circularity.

Avoid pseudo-replication: the 49 D127 compound instances are not 49 independent biological mechanisms; classes and overlapping compounds share labels. Report class-balanced macro averages and cluster bootstrap (compound, class, cohort pair). Pre-register any statistical test for a new holdout; use descriptive uncertainty for the already-seen D127 rerun.

### Proposed promotion policy (to freeze only after Stage 0 audit)

A candidate **must** (i) have sufficient coverage for at least the four A-wave mechanisms or transparently abstain, (ii) pass both anti-estrogen and BRAF genotype/context negative controls where testable, (iii) improve an independent biomarker-conditioned paired direction metric over layered_pct by a prespecified margin with CI not indicating meaningful regression, (iv) retain the current D127 known positives and not increase unexplained pan-cancer cytostatic dominance, and (v) show supporting direction in at least one external functional screen without requiring its labels for calibration. Mark "exploratory improvement" if only development pairs improve. If early controls fail, **retain production D122** and archive the candidate analysis instead of repeated threshold tuning.

The exact numerical improvements/uncertainty bands and cohort counts are **open pilot decisions** and should be frozen once availability is audited, not invented after results appear. Include a hard separation between **feasibility gate** and **biological validation gate**.

## 6. Edge-case register (minimum required coverage)

| Case | Why an aggregate score misleads | Required behavior |
|---|---|---|
| Drug target absent/low protein; gene expressed | Expression != engagement | Keep flag, score T separately, M unscorable without a footprint |
| Target mutation but no upregulated transcript (BRAF, EGFR) | Driver activity is not expression abundance | Use variant and downstream pathway, not target logFC alone |
| Oncogenic receptor fusion without structural data (ALK/ROS1/RET) | Missed by MC3 missense + GISTIC | Hold out/unscorable until fusion evidence |
| Direct TF target missed among top altered regulators (ESR1) | D126 failed | Force target-anchored footprint; report coverage |
| Inhibitor, agonist, partial agonist, degrader or mixed MoA | Single - sign is invalid | Explicit signed action or abstain/mixed |
| A drug inhibits a kinase but its mRNA rises | Protein inhibition != mRNA repression | Test proximal pathway response, not kinase mRNA sign |
| Shared genes in several regulons | Double count / correlated evidence | Unique-gene aggregation, separate explanatory TF scores |
| Conflicting target→TF paths | Spurious sign assignment | Ambiguous; retain as measured uncertainty |
| High-degree signaling hubs | Generic coverage advantage | Degree-matched null and hub-contribution diagnostic |
| Generic growth arrest / cell death | Inflated reversal, false mechanistic similarity | Matched-stress negatives, proliferation ablation, retain 79c specific profile |
| Dose/time reversal and toxic doses | Average masks biology | Stratum-level metrics and dose/time stability |
| Drug target selectivity / polypharmacology | Wrong mechanism attribution | Unit-level contributions; multi-target flags; potency sensitivity |
| Context mismatch / borrowed lines | Sign may invert | Separate direct/borrowed; genotype- and lineage-aware confidence |
| Mutant-specific inhibitor paradox (RAF) | Opposing effects in BRAF WT/RAS active | Use mechanism-specific contextual negative |
| Resistance bypass (RB1 loss, AR variant) | Target hit no longer blocks growth | Functional-state negatives; never infer response solely from target |
| Synthetic lethality (PARP/HRD) | Beneficial killing not normal-expression restoration | Separate R from viability outcome; do not treat low R as failed therapy |
| Immune, microenvironment, or antibody mechanisms | L1000 cell lines incomplete | Outside scope or separate modality arm |
| Cell-type composition / GTEx tumor-normal shifts | Apparent disease genes are stromal/lineage | Purity/single-cell or stable-subgroup sensitivity |
| Landmark vs inferred BING genes | Shared imputation error | Separate feature-space arms and reliability weights |
| Annotation/label leakage | Biomarkers used in both selection and evaluation | External holdout, frozen joins and tests |

## 7. Minimal implementation design and output files

Create a **new sandboxed experiment** (candidate D128), not an edit of 79e or the production lead scorer until accepted.

Proposed files:
- scripts/d128/a0_audit.py — cohort×drug×line×target coverage/eligibility and provenance;
- scripts/d128/a1_score.py — sparse, target-anchored gene/regulon M and R with per-gene attributions;
- scripts/d128/a2_controls.py — synthetic and experimental positive/negative controls;
- scripts/d128/a3_d127_compat.py — exactly same D127 instances; frozen development comparisons;
- scripts/d128/a4_genotype_pairs.py — genomic strata definition/joins and new holdouts;
- scripts/d128/a5_external.py — independent PRISM/GDSC phenotype tests;
- scripts/d128/config.py and tests/test_d128_*.py — thresholds, paths, fixtures and null invariants.

Output root: results/lincs/d128/ with:
- coverage.tsv, eligibility.tsv, source_line_context.tsv;
- gene_contributions.tsv.gz, regulator_contributions.tsv.gz, target_units.tsv;
- scores_by_line.tsv, scores_by_cohort.tsv (T/M/R, coverage, abstain reason);
- control_results.tsv, d127_seen_pairs.tsv, genotype_holdout_pairs.tsv;
- external_viability.tsv, nulls.tsv, sensitivity_arms.tsv, summary.json;
- run_manifest.json (git SHA, data checksums, eligibility criteria, score equations, seeds, environment, timestamp and frozen hypothesis/threshold version).

**Execution order:** A0 audit → freeze feasibility rules → A1 synthetic unit tests → A2 mechanistic controls → A3 D127 compatibility → freeze a single candidate formula → A4 genuinely independent patient molecular holdouts → A5 external functional validation → A6 robustness/adoption decision. Do not optimize to D127 results.

## 8. Literature anchors (not a substitute for source-specific data audit)

- Subramanian et al., *Cell* 2017, L1000 978 landmarks and best inferred expression: https://doi.org/10.1016/j.cell.2017.10.049
- Szalai et al., 2019, proliferation/death signatures confound perturbational MoA: https://pubmed.ncbi.nlm.nih.gov/31552418/
- ER-positive breast pharmacology, SERMs vs SERDs and context: https://pmc.ncbi.nlm.nih.gov/articles/PMC12018393/
- Lapatinib response depends on HER2 amplification: https://pubmed.ncbi.nlm.nih.gov/19620495/
- Chapman et al., *NEJM* 2011, vemurafenib in BRAF V600E melanoma: https://doi.org/10.1056/NEJMoa1103782
- Hatzivassiliou et al., *Nature* 2010, RAF inhibitor paradox in WT/RAS contexts: https://doi.org/10.1038/nature08833
- EGFR-mutant lung cancer gefitinib controlled trial: https://doi.org/10.1056/NEJMoa0909530
- TCGA integrated thyroid carcinoma, BRAF-like vs RAS-like: https://unclineberger.org/peroulab/wp-content/uploads/sites/1008/2015/04/TCGA-Thyroid-CELL-2014.pdf
- Palafox et al., *Nature Communications* 2022, RB1 loss/high p16 and CDK4/6 resistance: https://pubmed.ncbi.nlm.nih.gov/36071033/
- ASCO ovarian PARP guidance/HRD: https://pmc.ncbi.nlm.nih.gov/articles/PMC8942301/
- Liu et al., CARNIVAL causal mechanistic footprint inference: https://pmc.ncbi.nlm.nih.gov/articles/PMC6848167/
- GDSC and PRISM external cell-line sensitivity resources: https://doi.org/10.1093/nar/gks1111 and https://www.theprismlab.org/consortium-screens/data-analysis

## 9. Open decisions to resolve in the Stage 0 report (not by guessing)

1. Which drug classes have actual eligible LINCS small-molecule compounds and genotype-appropriate source lines? Inventory before committing specific controls.
2. Does patient-level MC3/GISTIC permit SKCM BRAF V600 and LUAD EGFR driver stratification with enough samples to rebuild disease gene signatures/MRs? The repository's pre-existing subgroups are not sufficient by themselves.
3. Is source-line genotype known for the relevant compounds (MCF7/ER, melanoma/BRAF, lung/EGFR)? What fraction is unknown?
4. What fraction of the new footprint is directly measured 978 landmarks vs BING inferred genes? Can controls be scored with >=10 independent targets?
5. Which matched-source data can support independent dose/time subprofiles versus one collapsed LINCS mean?
6. Which genuinely independent validation cohorts/sensitivity labels have **not** already influenced feature selection, cohort choice or parameter tuning?
7. After reviewing feasibility **without efficacy outcomes**, freeze the minimum-coverage rule, one primary integration formula (if any), effect-size acceptance bar, and the list of confirmatory molecular holdouts in a separate commit.

**Decision principle:** If the score cannot explain a target-linked mechanism in a high-information control without generic growth arrest, it should remain an exploratory diagnostic, not an extra ranking layer.