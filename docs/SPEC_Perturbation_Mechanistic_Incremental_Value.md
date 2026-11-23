# Research Specification — Does Perturbation-Derived Regulatory Information Add Disease-Relevant Mechanistic Evidence?

**Project:** `kirtihargunani-cmu/repurposing`  
**Status:** implementation-ready **proposed pre-registration**, not yet registered, executed, or adopted. **Amended 2026-11-23 (Amendment 1, directly below): it overrides any conflicting text in sections 0-15.**  
**Prepared:** 2026-11-23.  
**Primary scientific question:** *Can observed drug perturbations and their ARACNe-derived regulatory effects identify disease-relevant drug mechanisms that are not already explained by drug-target annotations and baseline disease biology?*  
**As amended:** *How much drug-mechanism signal can the production scorer resolve from noisy, incomplete perturbation data, step by step, and where is it resolved?*

> **Agent handoff:** Execute phases in order. Do not change production rankings, gates, or thresholds. Do not tune analysis choices after examining outcome results. Write completed results, negative controls, limitations, and interpretations for every test, including tests that cannot be executed. Commit a dated protocol and frozen configuration **before** running any new decision-bearing test, and record code SHA, input hashes, output manifests, data counts, and exact command lines. This document is a proposed protocol; committing it is not evidence that the experiments have run.

---

## Amendment 1 — 2026-11-23 (owner decisions; overrides any conflicting text below)

Owner, 2026-11-23: run it *"on the most recent results, use our new scoring function (filter with D/A/O and then gene
level regulon reversal), maintain our separated timepoint windows. let's not run PRISM/GDSC selective killing ...
report results optimistically, measure the amount of signal we can actually resolve via a processing step. take into
account the data are fundamentally noisy and we also have a limited amount/missing information so really any
information we can extract is a good sign for further analysis."* This amendment also resolves findings 1-9 of
`docs/AUDIT_SPEC_Perturbation_Mechanistic_Incremental_Value.md`.

### A1.1 Inputs: the most recent full run only

Every input comes from the 2026-11-23 build of `scripts/run_pipeline.sh` (13:45-15:17; reproducible per D210), or
from the disease-side outputs that run consumed unchanged (`results/step6/aracne_project_v2_tf[_strata]/group_mr.tsv`;
per-gene breakdown `regulon_contributions.tsv.gz`, built 2026-11-22). `inputs_manifest.json` records each input's
`.build_manifest.json` run id; Stage 0 **refuses** an input whose manifest predates that run (other than those
disease-side inputs). Exploratory outputs built outside the runner (D196-D206: programs, NES ambiguity, two-stage,
covariance, `interpretation/`) are not inputs; programs (Experiment E) are rebuilt through `runlog.sh` first. If any
production stage is rebuilt during the study, the study re-reads from the new run and records both run ids.

Reference counts from that run (to be re-verified in Stage 0): 3,721 compounds x 41 cancers (64,094 rows) plus 14
subgroups; timepoint data for 29 cancers (4,128 drug rows with both windows, 5,823 early only, 7,509 late only) and
8 of 14 subgroups (2,402 rows); 4,418 eligible (cancer, ingredient) pairs, 615 ranked early over 20 units and 882
late over 33; 3,382 eligible pairs with no timepoint data.

### A1.2 The scorer under study is production (D211)

- **Filter:** a (unit, compound) pair is eligible when `layered_T >= 50` on the D/A/O path (`79t3`; D120, D177).
  Chemotherapy / DNA-acting / essential-gene rows (no D/A/O) are analysed separately, never mixed in.
- **Score:** gene-level regulon reversal `rev_gene = - sum_i sum_j d_ij * C_ij / ||s||` (`docs/FORMULAS.md` s3),
  per (unit, window), percentile among eligible drugs (`79ar`; `results/lincs/gene_ranking_v2_tf.tsv`). Its exact
  per-regulator terms are `rev_gene_contrib` in `regulator_contributions_by_window_all[_strata].tsv.gz`.
- **Roles:** D, A, O and T are the *filter*, so none of them is ever an outcome. Regulator-level `rev` and `rev_mech`
  are reported as reference points only. Production is not changed by this study.
- The universe is FDA-only by construction (target evidence exists for FDA, gate-passing rows only).

### A1.3 Timepoints stay separated

Early (<= 6 h) and late (>= 24 h) are analysed separately at every step; no pooled score enters any analysis; a pair
with no timepoint data stays `unknown` with its D180 reason and is counted in coverage, not dropped. Subgroups are
included wherever their window run has data.

### A1.4 No PRISM/GDSC selective killing

The functional outcome of Experiment A is removed: section 2.2 `outcome.tsv` (PRISM/GDSC part), 2.3 (D206 rebuild),
4.2 `U_A_sens`, 5 A0-A2 (sensitivity outcome), A4, A6 (sensitivity coverage gate and +0.05 rule), and the D206 /
efficacy instructions in section 14. `scripts/79k_sensitivity_validation.py` is not run. Nothing in this study is an
efficacy claim.

### A1.5 Experiment A, redefined: how much signal each processing step resolves (PRIMARY)

Follow one drug-cancer pair through the production scorer and measure the signal that survives each step, per
window:

| Step | What it is | Source |
|---|---|---|
| S1 | Condition profile (LINCS level 5, activity floor 0.25 at 6 h, 0.40 otherwise) | `drug_condition_profiles_*` |
| S2 | Drug regulator activity on the cancer's ARACNe network | `drug_viper_specific_v2_tf_cond/` |
| S3 | Window aggregate over the cancer's matched in-key lines | `window_reversal[_strata].tsv.gz` |
| S4 | Regulator-level reversal `rev` (reference) | `regulator_contributions_by_window_all*.tsv.gz` |
| S5 | Gene-level reversal `rev_gene` | same, `rev_gene_contrib` |
| S6 | `rev_gene` within the T >= 50 eligible set (production) | `gene_ranking_v2_tf.tsv` |

At each step, three measurements, each with a cluster bootstrap 95% interval over independent drug / target-unit
clusters:

1. **Reproducibility.** Split-half reliability (Spearman-Brown) where >= 2 conditions exist, and the attainable
   ceiling `sqrt(reliability)`. This is the noise floor the later steps work against.
2. **Signal above a matched null.** Separation from permutations that preserve drug, line, window, condition count
   and target-unit clusters (section 3 rules).
3. **Agreement with mechanism labels that share no input with the filter or the score:** (a) the drug's own primary
   target ranks high in its regulator response (D158 style, enrichment over chance in the measured regulator
   universe); (b) on-target pathway moves in the expected direction (PROGENy, D99 style); (c) drugs sharing a primary
   target unit look more alike at the gene level than matched drugs that do not; (d) D127 biomarker pairs,
   descriptive only (already seen in D211).

The headline is a **resolvable-signal curve**: for each step, the signal as a fraction of the step's reliability
ceiling, with intervals, separately for early and late and for groups and subgroups. A step that keeps or sharpens
signal is the finding; a step where the curve falls to the null shows where the data stop resolving it.

### A1.6 Experiment B: an empirical noise model

Only 10.1% of conditions carry >= 2 signatures (median 1) and 400 drug x line pairs have >= 2 conditions in both
windows, so B's per-row confidence rule is replaced: estimate the noise distribution from split halves and
replicate signatures on those 400 pairs, apply it as a calibrated uncertainty to single-condition rows (including the
2,639 matched drug x line x dose pairs), and report the false sign-flip rate from the pseudo-replicate null. The
transition classes in B2 are then assigned under that model.

### A1.7 Reporting stance

The data are fundamentally noisy (median one signature per condition) and incomplete (12 of 41 cancers have no
window run; 3,382 of 4,418 eligible pairs have no timepoint data), so **any resolved signal is a good sign for
further analysis** and is reported as such:
- Lead every report with what was resolved: effect, interval, denominator, and its size relative to the
  reliability ceiling.
- **Decision rule:** a quantity counts as *resolved* when its 95% interval excludes the matched null. The practical
  thresholds elsewhere (+0.05 deltas, odds ratio >= 1.5, FDR < 0.05 enrichment bars) become descriptive benchmarks,
  not pass / fail gates.
- A result that does not separate from the null is reported as *not resolved at current data depth*, together with
  the noise and coverage that limit it, and the data that would resolve it; never as evidence that the mechanism is
  absent.
- What stays from the original protocol: rules are fixed before results, controls are run and reported beside every
  positive result, and no threshold moves after a number is seen.

### A1.8 Superseded text (for readers of the sections below)

Section 0.2 row A (now A1.5); section 1 row D202-D205 (`rev_gene` is production since D211); section 2.3 subgroup
bullet (8 of 14 subgroups have window data) and D206 bullet; section 4.2 `U_A_sens`; section 5 outcome parts and A6;
section 11.1 "+0.05 ... underpowered" language (A1.7 applies); section 13 rows on A ("do not add reversal to the lead
score": reversal is already production); section 14 (replaced below). Section 12.1 (feature branch and PR versus the
project's commit-to-`main` convention) is **still open** for the owner.

### A1.9 Handoff prompt (replaces section 14)

> Implement `docs/SPEC_Perturbation_Mechanistic_Incremental_Value.md` as amended by Amendment 1. Use only the outputs
> of the most recent full run of `scripts/run_pipeline.sh` (record its run ids; refuse older inputs). The scorer under
> study is production D211: T >= 50 on the D/A/O path, then gene-level regulon reversal `rev_gene`, per timepoint;
> D/A/O/T are never outcomes. Keep early and late separate everywhere and keep unknown windows as unknown. Do not run
> PRISM/GDSC or any efficacy endpoint. Run Stage 0, then the resolvable-signal analysis (A1.5) for steps S1-S6 with
> reproducibility, matched-null separation and mechanism-label agreement, then B with the empirical noise model
> (A1.6), then C, D, E and F. Report optimistically per A1.7: lead with resolved signal and its share of the
> reliability ceiling; describe unresolved results as limited by noise and coverage. Fix rules before results, run
> every control, change no production score.

---

## 0. Executive decisions and boundaries

### 0.1 What this project is **not**

- Not another undirected search for a drug score that boosts approval-recovery AUC.
- Not a rerun of D183 missing-window imputation; D183 was closed under its registered acceptance rules.
- Not an assertion that ARACNe network edges are causal drug-target signaling interactions.
- Not an attempt to infer target engagement solely from a VIPER activity value.
- Not a retraining experiment in which unapproved or unstudied disease–drug pairs are labelled ineffective.
- Not a reintroduction of the ARACNe–CollecTRI consensus, which was removed in D208. CollecTRI/DoRothEA are archived, not production arms.

### 0.2 Priority-ordered experiments

| Stage | Main deliverable | Distinct question | Priority |
|---|---|---|---|
| 0 | Frozen audit, data dictionary, and eligible analysis universes | Can we establish valid paired, independent observations? | Blocking |
| A | Incremental-value analysis and error taxonomy | Does the perturbation signal add information **conditional on** target/disease evidence? | 1 |
| B | Temporal mechanism decomposition | Which observed gene/regulator changes explain early/late divergence? | 2 |
| C | Primary/secondary target recoverability | Do confidently engaged secondary targets become detectable once exposure, expression, and downstream effects are accounted for? | 3 |
| D | Network representation sensitivity | When a regulator changes interpretation between contexts, which regulon edges cause it? | Follow-on |
| E | Disease-program accessibility | Which parts of disease-defined programs are reproducibly pharmacologically perturbable? | Follow-on |
| F | Integrated synthesis **without changing production ranking** | What evidence is present, missing, contradictory, or condition-specific for each candidate? | After A–E |

**Decision rule:** A negative result closes only the precise tested question on its eligible universe and comparator. Do not generalize a null from matched LINCS lines to untested tissues or from growth viability to all possible therapeutic phenotypes.

### 0.3 Required outputs of *each* stage

1. `protocol.md` (specific hypothesis, universe, locked choices, controls, criteria, amendments).
2. `inputs_manifest.json` (paths, SHA256, sizes, row counts, field names, units, cohort coverage, source SHA).
3. `analysis.tsv[.gz]` with *one independent observation per documented unit* or a clear parent identifier for clustering.
4. `summary.tsv` with effect sizes, sample sizes, bootstrap intervals, statistical tests, and missingness.
5. `controls.tsv` including deliberately adverse/negative controls.
6. `figures/` with at least a coverage plot, main paired effect plot, negative-control plot, and a selected mechanistic worked example when applicable.
7. `report.md` describing what was done, numbers found, what the result supports, what it does **not** support, and next conditional action.
8. `run.json` including git SHA, arguments, parameter registry check, environment, deterministic seed, timestamps, and resolved file hashes.

A stage with insufficient data still produces `coverage.tsv`, `not_evaluable.tsv`, a short `report.md`, and a `run.json` with `status=not_evaluable`; it does not silently disappear.

---

## 1. Established facts — do not re-test as though unknown

This is a **baseline inventory**, not a reproduction claim. The implementing agent must verify current input files and relevant decision records before using historical numbers.

| Record | What was tested or discovered | Established interpretation | Constraint on new work |
|---|---|---|---|
| C28/C29 | Approval recovery, PRISM/GDSC sensitivity, and DepMap target dependency | Reversal did not reliably track known use, selective killing, or dependency; approved drug targets did carry disease-specific dependency information | Distinguish mechanistic ranking from efficacy prediction; dependency is not a fully independent anchor if reused in a score |
| D99 | PROGENy pathway checks | Many expected pathway signs were recovered, but ARACNe ESR1 read anti-estrogens incorrectly | Use target/pathway perturbation controls; do not equate one ARACNe NES with actual target engagement |
| D128–D133 | Landmark/inferred gene dependence, permutation nulls, regulator correlations | Many features are correlated; naive independence-based significance can be inflated | Cluster/bootstrap at independent biological unit; nulls must respect dependence and group structure |
| D137 | Seven-link mechanistic tracing with negative controls | The procedure failed: negative controls could score more convincingly than real mechanisms | Do not recycle the tracing algorithm as proof of causal paths; require falsification controls |
| D154/D156/D157 | Proliferation slopes and generic/universal/line-selective decomposition | Proliferation-response components vary with time and can be drug-/line-selective | Do not automatically subtract all proliferation response or label all proliferation-linked effects artifacts |
| D158 | Exact-target recovery from regulator opposition | Primary targets enriched (top-1% 3.84× chance); secondary targets near chance | Analyze secondary targets after checking binding/exposure/expression and downstream effects; do not call anecdotal hits discovery |
| D169–D183 | Window separation and imputation | Early/late disagreement ~37%; D183 models failed registered imputation criteria despite measurable reliability | Work with **observed** windows only; keep unknown as unknown; never impute in this project |
| D179/D182/D193–D195 | Exact regulator→gene contribution tables | Regulator and target contributions decompose measured scores exactly; window-specific versions exist | Reuse exact contributions; do not reconstruct from unweighted z-scores |
| D192 | Per-cohort approval-recovery | Cohort AUC strongly heterogeneous and sometimes inverted; matched-row window resolution did not improve approval AUC | Do not use overall approval AUC as primary discovery test |
| D196–D201 | Program representation | Disease/drug program overlap low; the specific clean-program scorer failed; program coverage incomplete | Interpret program coverage; do not repeat the failed D201 ranker |
| D202–D205 | Gene-level score and covariance | `rev_gene` failed ranking; isolated covariance had target-role associations, but D205 combined test was invalid due to scale/file-clobber issues | Treat covariance as explanatory only; regenerate and independently verify provenance before any quantitative claim |
| D206 | Drug sensitivity table | Non-reproducible table: producer script/sign convention absent | Rebuild a signed, logged, unit-checked viability benchmark before using it for decisions |
| D208 | Production cleanup | `_v2_tf` ARACNe-only; stale consensus reads and missing evidence flags repaired; archived CollecTRI/DoRothEA and stand-ins off | Fail fast on consensus suffix; preserve `SPLIT_SFX=_v2_tf` |
| D209 | Target-deduplicated D127 | 16/20 independent target units correctly ordered (p=0.0059) on seven primary biomarker comparisons | Preserve unit-level independence; 20 units are not a broad clinical validation set |

**Important independence warning:** D120 `D` uses DepMap dependency and D120 `A` uses alteration. Neither can serve as an *independent* validation outcome for a predictor that includes the same component. A positive correlation with one of those components is a mechanistic alignment check, **not** out-of-sample validation of therapeutic efficacy.

---

## 2. Required source-of-truth and data contract

### 2.1 Repository files and scripts to inspect before implementation

The paths below are observed in repository source/docs as of 2026-11-23; the agent must resolve current versions and document any renames. **Never invent a missing file or column.**

| Source | Purpose | Minimum fields / checks |
|---|---|---|
| `docs/DECISIONS.md`, `docs/PROGRESS.md` | Latest method decisions and negative results | D120, D127, D137, D154–D158, D169–D184, D192–D209 |
| `config/adopted_parameters.yaml`, `scripts/check_params.py` | Adopted vs design-only settings | Verify before experimental builds |
| `scripts/fda_rows.py` | FDA ingredient grouping and current target-gated row selection | Confirm `SPLIT_SFX=_v2_tf`; log ingredient ↔ pert_id mapping |
| `results/lincs/dependency/gene_dep.tsv.gz`, `drug_dep.tsv.gz` | DepMap gene and drug target dependency | Cohort, gene/compound, dependency convention, screening coverage |
| `results/lincs/target_role/` and D120/D121 outputs | Disease-target dependency/alteration/overactivity and layered units | Discover actual filenames, exact component definitions, and target-role provenance |
| Current D122/79t/79t3 layered score outputs | Disease-target compatibility, **without reversal** | Discover actual filenames and exact raw/percentile fields; test that target unit scores are drug-independent when expected |
| `results/lincs/d127/{instances,pairs,summary,dedup_units,dedup_pairs,dedup_summary}.tsv` | Independent biomarker-unit reference | Verify source SHA and D209 controls |
| `results/systems_validation/integrated/window_reversal.tsv.gz` | Observed time-separated reversal | Disease/cohort, `pert_id`, `rev_early`, `rev_late`, window coverage, line/condition counts |
| `scripts/79k_window_long.py` and its outputs | Explicit early/late/unknown long representation | Check uniqueness of `(cohort,pert_id,window)` |
| `data/processed/lincs/drug_condition_profiles_index.tsv` | Conditions and coverage | `pert_id`, `ModelID`, `pert_dose`, `pert_time`, `n_sig`, `median_tas`, `generic_removed` |
| `data/processed/lincs/drug_condition_profiles_specific.f32` | Condition × gene signatures | Read with *documented* shape, dtype and feature map; do not assume shape |
| `data/processed/lincs/drug_viper_specific_v2_tf_cond/` | Condition × ARACNe regulator activity | Recover network-specific regulator ordering and binary shape |
| `results/lincs/regulon_contributions_v2_tf/regulator_contributions_by_window_all.tsv.gz` | Every regulator × drug × window contribution | `cohort`, `pert_id`, `window`, `regulator`, `disease_nes`, `drug_activity`, `rev_contrib`; optional `cos_contrib`, reliability fields |
| `results/lincs/regulon_contributions_v2_tf/target_contributions_by_window.tsv.gz` | Selected exact gene contributions | `regulator`,`target`, `drug_z`, `activity_contrib`, `disease_contrib_nes`; top-25 restriction **must be recorded** |
| `scripts/79e3_window_contributions.py` | Full exact reconstruction and targeted `QUERY` mode | `ALL=1` for complete regulator export; `TOP_REGS=0` or targeted query to avoid top-25 truncation when necessary |
| `results/lincs/programs/{programs,program_filter,program_drug,program_overlap}.tsv*` | Disease programs and flags | Determine members, program weights/aggregation, filtered coverage |
| `results/lincs/interpretation/` | D203–D206 prior outputs | Mark provisional if provenance not reproducible; never treat D206 as gold standard |
| `results/lincs/d134/`, `docs/SPEC_target_recoverability.md` | Prior target recoverability feasibility | Check S0 effective drug/line pairs and avoid re-counting annotation-only pairs |
| `docs/SPEC_dose_response_bmd.md` | Dose-response feasibility/plan | This document was marked design-only; do not assume BMD outputs exist |
| `results/lincs/drug_annotations.tsv`, ChEMBL cache, `results/lincs/d109/cmax.tsv` | Target/binding/clinical exposure annotations | Keep target class, type of potency, units and assay system; do not conflate Cmax with in-vitro concentration |

If a source is absent or stale, identify the producing script from the repository and reproduce it **only inside a namespaced sandbox** with build logs and manifests. No production file overwrite. Track `_v2_tf` consistently. Do not use archived CollecTRI/DoRothEA without separately registering a nonproduction comparison.

### 2.2 Normalized analysis tables

Construct the following **read-only canonical analysis views**, stored as Parquet (or compressed TSV with equivalent schema) with unique keys and strict joins.

**`drug_disease.tsv`** — one row per `(cohort, ingredient_id, pert_id, target_unit_id)` with target evidence and pooled production metadata. Unit: one scored LINCS compound; ingredient aggregation is a **separate** explicitly documented level. Fields: `cohort`, `tissue_key`, `pert_id`, `ingredient_id`, `target_unit_id`, `layered_raw`, `layered_pct`, `dependency_gene/target`, `alteration`, `overactivity`, `target_presence`, `source` (`matched`, `borrowed`, etc.), `line_match`, `transfer_conf`, `low_exposure`, `target_absent`, `local_route`, `n_lines`, `n_targets`, and feature-missing reason. Where multiple targets exist, retain child `(cohort,pert_id,target)` rows; do not collapse targets to max without an explicit sensitivity comparison.

**`drug_disease_window.tsv`** — one row per `(cohort,pert_id,window)` for `window∈{early,late}` **observed only**, plus a **separate** unknown-coverage inventory. Fields: reversal native units, window-relative percentile, condition and line counts, dose summary, TAS, gene/regulator coverage, and flags. No duplicating pooled scores as early or late.

**`drug_line_condition.tsv`** — one row per `(pert_id, ModelID, dose, time, replicate-or-condition-id)` with data provenance, response-quality diagnostics, and matched-line applicability. Normalize all dose units to µM in an *additional* field without destroying source units.

**`drug_target.tsv`** — one row per `(pert_id, target_gene, source_assay)` with primary/secondary status, measured Ki/Kd/IC50/EC50, potency value and unit, binding/inhibition modality, assay confidence, tissue/line target expression, and engagement plausibility. Missing potency = `unknown`, not zero.

**`regulator_contribution.tsv`** — one row per `(cohort,pert_id,window,regulator)` with `disease_nes`, drug activity, exact reversal contribution and optionally reliability/condition counts. **`edge_contribution.tsv`** — one row per `(cohort,pert_id,window,regulator,target_gene)` with signed disease and drug contributions and edge metadata; record whether full regulon or top-25-only.

**`outcome.tsv`** — one row per truly observed independent outcome `Y` with explicit assay ID, line identity, compound identifier, metric units, **beneficial-direction sign**, dose, treatment duration, sample size, and provenance. Cancer selective killing in PRISM/GDSC and CRISPR gene dependency **must not be pooled as if they were the same outcome**.

### 2.3 Data audit checks (blocking)

- All joins state key uniqueness, cardinality, and unexpected row gains/losses. Fail if any join silently multiplies observations. Assert `(cohort,pert_id,window)` unique in window table.
- Canonical drug mapping: ingredient, aliases, salt forms, multiple BRD IDs, all target unit IDs. Measure ambiguous IDs and count unresolvable drugs. Preserve 79aa/79ab ChEMBL resolution behavior; never convert absent indication to negative outcome.
- Numerics: confirm dose and time units; scores (`rev`, `rev_pct`, `layered_pct`, `gene_dep`) and their signs, and exact target scoring source; verify raw versus percentile are not inadvertently mixed.
- Reconstruct exact `rev = - d·a / ||d||` for a seeded representative sample, *using the existing window scorer's conventions*; verify sum of per-regulator contributions equals stored reversal. For per-gene tables verify `sum_j activity_contrib = drug_activity` and `sum_j rev_contrib = regulator rev_contrib` within documented float tolerances. Abort on failure.
- Confirm original condition and dose identifiers: a mean-of-doses or mean-of-lines is never labeled a replicate. Repeated conditions are not counted as independent experiments.
- Audit 41 disease groups, 14 subgroups, eligible 29 window-supported cohorts, and direct/borrowed/lineage-matched coverage. Counts here are **historical**; recompute from current files. Subgroups have historically lacked window data; mark `not_evaluable` until verified.
- Verify target role `D` truly uses held-out data if used as an outcome; if not, label it *in-sample mechanistic agreement*.
- Reproduce independent D209 dedup unit counts as a sanity check, without treating replication as a new discovery.
- Recompute D206 PRISM/GDSC from raw matrices with a committed producer, sign check using documented sensitive drugs, line IDs and common drug mappings; if unavailable, disable the sensitivity endpoint rather than reading its orphan table.

### 2.4 Immutable run layout

```
results/lincs/mechanistic_incremental_v1/
  00_audit/{protocol.md,coverage.tsv,join_audit.tsv,inputs_manifest.json,report.md,run.json}
  A_incremental/{...}
  B_temporal/{...}
  C_target_recovery/{...}
  D_network_sensitivity/{...}
  E_program_access/{...}
  F_synthesis/{...}
scripts/mechanistic_incremental_v1/
  common/{io.py,identifiers.py,validation.py,statistics.py,plots.py}
  00_audit.py
  A_incremental.py
  B_temporal.py
  C_target_recovery.py
  D_network_sensitivity.py
  E_program_access.py
  F_synthesis.py
  run_all.sh
config/mechanistic_incremental_v1.yaml
```

Use one deterministic, versioned `analysis_manifest.json` for the entire suite. Run production `scripts/check_params.py` in preflight; experimental values live under the new config and do **not** mutate `config/adopted_parameters.yaml` until specifically adopted through a separate decision.

---

## 3. Common mathematical definitions and statistical rules

Let `c` index disease cohorts, `x` drugs, `l` cell lines, `t∈{6h,24h}` or the project's early/late windows, `r` regulators, `g` target genes, and `k` disease-defined programs.

**Disease vector:** `d_cr` = ARACNe-derived disease-vs-normal VIPER NES on network selected by the production `net_for` rule.  
**Drug vector:** `a_cxtr` = matching network's perturbation-derived regulator response.  
**Reversal:** `R_cxt = - Σ_{r∈S_cxt} d_cr a_cxtr / sqrt(Σ_{r∈S_cxt} d_cr²)` following the exact production scorer; `S_cxt` is the actual shared finite regulator set, **not a newly chosen top-MR subset**.  
**Per-regulator contribution:** `q_cxtr = - d_cr a_cxtr / ||d_c,S||`, such that `R=Σ_r q`.  
**Temporal change:** `Δa_cxr = a_cx,late,r - a_cx,early,r`, and `Δq_cxr = q_cx,late,r - q_cx,early,r` on the **intersection of regulator/edge universes**; separately report changes in coverage if the universe differs.  
**Target role:** use the **actual D120/D121 definitions from source**, not a replacement formula. `D` is CRISPR gene-dependency percentile, `A` alteration percentile, `O` overactivity, `T` combination. Preserve which axes are independent of drug-response inputs.

**Ranking/scoring distinctions:**
- A within-cohort percentile is a rank within a changing candidate universe; paired tests should preferentially compare *raw* quantities on fixed rows and report percentiles secondarily.
- Different cohorts may use different ARACNe regulons for the same regulator. Do not compare NES values numerically across them without an explicitly defined mapping.
- Targeted hypothesis test ≠ outcome prediction. Rank-based enrichment of known targets does not show that a drug kills cells or benefits patients.

**Pre-specified inference rules:**
- Effect sizes and 95% uncertainty intervals **always**, including null outcomes. Primary confidence intervals should cluster bootstrap at the independent unit (drug/target-unit or disease as applicable), not at duplicate drug–disease rows.
- All model evaluation on held-out groups; never select features, thresholds, cohort subsets, or score weights on test folds. Drug-group CV primary; leave-one-tissue-key-out secondary. For target generalization, additionally hold out target family / independent target unit.
- Use paired tests on *identical row universes*; show selection/coverage differences separately. Two-sided p-values unless a test has a direction fixed before results. For multiple mechanistic hypotheses, report BH-FDR within declared families.
- Permutation nulls preserve drug identity/target-unit clusters, cell-line/tissue grouping, dose/time availability, annotation-count strata, and program size/degree as relevant. Never shuffle individual correlated regulon genes as if exchangeable.
- Specify effect-size scales before testing. Score gains cannot be manufactured by changing denominators or by mixing absolute and relative ranks.
- Pre-registration distinguishes **confirmatory**, **secondary**, **exploratory**, and **non-evaluable**. Exploratory discoveries can motivate later registration but do not pass confirmatory bars retroactively.

---

## 4. Phase 0 — Immutable audit and analysis-universe construction

### 4.1 Hypothesis and purpose

No biological hypothesis. Establish what exists, is independently measured, and can be joined without selection artefacts. This stage gates all later work.

### 4.2 Procedure

1. Checkout/record current default branch SHA; verify D208 `_v2_tf` production guard.
2. Read current scorer and parameter configs and produce a machine-readable DAG: source → transformation → scores → tests. Include whether a feature is derived from the target, disease, drug, cell line, or outcome.
3. Inventory each source table: shape, axis meaning, units, uniqueness key, n cohorts, n drugs, n matched lines, window coverage, target coverage, sample multiplicity, and null/missing codes.
4. Build **four distinct universes**, never forcing one universal intersection for everything:
   - `U_A_mech`: valid target-score + perturbation response + target/pathway mechanism labels.
   - `U_A_sens`: `U_A_mech` with independently reconstructible, direction-checked PRISM/GDSC outcomes.
   - `U_B`: same drug–cohort with **both observed** windows and exact contribution coverage.
   - `U_C`: compound–target–line–condition with target annotation, potency/exposure category, and a viable control panel.
5. Create `eligibility_reason` and `excluded_reason`; count exclusions by cohort, drug class, target class and tissue. Plot missingness.
6. Unit test actual path and feature resolution on a small set: BRCA-LumA fulvestrant, BRCA-Basal tamoxifen, LUAD erlotinib, SKCM vemurafenib, a top cytotoxic, and a negative/target-absent compound. **Examples are sanity checks, not statistical evaluation units**.

### 4.3 Stage acceptance

All exact decompositions pass, no silent many-to-many joins, `_v2_tf` verification passes, provenance of all chosen outcomes is reproducible, and all missing-data categories enumerated. Otherwise **stop downstream affected stage**, not the entire suite; report blockers.

---

## 5. Experiment A — Independent incremental value of perturbation information (PRIORITY 1)

### A0. Scientific design

**Primary question:** On comparable drug–disease pairs, does an observed ARACNe drug-response feature improve the ability to recover an *independent, disease-relevant functional readout* beyond target annotations and baseline disease biology?

**Important:** Multiple endpoints exist; none is a comprehensive definition of clinical efficacy. Test separate endpoints and label their roles:

| Outcome | Interpretation | Independence hazard |
|---|---|---|
| Signed PRISM/GDSC growth sensitivity, matched to drug and line | Functional drug response under a different assay | Drug/cell overlap and assay duration differ from LINCS; no shared scoring labels permitted in feature selection |
| D127-style same-drug paired disease/biomarker contrasts | Biological discrimination between relevant and irrelevant contexts | Layered score already built from target unit; only independent *increment* of perturbation can be credited |
| DepMap CRISPR target dependency | Genetic vulnerability | **Do not use as independent outcome if D120 `D` is part of baseline feature set**; perform D-held-out mechanistic alignment arm only |
| Known clinical indications/phase and trial outcomes | Weak, biased translational context | Secondary descriptive only; `no_record` ≠ negative |

### A1. Required inputs and row identity

- `U_A_sens`: independent compound–line–dose/assay sensitivity values, with line weights and matched coverage. Use a **drug × model** functional endpoint before any cohort-level aggregation where possible.
- Drug–cohort target-role components and current `layered_S` / `layered_pct`.
- Early/late observed ARACNe regulator matrices; pooled production scores recorded for reference **but not represented as independent timepoints**.
- Mechanism class, target-unit ID, effect direction, expression/mutation/alteration, line match, n conditions, dose, TAS, quality/route/exposure flags, gene/regulator coverage.

### A2. Reconstruct the functional outcome independently

For each drug `x` and screened line `l`, read assay outcome `y_xl` with raw units (e.g. log-fold change or AUC of viability vs dose). Determine **a priori from source documentation** whether larger or smaller values mean more killing. Define `s_xl = sign_source × standardized(y_xl)` with **higher = greater sensitivity** and write sign provenance. Standardization must be defined in training folds if used for prediction; for descriptive within-drug selectivity, z-score across all qualified lines for that drug only if sufficiently many lines exist. Keep raw outcome beside z.

For cohort `c`, let `w_cl ≥ 0` be *existing* Celligner matched, in-key weights over observed screened lines `L_cx`. Then

`Y_cx = Σ_(l∈L_cx) w_cl s_xl / Σ_(l∈L_cx) w_cl`.

Report `n_lines`, `weight_coverage = Σ_obs w / Σ_eligible w`, standard error or line-bootstrap interval, and whether every used line represents the cohort's lineage. Do **not** silently substitute unqualified borrowed lines or treat replicate assays as independent lines. Construct a separate direct-only and coverage-sensitive analysis. Never give `Y` a therapeutic interpretation beyond assay sensitivity.

**Positive controls:** gene/pathway-targeted drugs in known dependency-matched lines and expected sensitive/resistant comparisons where observed assay data exist.  
**Negative controls:** target-absent lines, unrelated target-unit/drug pairings, gene-expression-shuffled/mismatched disease contexts, and a deliberately generic drug-class baseline. Controls must be run before attributing a positive result to regulatory biology.

### A3. Freeze baseline and new perturbation feature families

Use a small **predeclared nested feature design** to prevent a kitchen-sink boost:

- **B0 annotation-only:** drug identity class, main/secondary target list, potency support, number of targets, assay/annotation coverage, FDA/known-use indicator for confounding analysis only; exclude indication labels when predicting indications.
- **B1 target + disease biology (PRIMARY BASELINE):** D120 independent target components `D`, `A` and eligible `O` if its origin is genuinely outside perturbation response, target expression/presence, target-unit `layered_S`, molecular subtype, driver features, matched-line relevance, and B0. If `layered_S` already encodes a feature, avoid adding its copies as separate predictors unless explicitly modeling its components; specify one canonical B1 set before test.
- **P1 global response (minimum perturbation increment):** observed `rev`, `rev_mech` where valid, early and late separately, raw response magnitude, sign consistency, `prolif_share`, number of conditions and line coverage. **P1 coverage-only placebo** includes quality/availability counts but excludes actual response values.
- **P2 mechanism-local response:** activity and exact `rev_contrib` only for target-connected regulators/pre-specified target/pathway footprints, with target↔regulator link sources recorded. Do not declare ARACNe associations causal. Include signed per-target gene opposition *only where full edge coverage is available*.
- **P3 temporal response:** early and late raw scores as separate features when both observed, `ΔR`, sign flip, persistence, early/late target-footprint change. This is a distinct subset and must be compared with B1/P1 restricted to **identical both-window rows**.

Avoid unbounded feature selection. Feature transformations (scaling, regularization, missing-value indicators) are learned only on training folds. Handle multi-target drugs by retaining mechanism units; no automatic max-dependency aggregation without a mean/top-target sensitivity arm.

### A4. Three complementary tests (do not conflate)

**A4.1 Descriptive conditional stratification.** Stratify disease–drug pairs into quintiles of B1 evidence *within cohort*, and within each bin compare high/low perturbation measures. Report outcomes, effect sizes, coverage and class composition. Within-bin Spearman or partial residual correlation: residualize `Y` and P1 on B1 and data-quality covariates in training subsets, then correlate residuals in held-out rows. This is descriptive and cannot establish causality.

**A4.2 Out-of-sample incremental model test (PRIMARY).** Fit a **simple regularized baseline** (e.g. ridge or elastic-net selected in nested CV) for `Y` with B1; fit identical learner with B1+P1; optionally B1+P2 and B1+P3 only after their predeclared coverage gates. Match rows across models, fit all preprocessing inside train, and compare out-of-fold predictions:

`Δmetric = metric(B1 + P, Y_test) − metric(B1, Y_test)`.

Primary metric: **within-cohort Spearman** where each cohort has at least 20 independent compound units and variation in outcome; summarize median across eligible cohorts. Also show pooled out-of-fold Spearman with cohort-stratified uncertainty, MAE or RMSE as appropriate, and calibration if predicting numeric sensitivity. For categorical biomarker contrasts use paired correct-order share and AUC with target-unit dedup; report coverage versus same-drug D127 reference. Keep per-cohort cells and do not hide sign inversions by pooling. Baseline versus augmented models on identical rows, with 95% paired, cluster-bootstrap CI for delta.

**A4.3 Added-information test within matched target units.** Restrict to `(cohort,target_unit)` groups with >=2 genuinely distinct drugs/conditions and varying perturbation response. Since layered target score can be identical within unit, test whether P1/P2 predicts independent observed sensitivity differences between compounds **after accounting for shared target evidence**. Use fixed effects for `(cohort,target_unit)` or equivalent within-unit contrasts. If no independent variations remain, report `not_evaluable`; do not treat repeated identical scores as separate observations.

### A5. Splits, confounders and controls

- **Primary split:** group by `ingredient_id` / compound family so salt/BRD versions of the same drug never leak across folds. 5 outer folds, nested tuning on 4 inner folds, deterministic seed frozen before outcomes.
- **Generalization:** leave-one-tissue-key-out and, where feasible, leave-one-target-family-out. Report poor support rather than pretending held-out groups are covered.
- **Confounders:** assay line and tissue, drug class, drug dose, n conditions, quality/TAS, timepoint, target annotation count, target evidence strength, LINCS matched-line status, proliferation magnitude, broad cytotoxicity, and whether the target unit was measured in source disease.
- **Placebo 1:** replace true perturbation values with shuffled values among drugs matched on tissue, window, target-count, assay coverage and broad mechanism class, without changing the metadata or missingness.
- **Placebo 2:** add only perturbation quality/availability variables to B1. If B1+quality accounts for the observed gain, the gain is not attributable to regulatory measurements.
- **Placebo 3:** drug-class-only baseline, and within-target-unit evaluation, to distinguish class familiarity from perturbation specificity.
- **Sensitivity:** matched/borrowed separated; single-line excluded; high-TAS; high-reliability; cytotoxics excluded; both windows subset; target annotation sources; direct phenotype-vs-target-dependency outcomes; pharmacologic plausibility subset.

### A6. Registration thresholds (proposal; freeze before running)

**Coverage gate:** primary test only if >=5 cohorts each have >=20 independent compounds with valid independent `Y`, with >=2 distinct target units among drugs, and positive controls have the expected sensitivity direction. If not, report coverage/feasibility, no imputation/re-ranking and no false positive/negative scientific conclusion.

**Incremental-value PASS (joint):** (i) median cohort-wise held-out Spearman delta ≥ +0.05; (ii) paired 95% bootstrap CI lower bound > 0; (iii) >=60% of eligible cohorts show positive delta; (iv) correct-direction signal survives at least one stringent drug-group and tissue/target-family holdout with adequate support; (v) not explained by quality-only placebo, drug-class baseline, or target-unit duplication. The +0.05 threshold is an **initial proposed practical improvement**, not an empirical discovery. It must be signed off/frozen before data inspection. If effect CI excludes 0 but is <0.05, report a small effect, not PASS.

**FAIL:** enough data and power, but criteria fail. **INCONCLUSIVE:** insufficient outcome coverage, high uncertainty, direction control fails, or test folds too sparse. **HETEROGENEOUS:** preregistered mechanism/disease strata show effects with FDR-controlled interactions despite overall FAIL; label as stratum-limited and replicate independently.

**A7. Required outputs:** `A_observations.tsv.gz`, `A_outcomes_signed.tsv.gz`, `A_feature_audit.tsv`, `A_cv_folds.tsv`, `A_models.tsv`, `A_paired_deltas.tsv`, `A_within_target_units.tsv`, `A_controls.tsv`, `A_per_cohort.tsv`, `A_error_taxonomy.tsv`, `A_shap_or_coefficients.tsv` (only interpretable held-out coefficients; no SHAP if inappropriate), coverage and residual plots, worked examples, full `report.md`.

**Interpretation required:** Distinguish “perturbation data predicts selective killing” from “perturbation data is mechanistically consistent with target biology.” Neither is equivalent to patient benefit.

---

## 6. Experiment B — Temporal, regulator and gene-level explanation (PRIORITY 2)

### B0. Hypothesis

The 6h↔24h changes are structured: persistent, delayed, reversing, and transient responses are associated with distinct gene-level contributions and pathway/target annotations. The experiment **explains observed transitions** and does not impute unobserved windows.

### B1. Eligibility and fixed comparisons

Use `U_B` with measured early and late windows. Analyzed units: `(cohort,pert_id,regulator)` and separately `(pert_id,ModelID,dose,regulator)` for conditions with matched dose/time where available. Do not assume the cohort-weighted early and late profiles used identical lines, doses or condition numbers: report `lines_common`, `doses_common`, n conditions, time-difference and weight coverage. Produce three nested arms:

1. **Primary controlled arm:** same drug, same line, same dose (or the nearest **exactly documented matched dose**, no hidden interpolation), both windows; quality filters applied consistently.
2. **Secondary production-like arm:** original window score based on all eligible matched conditions, explicitly reporting dose and line composition differences.
3. **Higher-confidence arm:** both windows with >=2 independent conditions and passing split-half reproducibility checks; flag loss of support.

### B2. Core calculations

Per observed `(x,l,t,r)`, compute the exact ARACNe activity `a_xltr`. Within identical network/edge universe:

`Δa_xlr = a_xl,late,r − a_xl,early,r`.

For cohort-level score:

`ΔR_cx = R_cx,late − R_cx,early`.

For exact regulator terms:

`Δq_cxr = q_cx,late,r − q_cx,early,r`, `Σ_r Δq_cxr = ΔR_cx` **only when common regulator norm/universe matches**; otherwise report the explicit coverage-composition remainder after recomputing both on the intersection.

For each fixed regulator edge `(r,g)`, compute:

`Δc_cxrg = c_cx,late,r,g − c_cx,early,r,g`.

Verify `Σ_g Δc = Δa_r` under the relevant full-target decomposition and `Σ_r Σ_g Δrev_contrib = ΔR` in matched-universe rows. Use full edges via `79e3` query/full-target mode, not top-25 truncated tables, for decision-bearing tests.

**Reproducible sign rule:** estimate early/late standard error or split-half interval using independent conditions; a response is `positive`/`negative` if its two-sided 95% uncertainty interval excludes zero, otherwise `uncertain`. Sign reversal only if confidently positive at one window and confidently negative at the other; other categories: `persistent_same_sign`, `early_only`, `late_only`, `both_uncertain`, `mixed/ambiguous`. Avoid fixed NES >0 without uncertainty, which would mistake small noisy effects for biology.

**Persistence measures:**

- `signed_delta = a_late − a_early`
- `abs_delta = |a_late − a_early|`
- `directional_similarity = cosine(a_early, a_late)` within common regulator space, if norms nonzero.
- `dominant_gene_concentration = (sum of largest 5 |Δc|) / sum_g |Δc|` and `n_eff_gene = (Σ|Δc|)^2 / Σ Δc²` (effective contributors, explicitly a participation-style measure, not number of independent genes).
- `gene_set_overlap = |top_edges_early ∩ top_edges_late| / |top_edges_early ∪ top_edges_late|`, with top edges chosen on predeclared fixed count 10 or all if fewer; report sensitivity to 5/20 **descriptively**, not tune by outcomes.
- `same_genes_flip`: fraction of high-contribution shared genes whose signed `c_early * c_late < 0`; `replacement`: fraction of late absolute mass on genes not in the top early set.

### B3. Attribution of early/late differences

For each high-confidence transition, partition the **absolute** change `Σ_g |Δc_g|` into annotated gene sets: Hallmark proliferation, stress/DNA damage, apoptosis, lineage/identity, and pre-specified known target pathways (make gene-set sources explicit and count overlaps; report overlapping and disjoint partition separately). Attribution is descriptive; do not claim unique cause when one gene belongs to several sets.

Assess if reversals are driven by (i) same genes flipping, (ii) different genes replacing one another, (iii) a few leverage genes, (iv) line/dose-composition differences. Use shared-model, matched-dose data to distinguish genuine time changes from compositional confounding.

Use D154 proliferation β by time only as a **diagnostic**. Compare observed β changes against gene-contribution changes and generic/universal/line-selective proliferation components from D157. Do not fit a pooled β and subtract it to manufacture a new score.

### B4. Nulls, controls and tests

- **Positive:** known early pathway-responsive drug mechanisms (on-target PROGENy D99 direction where independently supported); persistent pathway controls with replicated signatures.
- **Negative:** time labels permuted **within drug × line × dose replicate groups** where exchangeability is justified; matched pairs of unrelated drugs in same line/time with similar response magnitude; broad cytotoxic versus mechanism-specific controls.
- **Pseudo-replicate:** split an observed single window's conditions into halves, run the same transition classifier; estimate false sign-flip rate at the chosen confidence rule.
- **Stability:** 2+ / 3+ / 5+ conditions, high vs newly admitted TAS signatures, only same-dose matched pairs, excluding single-line responses.
- **Null for mechanism enrichment:** within `(cohort,target class,network-size,TAS,condition-count)` strata, permute mechanism labels by independent drug/target-unit clusters. Compare transition class × pathway enrichment; BH-FDR within each declared family.

**Proposed confirmatory criterion:** at least one predeclared known pathway/target class shows `transition class` enrichment FDR <0.05, effect odds ratio ≥1.5, with no comparable enrichment in the pseudo-replicate null, and same direction in the high-confidence matched-dose arm. Otherwise descriptive/negative; not a reason to re-rank. This threshold must be frozen before analyzing outcomes.

### B5. Specific worked cases and output

Required illustrative case files (only if coverage supports): fulvestrant/clomifene/tamoxifen in LumA versus Basal, erlotinib/LUAD, a MAPK inhibitor/SKCM, and one matched cytotoxic. For each: raw 6/24 profiles; top 20 regulators by |Δq|; top 25 edges by |Δc|; gene-level signs; target engagement evidence; dose/line comparability; data quality; alternative explanation.

Outputs: `B_matched_conditions.tsv`, `B_transitions.tsv.gz`, `B_gene_decomposition.tsv.gz`, `B_pathway_attribution.tsv`, `B_pseudoreplicate_null.tsv`, `B_enrichment.tsv`, `B_reliability.tsv`, `B_case_studies/`, paired spaghetti plots, time-window heatmaps, `report.md`.

**Interpretation:** even a reproducible 6→24 sign flip does not by itself establish feedback signaling or resistance; those are hypotheses requiring an independent perturbation or phosphoproteomic measurement.

---

## 7. Experiment C — Primary/secondary target recoverability and pharmacologic plausibility (PRIORITY 3)

### C0. Hypothesis

D158's secondary-target null could be partly due to mixing weakly bound, absent, or unengaged annotations and testing **exact target rank** when effects are downstream. This experiment tests both possibilities with independent controls.

### C1. Cohort/target universe and definition of truth

Unit `(drug,target,cell_line,window,dose)` when condition-level coverage exists; otherwise `(drug,target,cohort,window)` with the source mapped and flagged. Classify target role using an auditable hierarchy: (a) curated primary mechanism; (b) independently binding-confirmed secondary target; (c) annotated/unverified target; (d) no annotation. If primary/secondary designation is ambiguous, label unknown and exclude from primary-vs-secondary statistical claims.

Use D134 feasibility audit to find target units with >=2 genuinely distinct compounds sharing tested lines; do not count duplicate targets, multiple salts, or multi-BRD records as independent positives. Report coverage by role, target family, line, dose, potency source, and non-target controls.

### C2. Exposure and engagement plausibility (without abusing clinical PK)

Normalize potency measurements to **molar units** but retain assay endpoint type. `Ki`, `Kd`, `IC50`, and `EC50` are **not interchangeable**. Report separately for each assay type, and use a potency category based on a documented, fixed hierarchy of *cellular assay matching* rather than pooling all into a single physical binding constant.

For genuine equilibrium `Kd` (or `Ki` under justified binding assumptions) and **measured unbound** in-vitro concentration `C_free`, optional theoretical occupancy is:

`occupancy = C_free / (C_free + Kd)`.

When only nominal LINCS `C_nominal` is available, compute `log10(C_nominal / potency)` **as an exposure-to-potency proxy only**, clearly labeled not occupancy or target engagement. Account for unknown free fraction, uptake, metabolism, target residence time, and assay format. Plasma Cmax is a separate *clinical plausibility flag* and cannot be substituted for cell-culture free concentration.

Proposed predeclared descriptive bins for **matched-assay potency proxy**: `ratio <0.1`, `[0.1,1)`, `[1,10)`, `≥10`, plus `unknown`; do **not** call them target occupancy bins or select a best threshold from recovery. Always report the continuous log-ratio test as primary; bins are visual explanations only.

### C3. Two distinct recovery tasks

**C3.1 Exact-target recovery (reproduction and stratification of D158).** For each drug–line–window regulator list, calculate target percentile or rank in the expected direction and `hit@k`, with primary `top 1%` and secondary top 5%/10 and absolute top-10. The eligible universe is the regulators **actually measured in that same network**, so `chance@k = k/N` per observation. If annotated target not in the regulator universe, record `unscorable_target`; do not label it a recovered negative. Run primary and secondary separately, stratified on target expression, potency proxy, number of conditions, network regulator count, cell line matchedness, and time. D158 is a historical comparator, not assumed replicated.

**C3.2 Downstream mechanism recovery (new).** Build a **fixed, independent** target→pathway/regulator map from literature-curated directed/signed interactions where available (e.g. pinned OmniPath/SIGNOR), *not* from the same LINCS transcription profiles used to test it. For each target `g`, predefine a target-associated set `H_g` by 1-hop interactions, and a separate 2-hop sensitivity arm. Exclude unsupported directions and preserve ambiguous edges. Compute a set enrichment of response rank for `H_g` using a fixed weighted mean or rank enrichment statistic, with degree/size-matched permutations; avoid claiming direction of drug effect unless pharmacologic agonism/antagonism and interaction signs are documented. Compare target's known pathway to matched non-target sets with similar size, network degree, tissue expression and annotation density. **D137 negative controls must be repeated on any putative tracing rule before claiming mechanism**.

For a target pathway set with expected signed direction `s_gr` and observed standardized response `z_xtr`, one optional statistic (only for supported signed edges) is:

`M_xgt = Σ_(r∈H_g) w_gr s_gr z_xtr / sqrt(Σ_(r∈H_g) w_gr²)`.

Freeze weights from **external evidence** and standardize using matched nulls. If edge signs are incomplete, use unsigned enrichment and label interpretation accordingly. An ARACNe coexpression edge alone does not specify causal sign or hop.

### C4. Tests, controls, inference

- Recover known primary target mechanisms as a **positive yardstick**; if this fails, do not interpret secondary nulls as informative.
- Exact secondary target enrichment versus degree/expression/annotation-matched non-targets (not versus all random genes).
- Within secondary targets, test monotonic association between continuous `log10(C_nominal/potency)` and standardized retrieval rank, while adjusting for target expression, class, condition count and cell line, with grouped uncertainty by drug/target.
- Compare exact-target recovery versus downstream target-associated pathway recovery on the **same** pairs, with matched-control nulls and paired CI.
- Negative controls: target absent in line, well-documented non-target, wrong drug–target label with family/degree matched, irrelevant downstream pathway, target-free lines, and D137-style cohort/drug/line randomized controls. Controls must not close at the same rate as real mechanisms.
- Stratify kinase, nuclear receptor, GPCR, epigenetic enzyme, DNA/RNA effectors and other protein classes; explicitly mark where no OmniPath route exists (D137 TOP1/TOP2A etc.), rather than treating missing graph edges as a negative biological result.

**Proposed decision:** a secondary-target signal is provisionally supported only if (1) primary yardstick works; (2) secondary enrichment has FDR<0.05 with ≥1.5-fold lift over matched controls; (3) improvement persists after target-unit deduplication; (4) true-target recovery exceeds negative-control recovery with a held-out drug and target-family evaluation. If only downstream succeeds, call it **pathway-consistent effect**, not verified engagement. If coverage of strong binding secondary targets is <30 independent target units, call the stratum feasibility-limited, not a decisive null.

Outputs: `C_target_catalog.tsv`, `C_coverage.tsv`, `C_exact_recovery.tsv`, `C_downstream_recovery.tsv`, `C_potency_strata.tsv`, `C_controls.tsv`, `C_by_target_class.tsv`, `C_case_studies/`, `report.md`.

---

## 8. Experiment D — Why ARACNe regulator interpretations differ between contexts (FOLLOW-ON)

### D0. Hypothesis

Different regulon targets/signs and different baseline line environments, rather than only random measurement error, account for much of the same-drug/same-regulator variation across disease network contexts (D137/D197).

### D1. Controlled network-swap design

**Primary natural experiment:** Same drug–line conditions and same perturbation gene vector projected through BRCA LumA and BRCA Basal networks; evaluate the shared regulators. Separate a `same-network/similar-lines` reference (COAD/READ where feasible) and a `different-network/different-lines` descriptive reference (e.g. LumA/LUAD). Never compare mismatched samples without marking both network and line changes.

For each shared regulator `r`, calculate:

- `J_r = |E_A(r)∩E_B(r)| / |E_A(r)∪E_B(r)|` (target-set Jaccard).
- `sign_agreement_r = fraction_{g∈intersection}[sign(mode_A,rg)=sign(mode_B,rg)]`.
- `Δactivity_rx = a_rx^(A)−a_rx^(B)` on **identical input gene vector**.
- `Δshared`: activity difference reconstructed on common targets with documented per-network rank/normalization conventions; `Δunique`: difference attributable to unique targets. Because VIPER re-ranks and rescales target sets, these terms need not be an algebraic partition unless exact contribution code proves it. If not exact, report a counterfactual common-target re-score with the residual separately, and label it descriptive.
- `overlap_weighted` and gene-level contribution diagnostics; regulator-size and expression controls.

Use matched-null regulator pairs by size, original response magnitude, and expression. Bootstrap drug conditions by `pert_id` and line. Test whether target overlap/sign agreement explains held-out absolute activity disagreement, not only Spearman of full matrices. Do not choose the network that best reproduces a preferred disease label without independent evidence.

**Decision:** mechanistic interpretation supported only for regulators with reliable network-swap data, directional perturbation controls, and stable edge-level explanation on held-out drugs. Otherwise describe network dependence as unresolved.

Outputs: `D_network_pairs.tsv`, `D_edge_switches.tsv`, `D_activity_decomp.tsv`, `D_controls.tsv`, `D_worked_ESR1.md`, `report.md`.

---

## 9. Experiment E — Disease-program pharmacological accessibility (FOLLOW-ON)

### E0. Hypothesis

Although drug programs and disease programs poorly overlap (D200), drugs may reproducibly modulate *parts* of disease-derived programs; those accessible parts may be enriched for disease-relevant target vulnerabilities independently of global program reversal.

### E1. Define fixed disease programs and projection

Use D198 program memberships defined from untreated patient-level ARACNe activity, per cohort only; no cross-tissue program transfer. Use all program members, not only the ten summarized members. Keep D199 flags for plate, RNA quality, proliferation and generic response; **flag**, do not automatically delete all non-clean programs. Calculate coverage:

`coverage_kx = |P_ck ∩ S_cx| / |P_ck|`.

For `d_cr`, `a_cxtr`, define a descriptive module opposition:

`O_cxkt = -Σ_(r∈P_ck ∩ S_cxt) d_cr a_cxtr / sqrt(Σ_(r∈P_ck ∩ S_cxt) d_cr²)`.

Record `n_used`, `coverage`, `program identity`, `program flag`, and zero-norm `not_evaluable`. This is an aggregation of the same reversal contributions, not an independent discovery. Do not interpret high `O` as target engagement or invent a score improvement.

**New evidence:** For each program, ask which observed drugs repeatedly perturb its regulators across independent conditions/cell lines, and whether its independently known vulnerabilities and causal target pathways are disproportionately affected. Compute within-drug repeatability of program response (condition split-half ICC/correlation), cell-line selectivity (variance across lines after subtracting drug mean, with enough lines), and pathway target enrichment versus matched degree/size null. For mechanism-local accessibility, calculate response of regulators in a predeclared independently supported target-connected subset, not a post-hoc best-performing subset.

### E2. Controls

- Program-member shuffle within cohort preserving module size, regulator target-count/degree, baseline |NES|, and measurement reliability.
- Compare plate-flagged programs to clean programs: plate programs must not appear as stronger mechanistic evidence merely because they are stable.
- Drug-class- and amplitude-matched generic responses to test whether accessibility is nothing more than large response magnitude.
- Require repeated conditions/lines; single-line drugs cannot support cell-line-selective conclusions.

**Decision:** label a drug–program interaction as *reproducibly perturbable* only with sufficient repeat conditions and an uncertainty interval excluding zero; label as *mechanism-linked* only with externally supported target/pathway evidence and a passing matched-null test. The program arm is explanatory, not a new production ranker. No retrospective D201 AUC optimization.

Outputs: `E_program_coverage.tsv`, `E_program_response.tsv.gz`, `E_repeatability.tsv`, `E_target_enrichment.tsv`, `E_controls.tsv`, `E_accessibility_matrix.tsv`, `report.md`.

---

## 10. Experiment F — Integrated evidence and a mechanism report (NO AUTOMATIC RERANKING)

### F0. Objective

Create a reproducible explanatory report per drug–disease pair showing which claims are supported by which independent observations and what remains unknown. This is **not** one learned scalar or an efficacy claim.

### F1. Evidence layers and strict language

For each pair assemble:

| Layer | Evidence field | Allowed statement |
|---|---|---|
| Disease | Differential activity, mutation, target role, dependency | “Target/pathway is dysregulated or genetically vulnerable in these models” |
| Binding | Curated target, potency assay, expression | “Drug can bind/modulate the target in the stated assay; engagement in LINCS is uncertain/corroborated” |
| Perturbation | Direct observed drug gene/regulator response, exact signed contributions | “Observed in the specified drug/line/time/dose condition” |
| Temporal | Observed early/late transition and reliability | “Response persists/reverses/appears late in measured windows” |
| Network | ARACNe associations, pathway graph and edge support | “Network-consistent association”; no unjustified causal route |
| Functional | Independent CRISPR/PRISM/GDSC endpoint | “Genetic vulnerability” or “assay sensitivity”, not general efficacy |
| Translation | Clinically plausible exposure, route, trial data | “Supports/limits feasibility”; no-record is unknown |

Mark `confirmed`, `supported`, `conflicting`, `unknown`, `not_tested`, `failed_control` at each layer with explicit evidence/source IDs. Avoid calling target dependency, regulator reversal and approval three independent validations if shared annotations/measurements connect them.

### F2. Candidate selection for detailed explanations

Select **without outcome peeking**: (i) known positive controls, (ii) high target-role/low reversal, (iii) low target-role/high reversal, (iv) target-compatible with differential early/late response, (v) plausible secondary-target case, (vi) target-absent negative, (vii) one broad cytotoxic, and (viii) top novel (no recorded indication) with adequate matched-line evidence. Select one per category per represented cohort if possible; missing category should be explicitly empty.

Deliver:
- `F_mechanism_evidence.tsv.gz`: one row per `(cohort,drug,window,target_unit)` with raw evidence, direction, reliability, provenance and uncertainty; no fabricated timepoints.
- `F_mechanism_cases.md`: worked mechanisms with network/regulator/edge contribution chain, dose/time/model coverage, positive/negative controls, alternative explanations, and falsifiable follow-up experiment.
- `F_evidence_summary.tsv`: cohort-level counts of supported vs unknown vs contradictory claims, by evidence layer.
- `F_decision_matrix.md`: for A/B/C/D/E pass/fail/inconclusive and whether an additional controlled experiment or external data is required.
- `F_final_report.md`: executive answer to the primary question, numeric results, controls, limitations, reproducibility, and justified subsequent project decisions.

---

## 11. Statistical power, missingness and decision discipline

### 11.1 Power and uncertainty

Do not assume 4,128 both-window rows are 4,128 independent experiments: the historical set involves 939 drugs across 14 cohorts/8 tissue keys; line/model overlap and target duplication lower independence. For each endpoint, report independent drug count, target-unit count, tissue count, and source assay count. Use cluster-robust bootstrap, hierarchical analysis or group permutation consistent with the question.

Before model fitting, compute a **minimum detectable increment** from group-bootstrap variability on training folds only. If the proposed +0.05 correlation improvement is smaller than detectability, label the primary test *underpowered* rather than “negative.” Report attained confidence interval and observed power limitations without post-hoc power claims.

### 11.2 Multiple comparisons

Primary claims are A4.2 and B4/C4 if their data gates pass. FDR-control secondary per-regulator/per-pathway screens. If hundreds of programs/targets are examined, treat any standout examples as hypothesis generation unless they survive independent holdout/replication.

### 11.3 Interpretation of missing data

Four separate states: (1) not LINCS-tested in matched line; (2) measured but weak/inactive; (3) tested but timepoint missing/unknown; (4) annotation missing. Never collapse these to zero response, no target, or non-efficacy. Return the evidence status explicitly.

### 11.4 Audit of claimed independent outcomes

For each experiment, generate an `evidence_dependence_graph.tsv` with nodes = measured datasets / derived scores / labels and edges = “derived from” or “shares raw input”. If a predictor and outcome share input upstream, mark their evaluation as **consistency**, not independent validation. D120 target role and D127 layered target components, for example, are not independent of each other.

---

## 12. Engineering instructions for a coding agent

### 12.1 Workflow order

1. Record repo SHA, `git status`, current branch, adopted parameter registry, and all source file hashes; do not run expensive builds until audit.
2. Create `config/mechanistic_incremental_v1.yaml` and this spec in repo on a **new feature branch**. Freeze decisions and test IDs before computing outcomes. Commit spec+config first; log that commit SHA as protocol ID.
3. Build `00_audit.py` and validate the canonical tables with unit tests. In `report.md`, list skipped sources and why.
4. Run **A** using signed and independently reconstructed functional outcomes; if functional data are not reproducible or coverage/power gates fail, still generate A mechanistic controls and clearly separate non-evaluable primary outcome.
5. Run **B** strictly on observed early/late conditions. Do not run D183 imputation or write filled-in missing windows.
6. Run **C** with pre-exposure target annotation audit. Only evaluate secondary-target lift in the subset with adequate binding/line/dose evidence; record non-evaluable classes.
7. Run **D** and **E** only after A–C reports are written; these are explanatory, not opportunities to tune A/C criteria.
8. Run **F**, write a consolidated methods-and-results document, and generate a provenance index linking source files, run logs, figures, and decisions.
9. Open a PR for review. **Do not merge or modify production scoring.** Include a concise diff summary and commands to reproduce. No secrets or raw private identifiers in commits.

### 12.2 One-command runner skeleton

```bash
# Suggested interface, to be implemented; these commands do not exist yet.
bash scripts/mechanistic_incremental_v1/run_all.sh --preflight
bash scripts/mechanistic_incremental_v1/run_all.sh --stage audit
bash scripts/mechanistic_incremental_v1/run_all.sh --stage A
bash scripts/mechanistic_incremental_v1/run_all.sh --stage B
bash scripts/mechanistic_incremental_v1/run_all.sh --stage C
bash scripts/mechanistic_incremental_v1/run_all.sh --stage D
bash scripts/mechanistic_incremental_v1/run_all.sh --stage E
bash scripts/mechanistic_incremental_v1/run_all.sh --stage F
```

Every executable command should be wrapped in the repo's logging system (`scripts/runlog.sh` if available) and must verify the expected outputs/checksums before writing a stage success marker. A crashed partial run must be distinguishable from a successful one. Read-only inputs; write only namespaced `results/lincs/mechanistic_incremental_v1/` and its logs.

### 12.3 Required unit and integration tests

- Same drug different BRD IDs → consistent ingredient mapping; two different compounds sharing a target unit are not collapsed prematurely.
- No input row duplicates on `cohort×pert_id×window` and `drug×target×line×condition` keys.
- Time/dose units correct; unrecorded sentinel values (`-666`) always `unknown`, not early.
- Null/missing outcome labels never coded as negative efficacy.
- `rev` rebuild agrees with stored score at production rounding tolerance; full edge/regulator sums agree exactly within floating-point tolerance.
- Gene contribution normalization and sign preserved for early/late; same network and matched common feature universe.
- D209 independent unit count reproduced on frozen D127 data; no pseudo-replication in sign tests.
- Test fold leakage: no ingredient/pert_id or target-family collisions where held out; transformation fitted only inside train.
- Model comparison on exactly same rows; report how many rows remain after intersection.
- Raw PRISM/GDSC sign check, efficacy units, drug × cell-line mapping, and known positive controls.
- Shuffling response vectors destroys target-specific gain but preserves quality-only features.
- Repeated run with same seed produces same summary and table hashes (except timestamps/provenance fields).
- Confirm no modifications to production score or archived consensus files.

### 12.4 Resource-aware implementation

The full gene-level tables and condition matrices are large. Use memory-mapped `.f32` with source-provided shape, streaming Parquet/Arrow or chunked TSV, lazy/full extraction by selected cases, bounded thread counts, and staged materializations with checksums. Avoid loading >10 million gene-edge rows into a single wide pandas object if a streaming aggregation suffices. D179 and D193–D195 already provide exact decompositions; reuse them rather than repeating expensive VIPER unless stale/absent.

### 12.5 Reporting style for each test (mandatory)

Each `report.md` must answer, in plain language:

1. What did we hope to determine, and what would have falsified it?
2. Exactly which array was analyzed, what does each cell contain, its units, dimensions, and key?
3. What preprocessing/selection changed the population or scale?
4. What are the independent observations and the positive/negative controls?
5. What is the main numerical result, with confidence interval and denominator?
6. Does the result establish measurement reliability, mechanistic consistency, independent functional prediction, or none?
7. Which interpretation is ruled out, and which alternatives remain?
8. What additional raw data or intervention would decisively answer the remaining question?
9. Were any decision rules amended? If yes, which commit predated the new result?

Every figure needs a caption containing `n` independent units, axes/units, tested window, whether score is raw/percentile, and the relevant source file.

---

## 13. Decision tree after execution

| Outcome pattern | Proper next action |
|---|---|
| A passes broadly and negative controls pass | Register a **separate** study to integrate independent perturbation evidence into candidate ranking; hold out a new functional benchmark |
| A passes only for particular target classes/windows | Build mechanism-specific applicability/evidence flags; keep abstention for other classes |
| A fails with adequate independent functional outcomes | Do not add reversal to the lead score; prioritize target/perturbation interpretation or new measurements rather than score tuning |
| A inconclusive due to sensitivity coverage | Obtain matched functional drug-response assays or new data; do not call the perturbation arm biologically useless |
| B yields reproducible target-associated transitions | Develop a temporal mechanistic model or prospective targeted time-course experiment; do **not** revive missing-window imputation without new evidence |
| B mostly driven by dose/line composition | Improve matched-condition design and report temporal evidence only when comparable |
| C secondary recovery improves with exposure/plausibility | Register an engagement-aware secondary-target interpretation test on independent target/drug families |
| C exact-target fails but downstream pathway succeeds | Report pathway-consistent effects with proper uncertainty; do not call the secondary target experimentally engaged |
| C fails even in strong, expressed, dose-plausible cases | Stop exact secondary target recovery through disease-opposition rank; evaluate target-focused assays or different perturbation data |
| D shows a reproducible network-definition artifact | Add regulator-level interpretation flags and network sensitivity to reporting; don't choose network based on label recovery |
| E maps reproducibly accessible disease programs | Use them to design mechanistic follow-up experiments; don't automatically substitute program reversal for the lead score |

---

## 14. Final agent handoff prompt (copy and paste)

**Replaced by Amendment 1, A1.9 (2026-11-23). Do not use the prompt below: it asks for the PRISM/GDSC rebuild that A1.4 removes.**

> Implement the repository research specification in `SPEC_Perturbation_Mechanistic_Incremental_Value.md` for `kirtihargunani-cmu/repurposing`. Start by reading the current default branch's `docs/DECISIONS.md`, `docs/PROGRESS.md`, `config/adopted_parameters.yaml`, source-producing scripts, and D208/D209 production updates. Record the repo SHA. Treat this spec as a proposed registration: commit the protocol and locked config on a new feature branch **before** looking at any new outcome results. Do not modify production scores, thresholds, current rankings, or archived regulon data. Use current ARACNe-only `_v2_tf`, re-use exact contributions and existing D183/D158/D209 outputs with verified provenance, and respect unknown/unmeasured windows. Execute Stage 0 and experiments A, B, C, then D and E, followed by synthesis F, in that order. Do not skip stages that fail: write negative-control results, coverage, uncertainty, and plain-language interpretation. Rebuild D206 sensitivity benchmarking from identifiable raw inputs or declare independent efficacy evaluation non-evaluable. Match rows across baselines and augmented models; group split by drug and target unit; leave out tissue/target families when support exists; include placebo quality-only and drug-class baselines. For each test write all required TSV/Parquet results, plots, `report.md`, `run.json`, hashes and a consolidated final report. Do not adjust acceptance thresholds after results. On completion, open a PR with exact reproduction commands and an explanation of what was established, what failed, what remains unknown, and what future experiment is warranted. **No auto-merge.**

---

## 15. Provenance and caveat notes

This specification is based on the repository's documented methods and results inspected 2026-11-23, especially `docs/DECISIONS.md`, `docs/PROGRESS.md`, `docs/SPEC_window_imputation.md`, `docs/SPEC_covariance_correction.md`, `docs/SPEC_d127_target_dedup.md`, `docs/SPEC_target_recoverability.md`, `docs/SPEC_dose_response_bmd.md`, `docs/NOTE_proliferation_decomposition.md`, `scripts/79ac_window_imputation.py`, `scripts/79e3_window_contributions.py`, `scripts/79n_dependency.py`, `scripts/79k_window_long.py`, `scripts/79ai_programs.py`, `scripts/79ao_two_stage.py`, `scripts/79ap_covariance_correction.py`, and `scripts/d127/d127_biomarker.py`.

**Do not confuse planning with evidence:** no analyses specified in this file have been executed as part of its creation. Some paths in the historical documents may now be stale or removed after D208; Stage 0 must discover and verify current artifacts rather than trusting this text as an absolute file manifest. Several proposed decision thresholds are new design choices and must be accepted/frozen before the implementing agent examines outcomes. Positive findings from dependencies/target scores sharing inputs are consistency checks, not independent predictive validation. Inference from LINCS cell lines to patients is not established.