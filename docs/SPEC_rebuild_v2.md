# Rebuild spec v2: corrected signatures, networks and master regulators

Status: **agreed 2026-10-18, not yet executed.** Decisions D45-D47, D50-D51, D53-D60 in [`DECISIONS.md`](DECISIONS.md)
(D48, D49 and the lab-shift part of D52 were replaced by D57 on 10-18). The measurements behind them are in
[`PROGRESS.md`](PROGRESS.md) (2026-10-18) and `results/step4/reference_diagnostics/` (scripts `31`, `33`, `34`).

## Why a rebuild

1. **TCGA normal tissue and GTEx normal tissue differ a lot.** Same organ, a median of 3,844 genes (FDR < 0.05,
   |logFC| > 1; 17 groups); the difference correlates 0.40-0.86 (median 0.62) with the tumor-vs-GTEx signatures and
   its top genes are the recurring PILRB/CCNL2 module. It is **not a lab difference alone**: it mixes lab/protocol,
   surgery vs death, cancer patients vs non-cancer donors, and the state of tissue next to a tumor (field effect,
   inflammation, tumor cells). These cannot be separated with available data, so the gap is **not subtracted**.
2. **GTEx collection effects can be corrected from inside GTEx.** Ischemic time and death circumstance lower
   donor-to-donor variance in about half the tissues, and a reference shifted toward the best-preserved donors is
   closer to real best-preserved donors in all 20 tissues tested. This correction explains only 0-10% of the gap in
   point 1.
3. **GTEx colon and kidney contain samples from a different part of the organ** (muscle wall; medulla / failed
   samples), which made them far more variable than TCGA adjacent normals.
4. **Corrections never reached the drug-facing results.** Mappability and RNA-quality corrections were applied only
   to patient-level scores (Steps 8-9), not to the Step 4 signatures, the networks, or the group master regulators
   used by Step 10 and the drug side. A GC-only signature reproduces group master-regulator rankings (median rho
   0.20; KIRC 0.66).

Why GTEx stays the reference (D57, D58): it is healthy tissue, and tumor-vs-GTEx signatures reproduce across
disjoint patient halves as well as or better than tumor-vs-adjacent signatures in 8 of 11 organs (exceptions: kidney,
colon). TCGA adjacent normals are used as a check, and are being considered as a secondary drug target (open item).

## Principles

- **New version, old one kept.** Every v2 output goes to a `_v2` path; v1 stays reproducible for comparison.
- **Same methods, except where listed below.** limma (voomByGroup replaces voom, D69), ARACNe-AP settings, VIPER filters, nulls and thresholds
  are unchanged unless a step says otherwise.
- **Checks are fixed before running** (section "Acceptance checks"). They are reported whatever they show; nothing
  is tuned to make checks pass (D19).

## What changes, step by step

### Step 1 (disease groups) and Step 2 (expression data): no change
Same 49 groups, 10,211 patients, Toil GENCODE v23 data.

### Step 3 (covariates): additions and one switch
- **Add GTEx collection variables** to `results/step3/covariates_toil.tsv`:
  - `gtex_hardy` (Hardy scale 0-4, `DTHHRDY` in `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`) and
    `gtex_death` = ventilator (0) / fast (1-2) / slow (3-4);
  - `gtex_autolysis` (`SMATSSCR`, sample attributes; used only to define best-preserved donors);
  - `gtex_ischemic_min` is already present (missing for 21 of 6,662 GTEx normals; those donors are dropped).
  - **Done 2026-10-18** (`scripts/20_covariates.py`): over the 3,542 GTEx samples flagged one-per-donor, death
    circumstance ventilator 2,181 / fast 869 / slow 451 / missing 41; autolysis 0: 815, 1: 2,236, 2: 384, 3: 49,
    missing 58. All existing columns unchanged (checked).
- **Purity: ESTIMATE score everywhere (D45).** Step 4 already uses ESTIMATE; the patient-level corrections (Step 6)
  switch LGG back from ABSOLUTE to ESTIMATE.
- The RNA-quality index stays a within-role measure for patient-level correction only (D50).

### Gene set (applies to Steps 4, 5, 6, 10): drop hard-to-map genes (D46)
- Start: the protein-coding genes used now (19,797; 19,168 occur in a Step 4 signature).
- Drop genes with **gene mappability < 0.8** (Saha & Battle 2018, GRCh38/GENCODE v26 resource in
  `data/raw/crossmap/`): 1,511.
- Drop **readthrough loci** marked `drop == TRUE` in `data/interim/readthrough_genes_v23.tsv`: 344.
- Together 1,824 genes dropped. The 71 genes without a mappability value are **kept** and flagged (D54).
- Counted over all 19,797 protein-coding genes (not only those in a Step 4 signature): 2,013 dropped (1,664
  mappability < 0.8, 384 readthrough, 35 both), **17,784 kept**, 211 of them without a mappability value.
- The gene list is written once (`results/step3/genes_protein_coding_v23_v2.tsv`, dropped genes with reasons in
  `genes_v2_dropped.tsv`; `scripts/22_gene_set_v2.py`, runner step `gene_set_v2`) and every v2 step reads it.
  **Done 2026-10-18.**

### Step 4 (tumor-vs-normal signatures): prepared GTEx reference

**4a. Prepare the GTEx reference (new script `32_gtex_reference.R`; D57, D59).** One row per GTEx sample used by any
group, written to `results/step4_v2/gtex_reference.tsv`:
1. *One sample per donor per tissue*, as v1. Drop samples without ischemic time or death circumstance.
2. *Sub-region filter (D59)*, before anything else:
   - **Colon - Transverse:** keep samples whose mean log-CPM of 8 epithelial markers (CDX1, CDX2, VIL1, KRT20, CDH17,
     LGALS4, CEACAM5, EPCAM) is **>= 7**; the rest are muscle-wall samples. 165 -> 109.
   - **Kidney - Cortex:** keep samples whose mean log-CPM of 8 proximal-tubule markers (SLC34A1, LRP2, CUBN, SLC22A6,
     SLC22A8, SLC5A12, ALDOB, GATM) is **>= 5**; removed: 2 near-empty samples and 5 medulla-dominant samples.
     28 -> 21.
   - Log-CPM for the scores: TMM over the tissue's GTEx samples, prior count 2, protein-coding genes.
3. *Best-preserved donors*, per GTEx tissue after the filter: death circumstance = ventilator, ischemic time <= the
   tissue's 33rd percentile, and autolysis score <= 1 where scored.
4. *Anchor class per tissue*, from the number of best-preserved donors (n_bp):
   - **n_bp >= 45 -> "observed":** the reference is the average of the actual best-preserved donors (other donors
     still used, see 4b).
   - **10 <= n_bp < 45 -> "model":** the reference is the model's prediction at the best-preserved donors' median
     collection values.
   - **n_bp < 10 -> "average", flagged:** collection noise removed, reference = average donor (brain cortex 0,
     frontal cortex 5, kidney <= 7, cervix).
5. *Collection covariates per sample* (GTEx only; 0 for tumors), centered at the tissue's anchor:
   - `coll_isch` = log(ischemic minutes) minus the anchor value (median log ischemic time of the best-preserved
     donors; the GTEx-tissue mean for class "average"), divided by the tissue's SD of log ischemic time;
   - `coll_fast`, `coll_slow` = death-circumstance indicators minus the anchor value (0 for "observed" and "model",
     since best-preserved donors are ventilator cases; the tissue mean for "average");
   - `coll_nonbp` = 1 for GTEx donors that are not best-preserved, in tissues of class "observed" only; 0 otherwise.

**4b. Signature model (changes to `30_signatures.R`, run as v2; D57).**
- Samples: as v1 (one primary tumor per patient; GTEx samples from 4a; contrasts need >= 10 normals, CESC-AD 5).
- Genes: the v2 gene set, then `filterByExpr` + TMM as v1.
- **voomByGroup instead of voom (D69):** tumors and normals each get their own mean-variance trend, so each
  arm's values are weighted by that arm's own noise (`scripts/vendor/voomByGroup.R`, YOU-k/voomByGroup 38c8abf,
  `group = condition`). Tested on 13 groups (`scripts/36_voom_group_variance_test.R`): tumors 1.2-1.8x more variable
  than normals; slightly more reproducible signatures and regulators; controls unchanged; top-500 genes and FDR-cut
  regulator lists shift 10-30%. voom signatures are also written for check C12 (`logFC_voom` / `t_voom` / `z_voom`
  columns in the per-group files; `z_tumor_vs_<contrast>_voom.tsv`).
- **Model:** `~ condition + age + sex + ESTIMATE + tumor sample-type indicators + coll_isch + coll_fast + coll_slow
  (+ coll_nonbp)`.
  - Because the collection covariates are 0 at the anchor, the condition coefficient = tumor vs a GTEx donor at the
    anchor: the best-preserved average ("observed"; `coll_nonbp` absorbs the other donors' offset), the modelled
    best-preserved donor ("model"), or the average donor ("average").
  - All GTEx donors stay in the model, so sample size and the variance estimate are unchanged.
  - Groups with several GTEx tissues (SARC, SKCM, LGG/GBM, CESC, ESCA): covariates are centered per tissue; slopes are
    shared across the group's tissues.
  - No RNA-quality index (D50); no plate (collinear with condition, as v1). No TCGA-adjacent-based correction (D57).
- Signature value: moderated t -> signed z, as v1.
- **Adjacent contrast** (18 groups) rebuilt with the v2 gene set and the v1 model (voomByGroup as the GTEx contrast). It is a check (D58).
- **Done 2026-10-20** (`VERSION=v2 scripts/30_signatures.R`; checks `scripts/37_step4_v2_checks.R` -> `results/step4_v2/checks/`).
- Outputs: `results/step4_v2/{z_tumor_vs_gtex,z_tumor_vs_adjacent,signature_summary,gtex_reference}.tsv`,
  `data/processed/step4_v2/<group>__<contrast>.tsv.gz`.

### Step 5 (regulator networks): rebuild all 28 ARACNe networks on the v2 gene set (D46, D51)
- **Input per cohort as v1** (tumors only, subtypes pooled, COAD + READ -> COADREAD, log2(TPM+1), genes with
  TPM >= 1 in >= 10% of the cohort's tumors), **restricted to the v2 gene set**. The candidate regulator list loses
  the dropped genes.
- No GTEx-related correction (networks use TCGA tumors only). No purity/quality regression of the input (v1 check:
  < 3% of edges lose half their strength after purity, plate and quality; D21).
- ARACNe-AP settings unchanged: commit 30967ea, MI p 1e-8, 100 bootstraps, DPI, Bonferroni consolidation.
  Runtime about 2.5 days.
- **After the build:** remove edges between cross-mappable gene pairs (symmetric cross-mappability > 0; 12.7-14.2%
  of edges in v1), then `aracne2regulon` on the network + the v2 input matrix; no target cap (as v1).
- **Ordering:** the network inputs do not depend on Step 4, so the networks can be built while Step 4 is run and
  checked.
- **Two regulator sets (D64):** the full set (5,789 candidates; `<NET>.regulators.txt`, networks in `aracne_v2/`) and a
  TF-only set (Lambert TFs; `<NET>.regulators_tf.txt`, networks in `aracne_v2_tf/`, `REGSET=tf`). Same expression
  inputs and settings.
- **Extra network BRCA_LUM (D65):** LumA + LumB tumors only (769 tumors, 14,011 genes, 4,643 candidate regulators,
  1,241 TFs), built in both sets; luminal strata and groups can be scored with it next to the pooled BRCA network.
- **No proliferation regression of the input (D66).**
- Inputs built 2026-10-18 (`VERSION=v2 scripts/42_aracne_inputs.py`; `results/step5/aracne_cohorts_v2.tsv`): 29 networks,
  12,911-14,783 genes and 4,286-4,890 candidate regulators per network (v1: 13,838-15,925 and 4,522-5,163).
- Diagnostics rerun: bootstrap edge support (`45`), edge confounding (`46`).
- Outputs: `data/processed/step5/aracne_inputs_v2/`, `data/processed/step5/aracne_v2/<NET>/`,
  `data/processed/step5/regulons_aracne_project_v2.rds`, `results/step5/aracne_*_v2.tsv`.

### Step 6 (regulator activity, VIPER)

**6a. Group level (`50_viper.R` with v2 inputs).**
- msviper on each group's v2 Step 4 z signature with its cohort's v2 regulons.
- Regulator filters unchanged: expressed in the tumors (median logCPM >= 1), >= 10 targets and >= 50% of targets
  measured.
- **Null unchanged:** 1,000 random splits of the group's GTEx normals (the 4a set, so colon/kidney filtered) into two
  halves, Welch t on logCPM -> z, as v1 (v2 gene set). As in v1 the null is not covariate-adjusted while the signature
  is (assumption A3, kept to change one thing at a time; D55).
- CollecTRI and DoRothEA arms rerun on the v2 signatures (supporting evidence, D17).
- **Executed 2026-10-28 (D73).** The first v2 run built the null from all of a group's GTEx samples, not the 4a set;
  fixed in `50_viper.R` (the null now uses exactly the normals the signature used). A covariate-adjusted null (assumption
  A3) was tested (`64_null_adjust_test.R`): no false calls, more master regulators, but it lost one known-biology
  control, so by the pre-set rule it was not adopted. The test showed FDR-cut master-regulator lists are unstable
  between nulls while the NES ranking is not (rho >= 0.98); the drug side therefore weights every regulator by its NES
  instead of using an FDR list (D74).

**6b. Patient level (`61_quality_correct.R` with v2 inputs; D45, D57).**
- Tumors: gene-level regression on purity (**ESTIMATE in every group, including LGG**) + RNA-quality index, as v1.
- Normals (the 4a GTEx set): gene-level regression on RNA-quality index + the 4a collection covariates; each normal is
  set to its fitted value **at the anchor** (residual + prediction at coll_* = 0, quality index at its mean), so the
  normals each patient is compared with are centered on the same reference as Step 4.
- v2 gene set and v2 regulons; mappability edge filter already part of v2 regulons.
- Outputs: `data/processed/step6/aracne_project_v2/` (group), `data/processed/step6/aracne_project_v2_patient/`.

### Step 7 (master regulators): rules unchanged
- Exploratory: ARACNe FDR < 0.10. Confirmatory: ARACNe FDR < 0.05, or FDR < 0.10 with a same-sign FDR < 0.10 call in
  CollecTRI or DoRothEA (D18).
- Rerun: known-biology controls, three-way regulon comparison, regulator modules, recurring-module diagnostics,
  driver classes.

### Steps 8-10 and the drug side: rerun on v2 after the checks
**Status 2026-10-29 (TF arm, build `v2_tf_xmap`):** Step 8 clustering done with and without the recurring module; Step 10
subgroups chosen by a fixed keep rule (D71) on the clustering without the module, plus COAD as an exception (D79, 14
subgroups); subgroup signatures and group / patient VIPER done; drug ranking done for 41 groups and 14 subgroups with
consensus regulons (D75), NES weighting (D74), cell-line match and transfer-confidence flags (D76, D77), an evidence
status per disease (D78), specificity first in reports (D80). Step 9 replication and PRISM / GDSC validation not yet
(validation on hold, user). Details: PROGRESS 2026-10-28 / 29.

Original plan:
- **Step 8** clustering and breast axes on v2 patient scores; same pre-fixed rules (D30, D37).
- **Step 9** replication with v2 regulons (cohort-average reference; no GTEx, so no GTEx correction).
- **Step 10** stratum signatures with the v2 Step 4 model (4a GTEx reference, collection covariates, v2 gene set);
  stratum master regulators with v2 regulons.
- **Drug side** (celligner branch): rescore reversal/specificity with v2 master regulators; rerun the PRISM/GDSC
  validation. Celligner matching is unaffected. Possible secondary use of TCGA adjacent normals in drug ranking: open.
- **Open (PI):** the recurring module tracks the TCGA-normal vs GTEx axis (89% of its 555 members are higher in GTEx;
  all 17 core members in the top 0.6%). Keep it in clustering (D26) or remove it with a within-TCGA sample-level
  correction: to decide before Step 8 is rerun. Kidney reference (D60): to decide.

## Acceptance checks (fixed before running)

v1 values in brackets. Reported whatever they show.

| # | Check | Expected if the corrections work |
|---|---|---|
| C1 | GC tilt of (tumor-vs-GTEx minus tumor-vs-adjacent) z, Spearman with gene GC [median -0.24; KIRC -0.47, LumA -0.45] | Reported; the gene filter alone may not remove it, since the TCGA-vs-GTEx gap is no longer subtracted |
| C2 | GC-only signature vs group MR NES, aREA, Spearman [median 0.20; KIRC 0.66, WT 0.63] | Lower |
| C3 | Agreement of tumor-vs-GTEx with tumor-vs-adjacent signatures, 18 groups [median ~0.6; KIRC 0.33, PRAD 0.35] | Reported; expected to change little (the gap is kept) |
| C4 | Collection signal left in signatures: Spearman of tumor-vs-GTEx z with the within-GTEx ischemic-time effect [median -0.18; lung -0.43 to -0.48, LGG -0.38 to -0.42] | Near 0 |
| C5 | Known-biology controls [70 pass, 63 right direction below cutoff, 2 wrong; 6 of 16 tissue-specific pass] | Tissue-specific passes not fewer; proliferation passes not lost |
| C6 | Proliferation retained: Hallmark E2F / G2M enrichment in tumor-vs-GTEx z | Still strongly up |
| C7 | Master regulators at FDR < 0.10 per group [median ~600; CESC-AD 10, ESCA-SCC 4, UCS 2]; overlap with v1 (Jaccard) | Reported; large changes explained |
| C8 | Recurring module: present in how many groups at patient level [39 of 41]; association with plate | Reported (feeds the D26 decision) |
| C9 | Network stability: median bootstrap support, edges per network [10-12 / 100; 157k-416k edges] | Similar |
| C10 | Mappability: no gene with mappability < 0.8 among regulators or targets; no cross-mappable edges | Exactly 0 |
| C11 | Colon / kidney: split-half reproducibility with the filtered GTEx set vs TCGA adjacent [colon 0.956 vs 0.968; kidney 0.871 vs 0.920] | GTEx closer to or above adjacent |
| C12 | voomByGroup vs voom on the v2 signatures, all groups: whole-signature z rho, top-500 overlap, FDR < 0.05 / < 0.10 regulator-set Jaccard [test, 13 groups: rho 0.996-0.999, top-500 69-93%, Jaccard 0.68-0.92] | Reported; groups whose drug-side lists change most flagged |

Step 4 v2 results (2026-10-20; `results/step4_v2/checks/`):

| Check | v1 -> v2 | Verdict |
|---|---|---|
| C3 Agreement with tumor-vs-adjacent signatures | median 0.59 -> 0.65 (colon 0.56 -> 0.70, lung up) | Better |
| C4 Post-mortem signal left in signatures | median -0.18 -> -0.10; closer to 0 in 28 of 39 groups (lung -0.48 -> -0.30, colon -0.36 -> -0.16) | Better; brain unchanged (average-donor anchor, no shift) |
| C6 Proliferation retained | E2F / G2M up in 38 of 40 groups (v1 38) | Kept |
| C12 voomByGroup vs voom | top-500 genes shared: median 83%, lowest SKCM 51% | As in the test |
| C1 GC tilt (GTEx minus adjacent z) | median -0.24 -> -0.30 | Larger; accepted (D70) |
| C11 Split-half reproducibility, GTEx vs adjacent | READ, KIRC: GTEx ahead; COAD, KIRP: adjacent ahead; KICH similar | Mixed |

C1 larger GC tilt accepted (D70). Pediatric groups (TARGET-WT, NBL) flagged as unstable (A36), no action.

Step 6 v2 results (2026-10-25; `results/step6/checks_v2/`, `checks_v2_tf/`; regulons with cross-mappable edges removed):

| Check | v1 | v2 full | v2 TF-only | Verdict |
|---|---|---|---|---|
| C10 low-mappability genes / cross-mappable edges left | - | 0 / 0 | 0 / 0 | Met |
| C9 edges per network, median; bootstrap support, median | 342k; 11 | 266k (13% removed by the filter); 11 | 157k; 14 | Similar |
| C2 GC-only signature vs MR NES, median rho | 0.15 | 0.17 | 0.10 | Not lower in v2 full (KIRC 0.69, WT 0.65); fits C1 / D70 |
| C5 controls: proliferation pass, tissue-specific pass, wrong direction | 64, 6, 2 | 60, 6, 2 | 70, 8, 4 | Tissue-specific kept; 4 proliferation passes lost (direction kept) |
| C7 MRs FDR < 0.10 per group, median; v1-v2 Jaccard; NES rho | 587 | 567; 0.67; 0.977 | 172 | Reported |

## Assumptions carried by v2 (also in DECISIONS.md)

- A32: the 71 genes without a mappability value behave like normally mappable genes.
- A33: GTEx collection effects are additive per gene in log ischemic time and death circumstance; other agonal and
  handling factors are not modelled (autolysis, fixative time, RIN and batch added nothing in the held-out test).
- A34: best-preserved GTEx donors are the closest available stand-in for living tissue, but still post-mortem;
  surgical tissue has its own warm-ischemia changes that TCGA did not record. The remaining TCGA-vs-GTEx gap (lab,
  surgery vs death, donor population, tumor-adjacent biology) stays in the signatures.
- A35: the marker thresholds pick the right part of the organ in GTEx colon and kidney.
- Groups with several GTEx tissues share collection slopes across those tissues.

## Execution order and runtime (when approved)

| Order | Step | Script(s) | Runtime |
|---|---|---|---|
| 1 | GTEx collection columns; v2 gene list (**done 10-18**) | `20_covariates.py`, `22_gene_set_v2.py` | minutes |
| 2 | ARACNe inputs v2 (done 10-18); network build: full set then TF-only, 29 networks each (**done 10-24**, MAIL-15) | `42`, `43` (`VERSION=v2`, then `REGSET=tf`) | ~2.5 + ~0.7 days (background) |
| 3 | GTEx reference: sub-region filter, best-preserved donors, anchor classes, covariates (**done 10-18**) | `32_gtex_reference.R` | minutes |
| 4 | Step 4 v2 signatures; checks C1, C3, C4, C6, C11, C12 | `30_signatures.R` (v2) | < 1 h |
| 5 | Regulons v2 (+ cross-mappable edge removal; **done 10-25**); edge diagnostics (support done; confounding `46` done 10-24 for the full set, TF-only not run) | `44`, `45`, `46` | hours |
| 6 | Group VIPER v2 (3 regulon sets), MR comparison; checks C2, C5, C7, C9, C10 (**ARACNe full + TF-only and checks done 10-25; null fixed and all arms rerun 10-28 (D73)**: TF, CollecTRI, full done; DoRothEA subgroups running; `51` / `59` comparison not rerun) | `50`, `51`, `59`, `64` | hours |
| 7 | Patient VIPER v2, recurring module, check C8 (**done 10-25**: module still present in 36 of 41 groups, plate-linked; TF arm 61 done 10-26 on MAIL-15); modules (58) not yet | `61`, `58`, `62` | hours |
| 8 | Report; PI decisions on the recurring module (D26: subgroups now defined without it, D79; PI to confirm) and the kidney reference (D60, open) | | |
| 9 | Steps 8-10 and drug side on v2 (**10-28 / 29, TF arm**: Step 8 + no-module run, keep rule, 14 subgroups, subgroup signatures / VIPER / characterization, drug ranking 41 groups + 14 subgroups; Step 9 replication and PRISM / GDSC validation not yet) | `80`, `84`-`87`, `30` / `50` STRATA=; celligner `74b`, `79b`, `79d`-`79n` | ~1 day |
