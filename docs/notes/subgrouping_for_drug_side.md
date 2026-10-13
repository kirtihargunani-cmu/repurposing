# Subgrouping status: note for the drug-side session

From the disease-side session, 2026-10-12. Everything referenced is on branch `main` in `~/projects/repurposing`;
from the `celligner` worktree, read it with `git show main:<path>` (do not check out `main` there).
Background and numbers: `docs/PROGRESS.md` (Step 8, 9 and 10 entries).

## 1. What the strata are

Step 8 clustered patients on regulator activity in six groups (BRCA-Basal/Her2/LumA/LumB, TCGA-COAD, TCGA-READ).
Step 9 checked the results in independent cohorts. User decisions (2026-10-12) turned that into 18 strata, each to be
matched to drugs:

| Stratum | Parent group | Patients | What it is | Evidence |
|---|---|---|---|---|
| COAD-CMS1/3-like | TCGA-COAD | 92 | MSI-like; CMS1 (immune/MSI) + CMS3 (metabolic) | Replicated in GSE39582 (566) and CPTAC-COAD (105) |
| COAD-CMS2/4-like | TCGA-COAD | 193 | CMS2 (WNT/MYC) + CMS4 (mesenchymal) | Replicated, as above |
| READ-C1 / C2 / C3 | TCGA-READ | 25 / 39 / 27 | C1 proliferation-type modules, CMS mixed; C2 mostly CMS4/2; C3 CMS3-enriched | **Tentative**: plate-associated, no rectal cohort to replicate |
| 13 breast axes | BRCA subtypes | whole subtype | Continuous axes (principal components of regulator modules), see below | Checked in METABRIC and CPTAC-BRCA |

Breast had no robust clusters. Instead each subtype has continuous axes, and every patient has a score on each axis.
For drug matching each axis is cut into thirds: the **high end** (top third, positive score) and the **low end**
(bottom third) are treated as two strata. That gives 26 breast signatures.

| Axis | High end | Low end | Replication (METABRIC / CPTAC) |
|---|---|---|---|
| Basal PC1 | T cells, IFN-gamma, proliferation | Fibroblast stroma, EMT | partial / **replicated** |
| Basal PC3 | Myogenesis, junctions (FOXC1/SOX8) | Mitochondrial/housekeeping genes (possible technical) | partial / partial |
| Basal PC4 | Myeloid, complement, glycolysis | Ribosomal | partial / partial (low congruence) |
| Her2 PC2 | Glycolysis, hypoxia, EMT (FAT1/CAV1) | Lymphocyte genes | partial / partial |
| Her2 PC3 | Proliferation + immune | ESR1/XBP1 | partial / partial |
| Her2 PC4 | Interferon (STAT1) | MYC targets, ribosomal (possible technical) | **not replicated** / partial (n = 14); kept, flagged **weak** |
| LumA PC2 | Immune (IL10RA, CD74) | Myofibroblast stroma | partial / partial |
| LumA PC3 | Proliferation (AURKA, MYBL2) | ESR1/GATA3, cilia genes | **replicated** / partial |
| LumA PC4 | Cancer-associated fibroblasts (COL11A1, MMP11) | Endothelial / lymphocyte | **replicated** / partial |
| LumA PC5 | Glycolysis | Adipose (PPARG/FABP4), lymphocytes | partial / partial (low congruence) |
| LumB PC1 | EMT, glycolysis, hypoxia | Immune | partial / partial |
| LumB PC3 | Stroma, ESR1/GATA3 | T cells, checkpoints (LAG3, CTLA4) | partial / partial |
| LumB PC4 | Endothelial (CDH5, ERG, CD34) | Inflammatory / epithelial | **replicated** / partial |

The axis directions come from the 2026-10-12 PC sign. "High" and "low" are labels only; neither end is "worse".

## 2. Files

| File | Content |
|---|---|
| `results/step8/step10_strata.tsv` | One row per stratum or axis (kind, parent group, n, description, status) |
| `results/step8/step10_patient_strata.tsv` | COAD/READ patient -> stratum |
| `results/step8/breast_axes_scores.tsv` | Breast patient axis scores (thirds = high / low ends) |
| `results/step8/breast_axes_annotation.tsv`, `breast_axes_loadings.tsv` | What each axis is; module loadings |
| `results/step9/brca_axes_replication.tsv` | Breast axis replication tests and verdicts |
| `results/step10/z_strata_vs_gtex.tsv` | **Gene x stratum signed z, stratum tumors vs healthy GTEx tissue** (31 columns; Step 4 model) |
| `results/step10/signature_summary.tsv` | Per signature: n, changed genes, similarity to the parent Step 4 signature |
| `results/step10/z_axis_drivers.tsv`, `axis_driver_top.tsv` | Breast: genes that follow each axis (explanation only, not for matching) |

Stratum names, used as `group_id`: `COAD-CMS1/3-like`, `COAD-CMS2/4-like`, `READ-C1`..`C3`, `BRCA-<subtype>-PC<k>-high|low`.

## 3. The main finding for drug matching

Every stratum signature is almost the same as its parent group's Step 4 signature. Spearman over genes is 0.97-0.99,
and the two ends of a breast axis agree at 0.93-0.97. The shared tumor-vs-normal signal dominates, and the stratum-specific
part is a thin layer on top. That layer is real biology: the high-minus-low difference follows the axis driver genes
(rho 0.82-0.93). **Ranking drugs by reversal alone will give each stratum almost its parent's list.**

## 4. User decision: score by reversal AND specificity

The user wants strata ranked by both:

- **Reversal**: does the drug reverse the stratum's tumor-vs-normal state? This is your existing `rev`, computed on
  the stratum's master regulators.
- **Specificity**: does the drug reverse this stratum more than the others? Your current `spec` compares a disease
  with other tissues and deliberately excludes sibling subtypes. For strata the siblings are the comparison that
  matters, so a within-parent specificity is needed in addition. Proposal, to be agreed with the user:
  - `spec_sib = rev_z(stratum) - mean rev_z(sibling strata)`, using the same drug activity and the same robust z as
    `spec`. Siblings:
    - COAD-CMS1/3-like vs COAD-CMS2/4-like
    - READ-C1 vs C2 + C3 (and so on for C2 and C3)
    - the high end vs the low end of the same breast axis
  - Keep the existing cross-disease `spec` as well, with the strata's parent group's references.
  - A drug ranks well for a stratum if it reverses it (`rev` > 0, top tier) and is more specific to it than to its
    siblings (`spec_sib` high).

## 5. Inputs from the disease side and user decisions for the drug side (2026-10-12)

**Stratum master regulators (done, on `main`):** `results/step10/stratum_mr.tsv` is in the same columns as `group_mr.tsv`, plus
`parent_group`, with `group_id` = stratum name. Method as Step 6 / `50_viper.R` (aracne_project, the version your
`group_mr.tsv` comes from):
- msviper on the Step 10 stratum signature
- the parent's ARACNe network (BRCA / COADREAD)
- the 1,000 normal-vs-normal split null from the parent's GTEx normals
- regulator-expressed filter judged on the stratum's own tumors

Supporting files:
- `stratum_mr_summary.tsv`: MRs per stratum, similarity to the parent, MRs differing from the siblings
- `stratum_regulator_filter.tsv`
- `stratum_mr_parent_check.tsv`: the same code reproduces each parent's `group_mr.tsv`

Script: `scripts/95_stratum_mr.R`.
- Parent check: NES Spearman 0.9997 with `group_mr.tsv`, FDR < 0.10 sets Jaccard 0.85-0.94 (random null)
- MRs at FDR < 0.10: COAD-CMS1/3-like 216 (none at 0.05), CMS2/4-like 347, READ 231-292, breast 617-1,015
- NES Spearman with the parent: 0.91 for CMS1/3-like, 0.97-0.99 for the others

**Stratum membership:** `results/step10/stratum_patients.tsv` (`group_id` = parent, `stratum`, `sample`, `weak`) lists
every stratum, including the breast axis ends: exactly the tumors in each stratum signature. Sample IDs match
`weights_keyed.tsv.gz` (e.g. `TCGA-XX-XXXX-01`).

**User decisions for the drug side:**
1. **Cell-line matching per stratum = the mean of the stratum's own patients' Celligner weights** (per-patient hybrid
   weights, as `disease_line_weights` does for whole diseases). Not the parent group's weights.
2. **Target gate per stratum = judged on the stratum's own patients** (`expressed` = TPM >= 1 in >= 20% of the
   stratum's tumors, the same config rule).
3. **Scoring = reversal + specificity**, with within-parent sibling specificity as proposed in section 4 (formula to
   confirm with the user).
4. **Her2 PC4 is kept but flagged weak** (`weak = TRUE` in `stratum_patients.tsv` and `stratum_mr_summary.tsv`); carry
   the flag into drug outputs.

Network: unchanged, the parent's (BRCA for breast strata, COADREAD for colon/rectum).

## 6. Caveats to carry into the drug results

- READ-C1/C2/C3 are tentative: small (25-39), plate-associated, and not replicated.
- Breast axes:
  - Her2 PC4 did not replicate in METABRIC.
  - The low ends of Basal PC3 and Her2 PC4 look like housekeeping/ribosomal genes (possible residual technical
    signal). Rank these last.
  - Most axes are clear in METABRIC by gene content (top-100 driver genes 86-100% same direction) but don't come out
    as one clean principal component there. They are directions defined by TCGA loadings, not natural clusters.
- Breast axis thirds overlap across axes: one patient is in several strata.
- The recurring PILRB/CCNL2 regulator module (plate-linked, unexplained; kept by PI decision) loads on LumA PC2,
  Her2 PC3 and LumB PC1 (|rho| ~0.6).
- The sign of an axis is arbitrary; use the annotation, not the words high/low, to interpret it.
