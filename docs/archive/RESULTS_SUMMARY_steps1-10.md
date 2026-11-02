# Repurposing pipeline: results through Step 10

> **Status 2026-10-31: v1 results.** v2 has replaced Steps 4-10 and the drug side; for v2 see [`PIPELINE_V2.md`](../PIPELINE_V2.md) and `PROGRESS.md` 2026-10-25 to 10-29. Kept as the v1 record.

> **Update 2026-10-18:** Steps 4-10 are being rebuilt (v2) after a large TCGA-vs-GTEx lab difference was measured
> (it overlaps the tumor signatures at median rho 0.62 and matches the recurring PILRB/CCNL2 module). Numbers below
> are v1. Spec: [`SPEC_rebuild_v2.md`](../SPEC_rebuild_v2.md); decisions D45-D56 in [`DECISIONS.md`](../DECISIONS.md).

As of 2026-10-14. Written for the PI; shared version: https://claude.ai/artifact/BWMvLHK52MkqoBzy1bLegp.
Details and numbers: [`PROGRESS.md`](../PROGRESS.md), [`DECISIONS.md`](../DECISIONS.md),
[`notes/subgrouping_for_drug_side.md`](../notes/subgrouping_for_drug_side.md); drug side: `docs/notes/drug_side_decisions_and_caveats.md`
on branch `celligner`.

Steps 1-10 turned 19,131 Toil samples into 41 tumor-vs-normal signatures and 31 subgroup (stratum) signatures. Colon
splits robustly into an MSI-like and a non-MSI group (replicated in two cohorts); breast subtypes vary continuously,
captured as 13 axes. A first drug-side run on these results does not yet validate against measured drug sensitivity,
and the regulon step is the main suspect.

## At a glance

| Stage | Result | Dropped or found |
|---|---|---|
| Toil data | 19,131 samples, 60,498 genes | TCGA 10,535 · TARGET 734 · GTEx normals 7,862; one pipeline (Toil, GENCODE v23) |
| Disease groups (Steps 1-3) | 49 groups, 10,211 patients | Subtypes split (PAM50, IDH, histology); 30+ patients; one sample per patient |
| Tumor-vs-normal signatures (Step 4) | 41 signatures | 8 groups without usable normal tissue: LAML, TARGET-AML, ALL-B, DLBC, MESO, THYM, UVM, CHOL |
| Master regulators (Steps 5-7) | 28 ARACNe networks; median about 600 per group | Near-empty in CESC-AD, ESCA-SCC, UCS; only 6 of 16 tissue-specific checks pass |
| Subgroups (Steps 8-9) | 6 groups: 4 breast subtypes, COAD, READ | COAD MSI-like vs rest, replicated in 2 cohorts; READ tentative; breast continuous (13 axes) |
| Stratum signatures (Step 10) | 31 signatures | 5 colon / rectum subgroups + 26 breast axis ends; 0.97-0.99 similar to the parent |

## Steps 1-3: data and disease groups

The cohort is 49 disease groups with 10,211 patients, all processed by one pipeline (Toil recompute of TCGA, TARGET
and GTEx on GENCODE v23).

- **Data:** 19,131 samples: 10,535 TCGA, 734 TARGET, 7,862 GTEx normals; 60,498 genes, 19,797 protein-coding used for
  drug matching.
- **Disease groups:** 40 cohorts, split into subtypes where histology differs (BRCA by PAM50, LGG by IDH / 1p19q,
  ESCA, CESC, UCEC, PCPG, TGCT, SARC). Groups need at least 30 patients; BRCA normal-like dropped; breast is
  female-only.
- **Covariates:** one sample per patient (diagnosis sample first); sex, age, site, plate, ESTIMATE and ABSOLUTE
  purity, blast %.

## Step 4: tumor-vs-normal signatures

41 groups have a tumor-vs-normal signature; 8 have none because no usable normal tissue exists (LAML, TARGET-AML,
ALL-B, DLBC, MESO, THYM, UVM, CHOL).

- **Model:** limma-voom on protein-coding genes, adjusted for age, sex, ESTIMATE purity and sample type; no plate term
  (plate cannot be separated from tumor vs GTEx). Every gene keeps its signed z; no cut-off.
- **Reference:** 40 groups vs GTEx; BLCA vs adjacent normal only (9 GTEx bladders). 18 groups also have a
  tumor-vs-adjacent comparison.
- **GTEx vs adjacent agreement:** median Spearman about 0.6, from 0.33 (KIRC) and 0.35 (PRAD) to 0.80 (BRCA-Basal,
  UCEC-EEC). The two normal references disagree enough that the choice of reference matters.

## Steps 5-7: networks and master regulators

Our own ARACNe networks find master regulators in all 41 groups (median about 600 at FDR < 0.10), but the
tissue-specific known-biology checks mostly fail.

- **Networks:** 28 ARACNe networks, one per cohort (subtypes share their cohort's network; COAD and READ share one),
  full regulator set of about 6,200 (TFs, kinases, signaling proteins), no target cap. CollecTRI and DoRothEA are
  supporting evidence.
- **Edges:** median bootstrap support 10-12 of 100. Confounders barely move individual edges (< 3% lose half their
  strength after purity, plate and quality adjustment). Small cohorts (UCS, ACC, KICH) give the least stable networks.
- **Master regulators:** FDR < 0.10 in every group; near-empty in three (CESC-AD 10, ESCA-SCC 4, UCS 2). The curated
  sets find far fewer (CollecTRI: a regulator in 17 groups), and the three sets agree poorly (NES Spearman -0.02 to
  0.53).
- **Known-biology checks:** 70 pass, 63 right direction below the cutoff, 2 wrong. Proliferation regulators pass
  easily; only 6 of 16 tissue-specific checks pass (LumB ESR1 / FOXA1 / GATA3, PRAD AR / ERG, KIRC EPAS1, LIHC HNF4A
  fail).
- **Redundancy:** about 4,200 regulators per network collapse into 37-112 co-varying modules per group.

## Step 6: patient-level corrections

Two technical patterns dominated patient differences and were corrected before clustering; one unexplained module
remains.

- **Purity:** immune and stroma content was the main patient-to-patient pattern (6-35% of variance). Removed with
  ESTIMATE (ABSOLUTE for LGG).
- **Plate = RNA quality:** the plate pattern tracked a long-minus-short gene index (R^2 up to 0.94), read as RNA
  integrity differences between batches. A gene-level quality correction brought it to near chance everywhere except
  BRCA-Basal and STAD.
- **Mappability:** genes with mappability < 0.8 and cross-mappable edges removed (adopted as standard; no loss of
  biology). CQN normalization was tested and not adopted.
- **Recurring module (PILRB / CCNL2):** present in 39 of 41 groups, plate-linked above chance, not explained by RNA
  quality, GC content or multi-mapping. Kept in clustering by PI decision.

## Steps 8-9: subgroups and replication

Colon has one robust split (MSI-like vs the rest, replicated in two cohorts); rectum has tentative clusters; breast
has no clusters, so each subtype is described by continuous axes.

Clustering: M3C consensus clustering on regulator-module scores, each group separately, with the rule fixed before
running ("no clusters" accepted).

| Group | Step 8 result | Replication (Step 9) | Used as |
|---|---|---|---|
| TCGA-COAD | 2 clusters: MSI-like (92) vs rest (195) | GSE39582 (566) and CPTAC-COAD (105): MSI-like cluster holds 73 of 75 dMMR and 24 of 24 MSI-H tumors; = CMS1+3 vs CMS2+4 | COAD-CMS1/3-like, COAD-CMS2/4-like |
| TCGA-READ | 3 clusters (25 / 39 / 27), unstable, plate-associated | Not possible: no rectal cohort | READ-C1, C2, C3 (tentative) |
| BRCA-LumA | No robust clusters | 5 axes; PC3 (proliferation vs ESR1 / GATA3) and PC4 (cancer-associated fibroblasts) replicated in METABRIC | 4 axes kept |
| BRCA-LumB | No clusters | PC4 (endothelial) replicated in METABRIC; others partial | 3 axes kept |
| BRCA-Basal | No robust clusters | PC1 (immune vs stroma) replicated in CPTAC-BRCA; others partial | 3 axes kept |
| BRCA-Her2 | No clusters | Partial; PC4 (interferon) not replicated in METABRIC, flagged weak | 3 axes kept |

Breast axes were kept only if above a shuffled-data null, stable under resampling and not technical; 4 axes were
dropped as plate-linked.

## Step 10: stratum signatures

The 31 stratum signatures are 0.97-0.99 similar to their parent group's, so the stratum-specific signal is a thin
layer on top of the shared tumor-vs-normal signal.

- **What they are:** each stratum's tumors vs healthy GTEx tissue, same model as Step 4. 2 COAD + 3 READ subgroups,
  and the top and bottom thirds of the 13 breast axes (26).
- **The thin layer is real:** the high-minus-low difference follows each axis's driver genes (Spearman 0.82-0.93).
- **Master regulators per stratum:** computed with the parent's network; 216-1,015 at FDR < 0.10.
- **Possibly technical:** the low ends of Basal PC3 and Her2 PC4 look like housekeeping / ribosomal genes; ranked
  last.
- **Consequence for drugs:** ranking by reversal alone gives each stratum nearly its parent's list; strata are also
  scored by specificity against their siblings.

## Drug side: first run

The first, deliberately quick drug-side run works end to end, but its scores do not predict measured drug
sensitivity.

- **Built:** tumor-to-cell-line matching (Celligner, DepMap 24Q4), LINCS drug profiles, a target-expression check,
  and reversal and specificity scores for the 41 groups and 31 strata.
- **Known drugs:** some come out right (erlotinib in lung adenocarcinoma 96th percentile, 5-FU in colon 93rd,
  flutamide in prostate 88th); anti-estrogens come out wrong (top in basal breast cancer, middling in luminal A).
- **Generic drugs:** growth-signaling inhibitors (mTOR, PI3K, MEK, CDK4/6) rank near the top for almost every cancer.
- **Validation:** against PRISM and GDSC, no score predicts which drugs selectively kill a cancer's cell lines
  (median within-cancer Spearman about 0), although that data recovers known drugs (BRAF / MEK inhibitors in
  melanoma, HER2 drugs in Her2 breast).
- **Main suspect:** tumor-derived ARACNe regulons read drug responses wrongly. The breast ESR1 regulon is a
  luminal-identity module that lacks the estrogen-response genes (GREB1, PGR, TFF1), so fulvestrant reads as
  switching the ER program on. Splitting off proliferation did not fix it.

## Open issues to discuss

1. **Regulon assumptions:** co-expression across tumors is read as regulation; non-TF regulators (kinases, MKI67,
   TOP2A) get "targets"; target signs come from tumor correlation. Test consensus regulons (CollecTRI, DoRothEA) or
   gene-level scoring on the drug side?
2. **What the signatures measure:** tumor vs GTEx is dominated by proliferation and tissue composition, and the
   tissue-specific known-biology checks mostly fail.
3. **Master regulators are not independent:** hundreds at FDR < 0.10 move together in modules, so a few modules
   (proliferation, luminal identity) dominate any score.
4. **Shared networks:** breast subtypes, COAD / READ and all strata use their parent's network, which limits how
   different their master regulators can be.
5. **Unexplained recurring module:** kept in clustering; loads on several breast axes.
6. **Short drug signatures:** LINCS measures 6-24 h in cell lines; patient on-treatment biopsy data (for example
   aromatase inhibitors at 2 weeks and 3 months) could test whether regulons capture real drug response.
