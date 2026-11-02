> **v2 (2026-10-29):** this guide describes the v1 result files. For v2 outputs and where they are, see [`PIPELINE_V2.md`](../PIPELINE_V2.md) ("Where the main outputs are").

# Results and data guide

One- to two-line descriptions of every results table and data folder on the disease side (Steps 1-10).
Drug-side outputs (`results/celligner/`, `results/lincs/`, scripts 05-09 and 70-79b) are maintained in the
`~/projects/repurposing-celligner` worktree. Decisions and numbers: [`PROGRESS.md`](../PROGRESS.md).

## `results/` (small tables, in git)

### step1: disease groups
| File | What it is |
|---|---|
| `disease_grounding_toil.tsv` | Each of the 40 TCGA/TARGET cohorts mapped to standard disease, tissue and tumor-type IDs (MONDO, OncoTree, UBERON) and its matched GTEx normal tissue. |
| `mondo_descendants.tsv` | For each cohort's disease ID, all narrower disease terms under it. |
| `disease_groups_toil.tsv` | The 49 disease groups after subtype splits (e.g. BRCA -> LumA/LumB/Basal/Her2), with counts and inclusion/exclusion reasons. |
| `sample_groups_toil.tsv` | Every sample assigned to its group and role (tumor, adjacent normal, GTEx normal). |
| `target_all_lineage.tsv` | B- vs T-cell lineage calls for childhood leukemia (TARGET-ALL) samples from marker genes. |

### step3: covariates
| File | What it is |
|---|---|
| `covariates_toil.tsv` | One row per sample: sex, age, hospital, plate, RNA quality, ABSOLUTE purity, ESTIMATE scores, the one-sample-per-patient flag, and for GTEx the collection details (ischemic time, death circumstance, autolysis score). |
| `covariates_beataml.tsv` | The same for BeatAML (one sample per patient; drug-response flag). |
| `estimate_toil.tsv` | ESTIMATE stromal / immune / combined scores (non-tumor tissue mixed into each sample). |
| `genes_protein_coding_v23.tsv` | The 19,797 protein-coding genes used for signatures and drug matching. |
| `genes_extended_v23.tsv` | Protein-coding + lincRNA + antisense genes (mechanism work only). |
| `genes_protein_coding_v23_v2.tsv` | v2 gene set (rebuild spec): protein-coding genes minus hard-to-map (mappability < 0.8) and readthrough genes; 17,784 genes. |
| `genes_v2_dropped.tsv` | The 2,013 genes dropped from the v2 set, with the reason. |

### step4: tumor-vs-normal signatures (`step4_extended/` = same files, extended gene set)
| File | What it is |
|---|---|
| `signature_summary.tsv` | One row per group x comparison: sample counts, status, number of changed genes. |
| `z_tumor_vs_gtex.tsv` | Gene x group: how much each gene is up/down in tumor vs healthy GTEx tissue (signed z). |
| `z_tumor_vs_adjacent.tsv` | The same vs normal tissue next to the tumor (18 groups). |

`step4/reference_diagnostics/` (2026-10-18, `scripts/31_reference_diagnostics.R`; measurements behind the v2 rebuild spec, [`SPEC_rebuild_v2.md`](../SPEC_rebuild_v2.md))
| File | What it is |
|---|---|
| `overview.tsv` | One-line headline numbers for everything below. |
| `source_effect_summary.tsv` | Per group with both kinds of normal: how much TCGA normal tissue differs from GTEx normal tissue (the lab difference), its GC/length tilt, and how much it overlaps the tumor signature (incl. a split-half control). |
| `z_source.tsv` | Gene x group: TCGA adjacent normal vs GTEx normal (signed z). |
| `source_consensus_vs_signatures.tsv` | How strongly the lab difference shared across organs shows up in every group's tumor-vs-GTEx signature. |
| `source_consensus_hallmarks.tsv`, `source_consensus_top_genes.tsv` | What the shared lab difference consists of (pathways; top genes each way). |
| `postmortem_summary.tsv` | Within GTEx: genes changed by ischemic time and by ventilator death, and their overlap with the tumor signature and the lab difference. |
| `z_ischemic.tsv`, `z_ventilator.tsv` | Gene x group z for those two post-mortem effects. |
| `quality_index_joint.tsv`, `quality_by_role.tsv` | RNA-quality index computed on tumors and normals together; differences by role and agreement with GTEx RIN. |
| `mappability_genes.tsv` | How many Step 4 genes are hard to map or readthrough loci (dropped in v2). |
| `normal_variance_summary.tsv` (+ `_spread_reps`, `_splithalf_reps`, `_adjacent_drivers`, `_adjacent_samples`) | `scripts/33_normal_variance.R`: per organ, is TCGA adjacent normal tissue more variable than GTEx (equal n), what drives the extra spread (tumor-likeness, purity, plate), and which reference gives tumor signatures that reproduce across disjoint patient halves. |
| `gtex_collection_variance.tsv` | `scripts/34_gtex_collection_correction.R`: per GTEx tissue, how much each collection model (ischemic time, death circumstance, autolysis, fixative time, RIN, batch) lowers variance in held-out donors. |
| `gtex_collection_anchor.tsv` (+ `_anchor_folds`) | Which GTEx reference (average donor, model-shifted to best-preserved, plain best-preserved average) is closest to hidden best-preserved donors. |
| `gtex_collection_toward_adjacent.tsv` | Check only: whether shifting GTEx toward best-preserved donors moves it toward TCGA adjacent normals. |

### step4_v2: v2 rebuild (SPEC_rebuild_v2.md)
| File | What it is |
|---|---|
| `gtex_reference.tsv` | Every GTEx sample considered as healthy reference: groups using it, collection data, colon/kidney sub-region score, whether it is used (and why not), best-preserved flag, and the collection covariates Step 4 v2 adjusts for. |
| `gtex_reference_tissues.tsv` | Per GTEx tissue: samples used, best-preserved donors, how the reference is anchored (observed / model / average) and the anchor values. |

### step5: ID mapping, regulator networks, edge checks
| File | What it is |
|---|---|
| `id_map.tsv` | Our gene IDs -> current official names, Entrez and other IDs. |
| `tf_universe.tsv` | The 1,639 human transcription factors (Lambert 2018). |
| `regulons_collectri.tsv` / `regulons_dorothea_ABC.tsv` | CollecTRI / DoRothEA regulator -> target links in our gene IDs. |
| `regulons_core.rds` | CollecTRI + DoRothEA packaged for VIPER. |
| `regulon_coverage.tsv` | Per regulator and resource: number of targets and how many mapped. |
| `regulon_checksums.sha256` | Fingerprints of the regulon files used downstream. |
| `aracne_cohorts.tsv` | The 28 cancer types with our own ARACNe networks (patients, genes, regulators). |
| `aracne_regulators.tsv` | The ~6,200 regulators given to ARACNe, flagged as TF / kinase. |
| `aracne_edge_support.tsv` | Per network: how often links reappear across the 100 bootstraps. |
| `aracne_edge_confounding.tsv` | Per network: how much link strength survives removing purity, plate, subtype and RNA quality (almost all). |

### step6: regulator activity (VIPER)
`collectri/`, `dorothea_ABC/` and `aracne_project/` each contain the same core tables:

| File | What it is |
|---|---|
| `group_mr.tsv` | Every regulator's activity score (NES), p and FDR in every group: the main result. |
| `group_summary.tsv` | Per group: regulators tested, hits, patient match to group profile, number of independent patterns. |
| `regulator_filter.tsv` | Which regulators were kept/dropped per group and why. |
| `fdr10_group_summary.tsv` | Hits per group at FDR < 0.10 (up / down). |
| `fdr10_nes_matrix.tsv` / `fdr10_fdr_matrix.tsv` | Regulator x group grids of scores and FDRs for regulators that are hits somewhere. |
| `controls.tsv` | Known-biology checks (expected vs observed direction, pass/fail). |
| `patient_concordance.tsv` | Per patient: how well their regulator profile matches their group's. |
| `axes.tsv` / `axes_covariates.tsv` | Main patterns of variation between patients and how much purity, age, plate or hospital explain them. |

Additional step6 files:

| File | What it is |
|---|---|
| `collectri/subtype_test*.tsv` | Experiment: splitting groups by published subtype did not add CollecTRI hits. |
| `aracne_precomputed/summary.tsv`, `group_mr.tsv` | Results with the published precomputed ARACNe networks (0 biology checks passed). |
| `comparison/per_group.tsv` | Per group: hits in each of the three regulon methods, agreement and score similarity. |
| `comparison/consensus.tsv` | Regulators that are hits in >= 2 methods in the same direction. |
| `comparison/confirmatory.tsv` | Regulators meeting the confirmatory rule (ARACNe FDR < 0.05, or < 0.10 with curated support). |
| `comparison/controls.tsv` | Known-biology checks for all three methods side by side. |
| `aracne_project/purity_correction/` | Patient-score purity correction (ESTIMATE; ABSOLUTE for LGG): what was removed, per group and regulator. |
| `aracne_project/quality_correction{,_cqn,_map}/` | Gene-level purity + RNA-quality correction (plain / CQN test / mappability-filtered = adopted): before-after checks, patient concordance, per-sample quality score. |
| `aracne_project/diagnostics/pc_plate.tsv`, `feat_plate.tsv` | Plate signal in raw expression vs regulator scores, vs chance. |
| `aracne_project/diagnostics/estimate_pca.tsv`, `estimate_regulators.tsv` | First exploratory look at purity correction. |
| `aracne_project/diagnostics/redundancy.tsv` | How ~4,200 regulators collapse into ~600-770 co-varying groups. |
| `aracne_project/diagnostics/brca_heterogeneity_*.tsv` | Pooled BRCA: subtype vs plate vs purity variation and tests for hidden subgroups (none found). |
| `aracne_project/diagnostics/plate_quality.tsv` | Per group: does the main pattern track plate, hospital or RNA quality (answer: RNA quality differing between plates). |
| `aracne_project/diagnostics/regulator_modules_{purity,qc,map}.tsv` | Modules of co-varying regulators, hubs and what each tracks, per correction version. |
| `aracne_project/diagnostics/recurring_module_*_{qc,cqn,map}.tsv` | The module seen in almost every group: members, links to 50 biological programs and technical factors, member and target features. |

### step7
| File | What it is |
|---|---|
| `driver_classes.tsv` | Per group: share of patients with kinase-pathway vs TF-driver alterations (mutation / copy number). |

### step8: clustering
| File | What it is |
|---|---|
| `decision.tsv` | Per group and feature set: best number of clusters, p-values, verdict, early plate/purity/quality checks. |
| `m3c_scores.tsv` | M3C stability scores and p-values for 2-10 clusters. |
| `assignments.tsv` | Each patient's cluster for every number of clusters tried. |
| `module_features.tsv` | Which regulators make up each clustering feature (module). |
| `agreement.tsv` | Agreement between the main and check feature sets (ARI). |
| `known_subtypes.tsv` | Clusters vs TCGA molecular subtypes (COAD split = MSI vs CIN). |
| `defining_modules.tsv` | Modules that most distinguish each cluster, with top regulators. |
| `breast_axes_pca.tsv` | Breast: the first 5 principal components of the module scores per subtype, with null, stability and plate/quality/purity checks; `keep` = used as an axis. |
| `breast_axes_annotation.tsv` | What each breast axis is: top-loading modules and the Hallmark programs, immune and proliferation scores it tracks. |
| `breast_axes_loadings.tsv` / `breast_axes_scores.tsv` | Module weights of each axis (to project other cohorts) and each patient's axis scores. |
| `step10_strata.tsv` | The strata that get their own signature and drug matching: COAD CMS1/3-like and CMS2/4-like, READ-C1/C2/C3, 13 breast axes. |
| `step10_patient_strata.tsv` | Which COAD / READ patient is in which stratum. |

### step9: replication
| File | What it is |
|---|---|
| `gse39582_samples.tsv` | GSE39582 sample info: MMR/MSI, stage, KRAS/BRAF/TP53, CIMP, survival, ESTIMATE. |
| `cptac_coad_samples.tsv` | CPTAC-COAD sample info: MSI, stage, survival, ESTIMATE, RNA-quality score. |
| `processing_summary.tsv` | Per cohort: tumors, genes, regulators scored, MSI counts. |
| `metabric_samples.tsv` / `cptac_brca_samples.tsv` | Breast replication sample info: subtype, ER/HER2, survival (METABRIC), ESTIMATE, RNA-quality score (CPTAC). |
| `processing_summary_brca.tsv` | Per breast cohort and subtype: tumors, genes, regulators scored. |
| `brca_axes_replication.tsv` | Each breast axis in METABRIC and CPTAC-BRCA: does it emerge on its own, how much variation it captures, do the same genes follow it; verdict. |
| `brca_axes_features_check.tsv` | Module scores rebuilt as member averages vs the Step 8 scores, and projected vs original axis scores (agreement check). |
| `crc_features_check.tsv` | Module scores rebuilt as member averages vs the Step 8 scores (agreement check). |
| `crc_transfer.tsv` / `crc_assignments.tsv` | TCGA-COAD clusters transferred to GSE39582 / CPTAC (nearest centroid): MSI association; each patient's label. |
| `crc_igp.tsv` | In-group proportion: do transferred clusters hold together in the new cohort. |
| `crc_independent.tsv` / `crc_independent_scores.tsv` | Fresh M3C clustering in each new cohort: best K, p-values, agreement with transferred labels, MSI link. |
| `crc_module_effects.tsv` | Per module: MSI-vs-stable difference in each cohort vs the TCGA cluster difference. |
| `crc_cms.tsv` | CMS1-4 calls (CMScaller) crossed with clusters and MSI status in all cohorts, incl. TCGA-READ clusters. |
| `crc_summary.tsv` | One-line summary per cohort of the replication tests. |

### step10: signatures for the strata
| File | What it is |
|---|---|
| `signature_summary.tsv` | One row per stratum signature (COAD/READ subgroups, breast axis high / low ends): sample counts, changed genes, similarity to the parent group's Step 4 signature and between axis ends. |
| `z_strata_vs_gtex.tsv` | Gene x stratum: how much each gene is up/down in that stratum's tumors vs healthy GTEx tissue (signed z); the drug-matching input. |
| `z_axis_drivers.tsv` | Gene x breast axis: how strongly each gene rises toward the high end of the axis (explanation only). |
| `axis_driver_top.tsv` | The 25 genes most tied to each end of each breast axis. |
| `stratum_patients.tsv` | Which patients are in each of the 31 strata (incl. breast axis ends); `weak` flags Her2 PC4. |
| `stratum_mr.tsv` | Master regulators of each stratum (activity score, p, FDR), same format as Step 6 `group_mr.tsv`: the drug-side input. |
| `stratum_mr_summary.tsv` | Per stratum: number of master regulators, similarity to the parent group, regulators that differ most from the sibling strata. |
| `stratum_regulator_filter.tsv` | Which regulators were scored per stratum and why others were dropped. |
| `stratum_mr_parent_check.tsv` | Check that the same code reproduces the Step 6 parent-group results. |

## `data/` (local only, not in git)

| Location | What it is |
|---|---|
| `raw/solid_toil/` | Main expression data: counts and TPM for ~19,000 TCGA, TARGET and GTEx samples plus labels. |
| `raw/aml_beataml2/` | BeatAML expression, clinical and drug response (pinned commit). |
| `raw/tcga_pancan/` | TCGA pan-cancer clinical, subtypes, ABSOLUTE purity, MC3 mutations, GISTIC copy number. |
| `raw/gtex/`, `raw/target_clinical/`, `raw/tcga_gdc/` | GTEx donor/sample info; TARGET clinical; TCGA aliquot list (plates). |
| `raw/gencode/` | GENCODE v23 annotation and transcript sequences. |
| `raw/hgnc/`, `raw/mondo/`, `raw/uberon/`, `raw/oncotree/`, `raw/tf_lambert/` | Gene names, ontologies, TF list. |
| `raw/regulons/` | Saved CollecTRI / DoRothEA pulls (reused so results don't change). |
| `raw/msigdb/` | MSigDB Hallmark gene sets (50 programs). |
| `raw/crossmap/` | Saha & Battle gene mappability and cross-mappability (GRCh38). |
| `raw/replication/` | METABRIC, GSE39582, CPTAC-BRCA, CPTAC-COAD as downloaded. |
| `raw/checksums.sha256` | Fingerprints of every download. |
| `interim/gene_gc_v23.tsv`, `readthrough_genes_v23.tsv` | Per-gene GC / transcript length; readthrough genes dropped and why. |
| `processed/step4*/` | Full per-gene tumor-vs-normal tables and count matrices. |
| `processed/step5/aracne/`, `aracne_inputs/` | Per cancer type: ARACNe networks, bootstraps, edge support/confounding; network inputs. |
| `processed/step5/regulons_aracne_project.rds`, `regulons_aracne*.{rds,tsv}` | Our 28 networks as VIPER regulons; the published precomputed networks. |
| `processed/step6/<regulon>/` | Per-patient regulator scores per group (CollecTRI, DoRothEA, ARACNe). |
| `processed/step6/aracne_project_{purity,qc,cqn,map}/` | ARACNe patient scores after each correction; `_map` = current Step 8 input. |
| `processed/step8/features/` | Clustering inputs: module scores per patient (main and check sets). |
| `processed/step10/` | Full limma tables for each Step 10 stratum signature. |
| `processed/step9/<cohort>/` | Corrected expression and patient regulator scores for the replication cohorts. |
