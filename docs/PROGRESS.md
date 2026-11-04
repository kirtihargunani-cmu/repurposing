# Progress log

## 2026-10-01
- Project scaffolded at `projects/repurposing`.
- Conda env `repurposing` built (R 4.4.3; limma 3.62.1, edgeR 4.4.0, DESeq2 1.46.0, recount3 1.16.0, biomaRt 2.62.0, org.Hs.eg.db 3.20.0, ESTIMATE 1.0.13). Exact versions in `environment.lock.yml`.
- Raw data pulled (~2.9 GB) and pinned in `config/datasets.yaml`: Toil TCGA-TARGET-GTEx (GENCODE v23), BeatAML 2.0 @ `6e403e1`, MONDO v2026-09-01, UBERON 2026-10-01, OncoTree `oncotree_2025_10_03`, GENCODE v23 GTF.
- Note: Toil expected-count file is `TcgaTargetGtex_gene_expected_count.gz` (not `_rsem_`).

- Step 1 (Toil cohorts): `config/toil_disease_curation.tsv` (hand-curated) -> `scripts/10_ground_toil.py` -> `results/step1/disease_grounding_toil.tsv` + `mondo_descendants.tsv`. 40 cohorts (33 TCGA, 7 TARGET), all IDs validated against pinned MONDO/UBERON/OncoTree/GTEx labels; 10,530 tumor samples.
  - GTEx normal match: 13 good, 19 partial, 8 none (CHOL, DLBC, MESO, THYM, UVM, all leukemias).
  - GTEx tissue choices reviewed and approved by PI. Changes: HNSC -> Esophagus - Mucosa (was Minor Salivary Gland); SARC -> Muscle - Skeletal + Adipose (subcutaneous, visceral), to be narrowed per histology.
  - Mixed-histology cohorts to split before Step 4: LGG, BRCA, CESC, ESCA, PCPG, SARC, TGCT, UCEC, TARGET-ALL.
  - TARGET-RT has no tumor samples in Toil (5 normals only); TARGET-CCSK has 13 tumors.
  - DepMap lineage column empty until DepMap is pulled.
- Step 1b (subtype splits): `config/subtype_splits.tsv` -> `scripts/11_split_subtypes.py` -> `results/step1/disease_groups_toil.tsv`, `sample_groups_toil.tsv`, `target_all_lineage.tsv`. Labels from TCGA PanCanAtlas clinical/subtypes (now pinned in datasets.yaml).
  - Splits: ESCA (SCC/AD), CESC (SCC/AD), UCEC (endometrioid/serous), PCPG (pheo/PGL), TGCT (seminoma/NSGCT), SARC (LMS/DDLPS/UPS+MFS), LGG (IDH-mut non-codel / codel / IDH-wt), BRCA (PAM50 LumA/LumB/Basal/Her2).
  - TARGET-ALL lineage from marker expression: 192 B, 1 T (relapse), 1 unassigned -> group ALL-B.
  - Inclusion rule: >= 30 tumor *patients*. 49 groups included (10,383 tumors / 10,224 patients).
  - Excluded: BRCA normal-like (39, low tumor content), TARGET-AML-IF (32 samples but 24 patients), CCSK, RT, ALL-T, CESC adenosquamous, UCEC mixed, SARC MPNST/synovial/desmoid; 17 TGCT + 3 LGG tumors lack labels.
  - Step 4 note: TARGET cohorts have repeat samples per patient (relapse, blood + marrow); keep one primary sample per patient.
- Step 3 (covariates): `scripts/02_download_covariates.sh` (GTEx v8 annotations, TCGA plates via GDC API, LAML blasts, PanCanAtlas ABSOLUTE, TARGET clinical) -> `scripts/21_estimate.R` + `scripts/20_covariates.py` -> `results/step3/`.
  - `covariates_toil.tsv` (19,565 sample-group rows): source, sample type + class, `primary_per_patient`, `patient_has_untreated_tumor`, sex, age (GTEx = bin midpoint, `age_is_binned`), site (TCGA TSS / GTEx center / TARGET protocol), batch (TCGA plate / GTEx expression batch), RIN, ischemic time, ABSOLUTE purity, ESTIMATE stromal/immune/score, blast %, tumor content %.
  - primary_per_patient: diagnosis (01/09) > PB at dx (03) > new primary (05) > metastatic (06) > recurrent (02/04/40) > post-treatment (41); GTEx donors keep the first-listed tissue of the group. 18,175 rows flagged; chosen tumors: 9,948 primary, 365 metastatic (364 SKCM), 47 recurrent/post-treatment only.
  - ESTIMATE (illumina mode, scores only) on 15,324 samples; vs ABSOLUTE purity Spearman -0.59 overall, median -0.64 per cohort. ABSOLUTE covers 90% of TCGA primary tumors.
  - Blast %: TCGA-LAML 173/173, TARGET-AML 148/162, TARGET-ALL 7/151 (TARGET blast file mostly covers other patients).
  - `covariates_beataml.tsv`: 671 RNA samples / 613 subjects; one per subject (initial diagnosis > residual > relapse; marrow > blood); 472 flagged samples have drug response.
  - `genes_protein_coding_v23.tsv`: 19,797 genes.
  - The one ALL "T-lineage" relapse sample belongs to a clinically B-precursor patient (excluded anyway).
- BeatAML one-per-subject now prefers the drug-tested specimen within a disease stage (475/479 drug-tested subjects covered).
- Step 4 (signatures): `scripts/30_signatures.R` -> `results/step4/signature_summary.tsv`, `z_tumor_vs_gtex.tsv`, `z_tumor_vs_adjacent.tsv` (gene x group signed z); full topTables in `data/processed/step4/` (not in git).
  - limma-voom, protein-coding genes, filterByExpr + TMM; ~ condition + sex + age + ESTIMATE score + tumor sample-type indicators; one sample per patient/donor; contrasts need >= 10 normals.
  - 38 groups have a tumor-vs-GTEx signature; 18 also tumor-vs-adjacent (BLCA adjacent only: GTEx bladder n=9). Spearman agreement GTEx vs adjacent: median ~0.6 (0.33 KIRC, 0.35 PRAD to 0.80 BRCA-Basal/UCEC-EEC).
  - No signature (no usable normal): TCGA-LAML, TARGET-AML, ALL-B, DLBC, MESO, THYM, UVM, CHOL (9 adjacent).
- Step 4 revisions (PI decisions):
  - CESC-SCC normals = GTEx cervix + vagina (same non-keratinized squamous epithelium; 87 donors). CESC-AD = cervix only, min normals lowered to 5 (7 donors; low power). Config: `gtex_tissues`, `min_normals` in `config/subtype_splits.tsv`.
  - BRCA female-only throughout (`restrict_sex` in `config/subtype_splits.tsv`; PI: male breast cancer is biologically distinct): 79 female GTEx donors; 13 male tumors (5 LumA, 8 LumB) and 1 male adjacent normal dropped. No BRCA model has a sex term.
  - 40 groups now have a protein-coding signature (58 contrasts).
  - Extended signature (protein-coding + lincRNA + antisense; `GENESET=extended`) in `results/step4_extended/` for mechanism work only; drug matching uses protein-coding. Protein-coding z identical between runs (Spearman >= 0.997); ~5.9k lincRNA + 4.5k antisense genes pass filterByExpr.
- Scope decision (PI, 2026-10-01): Steps 5-7 for all groups with a signature; Steps 8-10 for BRCA and COAD/READ. De-prioritize hematologic tumors (AML/ALL/DLBC): limited and more variable data. Steps 5-10 data requirements in `docs/archive/PROJECT_PLAN_steps5-10.md`.
- Environment moved to R 4.5.3 / Bioconductor 3.22 when the Steps 5-10 packages were installed; Steps 3-4 rerun (ESTIMATE identical; z Pearson >= 0.9998 vs R 4.4).
- Hematologic normal references: recount3 / GDC search recorded in `docs/archive/normal_reference_search.md`; parked per PI.
- Step 5 (regulons):
  - `scripts/40_id_map.py` -> `results/step5/id_map.tsv`: GENCODE v23 Ensembl -> current HGNC (Ensembl match, then symbol, unique previous symbol, unique alias). 97.4% of protein-coding genes mapped (522 unmapped, mostly retired clone names); 5,113 v23 symbols have since been renamed, so symbol joins would have lost them. `primary_for_hgnc` picks one v23 gene per HGNC record (176 duplicates dropped).
  - `tf_universe.tsv`: Lambert 2018 'Is TF? = Yes' (1,639; 1,634 in v23).
  - `scripts/41_regulons.R`: CollecTRI and DoRothEA A-C pulled from the OmniPath REST API (decoupleR/OmnipathR 4.0 wrappers fail against the current server; raw pulls saved in `data/raw/regulons/`), NF-kB / AP-1 complexes collapsed as in decoupleR; aracne.networks (25 TCGA cohorts, Entrez) mapped to Ensembl.
  - Coverage (`regulon_coverage.tsv`): CollecTRI 1,178 regulators / 42,601 edges (99.6% mapped; 614 with < 10 targets); DoRothEA A-C 429 / 32,286 (99.6%; 158 < 10); ARACNe ~6,050 regulators per cohort, ~1,500 of them Lambert TFs, median 95-97% targets mapped.
  - `regulons_core.rds` (CollecTRI + DoRothEA) in git; ARACNe tables and `.rds` in `data/processed/step5/` (1.2 GB). sha256 of both in `regulon_checksums.sha256`.
- Step 5c/d (project ARACNe, PI decision 2026-10-02): own networks for all cohorts with a signature, full regulator set (Lambert TFs + every aracne.networks regulator: 6,201, incl. 1,633 TFs and 375 kinases) because drivers may be kinases / signaling proteins. `scripts/42_aracne_inputs.py` (per cohort, subtypes pooled, COAD+READ -> COADREAD; log2(TPM+1), protein-coding with TPM >= 1 in >= 10%) -> `scripts/43_run_aracne.sh` (ARACNe-AP commit 30967ea, OpenJDK 17 in conda env `aracne`; MI p 1e-8, 100 bootstraps, DPI, Bonferroni-consolidated). Benchmark on BRCA: AP 7.6 min / bootstrap (108k edges); ARACNe3 (commit 3d8791a) 5 min / subnetwork but only 34k edges. Running detached (`logs/aracne.log`), order: BRCA, COADREAD, small cohorts, then rest; ~2.5 days. Low-n cohorts: ACC 77, KICH 66, UCS 57.
- Step 6 (VIPER, CollecTRI arm) `scripts/50_viper.R`: per-group regulon filters (regulator expressed: median logCPM >= 1; >= 10 targets and >= 50% present), msviper on the Step 4 z signature, per-patient VIPER vs the group's normals, patient-vs-group concordance, PCA axes of mechanism variation.
  - Null model: label permutation leaked tumor signal into the null (large effect + unbalanced groups), so the null is now 1,000 random splits of the group's normals (keeps gene-gene correlation, no tumor signal). This null is conservative: |NES| peaks ~3-4.
  - PI decision 2026-10-02: permutation (normal-split) null only. The analytic null (msviper without null; decoupleR default) assumes independent genes; regulon targets are co-expressed, so it overstates significance (test groups: |NES| ~7, 3-4x more calls, same ranking rho ~0.99). Available only as an opt-in internal check (`ANALYTIC=1` in `scripts/51_compare_mr.R`), never reported.
  - Thresholds: exploratory selection FDR < 0.10 (PI); confirmatory claims FDR < 0.05 or agreement with the ARACNe arm. Patient concordance uses the group's top 50 regulators with p < 0.05.
  - `scripts/51_compare_mr.R`: per-group counts at FDR < 0.10, regulator x group NES matrix, known-biology controls.
- Step 6 CollecTRI results (all 41 groups): at FDR < 0.05 / 0.10 / 0.25, 12 / 17 / 33 groups have >= 1 regulator (64 distinct regulators at 0.10, 176 at 0.25). Passing count unrelated to group size (Spearman with n_tumor -0.25). Known-biology checks ~90% correct direction; 13 pass at 0.10, 26 at 0.25 (gain mostly E2F1 proliferation). Groups with nothing even at 0.25: THCA, SKCM, PAAD, LUAD, PRAD, READ, ACC, ESCA-AD.
- Mechanism heterogeneity (large groups): patients correlate with their group's MR profile (median rho 0.6-0.9) but show only 15-36% of the group's top MRs as clearly altered; PC1 of patient NES is largely tumor purity (ESTIMATE R^2 0.6-0.92 in most large groups) -> adjust patient NES for purity before Step 8 clustering. Effective dimensionality 10-25.
- Driver classes in our patients (`scripts/60_driver_classes.py` -> `results/step7/driver_classes.tsv`; MC3 PASS + GISTIC, functional events only): kinase-pathway alteration is similar in groups with and without regulators (median 68% vs 60%); kinase-driven null groups: THCA (BRAF 60%), SKCM (BRAF 48%, NRAS 30%), PAAD (KRAS 66%), LUAD, LumA (PIK3CA 36%). Literature summary: `docs/archive/driver_mechanisms_null_groups.md`.
- Subtype-split experiment (`scripts/52_subtype_test.R`; PanCanAtlas subtypes in 9 large groups, size-matched random controls): no gain (2 of 33 subtypes with any FDR < 0.10 regulator; median 0 for subtypes, random subsets and whole groups). Not pursued with CollecTRI.
- Diagnosis: CollecTRI regulons do not match tissue context. Regulators with large own-mRNA increases have flat or wrongly signed CollecTRI targets (ESR1 in LumA z +22 but activated targets flat; SOX2 in LUSC +15 with 'repressed' targets up; ERG in PRAD with 'repressed' targets up). Testing tissue-specific ARACNe regulons (`scripts/53_aracne_precomputed_test.R`).
- Step 6 (VIPER, project ARACNe arm; BRCA + COADREAD networks, 2026-10-04). `scripts/44_aracne_regulons.R` (aracne2regulon on network.txt + the cohort's ARACNe input matrix; full regulator set, no target cap - PI decision 2026-10-04) -> `data/processed/step5/regulons_aracne_project.rds` (BRCA 4,936 regulators / 284,818 edges; COADREAD 4,728 / 343,603). `scripts/50_viper.R` with `REGULON=aracne_project` (each group uses its cohort's network) -> `results/step6/aracne_project/`.
  - FDR < 0.10: Basal 910, LumB 813, Her2 731, LumA 690, READ 333, COAD 321 regulators (~20% Lambert TFs, ~7% kinases); CollecTRI had 9 / 13 / 4 / 0 / 0 / 0. COAD has only 1 at FDR < 0.05 (BH plateau).
  - Known-biology controls: all 28 tested in the correct direction, 20 pass FDR < 0.10 (ESR1/FOXA1/GATA3 in LumA; MYC in COAD/READ; proliferation E2F1/FOXM1/MYBL2). Precomputed aracne.networks passed 0 at 0.10.
  - NES agreement with CollecTRI on shared regulators is low (Spearman 0.11-0.19).
  - Patient level: median frac_hit 0.8-1.0 (CollecTRI ~0.3); effective dimensionality 5.5-10.5. BRCA PC1 is associated with TCGA plate (R^2 0.56-0.79; plate has many levels, so R^2 is inflated, but CollecTRI PC1 was ~0.2) and purity has moved to PC2 (R^2 0.63-0.87). COAD/READ PC1 not batch-dominated. Unresolved; to review before Step 8.
- Patient-NES diagnostics (`scripts/54_patient_nes_diagnostics.R` -> `results/step6/aracne_project/diagnostics/`; no pipeline outputs changed):
  - Plate is already in raw tumor expression in all BRCA subtypes: genes with plate p < 0.01 = 15-76% (chance 1%); raw-expression PC1 R^2 on plate 0.31-0.75 vs shuffled-plate chance 0.07-0.31. ARACNe NES amplify it (50-88% of regulators; PC1 = 28-38% of variance, R^2 plate 0.56-0.79, all p = 0.001). COAD/READ ARACNe PC1 is not plate-linked.
  - Redundancy: ~4,200 regulators collapse to ~600-770 correlation modules (|r| >= 0.7). Proliferation is a small module (7-15 regulators in BRCA, 78 in COAD; 1-4% of regulators track proliferation). The large modules are (i) the BRCA plate axis: 300-620 generic regulators (ZFX, KMT2E, ALS2, BAD, ARF5), = PC1; (ii) immune/purity: ~150-190 regulators (BTK, CD4, ITGAL), r 0.93-0.96 with ESTIMATE; (iii) in COAD/READ a stromal module (SNAI2, WWTR1, HGF), r 0.86-0.88 with ESTIMATE.
  - ESTIMATE regression of ARACNe NES removes 11-16% of variance; the purity PC disappears (ESTIMATE R^2 0.63-0.87 -> 0; ABSOLUTE R^2 after <= 0.11); 6-9% of regulators are > 50% explained (immune BTK, IKZF1, SPI1, PRDM1, FOXP3, CIITA; in COAD/READ also TWIST1, SNAI2, PRRX1, ZEB1). In BRCA it makes the plate PC relatively larger (e.g. Basal 37% -> 42%).
- BRCA patient heterogeneity (`scripts/55_brca_heterogeneity.R` -> `results/step6/aracne_project/diagnostics/brca_heterogeneity_*.tsv`; 4 PAM50 subtypes pooled, 1,040 patients x 4,018 shared ARACNe regulators; views: as-is, ESTIMATE removed, ESTIMATE + plate removed, subtype protected; the corrected views are what-ifs, not a chosen design):
  - Per regulator, plate explains more variance than PAM50 subtype as-is (median 19% vs 13%); after removing plate, subtype rises to 18% and pooled PC1 becomes 86% subtype. Top subtype regulators FOXA1, XBP1, GATA3, ESR1, FOXC1, CDKN2A.
  - PAM50 separation (silhouette, top 10 PCs): Basal distinct (0.25-0.32), Her2 0.10-0.18, LumA/LumB overlap (-0.01-0.08).
  - Within subtype: effective dimensionality 5-8 as-is, 10-14 after plate removal. k-means (k 2-6) silhouette vs a Gaussian null with the same spread (M3C-style): within noise everywhere except LumA k = 2, which is barely above it. As-is that LumA split is plate (R^2 0.54, chance 0.07); after plate removal it is long nuclear regulators (ATRX, ASXL2, PBRM1) vs short housekeeping genes (POLR2I, NDUFA13, ELOB), which looks like an RNA-quality pattern (unverified). No evidence of discrete groups inside subtypes yet; variation looks continuous.
- Project ARACNe arm, all cohorts (2026-10-07): all 28 networks finished (100 bootstraps each; 157k UCS - 416k LUAD edges); regulons for the 26 remaining networks built with `scripts/44_aracne_regulons.R` (BRCA/COADREAD untouched); `scripts/50_viper.R` + `51_compare_mr.R` rerun on all 41 groups (`REGULON=aracne_project`). Oct 4 BRCA/COAD/READ outputs kept in `data/processed/step6/aracne_project_run_2026-10-04/`; rerun reproduces them (group NES r 0.999, 98.7% of FDR < 0.10 calls agree; patient NES identical, so the diagnostics stand).
  - All 41 groups have >= 1 regulator at FDR < 0.10 (median ~600; 190-2,033 except CESC-AD 10, ESCA-SCC 4, UCS 2). Passing count unrelated to n_tumor (Spearman 0.00). The 8 groups with nothing in CollecTRI even at FDR < 0.25 now have 204-773. Several groups have 0 at FDR < 0.05 but hundreds at 0.10 (COAD, BLCA, HNSC, STAD, ESCA-SCC, UCS): broad moderate signal, few extreme p-values.
  - Known-biology controls: 70 pass, 63 correct direction below cutoff, 2 wrong (KIRC EPAS1 NES 0.00; UCS E2F1 -0.22), 4 filtered. Tissue-specific controls weaker than proliferation: 6 of 16 pass; KIRC EPAS1/HIF1A, LIHC HNF4A/HNF1A, PRAD AR/ERG, LumB ESR1/FOXA1/GATA3 do not reach FDR < 0.10. Proliferation controls fail mostly in slow-growing tumors (THCA, PRAD, KIRC, KICH, PCPG, LIHC).
- Purity correction, all 41 groups (`scripts/56_purity_correct.R`, adopted 2026-10-07): each regulator's ARACNe patient NES regressed on ESTIMATE score per group (residual + mean) -> `data/processed/step6/aracne_project_purity/`; summary in `results/step6/aracne_project/purity_correction/`. No patients dropped. Variance removed 6-35% (TGCT-SEM 35%, SARC-DDLPS 28%, PCPG-PGL 21%; most 9-19%); after correction no top PC tracks ESTIMATE (R^2 0). Top purity-driven regulators are immune everywhere (EVI2A, CD53, HCLS1, LCP2, CD4, CSF1R).
  - ABSOLUTE (DNA) purity is not fully removed where ESTIMATE and ABSOLUTE agree poorly: LGG-IDHwt (rho -0.33; ABSOLUTE PC R^2 0.26 -> 0.41), LGG-IDHmut-codel (-0.28; 0.28 -> 0.32), LGG-noncodel, TGCT-NSGCT (0.27 after), KICH (0.20). Weak ESTIMATE-ABSOLUTE agreement also in PAAD, ESCA-SCC, PRAD, THCA, CESC-AD.
  - Plate is not BRCA-only: after purity correction PC1 R^2 on plate far above chance in STAD (0.59 vs 0.04), PRAD (0.59 vs 0.04), UCEC-SER (0.76 vs 0.22), ESCA-AD, LGG, OV (0.28 vs 0.02), CESC-SCC, UCEC-EEC as well as BRCA.
  - TGCT-SEM: seminomas are lymphocyte-rich by definition, so the 35% removed is likely partly real biology.
- Edge bootstrap support (`scripts/45_edge_support.R` -> `data/processed/step5/aracne/<NET>/edge_support.tsv.gz`, summary `results/step5/aracne_edge_support.tsv`): per edge, number of the 100 bootstraps containing it. Median support 10-12 / 100 in every network (minimum 6-9 = consolidation threshold); only 0.3-7% of edges in >= 50 bootstraps (BRCA highest, 7%; UCS/ACC/KICH lowest); support correlates weakly with MI (Spearman 0.18-0.30). Small cohorts produce many more distinct edges across bootstraps (UCS 24M, ACC 19M vs BRCA 3.6M): less stable.
- Edge confounding diagnostic (`scripts/46_edge_confounding.R` -> `edge_confounding.tsv.gz` per network, summary `results/step5/aracne_edge_confounding.tsv`; partial Spearman of regulator-target in the network's input matrix after removing ESTIMATE, plate, both, or subtype; no edges removed). Edges are mostly strongly monotone (median |r| 0.47-0.71; 1-7% with |r| < 0.2). Confounders barely touch individual edges: median |r| retained 0.92-0.99 after ESTIMATE + plate; edges losing > 50%: 0.07-0.8% (shuffled-plate control 0%). Between-subtype-driven edges: ESCA 3.6%, TGCT 2.0%, BRCA 1.7%, LGG 1.3%. Low-support edges are ~5-10x more often confounded than the top-10%-support edges, but both are < 1.5%. Implication: the plate / purity patterns in patient NES are not caused by confounded edges; they come from small, coherent shifts across many genes that VIPER adds up, so filtering edges would not fix them.
- Pipeline runners (2026-10-08): `scripts/run_steps1-4.py` (00-30) and `scripts/run_steps5-7.py` (03, 40-46, 50/51 per regulon, 56, exploratory 52-55 + 60). Python orchestrators over the existing step scripts (PI choice: no Python re-implementation, results unchanged); skip steps whose outputs exist, `--from/--only/--to/--force/--dry-run/--list`, per-step logs in `logs/pipeline/`, stop on first failure. Download steps only re-run with `--redownload` (some sources are fetched as "latest"). ARACNe build is opt-in (`--build-aracne`). Verified: re-running `ground` and `edge_support` through the runners reproduces byte-identical outputs. DoRothEA VIPER has never been run; the default regulon set includes it.
- 2026-10-08 decisions (PI): LGG purity correction uses ABSOLUTE (56_purity_correct.R ABSOLUTE_GROUPS; 16 of 508 LGG patients without ABSOLUTE dropped; ABSOLUTE-linked PC R^2 0.10-0.28 -> 0; other groups unchanged). TGCT-SEM immune question tabled (immune-defined tumors later). Step 8: "no clusters" is an accepted outcome when variation is continuous.
- DoRothEA A-C VIPER, all 41 groups (via run_steps5-7.py): 21 groups with >= 1 regulator at FDR < 0.10; controls 18 pass, 111 correct direction below cutoff, 6 wrong (FOXA1/GATA3 in LumA/LumB, ERG in PRAD). Three-way comparison (`scripts/59_compare_regulon_sets.R` -> `results/step6/comparison/`): consensus (FDR < 0.10 in >= 2 sets, same sign) = 128 group-regulator pairs in 20 groups, mostly proliferation (E2F1-4, TFDP1, FOXM1, MYBL2); 21 groups have no consensus regulator (incl. LumA, COAD, READ, LUAD, LUSC, PRAD, LIHC). NES agreement: CollecTRI-DoRothEA rho 0.15-0.53, DoRothEA-ARACNe 0.13-0.45, CollecTRI-ARACNe -0.02-0.30.
- Plate investigation (`scripts/57_plate_quality.R` -> `diagnostics/plate_quality.tsv`; purity-corrected ARACNe NES): plate is not mainly hospital. Plate and TSS overlap (Cramer's V 0.4-0.7), but within a hospital plate still explains PC1 (extra R^2 0.36-0.53 in BRCA, PRAD, STAD) while hospital adds little within plate (0.04-0.08). An RNA-quality proxy (per-sample long-gene minus short-gene expression index) explains PC1 at R^2 0.89-0.94 in all BRCA subtypes, PRAD, UCEC-SER, LGG-codel, 0.74-0.87 in STAD, LIHC, OV, SKCM; plate explains 38-76% of that index. Interpretation: plate effect = RNA integrity differences between processing batches (not verified with RIN/TIN; Toil has no BAMs).
- Regulator modules (`scripts/58_regulator_modules.R` -> `diagnostics/regulator_modules.tsv`; purity-corrected NES, |r| >= 0.7): 37-112 modules of >= 10 regulators per group. Members share few targets (median Jaccard 0.06), so modules are co-varying regulators, not duplicated regulons. Hub's regulon contains a median 42% of module members (chance 0.7%), but ARACNe edges are undirected, so hub = best representative, not proven upstream. Recognizable modules: proliferation (CCNA2/FOXM1/MKI67, all BRCA, OV, PAAD), endothelial (TIE1/DLL4/FLT4/SOX7/18; KIRC ERG-EPAS1 module), fibroblast (PDGFRB/DDR2/THY1), interferon (STAT1/PARP14/CD274), myeloid (MS4A6A/TLRs), melanocyte (MITF), chr19 KRAB-ZNF cluster (LGG). ~29% of module members sit in modules tracking the RNA-quality index (e.g. LumA largest module, 927 regulators, r = -0.94); 7% in proliferation modules.
- 2026-10-08 decisions (PI), adopted from the standards table: ARACNe is the primary regulon set, CollecTRI/DoRothEA as supporting evidence; failing known-biology controls are reported as limitations (no tuning); confirmatory = ARACNe FDR < 0.05, or ARACNe FDR < 0.10 with a same-sign FDR < 0.10 call in CollecTRI or DoRothEA (`results/step6/comparison/confirmatory.tsv`); regulator modules rebuilt after quality correction, one score per module for Step 8; ARACNe consolidation threshold kept; DepMap needed for drug matching; sensitivity checks after plate/quality correction. Small/weak groups (CESC-AD, ESCA-SCC, UCS) are KEPT for drug matching (PI: see what matches). Replication: use both cohorts per cancer (breast METABRIC + CPTAC-BRCA; colorectal GSE39582 + CPTAC-COAD); not yet downloaded.
  - FDR < 0.05 counts: CollecTRI 25 distinct regulators in 12 groups; DoRothEA 30 in 16 groups; ARACNe 4,483 in 35 groups (0 in BLCA, COAD, HNSC, STAD, ESCA-SCC, UCS). Confirmatory rule adds 17 group-regulator calls via curated support (e.g. E2F1-3 in WT, TFDP1 in DDLPS, MEIS2/NR4A1 in BLCA); 36 groups have >= 1 confirmatory regulator.
- RNA-quality correction (`scripts/61_quality_correct.R`, PI decision 2026-10-08; new Step 8 input `data/processed/step6/aracne_project_qc/`): per-sample long-minus-short-gene index; gene-level regression (tumors: purity [ESTIMATE; ABSOLUTE for LGG] + quality; normals: quality), then patient VIPER with the same ARACNe regulators. Results (vs NES-level purity correction): quality-linked PC R^2 up to 0.94 -> <= 0.014 in every group; BRCA PC1 plate R^2 0.56-0.79 -> 0.14-0.36 (chance 0.07-0.31), PRAD 0.59 -> 0.06, OV 0.28 -> 0.01, UCEC-SER 0.76 -> 0.26; residual plate above chance in BRCA-Basal (0.36 vs 0.19) and STAD (0.14 vs 0.04). Effective dimensionality BRCA 5-8 -> 17-19. Patient concordance with group MRs unchanged (median change +0.007; largest drop CESC-SCC 0.62 -> 0.54). Caveat: the index agrees only partly and inconsistently with GTEx RIN (r -0.35 to +0.48 by tissue) and ischemic time (rho -0.37 to +0.50), so it is a technical axis tied to plate, not a validated RIN substitute.
- Replication cohorts downloaded (`scripts/04_download_replication.sh`, registered in config/datasets.yaml): METABRIC 1,980 x 20,603 (PAM50/claudin-low, ER IHC); GSE39582 585 samples (CIT; KRAS/BRAF/TP53, MMR, CIMP, stage, survival); CPTAC-BRCA 133 and CPTAC-COAD 106 STAR-count tumors (GENCODE v36; md5-verified) + cBioPortal clinical (PAM50 / MSI), proteome and phosphoproteome. CPTAC has no normal RNA-seq.
- Celligner (2026-10-08; drug side, patient -> cell-line matching). DepMap 24Q4 from Figshare (`scripts/05_download_depmap.sh`; newest full release there, portal API is behind a browser challenge; md5-verified). Env `celligner` (`scripts/06_setup_celligner_env.sh`, `environment.celligner.lock.yml`): Python port broadinstitute/celligner d9c9246 + DeKegel/mnnpy bf1dd73 (PyPI 1.1.0 is the older API).
  - `scripts/70_celligner.py`: full fit, package defaults (PCA 70, cPCA 4, top 1,000 DE genes, MNN k1 5 / k2 50). Reference: all 1,673 DepMap models, log2(TPM+1). Target: all 10,211 primary_per_patient tumors of the 49 included groups (Toil log2(TPM+0.001) -> log2(TPM+1)), no purity correction (ESTIMATE / ABSOLUTE used afterwards as diagnostics). Genes: DepMap Entrez -> v23 Ensembl via id_map: 18,946 of 19,193 mapped, all shared. 37 line / 51 tumor clusters, 1,698 DE genes, cPCA eigenvalues 12.4 / 5.6 / 5.3 / 1.6. 12 min. Full refit, not projection into the published alignment: our tumors are the TCGA/TARGET-type target the published run uses. BeatAML not projected (not TPM scale; heme de-prioritized). Package quirk kept as published: its limma wrapper casts expression to int before the cluster DE step.
  - `scripts/71_celligner_outputs.py` -> `results/celligner/`: per-patient top-10 lines, patient QC, group summary, line coverage, purity diagnostic, lesion check; softmax weights in `data/processed/celligner/weights_softmax.tsv.gz` (tau 7.8 = median (d10 - d1) / ln 10, provisional). Candidates exclude 85 non-cancerous models (1,588 remain); 189 never appear in any patient's top 10.
  - Lineage sanity check: nearest line in the group's own tissue for >= 85% of patients in LAML/TARGET-AML (99-100%, heme lines), NBL, UVM, GBM/LGG, KIRC, LIHC, PRAD, THCA, COAD/READ, OV, BRCA LumA/LumB/Her2 (88-94%; Basal 71%), SKCM, TGCT. Low but biologically sensible: CESC-SCC / ESCA-SCC / LUSC -> head-and-neck and esophageal squamous lines; THYM -> T-ALL lines (no thymic lines; thymocytes); TARGET-WT -> rhabdomyosarcoma (DepMap has 5 Wilms lines, CCLF_PEDS_*, but unrestricted matching does not pick them); PCPG -> neuroblastoma (2 adrenal lines in DepMap).
  - Hub lines (one line takes most of a group's weight): KNS-42 for all GBM/LGG (mean weight 0.92-0.96), SK-N-MM for PHEO (0.87), Hep G2 for LIHC (0.75), NCI-H295R for ACC (0.75; the only ACC line), A-704 for KIRC (0.45). LGG sits far from every line (median nearest distance 61-116 vs ~35-70 elsewhere): no IDH-mutant glioma lines, so LGG has no real model. KICH matches nothing (median distance 126, weight spread over ~360 lines; no chromophobe lines).
  - Purity: the runbook's concern (low purity -> far from lines) is mostly absent; where distance tracks ABSOLUTE purity it runs the other way (purer tumors farther): FDR < 0.05 in LUAD (rho 0.38), OV 0.28, LIHC 0.25, PHEO 0.49, KICH 0.62, THYM 0.57, LumA 0.12, SKCM 0.13; negative in TGCT-SEM (-0.60), SARC-LMS (-0.35), KIRP (-0.17). So low-purity tumors are not flagged by distance and can get confident matches; `patient_qc.tsv` flags them (ABSOLUTE < 0.3 or blasts < 60%; 875 patients, incl. 127 TCGA-LAML, 104 LUAD) instead.
  - Lesion check (TCGA, MC3 hotspots + GISTIC amp vs DepMap hotspots / copy ratio >= 3): expression matching follows lesions only where they leave a transcriptional imprint. Fraction of top-10 lines carrying the lesion, carriers vs non-carriers in the same group: THCA BRAF 0.34 vs 0.15, PAAD KRAS 0.76 vs 0.47, Her2 ERBB2 amp 0.52 vs 0.38, LUSC NFE2L2 0.27 vs 0.16, BLCA FGFR3 0.14 vs 0.07, COAD BRAF 0.25 vs 0.13; no difference for LUAD/COAD KRAS, LumA PIK3CA, most amplifications; IDH1 carriers' lines are 10% IDH1-mutant (= non-carriers). Matched pairs need an explicit lesion filter before use as evidence.
  - Normal-tissue negative controls (`scripts/72_celligner_normal_control.py`, `CONTROL=adjacent|gtex|none`): normals projected with all tumors as one batch (Celligner centers each batch on its own mean) through the main run's cPCs, DE genes and MNN, against the main run's cPC-regressed lines (the package's compute_cPCs=False path expects makeNewReference() output; with the raw lines tumors kept their nearest line only 51% of the time, so that first run was discarded). Question: is a tumor's closeness to a line cancer-specific, or would normal tissue of the same organ match the same lines as closely?
    - Reproducibility: tumors alone re-projected = main run (d1 rho 0.999, same nearest line 94%, nearest in main top 3 99.8%). Adding normals moves tumors more (adjacent 6.7% of batch: rho 0.98, same nearest 83%; GTEx 9%: 0.97, 78%): top-1 line identity depends on batch composition; use top-k / weights, not the single nearest line.
    - Adjacent normals (18 groups with >= 8 same-patient pairs; TCGA sequenced adjacent normals for ~10% of patients, none for heme, pediatric, OV, GBM/LGG, TGCT, UVM, MESO, ACC, UCS; SKCM mostly metastases): normal farther from its nearest line than its own tumor in 73% of 683 pairs (p 5e-33); tumor and own normal share the nearest line in 13%. Clear separation (AUC >= 0.83): LUAD, THCA, LUSC, HNSC, KIRC, KIRP, COAD, LIHC, CHOL. Weak (0.56-0.72): BRCA Basal/LumB/Her2, STAD, PRAD, ESCA-AD, BLCA. Reversed: LumA (0.36; normal breast closer to lines than LumA tumors), KICH (0.08).
    - GTEx (40 groups, <= 50 donors per group; unpaired): median AUC 0.88. Reversed (normal tissue closer than tumors): KICH, UCEC-EEC/SER and UCS (GTEx uterus = myometrium, lands on sarcoma lines), SARC-LMS/UPS (muscle/adipose), STAD (0.21). Weak: DDLPS 0.57, LumA/LumB 0.65, Her2 0.71, ACC 0.72.
    - Tissue-driven matches: normal tissue hits the same lines as the group's tumors (share of normals whose nearest line is among the tumors' 5 most frequent): LIHC 100% (Hep G2 for every normal liver), THCA 74%, PRAD 72% (adjacent: same nearest line 29%, top-10 overlap 0.49), ACC 68%, DDLPS 64%, STAD 62%, LGG 44-60% (normal brain -> KNS-42). In these groups line choice reflects organ identity; normal tissue is still much farther from the lines than tumors in LIHC and THCA (distance separates, line identity does not); in PRAD only weakly (adjacent AUC 0.63, GTEx 0.84).
    - Decision (user, 2026-10-08): matching is for drug applicability to the disease, not for separating patients within a disease. Tissue-level / organ-driven matches are acceptable (a drug perturbing liver lines applies to liver cancer); patients are not filtered on the normal baseline. What still matters: groups whose tumors do not match their own tissue's lines (KICH, THYM, TARGET-WT, PCPG, CESC-AD; LGG far from all lines), single-line hubs (KNS-42, Hep G2, SK-N-MM, NCI-H295R), and LINCS coverage of the matched lines.
  - Tissue keys (`config/celligner_tissue_keys.tsv`, `scripts/73_celligner_tissue_keys.py`; user 2026-10-08): each patient matched only to lines of the group's key. Default organ (DepMap OncotreeLineage = OncoTree tissue of the group); substitutes where the organ has no usable lines: THYM -> T-cell lines (no thymus lines), PCPG -> neuroblastoma (adrenal lines are adrenocortical carcinoma only); hematologic groups (LAML, TARGET-AML, ALL-B, DLBC) -> Blood (myeloid + lymphoid, 259 lines). Outputs `keyed_top_lines.tsv`, `keyed_group_summary.tsv`, `keyed_line_coverage.tsv`, `data/processed/celligner/weights_keyed.tsv.gz` (same tau 7.8).
    - Cost of forcing the key (median extra distance to the nearest in-key line vs unrestricted): WT +29.7, CESC-AD +24.6, UCS +11.1, CESC-SCC +10.7, LUSC +7.3, ESCA-SCC +7.2, CHOL +6.9, ESCA-AD +5.5, KICH +4.5, STAD +4.2; 0 to +1.2 elsewhere. WT within kidney: nearest line Wilms (CCLF_PEDS_0048_T) for 82 / 121 patients, rhabdoid kidney lines for 39.
    - Blood key: AML patients stay on myeloid lines (99.6-99.9% of top 10), but 30% of ALL-B and 36% of DLBC patients have a myeloid line as nearest (top 10: 18% / 8% myeloid). An umbrella Blood key would let myeloid-line drug evidence count for B-cell disease.
    - Single-line keys: ACC (2 lines; NCI-H295R weight 0.90), PRAD (7), Testis (10); hubs persist within key (KNS-42 0.94-0.98 for GBM/LGG, Hep G2 0.77, SK-N-MM 0.64-0.88).
    - Hybrid mode (`MODE=hybrid`, user 2026-10-08): organ-key lines plus, per patient, any line closer than the patient's nearest organ line. Weights are distance-based, so closer non-organ lines dominate where they exist: organ share of the disease weight 3% WT (rhabdomyosarcoma lines), 8% CESC-AD, 23% CESC-SCC, 37% UCS, 45-49% ESCA-SCC / KICH / LUSC, 52-64% CHOL, ESCA-AD, STAD, BLCA; >= 73% elsewhere, 100% for blood, NBL, UVM, READ, PCPG-PGL, LGG-codel.
  - Matrices (`scripts/74_celligner_matrices.py`, `celligner` env): `disease_line_weights.tsv` (49 groups x 1,588 lines, rows sum to 1; mean of patients' hybrid weights, i.e. share of the disease's patients resembling each line) + long form with `in_key`, patient counts, cumulative weight; lines covering 80% of a disease's weight: median 20 (1 for GBM/LGG, ACC, PHEO; 111 CESC-AD). Line x line Euclidean distances in a seeded 70-PC space of the aligned matrix (`data/processed/celligner/pcs_70.npy`, `dist_line_line.tsv.gz`; package PCA was unseeded: tumor-line distances r 0.99987, same nearest 93% vs main run), `line_neighbors.tsv` (10 nearest lines per line), `line_distance_summary.tsv`. Median nearest-neighbor line distance 48.9 (tumor -> nearest line 55.5); within lineage 102, across 143; nearest neighbor same lineage 65% (lymphoid 96%, myeloid 89%; thyroid 18%, biliary 19%, esophagus/stomach 23%). Basis for extrapolating LINCS drug data to unprofiled lines; assumes baseline-expression neighbors respond alike (test with LINCS leave-one-line-out).
- LINCS 2020 metadata (2026-10-08; `scripts/07_download_lincs.sh`, public S3 build, cellinfo / siginfo / compoundinfo / geneinfo; signature matrices not pulled). `scripts/75_lincs_coverage.py` -> `results/lincs/`, `data/processed/lincs/`.
  - 240 LINCS cell types, 230 with compound signatures; 34,418 compounds with signatures, 12,206 with a high-quality one (is_hiq). Lines with >= 100 high-quality compounds: 79, >= 1,000: 12 (PC-3, A-375, MCF7, VCaP, HA1E, U-2 OS, A549, HT-29, HCC515, Hep G2, MCF 10A, HeLa); HA1E and MCF 10A are non-cancerous. 143 of 240 map to DepMap (CCLE name 137, RRID 5, name 1); unmapped are mostly primary / normal / stem cells and XC.* patient-derived samples, plus DLBCL TMD8 and HBL1 (not in DepMap 24Q4). Blood lines with >= 100 compounds: THP-1 790, Jurkat 433, U-937, SKM-1, PL-21, K-562, NOMO-1, OCI-LY3/LY19/Ly10, Mino, BJAB, NALM-6, KMS-34 (~100-200 each).
  - Disease weight on lines with >= 100 high-quality compounds (Celligner hybrid weights): LIHC 0.77 (Hep G2), PRAD 0.17 (VCaP), LumB 0.11, every other disease <= 0.09; 0 for GBM/LGG, sarcomas, PCPG, NBL, UVM, THCA, SEM. No compound covers >= 25% of any disease's weight except in LIHC. Weighted mean distance from a disease's unprofiled lines to the nearest profiled line: 50-75 for blood, SKCM, COAD/READ, PAAD, LUAD, THYM, ESCA-AD, BRCA; 80-100 for most carcinomas; 117-156 for KIRC, THCA, NSGCT/SEM, PRAD, NBL, PCPG, LGG/GBM, ACC (line nearest-neighbor median 49, within-lineage median 102). Extrapolation from profiled neighbors is therefore required for nearly every disease and is long-range for brain, neural crest, germ cell, thyroid, adrenal, kidney.
- LINCS level 5 compound signatures (2026-10-09): `scripts/08_download_lincs_level5.sh` (35.5 GB GCTX, MODZ, 720,216 x 12,328; size-checked, sha256 recorded) -> `scripts/76_lincs_level5_subset.py` (`celligner` env, h5py; 6 h, mostly slow early disk reads) -> `data/processed/lincs/level5_landmark.npy` (554,650 signatures from the 135 LINCS lines mapped to cancer DepMap models x 978 landmark genes, float32, 2.2 GB, no NaN) + `level5_landmark_sigs.tsv.gz` (drug, line, ModelID, dose, time, is_hiq, TAS) + genes. 90,988 high-quality signatures, 33,603 compounds.
  - Cross-line agreement (high-quality signatures, median Pearson over landmark genes): same drug / different line 0.13 vs random pairs 0.03 overall. By signature strength (TAS; median 0.28): TAS < 0.2: same drug 0.07, same line / different drug 0.06, random 0.02; 0.2-0.4: 0.18 / 0.07 / 0.03; 0.4-0.6: 0.37 / 0.18 / 0.09; >= 0.6: 0.52 / 0.33 / 0.21. Drug effects carry across lines only for active signatures, and strong signatures share a generic response (random strong pairs r 0.21). Extrapolation should use active signatures (TAS >= ~0.4) and account for the shared generic component; to be tested leave-one-line-out.
  - Direct coverage, 41 signature groups (`results/lincs/direct_drugs_by_group.tsv`; active = high-quality and TAS >= 0.4, measured in a LINCS line of the group's own strict tissue key, no neighbor borrowing): 2,898 distinct drugs over the 41 groups. PRAD 1,587 (PC-3 1,201, VCaP 601; PC-3 has ~0 Celligner weight for PRAD), BRCA x4 1,295 (MCF7), SKCM 1,099 (A-375), LUAD/LUSC 1,092 (A549, HCC515), COAD/READ 933 (HT-29), LIHC 609 (Hep G2), CESC x2 361 (HeLa), PAAD 284 (YAPC), STAD/ESCA 35, NBL/PCPG 21, OV 12, UCEC/UCS 9, GBM/LGG and BLCA 1; 0 for kidney (KIRC, KIRP, KICH, WT), THCA, TGCT x2, ACC, SARC x3, HNSC (12 groups; kidney / HNSC / brain / bladder lines exist in LINCS but their signatures are weak or not high-quality).
  - Extrapolation test (`scripts/77_lincs_extrapolation_test.py`, user 2026-10-09: active signatures only, drug-specific effect, leave-one-line-out). 29,754 active signatures -> 3,056 drugs x 72 lines (mean over doses / times per drug x line). Generic response = per-line mean active profile, projected out: random-pair r 0.127 -> -0.007, same drug / different line 0.435 -> 0.330 (gap kept, 0.31 -> 0.34). Leave-one-line-out on 23 lines with >= 50 active drugs (7,687 drug x line pairs), predicting the held-out specific profile: mean over all other lines with the drug r 0.485 (80% > 0.3); distance-weighted tau 50 0.483, tau 25 0.467, tau 10 0.425; nearest line only 0.360; random other drug in the nearest line -0.011. Ceiling (two single signatures, same drug and line): 0.553. Accuracy of the nearest line barely depends on its distance (r 0.39 at < 75, 0.35 at 75-150, 0.29 at 150-200).
  - Implication: at landmark-gene level a drug's specific effect is mostly shared across lines; the cross-line consensus recovers most of the reproducible signal (0.49 of a 0.55 ceiling), and Celligner expression distance adds nothing beyond it. Caveats: held-out lines are the core lines (mostly different tissues; few close neighbors exist), averaging also reduces noise, doses / times pooled, landmark genes only. Line-specific responses (e.g. target present vs absent) are not predicted by expression proximity in this test.
  - Leave-one-line-out accuracy of the cross-line mean rises with the number of lines the drug is active in: 1 line r 0.30 (50% > 0.3), 2: 0.38, 3-4: 0.46, 5-8: 0.53, 9+: 0.57 (88-90% > 0.3).
  - Drugs per tissue (`scripts/78_lincs_tissue_drug_counts.py` -> `results/lincs/tissue_drug_counts.tsv`, `disease_drug_counts.tsv`; DepMap lineages of LINCS-mapped cancer lines; counts from siginfo). 28 DepMap lineages with cancer lines (+ 'Other'); 17 have LINCS compound signatures, 15 have >= 1 active drug, 11 have >= 100 (Prostate 1,587, Breast 1,295, Skin 1,099, Lung 1,092, Bowel 933, Liver 609, Cervix 361, Bone 358, Myeloid 305, Pancreas 284, Lymphoid 187). None at all for Biliary, Adrenal, Eye, Ampulla, Head and Neck (1 line, 119 compounds, none high-quality), Pleura, Soft Tissue, Testis, Thyroid, Vulva/Vagina; Kidney (3 lines, 725 tested) and CNS / Bladder have none or 1 active. 3,058 drugs active in some cancer line; 1,144 in >= 3 lines, 694 in >= 5. With borrowing every tissue / disease reaches 3,058 drugs (>= 1 line), 1,144-1,871 (>= 3 lines) or 694-1,671 (>= 5 lines); tissues without direct data get exactly the borrowed set (1,144 at >= 3 lines). Non-cancerous LINCS lines not counted: HA1E (kidney, 877 active drugs), MCF 10A (breast, 393).
  - Tissue dependence of drug response (ad hoc check, 2,418 drug x held-out line comparisons, 563 drugs): one same-tissue line predicts the held-out line slightly better than one other-tissue line (median r 0.35 vs 0.32; wins 57%). Carcinomas barely differ (lung 0.40 vs 0.36, breast 0.36 vs 0.36, prostate 0.34 vs 0.33, bowel reversed); larger gaps for lymphoid (0.28 vs 0.21), myeloid (0.35 vs 0.28), skin (0.38 vs 0.28), liver (0.41 vs 0.31, few drugs). Tissue-tiered borrowing proposed, not built (user chose target gate + pharmacology first).
- Drug annotations and target gate (user 2026-10-09). `scripts/09_download_drug_annotations.sh` (Repurposing Hub 2020-03-24, md5-verified) + `scripts/79a_drug_annotations.py` (ChEMBL 37 REST, raw JSON cached) -> `results/lincs/drug_annotations.tsv` for the 3,058 drugs active in a cancer line: Hub match 1,115 (419 launched), ChEMBL 1,866 (416 with ATC codes, 404 max phase 4); targets annotated for 1,404. Locally acting flag (topical-only route, or every ATC code in a local group: A01, A06, A07, D except oral D01B / D05B / D10B, G01, R01, R02, S01-S03, V04, V08, V09): 50 drugs (e.g. bisacodyl, tegaserod, topical antifungals / steroids / antiseptics). Most LINCS compounds are tool compounds without route / ATC, so the pharmacology check covers ~400 approved drugs.
  - `scripts/79b_target_gate.py` -> `results/lincs/disease_drugs.tsv` (disease x drug; direct or borrowed from >= 3 other lines), `disease_drugs_summary.tsv`, `data/processed/lincs/target_expression.tsv.gz`. Target symbols -> v23 genes via HGNC current / previous / alias (99.9% mapped); expressed = median TPM >= 1 in the group's tumors. Per disease 1,144-1,871 drugs: 575-818 with an expressed target, 43-108 removed (all targets not expressed; 4-9%), 466-945 without target annotation (pass, flagged); 1,041-1,763 pass, 141-227 of them approved (ChEMBL phase 4), 12-22 flagged locally acting. Removal is tissue-specific (most in KICH 103, PRAD 108, SKCM 116, ACC 98) but crude: removals include metabolizing enzymes listed as targets (simeprevir: CYP1A2 / CYP3A4), isoform-level annotations of pan-isoform drugs (thapsigargin: ATP2A1, muscle SERCA), and functional low-mRNA targets (BIBR-1532: TERT). Refinements pending.
  - Gate rules made config-driven (user 2026-10-09: reproducible, consistent): `config/drug_gate.yaml`, gene classes from HGNC gene groups of the pinned HGNC snapshot; each disease x drug row records `gate_reason`; run record `results/lincs/disease_drugs_run.json` (config + HGNC sha256, thresholds). Rules: expressed = TPM >= 1 in >= 20% of the disease's tumors (subset-level; replaces median TPM >= 1); ADME genes never used as targets (CYP families 1-3, UGT, cytosolic SULT, arylamine NAT, FMO, GSTs, carboxylesterases, SLC22, SLCO + ABCB1/ABCC1-3/ABCG2: 133 genes, 26 occur as drug targets; CYP17A1 / CYP19A1 / CYP11 stay targets); family rescue through HGNC groups of 2-15 genes. Rerun is byte-identical.
    - Per disease now 4-44 removed (was 43-108), 1,119-1,827 pass (160-245 approved). simeprevir -> adme_only, thapsigargin -> target_expressed (20% rule), BIBR-1532 -> family rescue. Removed in KIRC include lomitapide (MTTP, liver / gut), GABA-A, cannabinoid, CCK and tachykinin receptor ligands, an ALK inhibitor.
    - Family rescue is too permissive (144 drugs rescued somewhere): HGNC groups mix receptor subtypes with different pharmacology (Adrenoceptors: beta1-selective metoprolol rescued by alpha receptors) and complexes (Telomerase holoenzyme: TERT rescued by DKC1 etc.). Options: restrict to close paralogs (Ensembl Compara, identity cutoff) or disable.
    - Family rescue disabled (user 2026-10-09: drugs should be specific to their annotated receptor / target, limit off-target binding; `family.enabled: false`). Per disease 17-76 removed, 1,093-1,795 pass (157-243 approved, 14-24 flagged locally acting); BIBR-1532 (TERT) and lomitapide (MTTP) removed in KIRC.
  - Not done: distance gap to LINCS core lines (no LINCS data yet); comparison with the published DepMap Celligner alignment.
- Quality edge check (46_edge_confounding.R now also adjusts for the RNA-quality index, alone and with ESTIMATE + plate): edges losing > 50% of strength after quality alone 0.2-1.1% (SARC 2.3%, UCS 7.3%); after ESTIMATE + quality + plate 0.6-3.1% (SARC 5.7%, UCS 8.7%). Networks kept as built; UCS (57 patients, single plate) is the least robust network.
- Regulator modules rebuilt on the quality-corrected scores (58_regulator_modules.R, NES_DIR default aracne_project_qc -> `diagnostics/regulator_modules_qc.tsv`; previous version renamed `regulator_modules_purity.tsv`). Regulators in quality-tracking modules 29% -> 0%; in plate-tracking modules (R^2 > 0.5) 7.6% -> 0.8%; largest module BRCA-LumA 927 -> 141, PRAD 1,120 -> 206. Modules now read as biology: BRCA fibroblast/CAF (DPYSL3, PRRX1, SNAI2, COL1A1), proliferation (CCNA2, FOXM1), T cells (LCK, CD3E), endothelial (ERG, SOX7), interferon (STAT1); COAD WNT/stem (ASCL2, PLAGL2), fibroblast (TWIST1, PDGFRB), secretory (SPDEF); THCA largest module MET/TACSTD2/PTPRE (BRAF-like axis); TGCT-NSGCT GDF3/FGF4 (embryonal histology). Residual: a recurring PILRB/HSF4/CCNL2/AKAP8L module in 27 groups (plate R^2 median 0.24); possible remaining technical axis, cause not identified.
- Recurring module investigation (`scripts/62_recurring_module.R` -> `diagnostics/recurring_module_*.tsv`; Hallmarks from MSigDB 2025.1 in data/raw/msigdb/; transcript GC computed from GENCODE v23 transcripts -> data/interim/gene_gc_v23.tsv). Found in 39 of 41 groups; 74 core members incl. paralog family AGAP4/5/6/9, readthrough JMJD7-PLA2G4B, duplicated CDK11A, CCNL1/2, many ZNFs; members GC-rich (median 0.52 vs 0.48) and enriched on chr16/17/19/22 (38% vs 26%); their ARACNe targets GC-rich in 25 of 28 networks. Not immune (rho with cytolytic activity -0.01, ESTIMATE immune -0.02, interferon +0.02-0.04). Apoptosis -0.38, p53 -0.25, DNA repair -0.22, but 36 of 50 Hallmarks are negative (none positive): a global, non-specific axis, not a pathway. Plate R^2 0.23 vs chance 0.06 (above chance in all 39 groups); GC index +0.10; RNA-quality index ~0. No Hallmark enriched among targets (min FDR 0.05). Interpretation: most consistent with a technical axis (batch-dependent GC bias and/or multi-mapping paralog quantification), not verified.
- CQN + readthrough cleanup (2026-10-08/09; `61_quality_correct.R` NORM=cqn -> `data/processed/step6/aracne_project_cqn/`, checks in `quality_correction_cqn/`; cqn 1.56.0 added to the env; 384 protein-coding readthrough loci dropped, list in data/interim/readthrough_genes_v23.tsv). Result: the recurring module is NOT removed (38 of 41 groups, 66 core members, still anti-correlated with 34 of 50 Hallmarks, plate R^2 0.27 vs chance 0.06) although its GC link dropped (GC index +0.10 -> +0.02): GC bias is not the cause. Top members still the AGAP4/5/6/9 paralog family, PILRB, ZFC3H1, CCNL2, DMTF1. CQN otherwise mixed: plate improved in STAD (0.14 -> 0.04) and BRCA, worsened in PRAD (0.06 -> 0.13) and UCEC-SER; patient concordance with group MRs dropped in GBM (0.52 -> 0.32), SARC-UPS, KICH, SARC-DDLPS. TMM version (aracne_project_qc) stays the Step 8 input pending PI decision.
- Mappability filter (2026-10-09; Saha & Battle 2018, as used for GTEx v8 trans-eQTLs; GRCh38/GENCODE v26 resources from figshare collection 4297352 in data/raw/crossmap/, md5 verified; `61_quality_correct.R` MAPPABILITY=1 -> `data/processed/step6/aracne_project_map/`, `quality_correction_map/`): genes with mappability < 0.8 removed (expression, regulators, targets) and edges between cross-mappable genes removed (12.7-14.2% of edges per network). AGAP4/5/6/9 have mappability 0.13-0.14 (true multi-mapping artifact) and leave the module; other core members are normally mappable (median 0.99). The recurring module persists (39 of 41 groups; median size 56 -> 33; 55 core members: CCNL2, PILRB, ZNF692, ANKZF1, HSF4, WDR90, PLA2G4B, MAPK8IP3, AKAP8L, GIGYF1; still anti-correlated with 35 of 50 Hallmarks; plate R^2 0.23 vs chance 0.06). Mappability filter costs no biology (patient concordance change median -0.006, range -0.056 to +0.04). Tested and rejected as the module's cause: RNA quality, GC bias (CQN), multi-mapping (mappability filter). Cause unresolved; plate-linked.
- 2026-10-09 decisions (PI): mappability filter adopted as standard (61_quality_correct.R default; Step 8 input = `data/processed/step6/aracne_project_map/`; 58/62 default to it). Recurring module kept in clustering; PI interpretation: growth / mitochondrial / immune-rescue regulators altered in cancer to support proliferation. Data note: in 39 groups the module's score is not associated with proliferation (E2F targets +0.02, G2M +0.04) or immune activity (~0) and is negatively associated with oxidative phosphorylation (-0.50), ROS (-0.42), mTORC1 (-0.40), MYC targets (-0.32).
- Step 8 (`scripts/80_consensus_clustering.R`, `81_cluster_characterize.R` -> `results/step8/`; PI: each of the 6 groups separately; main features = one score per module of >= 10 regulators, sensitivity = all modules + single regulators). M3C (PAM, 100 resamples, 100 Monte Carlo references, K 2-10); rule fixed before running: best K = max RCSI, clusters only if Monte Carlo and beta p < 0.05.
  - Main: BRCA-Basal, Her2, LumB no significant structure. BRCA-LumA K = 3 borderline (p 0.050 / 0.036, silhouette 0.07), associated with plate (p 5e-4), not significant with the sensitivity features, ARI 0.33 between feature sets -> not robust. TCGA-COAD K = 2 (RCSI 0.65, p 0.020 / 0.006, silhouette 0.19, no plate/purity/quality association), ARI 0.75 with the sensitivity clustering (whose own p = 0.12); cluster 2 (n = 92) holds 43 of 48 hypermutated-indel (MSI) tumors -> rediscovers the known MSI vs CIN split. TCGA-READ K = 3 (p 0.020 / 0.018) and significant in sensitivity too, but ARI 0.41 and plate-associated (p 0.031; 0.0015 sensitivity), not explained by PanCanAtlas subtypes -> tentative.
  - Sensitivity only: BRCA-Basal K = 2 significant (p 0.010 / 0.0007) but ARI 0.11-0.19 with the main run -> not robust.
- Step 9a colorectal replication processing (`scripts/90_process_replication_crc.R` -> `data/processed/step9/{gse39582,cptac_coad}/`, `results/step9/`; PI 2026-10-09: reference = each cohort's average tumor, viper method "scale"). GSE39582: 566 tumors (single-gene GPL570 probes, Entrez -> v23, MaxMean collapse -> 16,212 genes; ESTIMATE affymetrix; purity-corrected). CPTAC-COAD: 105 colon tumors (STAR counts v36 -> v23 protein-coding, filterByExpr + TMM, 13,131 genes after the mappability filter; ESTIMATE illumina; purity + RNA-quality corrected). TCGA COADREAD regulons with the mappability filter and the TCGA-COAD regulator set: 3,839 / 3,835 regulators scored (98%). Sanity: ESTIMATE immune higher in dMMR (959 vs 443, p 2e-9) and MSI-H (375 vs -36, p 3e-3), as expected. READ clusters cannot be replicated with these cohorts (colon only; CPTAC's single rectal tumor dropped).
- Step 9b colorectal replication (`scripts/91_replicate_crc.R` -> `results/step9/crc_*.tsv`; PI 2026-10-12: all three methods + CMS; clusterRepro 0.9 (CRAN) and CMScaller 2.0.1 (GitHub Lothelab/CMScaller @f1541730) installed into the env). All cohorts scored on the 78 TCGA-COAD module features (mean of sign-aligned z-scored member NES; r 0.998 with the Step 8 PC1 scores).
  - Transfer (nearest TCGA centroid, Spearman): the MSI-like cluster captures 73 of 75 dMMR tumors in GSE39582 (OR 60, p 5e-25) and 24 of 24 MSI-H in CPTAC (p 3e-11); TCGA reference 43 of 48. MSI rate in the MSI-like cluster 30-52% vs 0-3% in the other. In-group proportion 0.92-0.95 (TCGA, GSE39582 p < 0.001; CPTAC p 0.035 / 0.053).
  - Independent M3C: GSE39582 best K = 3, not significant (p 0.099 / 0.059); CPTAC best K = 4, significant (p 0.030 / 0.010). Their K = 2 splits agree with the transferred labels (ARI 0.59 / 0.45) and are strongly MSI-associated (p 2e-23 / 1e-11). The split is a strong axis rather than a sharply separated 2-group structure.
  - Defining modules: MSI-vs-MSS effect sizes across the 78 modules correlate with the TCGA cluster effects (rho 0.89 GSE39582, 0.94 CPTAC); TCGA top-10 modules same direction 10/10 in both.
  - CMS (CMScaller): MSI-like cluster = CMS1 + CMS3 (TCGA 45 + 34 of 92; GSE39582 83 + 80; CPTAC 18 + 16); other cluster = CMS2 + CMS4 (TCGA 82 + 77). So the split is CMS1/3 vs CMS2/4, broader than MSI. TCGA-READ clusters vs CMS: cluster 2 mostly CMS4/CMS2, cluster 3 CMS3-enriched, cluster 1 mixed (READ itself cannot be replicated in these colon cohorts).
  - Bug caught before reporting: a data.table scoping error scrambled TCGA labels in the first run (reference row showed no MSI association); fixed with an explicit lookup + check.

- 2026-10-12 decisions (user): colon groups named **COAD-CMS1/3-like** (Step 8 cluster 2, MSI-like, n = 92) and **COAD-CMS2/4-like** (cluster 1, n = 195); this is the colon stratification for Step 10 signatures and drug matching. The three TCGA-READ clusters (25 / 40 / 27) are kept as separate subgroups (READ-C1/C2/C3), each matched to drugs independently, flagged tentative. Breast: continuous axes, choice of method delegated to Claude.
- Breast continuous axes (`scripts/82_breast_axes.R` -> `results/step8/breast_axes_{pca,annotation,loadings,scores}.tsv`). Method: principal components of the Step 8 modules10 features per subtype (few independent axes, no hand-picked biology, fixed loadings can be projected onto METABRIC / CPTAC-BRCA). Kept if (1) variance above the 95th percentile of 500 column-permuted datasets (parallel analysis), max 5; (2) bootstrap-stable (median |cosine| of loadings over 200 patient bootstraps >= 0.8); (3) not technical: plate R^2 excess over chance, (R^2 - chance) / (1 - chance), < 0.2, and |rho| < 0.5 with RNA-quality index and ESTIMATE. The first run used a fixed raw plate R^2 < 0.25, which ignores that chance R^2 is 0.07-0.31 depending on group size (it dropped Basal PC1, R^2 0.33 vs chance 0.19); replaced before use.
  - Kept 13 axes. Basal PC1 (23%; stroma/EMT low vs immune + proliferation high), PC3 (6%; myogenesis/apical junction, FOXC1/SOX8), PC4 (5%; complement/glycolysis, myeloid NLRC4/SPI1 vs ribosomal). Her2 PC2 (15%; glycolysis/hypoxia/EMT, FAT1/CAV1), PC3 (13%; mitotic spindle/immune vs ESR1/XBP1), PC4 (6%; interferon/STAT1 vs MYC targets). LumA PC2 (14%; stroma/EMT low vs immune IL10RA/CD74), PC3 (12%; proliferation vs ESR1/GATA3), PC4 (6%; EMT/fibroblast DPYSL3), PC5 (5%; glycolysis, adipose PPARG/FABP4 low). LumB PC1 (16%; EMT/glycolysis/hypoxia), PC3 (11%; ESR1/GATA3 + low immune/proliferation), PC4 (8%; angiogenesis/endothelial).
  - Dropped as technical (plate excess >= 0.2): Basal PC2 (0.24; recurring module rho -0.62), Her2 PC1 (0.204, borderline; interferon/proliferation), LumA PC1 (0.38; anti-MYC/mTORC1/OXPHOS, recurring module +0.43), LumB PC2 (0.22). Unstable: PC5 of Basal, Her2, LumB.
  - Recurring module still loads on some kept axes (LumA PC2 +0.59, Her2 PC3 +0.57, LumB PC1 -0.62); kept per the 2026-10-09 decision.
- Step 10 signature reference (user 2026-10-12): each subgroup (COAD-CMS1/3-like, COAD-CMS2/4-like, READ-C1/C2/C3) is compared with healthy GTEx tissue of the organ, as in Step 4. Breast axes (user 2026-10-12, option C + A): patients split into thirds by axis score; top third ("high", positive PC score) and bottom third ("low") each vs GTEx breast, as Step 4 -> 26 signatures used for drug matching; companion per axis = genes that follow the axis across all of the subtype's tumors (explanation only).
- Step 10 signatures (`scripts/92_step10_signatures.R` -> `results/step10/`, full tables `data/processed/step10/`): 31 signatures (5 COAD/READ subgroups + 13 axes x 2 ends), Step 4 model (~ condition + age + ESTIMATE (+ sex)), protein-coding, filterByExpr/TMM/voom, moderated t -> z. Companion: limma-voom ~ axis score + age + ESTIMATE across the subtype's tumors (`z_axis_drivers.tsv`, top 25 genes each way in `axis_driver_top.tsv`).
  - Every stratum signature is close to its parent group's Step 4 signature (Spearman 0.967 COAD-CMS1/3-like, 0.977-0.992 all others); the two ends of a breast axis agree at 0.93-0.97. The shared cancer-vs-normal signal dominates; the stratum-specific part is a thin layer (expected, option C trade-off). The high-minus-low difference follows the companion axis genes (rho 0.82-0.93), so that layer is the axis biology.
  - Companion genes read as biology: Basal PC1 T cells/IFNG (CD8A, CXCL10, TIGIT) vs fibroblast collagen (COL1A1, PDGFRB); LumA PC2 immune (IL10RA, DOCK2) vs myofibroblast (ACTA2, TAGLN); LumA PC3 proliferation (AURKA, MYBL2) vs ciliary genes; LumA PC4 cancer-associated fibroblast (COL11A1, COL10A1, MMP11) vs endothelial/lymphocyte GIMAPs; LumB PC4 endothelial (CDH5, ERG, CD34); LumB PC3 stroma vs T cells/checkpoints (LAG3, CTLA4). Low ends of Basal PC3 and Her2 PC4 are mitochondrial/ribosomal housekeeping genes (possible residual technical component).
  - Implication for drug matching: drug rankings of subgroups / axis ends will largely follow their parent group; differences between strata should be read from the stratum-specific layer (e.g. ranking differences between ends, or companion lists), to be decided with the drug side.
- Step 10 strata (`scripts/83_step10_strata.R` -> `results/step8/step10_strata.tsv`, `step10_patient_strata.tsv`): 2 COAD + 3 READ subgroups and 13 breast axes. Breast axes are TCGA-only until replicated.

- Step 9c/9d breast axis replication (user 2026-10-12). Processing `scripts/93_process_replication_brca.R` -> `data/processed/step9/{metabric,cptac_brca}/<group>/`: METABRIC (Illumina HT-12, cBioPortal log2; Entrez -> v23, MaxMean, 16,983 genes; CLAUDIN_SUBTYPE Basal 209 / Her2 224 / LumA 700 / LumB 475) and CPTAC-BRCA (STAR counts, per-subtype filterByExpr + TMM, ~12,700-13,300 genes; PAM50 Basal 27 / Her2 14 / LumA 56 / LumB 17); ESTIMATE (+ RNA-quality index for CPTAC) gene-level correction within subtype; TCGA BRCA ARACNe regulons + mappability filter, each subtype's TCGA regulator set (94-95% scored); reference = cohort-subtype average tumor. Tests `scripts/94_replicate_brca_axes.R` -> `results/step9/brca_axes_{replication,features_check}.tsv` (rule fixed before running): (1) congruence of the TCGA loading vector with the cohort's own PCs 1-10 vs loadings shuffled across modules; (2) variance captured along the TCGA direction vs the same null; (3) genes following the projected axis vs the TCGA companion driver genes. Replicated = congruence >= 0.85 + p < 0.05 + gene rho >= 0.5; partial = p < 0.05 + gene rho >= 0.3. CPTAC subtypes < 30 patients flagged underpowered.
  - All 26 axis x cohort tests are above chance (congruence p <= 0.05; top-100 TCGA driver genes same direction 86-100%; gene rho 0.29-0.76). Replicated (strict): Basal PC1 (CPTAC 0.87), LumA PC3 (METABRIC 0.90), LumA PC4 (0.88), LumB PC4 (0.89); each is partial in the other cohort. Partial in both: Basal PC3, Basal PC4 (congruence 0.55), Her2 PC2, Her2 PC3, LumA PC2, LumA PC5 (0.66 / 0.58), LumB PC1, LumB PC3. Weak: Her2 PC4 (METABRIC not replicated, gene rho 0.29, congruence 0.54; CPTAC partial with congruence at the null, n = 14).
  - Reading: the axes are real and their gene programs carry over to an independent array cohort; most do not come out as a single clean PC in the new cohort (congruence 0.55-0.80; they mix with neighboring PCs of similar size), so they are best used as projected directions (fixed TCGA loadings), not re-derived per cohort.
- Drug ranking for strata (user 2026-10-12): score by reversal and specificity. Note for the drug-side session: `docs/notes/subgrouping_for_drug_side.md`.

- Step 10b stratum master regulators (`scripts/95_stratum_mr.R` -> `results/step10/stratum_mr.tsv` + `stratum_mr_summary`, `stratum_regulator_filter`, `stratum_patients`, `stratum_mr_parent_check`; user 2026-10-12). Method as 50_viper.R (aracne_project, the version behind group_mr.tsv): msviper on each Step 10 stratum signature, parent's network, 1,000 normal-split null from the parent's GTEx normals; regulator-expressed filter on the stratum's own tumors (user). Check: the same code on each parent's Step 4 signature reproduces group_mr.tsv (NES Spearman 0.9996-0.9997; FDR < 0.10 sets Jaccard 0.85-0.94, differences from the random null).
  - MRs at FDR < 0.10: COAD-CMS1/3-like 216 (0 at FDR < 0.05; parent 339), CMS2/4-like 347; READ-C1/C2/C3 292 / 292 / 231; breast strata 617-1,015. NES Spearman with the parent 0.91 (CMS1/3-like) and 0.97-0.99 elsewhere. MRs not called in the sibling stratum: 18 (READ-C3) to 311 (LumA PC2 high).
  - Regulators that most separate siblings read as expected biology: CMS2/4-like ASCL2, PLAGL2 (WNT/stem), EREG up vs CMS1/3-like; CMS1/3-like DUSP4 up; Basal PC1 low (stroma end) ITGA11, INHBA, ANTXR1 up and high (immune end) PDGFRB down; LumA PC3 high MYBL2, AURKB up.
- Running decisions/assumptions list started: `docs/DECISIONS.md` (user 2026-10-12: ongoing list; update with every new decision).

## 2026-10-18
- Retrospective: the RNA-quality and mappability corrections (Step 6e) reached only patient-level scores, never the Step 4 signatures, group master regulators or networks used by Step 10 and the drug side. Ad hoc check: the Step 4 tumor-vs-GTEx z carry a GC tilt (tumor-vs-GTEx minus tumor-vs-adjacent vs GC: rho -0.03 to -0.47 in 16 of 17 groups), and a GC-only signature reproduces group MR rankings (aREA, median rho 0.20; KIRC 0.66, TARGET-WT 0.63; 93% of KIRC FDR < 0.10 MRs have the GC-predicted sign).
- Decisions D45-D47 (DECISIONS.md): ESTIMATE everywhere; mappability < 0.8 + readthrough genes removed before Step 4 and the network rebuild; measure lab / post-mortem / quality first.
- Reference diagnostics (`scripts/31_reference_diagnostics.R` -> `results/step4/reference_diagnostics/`; no pipeline output changed):
  - Source effect, adjacent vs GTEx normals (both healthy; 17 groups with >= 10 of each; ~ source + age + sex + ESTIMATE): median 3,844 genes FDR < 0.05 and |logFC| > 1 (2,076 KICH - 6,148 COAD), i.e. as large as many tumor-vs-normal contrasts. GC tilt median -0.28 (KIRC -0.52, LumA -0.47); GC + length explain a median 10% of it (2-27%), so CQN alone would not remove it. Shared across tissues: median rho 0.61 between tissues, 0.41-0.86 with the leave-one-tissue-out consensus (groups sharing GTEx donors excluded).
  - Overlap with Step 4: source effect vs tumor-vs-GTEx rho 0.40-0.86 (median 0.62; PRAD 0.86, THCA 0.78), R^2 median 0.38. Split-half GTEx control (source vs half 1, tumor vs half 2) gives the same rho (0.42-0.83), so this is not shared-sample noise. The cross-tissue consensus correlates with the tumor-vs-GTEx signature of groups without adjacent normals too (median 0.39; -0.15 NBL to 0.62 Her2).
  - What the shared component is: higher in TCGA normals = OXPHOS, protein secretion, MYC targets, mTORC1, E2F (cameraPR FDR <= 0.008), histone HIST1H4C; higher in GTEx = CLASRP, CCNL2, LENG8, PILRB, HSF4, AGAP6, NPIPB5. These are the recurring-module genes: 89% of the 555 recurring-module members (Step 6, mappability version) are higher in GTEx (median z -3.9 vs +0.7 for other genes, p 3e-150); the 17 core members (>= 20 groups) all are, median rank in the top 0.6%. The recurring module (plate-linked within TCGA, anti-correlated with OXPHOS / mTORC1 / MYC) is the same axis as the TCGA-vs-GTEx lab difference. Cause still unverified (candidates: nuclear / intron-retained RNA capture, library protocol).
  - Post-mortem within GTEx (~ ischemic time + ventilator death + age + sex): ischemic time changes few genes in most tissues (median 68 FDR < 0.05) but thousands in lung (5,221), esophagus/HNSC (2,981 / 4,718), OV (2,194); where strong it points against the tumor signature and the source effect (lung rho -0.43 to -0.48 with tumor-vs-GTEx, -0.47 to -0.52 with source; COAD/READ -0.36 / -0.53 to -0.61). Ventilator (Hardy 0) donors look more like surgical TCGA normals (rho with source up to 0.59 COAD). Cross-tissue agreement of the ischemic effect is low (0.13).
  - RNA-quality index computed jointly (tumors + GTEx + adjacent): tumors lower than GTEx (median d -0.70) but TCGA adjacent normals higher than GTEx (median d +0.22; kidney +1.2 to +3.2, breast +0.6 to +1.4); index tracks GTEx RIN weakly (median rho 0.13). The tumor-vs-GTEx index difference is therefore not a lab/RNA-integrity difference (same-lab normals go the other way); adding the index to Step 4 would remove tumor biology.
  - Mappability: 1,511 of 19,168 Step 4 genes have mappability < 0.8, 344 readthrough loci to drop; 1,824 either.
- User (2026-10-18): the adjacent-vs-GTEx difference is not a lab difference alone (surgery vs post-mortem, cancer patients vs non-cancer donors, tumor-adjacent tissue state); D48 put under review. Goal restated: a reference that generalizes to new patients.
- Normal-tissue variance and generalization (`scripts/33_normal_variance.R` -> `results/step4/reference_diagnostics/normal_variance_*.tsv`; 11 organ units, equal n per arm, 30 resamples): TCGA adjacent normals are not consistently more variable than GTEx. Multivariate spread (top-20-PC distance to centroid, adjacent / GTEx): more variable in head-and-neck 1.72 (several sites vs one GTEx tissue), stomach 1.46 (tracks ESTIMATE, rho -0.75), esophagus 1.27 (n = 10); similar in breast, prostate, thyroid, uterus (1.05-1.09); less variable in colorectum 0.55, liver 0.72, kidney 0.79, lung 0.82 (GTEx colon / kidney contain very different sample subsets: per-gene variance ratios are bimodal across resamples, so the PC spread and BCV are the robust measures; BCV GTEx colon 0.71 vs adjacent 0.44, kidney 0.79 vs 0.48). Tumor-likeness (same-lab score: deviation from other adjacents along the tumor-vs-adjacent signature) median ~0 and higher in the more unusual adjacent samples in kidney (rho 0.67), esophagus and uterus (0.50), lung (0.41), thyroid (0.40). Plate above chance in breast (0.41 vs 0.15). Split-half signature reproducibility (patients split into disjoint halves, same tumor halves and normal n for both references): GTEx higher in 8 of 11 units (e.g. prostate 0.969 vs 0.931, head-and-neck 0.965 vs 0.928, stomach 0.949 vs 0.897), adjacent higher in kidney (0.920 vs 0.871) and colorectum (0.968 vs 0.956). Caveat: a consistent collection offset reproduces perfectly, so stability does not show freedom from bias. Half-signatures agree with the other reference's full signature at 0.33-0.80.
- GTEx collection correction test (`scripts/34_gtex_collection_correction.R` -> `results/step4/reference_diagnostics/gtex_collection_*.tsv`; per GTEx tissue, v2 gene set, one sample per donor):
  - Held-out variance (5-fold CV over donors; MSE vs M0 = age + sex + ESTIMATE): M1 (+ log ischemic time + death circumstance ventilator / fast / slow) lowers held-out variance in muscle, lung, liver, visceral adipose (~6%; 69-89% of genes), adipose, adrenal, esophagus mucosa, skin, thyroid (2-3%); no gain in breast, prostate, colon, stomach, pancreas, testis, ovary, uterus, vagina; worse in small sets (kidney 27 donors, brain BA9 49). Adding autolysis + fixative time (M2), RIN (M3) or isolation batch (M4) does not improve on M1 overall (median ratio 0.997 / 1.009 / 0.998 / 1.002).
  - Shifting all samples by a constant does not change variance, so variance cannot choose the target. Target test: hide a fifth of the best-preserved donors (ventilator, ischemic time in the tissue's shortest third, autolysis <= 1; 16-120 per tissue, 0 in brain cortex, 7 kidney), fit on the rest. Model-shifted reference (L2, M1) is closer to the hidden best-preserved donors than the average-donor reference (L1) in 20 of 20 tissues (gap closed median 0.13; muscle 0.45, adipose 0.36, lung 0.33, liver 0.25; stomach 0.01) and points the right way (rho 0.19-0.84). The plain average of the training best-preserved donors (L3) is closer still (median 0.18; lung 0.51, thyroid 0.43), except where few such donors exist (prostate, uterus, vagina, breast, testis, ovary: 25-53).
  - Check only: the L1 -> L2 shift moves GTEx toward TCGA adjacent normals in 7 of 9 organs (rho with the gap: colon 0.65, lung 0.36, uterus 0.33, prostate 0.28, liver 0.25) but closes only 0-10% of that gap; most of the TCGA-normal vs GTEx difference is not explained by GTEx-internal collection variation.
- GTEx sub-regions (colon, kidney): marker scores (mean log-CPM) are bimodal in GTEx Colon - Transverse (epithelial markers CDX1/2, VIL1, KRT20, CDH17, LGALS4, CEACAM5, EPCAM: a mucosa mode at 8-10 and a muscle-wall mode at -2 to 6; TCGA adjacent colon all >= 7) and Kidney - Cortex (proximal-tubule markers SLC34A1, LRP2, CUBN, SLC22A6/8, SLC5A12, ALDOB, GATM: 2 near-empty samples, 5 medulla-dominant samples below 5). Keeping colon >= 7 (165 -> 109) and kidney >= 5 (28 -> 21) brings GTEx variability to the TCGA adjacent level: colon BCV 0.67 -> 0.41 (adjacent 0.44), PC spread 92 -> 54 (56); kidney BCV 0.78 -> 0.47 (0.50), spread 87 -> 71 (70).
- Decisions (PI, 2026-10-18), D57-D60: GTEx stays the reference, corrected from within GTEx (log ischemic time + death circumstance, all donors kept, centered on the best-preserved donors: observed average where >= 45 such donors, model prediction where 10-44, average donor where < 10); no subtraction of the TCGA-adjacent vs GTEx gap (D48, D49 and the lab-shift part of D52 struck); TCGA adjacent normals = check only, secondary drug-ranking use under discussion; colon and kidney GTEx sub-region filter; kidney reference flagged for the PI. Spec rewritten: [`SPEC_rebuild_v2.md`](SPEC_rebuild_v2.md) (Step 4a = GTEx reference preparation; checks C4 and C11 replaced/added).
- Specs written (user 2026-10-18; not executed): tissue composition per cohort ([`SPEC_tissue_composition.md`](SPEC_tissue_composition.md), D61) and Celligner projection of all v2 GTEx reference donors per cohort (`celligner` branch `docs/notes/SPEC_celligner_normals.md`, D62). TCGA adjacent normals in drug ranking: the user picked a safety filter and the inflammation around tumors as a secondary target (sketch in the drug-side notes).
- v2 execution, step 1 (2026-10-18): `scripts/20_covariates.py` adds `gtex_hardy`, `gtex_death` (ventilator / fast / slow) and `gtex_autolysis` to `results/step3/covariates_toil.tsv` (GTEx rows only; all existing columns identical to before). Over 3,542 one-per-donor GTEx samples: ventilator 2,181, fast 869, slow 451, missing 41; autolysis 0 / 1 / 2 / 3 = 815 / 2,236 / 384 / 49, missing 58. `scripts/22_gene_set_v2.py` (runner step `gene_set_v2`) -> `results/step3/genes_protein_coding_v23_v2.tsv`: of 19,797 protein-coding genes, 2,013 dropped (1,664 mappability < 0.8, 384 readthrough, 35 both; list with reasons in `genes_v2_dropped.tsv`), 17,784 kept (211 without a mappability value, flagged). Rerun via the runner is byte-identical.
- v2 execution, step 2 (2026-10-18): network questions settled first (D63-D66): full regulator set plus a TF-only build; extra luminal-only breast network BRCA_LUM (LumA + LumB, 769 tumors); no proliferation regression; gene filter kept (removes 28 of 978 L1000 landmark genes and 412 candidate regulators, 100 TFs, no known-biology control). Inputs (`VERSION=v2 scripts/42_aracne_inputs.py` -> `data/processed/step5/aracne_inputs_v2/`, `results/step5/aracne_{cohorts,regulators}_v2.tsv`): 29 networks, 12,911-14,783 genes, 4,286-4,890 candidate regulators (TF-only lists 1,114-1,360). The 28 cohort inputs are byte-identical to the first v2 input build. ARACNe build started detached: full set (`VERSION=v2 scripts/43_run_aracne.sh`, log `logs/aracne_v2.log`), then TF-only (`REGSET=tf`, `logs/aracne_v2_tf.log`).
- v2 execution, Step 4a (2026-10-18): `scripts/32_gtex_reference.R` -> `results/step4_v2/gtex_reference.tsv` (sample level) and `gtex_reference_tissues.tsv`. 3,542 GTEx samples (one per donor per tissue), 3,425 used; excluded: 53 colon muscle-wall samples (epithelial score < 7), 6 kidney samples (proximal-tubule score < 5), 50 missing collection data, 8 absent from the Toil count matrix. Colon 166 -> 109, kidney 28 -> 21 (as in the exploration). Anchor classes: observed (>= 45 best-preserved) 10 tissues (subcutaneous adipose, esophagus mucosa, lung, muscle, pancreas, both skins, stomach, testis, thyroid), model (10-44) 10 (visceral adipose, adrenal, breast, colon, GE junction, liver, ovary, prostate, uterus, vagina), average (< 10, flagged) 6 (bladder, brain cortex, frontal cortex, both cervix tissues, kidney). Anchor ischemic times 80-222 min where anchored.
- Tissue composition (D61; `scripts/04b_download_composition_refs.sh` pins HPA v24 single-cell types and MCP-counter commit b6eac73 / package 1.2.0; `scripts/35_tissue_composition.R` -> `results/step4/composition/`): 40 groups, 56 scores per sample (10 MCP-counter populations, 46 HPA cell-type panels, G2M proliferation), tumors vs the v2 GTEx reference vs TCGA adjacent. Tumor vs GTEx: 1,958 of 2,240 cell-type tests at FDR < 0.05, 698 with |Cliff's delta| > 0.8. Composition share of the signature (out-of-fold ridge R^2 of the tumor-vs-GTEx z on HPA cell-type specificity): median 0.19 (0.11-0.38) with the ESTIMATE term, 0.20 without, so ESTIMATE removes little of it (it targets immune / stroma, not parenchymal cell types); highest in cervix 0.38 / 0.37, endometrium 0.34 / 0.33, ovary 0.28, BRCA-Basal 0.27, LUSC 0.27 (normal side: smooth muscle, ovarian stroma, adipocytes, alveolar type 1). Proliferating HPA types (undifferentiated cells, cytotrophoblasts) carry the tumor side, i.e. proliferation. Positive controls: 9 of 9 evaluable pass (breast adipocytes lower, Cliff -0.97 to -1.00; proximal tubule lower in KIRC -0.70 / KIRP -0.94; melanocytes higher in SKCM +0.57; T cells / B lineage higher in seminoma +0.997 / +0.97); 9 not evaluable because the fixed marker rule (>= 4x every other HPA type, >= 8 markers) gives no panel for cell types that share genes with close relatives (smooth muscle 1 marker, alveolar type 1 5, distal enterocytes 6, basal keratinocytes 0). Fix proposed, not applied: score related types as families (smooth muscle with fibroblasts excluded from the comparison set, proximal + distal enterocytes, basal + suprabasal keratinocytes, alveolar 1 + 2).
- Celligner projection of all v2 GTEx reference donors (D62, `celligner` branch `scripts/72b`): chunks of ~390 failed the pre-set drift check (repeated anchor samples keep their nearest line in 66-78% of chunks, need >= 90%; nearest-distance rho 0.87-0.97, need >= 0.98); results not used. Rerunning with chunks of 100 as the spec prescribes, with weight-level drift added.
- Composition rerun with cell-type families (D67): 51 of 56 HPA units scored (was 46 single types). Positive controls 11 of 18 pass, 2 fail, 5 not evaluable: newly passing LUAD alveolar cells lower (Cliff -0.32) and SKCM keratinocytes / squamous lower (-0.98); failing COAD / READ enterocytes + Paneth "similar" (tumors far lower, -0.96 / -0.95): the panel is differentiated-enterocyte genes, which colorectal tumors lose, so this measures dedifferentiation as well as cell numbers (assumption in the spec); not evaluable: smooth muscle (+ myoid) has 3 markers, shared with HPA Sertoli, skeletal muscle and myoepithelial cells. Composition share unchanged (median R^2 0.19 / 0.20).
- Rebuild spec v2 agreed (PI, 2026-10-18; not yet executed): [`SPEC_rebuild_v2.md`](SPEC_rebuild_v2.md), decisions D48-D56. v2 gene set (1,824 hard-to-map / readthrough genes dropped); Step 4 lab correction = leave-one-unit-out median adjacent-vs-GTEx logFC over 10 organ units added to GTEx normals; GTEx-only ischemic time + ventilator covariates; no RNA-quality index in Step 4; all 28 ARACNe networks rebuilt on the v2 genes with cross-mappable edges removed; ESTIMATE purity everywhere (LGG included); group and patient VIPER rerun; acceptance checks C1-C10 fixed in the spec. Steps 8-10 and the drug side follow after the checks and a PI decision on the recurring module (D26 vs its lab-axis link).

## 2026-10-19
- ARACNe v2 full build paused at 11:27 (user: machine needed by a colleague). State: COADREAD done (323,571 edges vs 343,608 in v1); BRCA 68 / 100 bootstraps, BRCA_LUM 84 / 100 (the two bootstraps in progress were discarded; all saved bootstrap files complete); the other 26 networks and the TF-only build not started. Pace before the pause: BRCA ~11.5 min, BRCA_LUM ~7 min per bootstrap.
- Resume (same or another machine): the script is resumable (keeps finished bootstraps, continues with seed = n saved + 1, skips networks with `network.txt`). Needs `data/processed/step5/aracne_inputs_v2/` and `data/processed/step5/aracne_v2/` (gitignored, ~2.6 GB; copy them over, or rebuild the inputs with `VERSION=v2 scripts/42_aracne_inputs.py`, byte-identical given the same Toil data), the conda env `aracne` (OpenJDK 17) and ARACNe-AP commit 30967ea at `~/tools/ARACNe-AP/dist/aracne.jar`, both set up by `scripts/07_setup_aracne.sh` (checked: a fresh build gives byte-identical threshold and bootstrap output to this machine's jar on KICH, seed 1). On a machine with a different core count set `THREADS` / `PARALLEL` (and `ARACNE_JAVA` / `ARACNE_JAR` for other paths); bootstrap output does not depend on the thread count (KICH seed 1, 8 vs 2 threads identical). Restart:
  `nohup setsid bash -c 'VERSION=v2 scripts/43_run_aracne.sh >> logs/aracne_v2.log 2>&1; VERSION=v2 REGSET=tf scripts/43_run_aracne.sh > logs/aracne_v2_tf.log 2>&1' > /dev/null 2>&1 < /dev/null &`
  To pause again, stop the whole process group (`kill -TERM -- -<pgid of the setsid bash>`); killing only the java processes would let the loop log the bootstrap as finished and start the TF-only build.

## 2026-10-20
- Subtype networks (D68): only the breast luminal test (BRCA_LUM vs BRCA); other subtypes stay on the shared cohort network.
- Test of A2, tumors and normals equally variable in limma (`scripts/36_voom_group_variance_test.R` -> `results/step4/voom_group_variance/`; v2 Step 4 model, 13 groups, voom vs voomByGroup (vbg, vendored `scripts/vendor/voomByGroup.R`, YOU-k/voomByGroup 38c8abf) vs voomWithQualityWeights by arm (vqw); master regulators with the v1 project ARACNe regulons). Decision rule fixed before running (script header).
  - The assumption is wrong: tumor BCV / normal BCV 1.19-1.77 (>= 1.2 in 12 of 13; COAD 1.77, Basal 1.74, LIHC 1.70; PRAD 1.19).
  - Split-half reproducibility (10 random halves of tumors and normals): z rho vbg - voom +0.0000 to +0.0061 in 12 of 13 (PCPG-PGL -0.0008), median +0.0013; NES rho (3 halves) higher in 10 of 13, median +0.0013 (PCPG-PGL +0.012, STAD +0.0065). vqw similar, slightly smaller gains.
  - Known-biology controls unchanged (10 pass, 42 right direction below FDR 0.05, 1 wrong: LIHC HNF1A, under all three methods).
  - Size of the change: whole-signature z rho with voom 0.996-0.999 and NES rho 0.996-0.9996, i.e. "no practical effect" under the pre-set rule (>= 0.99 / >= 0.98). The tops move more: 500 largest-|z| genes overlap 69-93% (median 85%); FDR < 0.05 master-regulator sets Jaccard 0.68-0.92 (READ 157 -> 230 regulators, THCA 1 -> 55).
- Decision (user, D69): adopt voomByGroup for the v2 Step 4 signatures; the top-of-list shifts are judged significant. voom signatures are kept alongside as check C12; the test results stay in `results/step4/voom_group_variance/`. A2 resolved. Not yet implemented in `30_signatures.R` (v2 Step 4 not run).
- v2 Step 4 signatures built (`VERSION=v2 scripts/30_signatures.R` -> `results/step4_v2/`, `data/processed/step4_v2/`): v2 gene set, Step 4a GTEx reference (sub-region filter, collection covariates centered at each tissue's anchor), voomByGroup (voom kept for C12). 40 tumor-vs-GTEx + 18 tumor-vs-adjacent signatures, as v1. GTEx normals per group by anchor class: observed (lung, skin, stomach, thyroid, testis, pancreas, muscle / adipose for sarcoma, esophagus mucosa), model (breast, colon, prostate, liver, adrenal for PCPG / NBL / ACC, uterus, ovary), average (kidney 21, brain 151, CESC-AD 7). Genes changed (FDR < 0.05, |logFC| > 1) mostly 5-15% fewer than v1 (e.g. KIRC 6,346 -> 5,081, COAD 6,700 -> 5,656).
- Step 4 v2 checks (`scripts/37_step4_v2_checks.R` -> `results/step4_v2/checks/`; v1 computed the same way for comparison):
  - C1 GC tilt of (GTEx minus adjacent) z: median -0.239 -> -0.298, larger in 14 of 17 groups (expected not removed; the gap is no longer subtracted). Decomposition: v2 model with voom -0.276, voomByGroup adds -0.02; the tumor-vs-GTEx z itself tilts slightly more with GC (median -0.126 -> -0.146; COAD -0.09 -> -0.14, READ -0.11 -> -0.17), the adjacent z not (0.119 -> 0.122). The C2 check at Step 6 (GC-only signature vs MR NES) measures whether this reaches the regulators.
  - C3 agreement of GTEx and adjacent signatures: median 0.593 -> 0.653 (COAD 0.555 -> 0.698, LUAD 0.593 -> 0.653, LUSC 0.726 -> 0.783, READ 0.572 -> 0.627; KIRC 0.33 -> 0.32).
  - C4 ischemic-time signal left in the GTEx signatures: median -0.176 -> -0.098, closer to 0 in 28 of 39 groups (lung -0.48 / -0.43 -> -0.30 / -0.24, COAD -0.36 -> -0.16, READ -0.33 -> -0.12); unchanged in brain (LGG -0.38 to -0.42, GBM -0.32): brain is an "average" tissue, so no anchor shift.
  - C6 proliferation retained: Hallmark E2F and G2M up at FDR < 0.05 in 38 of 40 groups (v1 38); not in PCPG-PGL / PHEO (as v1).
  - C11 split-half, GTEx (v2) vs adjacent, equal normal n per half: COAD 0.953 vs 0.966 (GTEx better in 3 of 10 splits), READ 0.887 vs 0.856, KIRC 0.927 vs 0.911, KIRP 0.922 vs 0.960, KICH 0.930 vs 0.938. Target (GTEx closer to or above adjacent) met for READ and KIRC, not for COAD and KIRP; halves are small (10 normals; READ 5).
  - C12 voomByGroup vs voom: z rho median 0.998 (TGCT-NSGCT 0.958 to 0.9996), top-500 overlap median 0.83 (SKCM 0.51, TGCT-NSGCT 0.64, UCEC-EEC 0.65, PAAD 0.68 to 0.95).
  - Flag: TARGET-WT genes changed 2,801 -> 6,916 (v1 -> v2; 28 -> 21 kidney normals). Pediatric tumors vs adult donors: age is almost collinear with tumor / normal, so the WT (and NBL) condition estimate is unstable to small changes in the normal set.
- Step 4 v2 result summary (recorded in the spec doc):

| Check | v1 -> v2 | Verdict |
|---|---|---|
| C3 Agreement with tumor-vs-adjacent signatures | median 0.59 -> 0.65 (colon 0.56 -> 0.70, lung up) | Better |
| C4 Post-mortem signal left in signatures | median -0.18 -> -0.10; closer to 0 in 28 of 39 groups (lung -0.48 -> -0.30, colon -0.36 -> -0.16) | Better; brain unchanged (average-donor anchor, no shift) |
| C6 Proliferation retained | E2F / G2M up in 38 of 40 groups (v1 38) | Kept |
| C12 voomByGroup vs voom | top-500 genes shared: median 83%, lowest SKCM 51% | As in the test |
| C1 GC tilt (GTEx minus adjacent z) | median -0.24 -> -0.30 | Larger; accepted (D70) |
| C11 Split-half reproducibility, GTEx vs adjacent | READ, KIRC: GTEx ahead; COAD, KIRP: adjacent ahead; KICH similar | Mixed |

- Decision (user, D70): C1 accepted. A GC difference between healthy donors and tumors may be partly real and is not by itself a problem; the aim is an accurate, reproducible workflow, not removing every edge-case signal. C2 at Step 6 is still reported. Pediatric groups (TARGET-WT, NBL) recorded as more variable (children's tumors vs adult donors; A36); no action unless a correction is needed for all cohorts.

## 2026-10-25
- ARACNe v2 build finished on MAIL-15 (resumed 2026-10-20 from the MAIL-10 state; last full-set network KIRC 2026-10-24 02:51), then the TF-only set; networks + all 100 bootstraps per network copied back to MAIL-10 (`data/processed/step5/aracne_v2/`, `aracne_v2_tf/`).
- Fix: the v2 regulons built on 2026-10-24 (and the group VIPER run on them) still contained the cross-mappable edges (~13%); `44_aracne_regulons.R` never had the filter (v1 applied it only at patient level in 61). Both are kept as `*_unfiltered` (regulons, `results/step6/aracne_project_v2{,_tf}_unfiltered/`), not used.
  - `44_aracne_regulons.R` (VERSION=v2, REGSET=tf): edges between cross-mappable pairs (symmetric cross-mappability > 0) removed from `network.txt` before `aracne2regulon`, as the spec says (edge likelihood is scaled by the network's largest MI, and paralog edges are among the strongest). Removed per network: full set 11.8-14.6% (median 13.1%; v1 12.7-14.2%), TF-only 10.7-14.9% (median 11.6%); counts in `results/step5/crossmap_edges_removed_v2{,_tf}.tsv`. Networks after the filter: full 128k-342k edges, 4,266-4,889 regulators; TF-only 64k-190k edges, 1,105-1,360 regulators.
  - Pseudoautosomal genes: the mappability file lists each PAR gene twice (X copy ~1.0, `_PAR_Y` copy 0). The v2 gene set (22) keeps the X copy; the edge filter gives the same edges with or without the `_PAR_Y` rows (0 differences over 58 networks). `61_quality_correct.R` and the checks now skip `_PAR_Y` rows (61 would have dropped the 13 PAR genes).
- Group VIPER v2 rerun on the filtered regulons (`50_viper.R`, VERSION=v2, full and REGSET=tf; 41 groups each); controls (`51_compare_mr.R`, REGULON=aracne_project_v2{,_tf}); bootstrap support (`45_edge_support.R`, VERSION=v2 [REGSET=tf] -> `results/step5/aracne_edge_support_v2{,_tf}.tsv`).
- Step 6 v2 checks (`scripts/63_step6_v2_checks.R`; SET=full: v1 vs v2 full -> `results/step6/checks_v2/`; SET=tf: v2 full vs v2 TF-only -> `checks_v2_tf/`):

| Check | v1 | v2 full | v2 TF-only |
|---|---|---|---|
| C10 low-mappability genes / cross-mappable edges left | - | 0 / 0 | 0 / 0 |
| C9 edges per network, median | 342k | 266k | 157k |
| C9 bootstrap support, median of network medians (range) | 11 (10-12) | 11 (10-12) | 14 (13-16) |
| C2 GC-only signature vs MR NES, median rho (groups > 0.4) | 0.15 (2) | 0.17 (2) | 0.10 (3) |
| C5 controls pass (proliferation + tissue-specific) / wrong direction | 70 (64 + 6) / 2 | 66 (60 + 6) / 2 | 78 (70 + 8) / 4 |
| C7 MRs FDR < 0.10 per group, median | 587 | 567 | 172 |

  - C2: computed the same way for every column (aREA of each group's kept regulons on minus z-scored gene GC; rho > 0 = MR ranking follows GC). Not lower in v2 (expected "lower"), in line with the larger C1 tilt accepted in D70; concentrated in KIRC (0.65 -> 0.69) and TARGET-WT (0.63 -> 0.65; ~90% of their FDR < 0.10 MRs have the GC-predicted sign). TF-only third group: UCS -0.48 (opposite sign).
  - C5 full set: tissue-specific passes unchanged (6); proliferation passes 64 -> 60 (spec: not lost): TARGET-WT E2F1 / MYBL2, BLCA MYBL2, COAD FOXM1, STAD FOXM1 drop below the cutoff (direction still correct); TGCT-NSGCT E2F1 now passes; KIRC EPAS1 wrong -> correct direction; LIHC HNF1A correct -> wrong direction. TF-only (vs v2 full): newly passing PRAD AR, LumB FOXA1, proliferation in BLCA, LUAD, READ, SKCM, UCS, LIHC MYBL2; wrong direction PCPG-PGL / PHEO E2F1, KIRC MYBL2 (and LIHC HNF1A).
  - C7 v1 vs v2 full: NES rho median 0.977, FDR < 0.10 Jaccard median 0.67 (voomByGroup, GTEx 4a reference and the edge filter together). v2 full vs TF-only: 90% of full-set TF MRs are also TF-only MRs, 72% of TF-only MRs are full-set MRs (median per group), sign agreement 100%, NES rho on shared regulators 0.97 (Jaccard 0.19 is low by construction: ~1,000 vs ~4,000 regulators scored).
- Patient VIPER v2 (`61_quality_correct.R` VERSION=v2; spec Step 6b): ESTIMATE purity in every group (LGG included); GTEx normals = Step 4a reference, regressed on quality index + collection covariates and set to the anchor (residual + intercept); filtered v2 regulons; checks vs v1 `aracne_project_map`. Outputs `data/processed/step6/aracne_project_v2_patient/`, `results/step6/aracne_project_v2/quality_correction/`. Done for 41 groups:
  - Patient concordance with group MRs unchanged (median rho 0.676 -> 0.675, frac_hit 0.98 both); drops > 0.05 in ESCA-SCC (0.56 -> 0.50) and READ (0.66 -> 0.59). Effective dimensionality median 9.2 -> 8.8. Quality index on PC1 R^2 <= 0.01 everywhere. PC1 plate excess over chance > 0.2 in BRCA-Basal only (0.36 vs chance 0.19); STAD 0.21 vs 0.04 (excess 0.18).
  - LGG (ESTIMATE replaces ABSOLUTE, D45): the purity signal moves rather than disappears. v1 (ABSOLUTE-corrected) left ESTIMATE on PC1 (R^2 0.31 codel, 0.21 noncodel); v2 (ESTIMATE-corrected) leaves ABSOLUTE on PC1 (R^2 0.46 codel, 0.24 IDHwt, 0.11 noncodel). In glioma the two measure different things (immune / microglia vs copy-number tumor fraction). For the user / PI; no change made.
- C8 recurring module (`62_recurring_module.R` NES_DIR=aracne_project_v2_patient -> `results/step6/aracne_project_v2/diagnostics/recurring_module_*_v2.tsv`; v1 = `_map`): still present after all v2 corrections. Found in 36 of 41 groups (v1 39); core (>= 30% of groups) 77 members (v1 55; 54 of the 55 kept, the one lost is JMJD7-PLA2G4B, dropped by the v2 gene filter as readthrough); median module size 52 (v1 33). Top members unchanged (PILRB, ZNF692, CCNL2, PLA2G4B, WDR90, HSF4, ANKZF1, MAPK8IP3, STK36). Plate: median R^2 0.24 (v1 0.23), excess over chance median 0.17 (v1 0.16), > 0.2 in 12 of 34 groups (v1 13 of 36). Still anti-correlated with metabolic Hallmarks (adipogenesis -0.52, xenobiotic -0.51, fatty acid -0.49, androgen -0.49, OXPHOS -0.45) and not with immune (ESTIMATE immune +0.01, cytolytic +0.02), purity (-0.08), quality index (-0.06) or GC index (+0.11). Feeds the D26 decision (keep in clustering vs remove), which is needed before Steps 8-10.

## 2026-10-28
- Edge confounding v2 (`46_edge_confounding.R` VERSION=v2, full set, all consolidated edges before the cross-mappable filter; run on MAIL-15 2026-10-24, committed now -> `results/step5/aracne_edge_confounding_v2.tsv`, 29 networks). Same picture as v1 (A9): edges retaining < 50% of |r_raw| after ESTIMATE + quality + plate median 1.9% (v1 2.0%), max UCS 8.5% (v1 8.7%); ESTIMATE alone 0.24% (0.23%), quality 0.44% (0.46%), plate 0.01% (0.01%); subtype (pooled cohorts) 0.64% (0.80%), max ESCA 3.8%. Median retained |r| after full adjustment 0.887 (0.890). TF-only set not run.
- The TF regulons built on MAIL-15 on 2026-10-24 (before the edge filter) are renamed there to `regulons_aracne_project_v2_tf_unfiltered.rds`; MAIL-10 holds the filtered version (checksum in `results/step5/regulon_checksums.sha256`).

## 2026-10-28 / 29: v2 TF arm through to drug ranking (decisions D71-D80)

Summary in one paragraph: the v2 pipeline now runs end to end on the TF-only ARACNe networks (the active build, `v2_tf_xmap`). Along the way we fixed a bug in the group-level null, tested whether a covariate-adjusted null helps, chose patient subgroups with a fixed rule, and rebuilt how drugs are scored so that the results say honestly how much evidence lies behind each drug list. Decisions are D71-D80 in `docs/DECISIONS.md`; the drug-side details are in `docs/notes/drug_side_decisions_and_caveats.md` on the `celligner` branch (C18-C22). Requirements for the future LLM reporting layer are collected in `docs/LLM_COORDINATION_NOTES.md`.

### Work from MAIL-15 committed
- Build registry (`config/builds.yaml`, `R/paths.R`, `scripts/00_status.R`): each build (v1, v2_full_xmap, v2_tf_xmap) is named once with all its paths, so scripts cannot silently mix versions. Active build: `v2_tf_xmap`.
- Step 6e on the TF arm (purity and quality corrections, `61` REGSET=tf) and Step 8 consensus clustering on the TF arm (`results/step8_v2_tf/`).

### 1. Group-level null: a bug fixed, and the adjusted null tested (D73)
- **What the null is.** To decide whether a regulator's tumor-vs-normal difference is real, we compare it with differences that arise by chance: the healthy samples are split randomly into two halves 1,000 times and compared with each other.
- **Bug.** For v2, `50_viper.R` built this null from all of a group's GTEx samples, including samples Step 4a had removed from the signature (colon muscle wall, kidney medulla, samples without collection data). The null should use exactly the samples the signature used. Affected 33 of 40 groups (colon: 166 instead of 109 normals). Fixed (fix "B").
- **Test of a covariate-adjusted null ("C")** (`scripts/64_null_adjust_test.R`, 11 groups): also removing age, purity, sex and collection effects from the healthy samples before splitting, as the signature model does. A safety check (null from one half of the normals, fake signature from the other half) found 0 false calls for both B and C. C raised the median number of FDR < 0.05 master regulators (19 -> 74) and gained 4 known-biology controls but lost one (MYBL2 in luminal A breast). The rule fixed before running said "lose no control", so **B was adopted, C not**.
- **The key finding:** FDR-cut lists are unstable. The number of FDR < 0.05 regulators jumped between the three nulls (colon 425 / 0 / 1; stomach FDR < 0.10 210 / 0 / 220) while the regulator ranking stayed the same (Spearman >= 0.98). Hundreds of regulators sit right at the threshold. This led to D74.
- The test got stuck for hours on CESC-AD: with 7 normals, the safety check made a 2-vs-1 comparison, which produced an all-missing signature on which msviper never returns. A bug in the test only; the other 11 groups were fine.
- **Reruns on the corrected null (B):** TF arm (41 groups + 14 subgroups), CollecTRI (41 + 14), full regulator set (41), DoRothEA (41; subgroups running). TF arm: old vs corrected null NES Spearman median 0.999; FDR < 0.10 master regulators median 172 -> 151. CollecTRI reaches FDR < 0.10 in few groups (median 0), which no longer matters for the drug side (D74).

### 2. Patient subgroups (D71, D79)
- **Keep rule (D71, fixed before deciding):** a Step 8 clustering becomes subgroups only if the clusters are statistically supported (M3C), reproducible on held-out halves of the patients (prediction strength >= 0.75), large enough (>= 30 patients), found by both feature sets (agreement >= 0.3), and not driven by plate, purity or RNA quality. Why: subgroups get their own signatures and drug lists, so an unstable or technical split would carry through to the drugs. Script `84_subgroup_consistency.R`.
- **Correction:** COAD's prediction strength is 0.748. It was first passed by rounding to 0.75; strictly it fails. Recorded in D71.
- **Without the recurring module.** The plate-linked recurring module (D26) was removed and the clustering repeated for all 41 groups (`80` EXCLUDE_MODULE=recurring). The decision changed in only 5 of 41 groups, mostly borderline ones, so the module does not create or hide structure broadly. BLCA and HNSC keep the same splits (agreement 0.87 / 0.91 with the original).
- **Decision (D79, option B):** define all subgroups from the clustering without the module (treating it as a technical nuisance), plus COAD as an explicit exception because its MSI-like vs CIN-like split matches the TCGA subtypes and was replicated in two external cohorts in v1. ACC is dropped (its split is not significant without the module). 14 subgroups, named from the TCGA subtypes of their members:
  - BLCA basal-squamous-like / luminal-like; HNSC basal-mesenchymal-like / atypical-classical-like; LGG-IDHmut-noncodel LGr2 / LGr3; PCPG-PHEO pseudohypoxia-Wnt / kinase-signaling; THCA BRAF-like-classical / RAS-like-follicular; TGCT-SEM C1 / C2; COAD CIN-like / MSI-like.
  - Table: `scripts/85_step10_strata_v2.R` -> `results/step10_v2/strata.tsv`, `stratum_patients.tsv`.
- **Subgroup signatures and VIPER:** `30_signatures.R` STRATA= (each subgroup's tumors vs the parent's normals, same v2 model) and `50_viper.R` STRATA= (group and patient level, parent's network and normals). Group level gives each subgroup's master regulators; patient level lets us measure variability inside each subgroup.
- **What defines each subgroup** (`87_subgroup_characterize.R` -> `results/step10_v2/defining_*.tsv`): BLCA basal-squamous-like has the PPARG/GATA3 luminal program down and FOXM1/E2F proliferation up; COAD MSI-like has CDX2/SATB2 intestinal differentiation down and NF-kB (RELB, NFKB2) up; THCA BRAF-like has the THRB differentiation program down and an inflammatory program up; HNSC basal-mesenchymal has the SOX2 program down. Many smaller programs are groups of chromosome-19 zinc-finger TFs, which tend to move together for reasons that may be technical; check before reporting any of them.
- **Variability within subgroups** (`86_subgroup_variability.R`, earlier D72 subgroups): patients match their subgroup's profile better than their parent group's in most subgroups. Flags: in COAD CIN-like the two largest axes of variation track plate (R^2 0.31 vs 0.07 by chance), so patient-level differences inside it need care; ACC proliferative was the least coherent.
- Earlier subgroup list D72 (BLCA, HNSC, COAD, ACC, LGG-noncodel; 10 subgroups) is replaced by D79.

### 3. Drug side (celligner branch; details in its notes C18-C22)
- **Ported to the v2 TF arm:** `79e`, `79h`, `79i` (VERSION / REGSET switches); `79b` needs no v2 version (same groups and tumors). Subgroups: line weights (`74b`), target gate (`79b` STRATA=1), and a subgroup mode in `79e` (STRATA=1, replaces `79g` for v2) with a sibling-specificity score.
- **Anti-estrogens misread (ESR1).** In MCF7, anti-estrogens lower the classic ER targets but raise GATA3 and FOXA1, which are ESR1's strongest targets in the tumor network, so VIPER read ER blockade as ER activation. Shadow (pleiotropy) correction did not fix it; the luminal-only network halved the error; the CollecTRI regulon (direct ER targets) read it correctly. A LINCS knockdown benchmark showed the ARACNe ESR1 regulon itself is fine (it reads ESR1 knockdown correctly): the problem is specific to these drugs.
- **Knockdown benchmark** (`79m`, `79n`; LINCS shRNA consensus + CRISPR, 1,803 usable signatures, 336 TFs): all regulon sets detect knockdowns only weakly and about equally (ARACNe vs CollecTRI: better in 51% of TFs). Too little data to choose a regulon set per TF.
- **Decisions:** D74 every regulator weighted by its disease NES, no master-regulator cutoff (because FDR lists are unstable, D73); D75 drugs scored with a consensus of ARACNe TF and CollecTRI (average per regulator over both sets, the same on the disease and drug side); D76 line-match flag for drugs tested only in cell lines that do not resemble the disease; D77 transfer-confidence flag for borrowed drugs; D78 evidence status per disease; D80 reports lead with disease specificity.
- **Results, 41 groups:** 14 groups have drugs tested in matching cell lines ("matched"); 27 rest only on borrowed drugs ("borrowed only"), including basal breast (no basal-like line has LINCS compound data), all brain tumors, kidney, sarcoma, testis, thyroid, adrenal. Known drugs: lung adenocarcinoma erlotinib 96th percentile (osimertinib #1), prostate flutamide 85th, colon irinotecan 98.5th on specificity. Basal-breast anti-estrogens still score high but are flagged and kept out of the recommended tiers. By overall reversal the same generic drugs top most diseases (statins, palbociclib, MEK inhibitors, anthracyclines); by specificity none does, hence D80.
- **Results, 14 subgroups:** only the COAD subgroups are "matched". COAD MSI-like: top drugs are BRAF / MEK inhibitors, and the drugs most specific to MSI-like over CIN-like are HSP90 inhibitors (MSI colon cancer is enriched for BRAF V600E, an HSP90 client). THCA BRAF-like vs RAS-like: HDAC inhibitors (studied for re-differentiating BRAF thyroid cancer). LGG and TGCT-SEM subgroups barely differ in predicted drug response.
- PRISM / GDSC validation on the v2 scores: on hold (user).

### Open
- D26 (recurring module): the evidence now favors treating it as a nuisance axis (subgroups are defined without it, D79); PI to confirm.
- D60 kidney reference; LGG purity covariate (A37, tabled); BRCA_LUM material-change cutoff (D68).
- Replication of the new subgroups in external cohorts (CPTAC-HNSCC, a bladder cohort; COAD GSE39582 / CPTAC-COAD as in v1).
- Small chromosome-19 zinc-finger programs in the subgroup characterization: check before reporting.

## 2026-10-29: Step 7 diagnostics on the corrected null
- **Known-biology controls** (`51`, every v2 set): TF arm 77 pass / 55 correct direction below FDR 0.10 / 4 wrong
  direction (PCPG E2F1, KIRC MYBL2, LIHC HNF1A; all NES near 0, FDR >= 0.83, so non-calls rather than contradictions).
  The null fix (D73) moved 5 of 139 verdicts across the cutoff (COAD / READ E2F1 and PRAD AR now just below; STAD
  FOXM1 / MYBL2 now pass); none changed direction. Full set: 67 / 67 / 2. CollecTRI and DoRothEA pass few controls
  (13 and 20), and their wrong-direction calls are also non-significant.
- **Three regulon sets and the D18 confirmed rule** (`59`, now with v2 switches and `confirmed.tsv`): median 117
  confirmed regulators per group on the TF arm (5,215 total; subgroups median 80). CollecTRI (median 0 per group) and
  DoRothEA (1) rarely reach FDR < 0.10, so the confirmed list is in practice ARACNe FDR < 0.05. Regulator rankings agree
  weakly between sets: NES Spearman ARACNe vs CollecTRI 0.18, vs DoRothEA 0.27, CollecTRI vs DoRothEA 0.37 (median per
  group). Groups with few confirmed regulators: CESC-AD 0, BLCA / UCS / STAD 1, HNSC 4.
- **Regulator modules** (`58`, v2 mode): 714 modules of co-varying regulators (median 18 per group, 15 members); 59 track
  proliferation (|r| >= 0.7), 62 plate (R^2 >= 0.3). The largest modules are often KRAB zinc-finger groups (ZNF333,
  ZNF224, ZNF141, DMTF1, CREBZF), overlapping the recurring module.
- **TF-only edge confounding** (`46` REGSET=tf): edges that lose half their strength after ESTIMATE + quality + plate,
  median 4.0% per network (full set 1.9%), max UCS 12%; plate alone 0.01%, ESTIMATE alone 0.9%.
- Questions raised for the user: keep D18 as is; implications of the weak ARACNe-CollecTRI agreement for the
  drug-side consensus (D75); whether 4% (UCS 12%) TF-edge confounding needs action.

## 2026-10-31: advisor code review (5 points), all confirmed
- **BRCA_LUM test not implemented (D65, D68).** The luminal-only networks exist, but nothing compares them with the
  pooled BRCA network, no cutoff for a "material change" was ever recorded, and `50_viper.R` maps LumA / LumB to the
  pooled BRCA network by cohort. Open: the user sets the cutoff, then the comparison is written and run.
- **Hidden v1 dependency in `42_aracne_inputs.py`:** chose cohorts from the v1 `results/step4/signature_summary.tsv`
  even with VERSION=v2. Fixed (reads `results/step4_v2/`); the v1 and v2 ok-group sets are identical (41), so the v2
  inputs of 2026-10-18 are unaffected.
- **Sex term:** `30_signatures.R` also drops sex when the smaller sex has < 3 samples, while the docs said only "when one
  sex". Documented (script header, PIPELINE_V2); checked that in v2 this second case never occurs, so no result changes.
- **Count cache:** `30_signatures.R` rebuilt its cached count matrix only when samples were missing. It now also checks
  the gene set (the current cache holds all 19,797 protein-coding genes, so nothing is rebuilt).
- **Uncorrected patient NES:** `50_viper.R` writes uncorrected patient NES to `data/processed/step6/<set>/`, used only for
  its own diagnostics; downstream steps use the corrected NES of `61` in `<set>_patient/`. Labelled in the script header
  and PIPELINE_V2.

## 2026-11-01: BRCA_LUM test set up (D81)
- **Cutoff fixed before any comparison (D81).** Material if, in LumA or LumB, (a) the weighted Spearman of NES on BRCA_LUM
  vs pooled BRCA is below that of every size-matched control and at least 0.10 below their median, or (c) ESR1, FOXA1,
  GATA3 or XBP1 flips sign (|NES| >= 1.96 on one side) in a way no control does. Top-list overlap was rejected as a
  rule (cut lists are unstable).
- **Why controls:** BRCA_LUM has fewer tumors (769 vs 1,040) and only luminal ones. Pooled BRCA can't separate the two;
  three networks of 769 random BRCA tumors (`47_brca_size_controls.py`) change only the size.
- Controls BRCA_SUB1-3 built on MAIL-10 (11:22-13:51); comparison script `66_brca_lum_test.R` written and dry-run
  checked (every network replaced by pooled BRCA). `44_aracne_regulons.R` gained `OUT_TAG` so the controls stay out of
  the production regulon file.
- Added `67_esr1_regulon_content.R` (descriptive, not in the rule): do the luminal regulons contain more
  estrogen-response genes on BRCA_LUM than on pooled BRCA and the controls? A test run on pooled BRCA and BRCA_LUM
  showed regulon sizes (ESR1: 450 targets pooled, 154 on BRCA_LUM) before the controls exist.
- Known before the test (pooled BRCA, already reported): ESR1, FOXA1, GATA3, XBP1 all positive in LumA and LumB.

- **BRCA_LUM result (D81): not material.** (a) weighted Spearman LumA 0.985 vs controls 0.994 (lower than all three,
  but by 0.01; the rule needs 0.10); LumB 0.989 vs 0.989-0.995; (c) no sign flips. LumA / LumB stay on the pooled
  network and no other subtype networks are built (D68). Pooled-BRCA sanity vs `50`: r = 0.999.
- Context, not in the rule: on BRCA_LUM, ESR1 and GATA3 read stronger in LumB (1.76 -> 3.06, 1.48 -> 3.13; controls
  1.5-1.8). Regulon content (`67`): BRCA_LUM regulons are a third the size (ESR1 450 -> 154) and not more
  estrogen-specific (9.1% vs 11.1%; controls 10.4-11.9%). Cross-check with D82: the ESR1 targets BRCA_LUM drops are 59%
  identity-only, the ones it keeps 85% signaling.

## 2026-11-01: identity vs signaling edges (D82); lineage-matched data (C24, R3)
- **Why:** a network built from mixed subtypes can call a gene a target only because both genes mark one subtype
  ("identity"). Drugs change signaling within hours, not identity, so identity targets can make drugs misread.
- **Rule (D82, fixed before running):** per edge, pooled Spearman vs the correlation left after removing subtype
  averages; identity-only < 50% kept (the `46` cutoff). Labels: histology (our subtype groups) and established
  molecular subtypes (PanCanAtlas mRNA, else selected), reported separately.
- **Result (`48_identity_edges.R`):** 0-7% of edges per network are identity-only (ESCA 6-7%, LGG / TGCT 4%, SARC 3.5%,
  BRCA 3.3%, COADREAD 0%), but they concentrate in lineage regulators: BRCA ESR1 51%, XBP1 46%, GATA3 45%, FOXA1 35%;
  ESCA HNF1A / B, FOXA2 / 3, PDX1 82-89%, HNF4A 73%, GATA6 58%, TP63 47%; BLCA (molecular) GATA3 75%, PPARG 74%.
  Proliferation regulators are signaling (BRCA MYC 1%, E2F1 7%). Subtype-specific edges (rewiring evidence, A9): up to
  17% (SARC), 15% LGG molecular; BRCA GATA3 59% of its targets.
- **Open:** whether drugs are scored with signaling-only regulons (separate decision, rule first).
- **Lineage-matched cell lines (C24, celligner):** the cell type must match, not only the organ. Squamous groups have
  no lineage-matched LINCS data (only HCC-15, 3 compounds, resembling none). New reporting rule R3: report "no
  good-quality data at <step>" instead of a list (`79r_data_gaps.py`; minimum 10 compounds on matched lines, user):
  14 of 41 groups and 2 of 14 subgroups have good-quality lineage-matched data.

## 2026-11-01: signaling-only target lists, comparison (D83)
- Three arms with a shared null (one seed per group): base (current lists), sig1 (identity-only edges under histology
  labels removed: 0.8% of all edges, 12 regulators fall below 10 targets), sig2 (histology + molecular: 2.3%, 34).
- **Rule met for both, effect small.** Anti-estrogen median `spec_mech` percentile LumA 9.7 -> 12.3, LumB 23.4 -> 26.0,
  Basal 95.6 -> 94.6; knockdown share read as lower activity 0.597 -> 0.599; controls 77 pass / 4 wrong in all arms;
  drug ranking agreement with base median 1.00 (lowest ESCA-AD 0.986 / 0.978). Base vs production: rho >= 0.999.
- Disease side moves in the expected direction: ESR1 NES LumA 2.71 -> 3.01, LumB 1.73 -> 2.11, Basal -2.43 -> -1.53.
- **Why the drug side barely moves:** the ARACNe ESR1 list (also sig1) reads anti-estrogens in MCF7 as ESR1 activation
  (+2.7 / +2.2, all 10 profiles; CollecTRI -2.4). Under anti-estrogens GREB1 falls (z -1.3) and the Hallmark estrogen
  genes fall (-0.2), but the genes ARACNe calls ESR1's activated targets rise (+0.36): co-expression across tumors is
  not the estrogen response of cells. The consensus averages +2.2 and -2.4 to about zero. A proliferation explanation
  was checked and ruled out (1 of 122 negative targets is a proliferation gene).
- Script fix during the run: `79s` first counted controls with the wrong label (0 passing); corrected before reading.

## 2026-11-02: docs audit, celligner merged into main, plain-language versions
- `celligner` merged into `main` (1e7d31b): drug-side code, results and docs are now on `main`, which holds all docs
  from here on (user). Main kept its own versions of the docs it maintains.
- Both technical docs checked against the code and results (paths, scripts, thresholds, counts): correct apart from
  status lines, branch pointers and open items, now updated. `run_v2_drug.sh sync_disease` is a no-op on `main`.
- New `DISEASE_SIDE_PLAIN_LANGUAGE.md` beside `DRUG_SIDE_PLAIN_LANGUAGE.md`; v1-era and finished docs moved to
  `docs/archive/` (index there); README rewritten as the entry point.

## 2026-11-02: variability inside the D79 subgroups, compared with the 41 cohorts
- `86` rerun on the 14 D79 subgroups, plus every cohort with all its patients (user). Size-matched measures (30
  patients, 20 subsamples) fixed before looking; seeding made reproducible across workers (L'Ecuyer).
- Total spread (summed NES variance): 13 of 14 subgroups below their cohort (ratio 0.53-0.94; LGG LGr2 1.07); 11 of 14
  in the lowest quarter of cohorts (cohort median 738, range 377-1,933). Coherence at n = 30 up in 10 of 14; own vs
  sibling 0.85-1.00 (median 0.96).
- Pre-set rule (coherence up and effective dimension down): 0 of 14, because effective dimension rises in 13 of 14:
  removing the parent's main axis (the split) spreads the rest over more axes. Total spread was added after the first
  run and is labelled so in the script.
- PC1 plate excess > 0.15 in TGCT-SEM C1 / C2 and COAD-CIN.

## 2026-11-02: ARACNe vs CollecTRI agreement on drug profiles (C25 follow-up, `79t`)
- **Question (user):** how many regulators behave like ESR1 (read in opposite directions by the two regulon sets)?
  Diagnostic only. Classes fixed before looking: Spearman over all 9,678 drug x line profiles, disagree < 0, weak
  0-0.3, agree >= 0.3; ESR1-type = opposite sign in > 50% of profiles where CollecTRI |NES| >= 2 (>= 30 profiles).
- **Result:** of 8,148 network x regulator pairs (375 regulators, 28 networks) 2,658 disagree (33%), 4,104 weak,
  1,386 agree; 2,638 ESR1-type. Median rho 0.09; curated-vs-curated baseline (DoRothEA vs CollecTRI) 0.23. 287
  regulators disagree in >= 1 network, 108 in at least half. Strongest: SP3, SP1, NR2F1, ETV3, BACH1 (rho -0.3 to -0.4
  in all 28 networks), also AR, NR3C2, TP73. Identity share (D82) adds a little: disagree 32% (< 10% identity) vs 46%
  (>= 30%). ESR1 disagrees in BRCA, CESC, SKCM, UCEC, OV; weak elsewhere.
- **How to read it (corrected the same day):** low agreement over all profiles does not show that either set is wrong:
  most drugs do not act on most TFs, so the correlation mostly reflects which general programs each target list picks
  up (the two curated sets, which share sources, agree at only 0.23). On knockdowns the two sets read TF loss equally
  (AUC shRNA 0.55 vs 0.54, CRISPR 0.63 vs 0.62; both read ESR1 knockdown as lower, -1.6 / -1.8). The ESR1 inversion is
  specific to the drugs, not a wrong ESR1 regulon.
- ESR1 in MCF7 (scratch check): ARACNe positive targets that are Hallmark estrogen genes (40) fall slightly under
  anti-estrogens (-0.09); the other 132 rise (+0.37; GATA3 +4.1, FOXA1 +1.1), on measured landmark genes (+0.71) more
  than inferred ones (+0.18). The line-specific projection (`79c`) roughly doubles the misread (raw +1.1, specific +2.4).
  One estriol profile only, too weak (GREB1 +0.45) to check the agonist direction.
- Ground truth available for a per-regulator choice: drugs acting directly on a TF with a known direction, mainly
  nuclear-receptor ligands (NR3C1 29 agonists, ESR1 13 antagonists / 9 agonists, RARA/B/G, AR, PGR, PPARG, VDR, NR3C2,
  RXR): about 12 TFs. Open: the per-regulator rule (fix before running).

## 2026-11-02: perturbation benchmark with overexpression (D84)
- **Rule (D84, fixed before running):** LINCS overexpression (wild-type ORF projects OEB / OEC; allele screens excluded)
  added to the knockdown benchmark, expected direction up; a set replaces the consensus only if better on both kinds,
  p < 0.05 and >= 5 points; per-TF switch only with >= 6 signatures in >= 2 lines; flags, not removal.
- **Data:** 614 overexpression signatures (197 TFs) subset (`08b`, `79u`); 533 scored; own mRNA up in only 94 of the
  365 where measured (L1000 probes may miss ORF transcripts), so the run without that check is reported too.
- **Result: keep the consensus.** 57 TFs testable on both sets: ARACNe better 25, CollecTRI 18, tied 14 (p = 0.15;
  without the own-mRNA check 90 TFs, p = 0.35). Both read knockdowns correctly 60% and overexpression ~70%.
  Per TF: ETV1, ATF1, SPDEF -> ARACNe, NCOA3 -> CollecTRI. Flags: 32 tested, 99 not supported, 197 insufficient,
  262 untested.
- **Ligands (report only):** CollecTRI ahead on line-specific profiles (0.64 vs 0.48, mainly nuclear receptors: ESR1,
  RARA, RXRB, PPARG agonists; ARACNe ahead on TP53, VDR), tied on raw profiles (0.58 vs 0.57). The line-specific
  projection makes ARACNe's ESR1 misread in MCF7 worse (+0.9 raw -> +2.0).
- **Found:** production drug activity predates the current regulator filter (30 rows more, 8 fewer on a rerun);
  `79d` raw run rewrites the regulator table, restored from a copy.

## 2026-11-02: exploration, can ARACNe edges be corrected? (not a decision; scratch scripts)
- **Edge-level test:** every ARACNe edge vs the D84 knockdown / overexpression signatures of its regulator
  (lineage-matched networks, 116k edge x signature pairs): overall 55% move as predicted. Best predictor: recurrence
  (share of the regulator's networks with the same signed edge): 52% at <= 10% -> 68% at > 75%; holds within the same
  signature (+11 points, 73% of signatures, p ~ 1e-24) and for shRNA, CRISPR and ORF. Negative-mode and D82 identity
  edges ~ chance; hubs, targets that are regulators, generic responders: no effect.
- **Hidden factors (BRCA, LUAD, COADREAD, SKCM):** removing the top 5 / 10 / 20 principal components of the tumor
  expression removes >= half the correlation of 52-64% / 72-80% / 84-87% of edges (D82 identity edges: 98% removed;
  D82 signaling: 76% at k = 10). Surviving edges agree more (k = 10: 61% vs 55%; +6 points within signature), partly
  the same signal as recurrence. Removed edges stay above chance (54-55%).
- **Regulon level: no gain.** Down-weighting by recurrence: +0.007 AUC on held-out TFs (5-95% -0.005 to +0.016); by
  factor survival or both: -0.0002. Per regulator, the share of factor-removed edges predicts neither benchmark
  failure (rho +0.10) nor drug-side disagreement with CollecTRI (-0.07).
- **Why:** major regulators drive the major axes (BRCA ESR1 / FOXA1 / GATA3 / XBP1 92-93% removed at k = 5, SKCM
  MITF / SOX10 75%), so removing the factors removes their programs too. ESR1: 167 of 172 positive targets removed;
  the 5 left still rise under anti-estrogens. From tumor co-expression, a lineage regulator's targets cannot be told
  apart from its lineage program; the ESR1 drug misread is structural, not an edge-quality problem.
- **CollecTRI edges, same test** (`scripts/explore/collectri_edges.R`): 52% move as predicted (ARACNe 55%; not the
  same signature mix: CollecTRI over all lines, ARACNe lineage-matched); repression edges 51%; CollecTRI edges also in
  ARACNe 56% vs 51% for the rest; 32 of 80 TFs below 50%. Per edge, neither set is validated; edges both sets support
  are the most reliable (ARACNe edges with the same sign in CollecTRI: 63%).

## 2026-11-02: perturbation-network pilot (D86): stop
- **What:** ARACNe on LINCS drug profiles within one cell line (MCF7, PC3, A375; 1,081-1,221 drugs each; 978
  measured genes, 81 measured TFs), raw and line-specific, 36 networks (all drugs + 5 mechanism-grouped held-out
  builds). Files in `scripts/pnet/`, `results/lincs/pnet/`; plain language in
  `PERTURBATION_NETWORK_PILOT_PLAIN_LANGUAGE.md`.
- **Result: no advantage, rule not met.** Pre-set check (knockdown / overexpression in the same lines, regulators all
  three sets score): 45 signatures, 13 TFs; perturbation network 62%, tumor ARACNe 69%, CollecTRI 62% (p = 0.63).
  Wider exploratory comparison (133 signatures, 37 TFs): 57% vs 59%; per line inconsistent. Networks moderately
  stable (61-65% of edges in >= 4 of 5 held-out builds).
- **Also found:** L1000 inferred genes (ESR1, AR, FOXA1, E2F1, SP1, RELA ...) are ~88% a fixed mix of landmark genes
  on held-out profiles (measured genes from the other 977: 62%), so they cannot serve as regulators in such networks.
- **Caveats:** the test is small (planned ~200 signatures, got 45 under the rule); the drug folds group by mechanism
  wording, which split tamoxifen from fulvestrant (affects only the ligand report).

## 2026-11-02: hybrid regulons (D87): keep CollecTRI
- **What:** CollecTRI's target lists with each link's direction relearned from LINCS drug profiles: per target, the
  correlation across ~1,000 drugs (per line; pooled for lines with few drugs) with the regulator's activity read from
  its other targets. Ligand drugs were left out of learning. Files in `scripts/hybrid/`, `results/lincs/hybrid/`.
- **Result: slightly worse, rule not met.** On the 828 D84 knockdown / overexpression signatures: 59.8% vs CollecTRI
  62.2% (p = 0.18), same direction in every split (knockdown, overexpression, measured and unmeasured TFs, own and
  pooled directions, raw profiles, all three variants) and in the consensus with ARACNe (65.6% vs 67.0%).
- **Also found:** learning flips about a quarter of CollecTRI's links in every line, so the links' co-movement across
  drugs often disagrees with the literature sign. Those flips do not make the perturbation reading better. The ESR1
  anti-estrogen reading in MCF7 got stronger (7 of 7 negative, median -3.5 vs -2.5), but one regulator does not
  carry the rule.

## 2026-11-03: raw vs line-specific drug profiles (D85): keep line-specific
- **What:** the drug-side ranking (79d, 79e, 79i, consensus regulons) rerun on raw and on line-specific profiles
  (`PROFILE`, `RUN_TAG=_d85`; production outputs untouched), plus the D84 benchmark with each line's generic drug
  response projected out (79n `PROJECT=1`). Decision script `scripts/d85/d85_decide.py`, outputs `results/lincs/d85/`.
- **Result: rule not met.** Raw profiles read TF ground truth no better (+0.4 points, needed +2). The checks pull in
  opposite directions: raw helps ARACNe (+1.8) and hurts CollecTRI (-2.1). On raw profiles 82% of broadly tested
  drugs reverse at least 80% of the diseases, vs 13% on line-specific ones: growth arrest dominates again, which is
  what the projection was added to stop. Known-drug controls are similar.
- **Also found:** the final specificity ranking (`spec_mech`) hardly depends on the choice (median Spearman 0.96).
  ESR1 in MCF7: on line-specific profiles the consensus does not read anti-estrogens as inhibition; ARACNe
  (+2 to +4) and CollecTRI (-2 to -3) cancel to about 0. On raw profiles it leans negative (12 of 15 profiles).

## 2026-11-03: exploration: ESR1 targets within MCF7 under drugs
- `scripts/explore/esr1_mcf7_modules.R`: ESR1's BRCA ARACNe targets across 1,145 MCF7 drug profiles (ligand drugs
  left out). 82% of targets move with the rest of the regulon in ARACNe's direction; luminal-identity genes (GATA3,
  FOXA1, ERBB4, ANXA9) cohere as strongly as estrogen-response genes (GREB1, MYB, XBP1); the two modules correlate
  +0.57 (raw +0.61). Under most drugs the luminal program and ER signalling move together; anti-estrogens are the
  exception that separates them. Co-variation within a line therefore cannot by itself sort identity from signalling
  targets; that needs the regulator's own perturbations (knockdown, overexpression, ligands).

## 2026-11-03: exploration: anti-estrogen time course in MCF7
- `scripts/explore/antier_timecourse.py` (level-5 hiq signatures of 8 anti-estrogens in MCF7): at 6 h the
  estrogen-response genes are already down (mean z -0.88: XBP1, NRIP1, CCND1, GREB1) while the luminal-identity
  genes have not moved (-0.06); by 24 h the luminal genes are up (+1.09: GATA3 +2.0, HEBP1, MINDY1 +2.4) and stay up
  at 48 h (+0.83). Timing separates them: direct ER signalling responds early, and the rise of the luminal genes that
  ARACNe's ESR1 regulon reads as ESR1 activity comes later. Drug profiles (`79c`) currently average all time points.

## 2026-11-03: exploration: drug timing and dose; identity vs drug-response links
- **Time and dose** (`scripts/explore/drug_time_dose.py`, `results/lincs/time_dose/`): 6 h data for 32% of active
  drug x line pairs (3,121 of 9,678); 1,292 pairs (504 drugs) have 6 h and 24 h. Of 2,728 pairs active at either time,
  46% are active at 6 h and 24 h, 45% only at 24 h (late) and 9% only at 6 h. Proteasome, CHK, topoisomerase, HDAC and
  CDK inhibitors respond early; glucocorticoids, MEK, EGFR, tubulin and GSK3 inhibitors mostly late. Dose (24 h, >= 3
  doses, 6,649 pairs): minimum active dose plausible as a potency proxy (bortezomib 0.03 uM, fulvestrant 0.04,
  trametinib 0.12, palbociclib 0.37, tamoxifen ~15); in 918 pairs the high-dose profile is a different program
  (r < 0.5 vs the lowest active dose; e.g. vorinostat, gefitinib). Anti-estrogens in MCF7: the luminal-gene rise is
  absent at 6 h and present at every dose at 24 h (+1.0 from <= 0.1 uM), so it is a slow on-target consequence of ER
  blockade, not an off-target high-dose effect. Time separates it; dose does not.
- **Baseline networks** (`scripts/pnet/depmap_baseline_inputs.py`, `compare_baseline_networks.R`): ARACNe on DepMap
  24Q4 expression (1,588 cancer lines; 69 breast lines), same 978 genes, 81 regulators and settings as the D86 drug
  networks. Overlap over chance: tumor networks with each other 13-17x, with the DepMap baseline 4-6x, with the drug
  networks ~2x; drug networks with each other 3x, with the baseline 1.4x. Tumor co-expression resembles cross-line
  (identity) co-expression; drug networks capture a largely separate structure. Tumor links (BRCA / PRAD / SKCM, on
  landmark genes) also in the baseline: 23-32%; in the line's drug network only: 7-9%; in neither: 54-62%. Against
  the line's knockdown / overexpression signatures (pooled): links in both 64% agree, drug-only 56%, baseline-only
  57%, neither 50%; baseline-only beats neither within signature (p = 0.04). Identity-type links are not the
  unreliable ones; tumor links that no other network supports are. Limits: landmark genes and measured regulators only
  (no ESR1), few signatures per line.

## 2026-11-03: time- and dose-resolved drug profiles (D88): descriptive only
- **What:** drug x line profiles from 6 h vs 24 h signatures (2,854 pairs) and from the lowest active vs the highest
  dose at 24 h (2,759 pairs), scored with tumor ARACNe as the drug side does; checked on the D84 ligand profiles.
  `scripts/timedose/`, `results/lincs/timedose/`.
- **Result:** 6 h profiles read ligands better (+6 points vs 24 h and vs the current profiles) and are less coupled to
  growth arrest, but with 171 rows the gain is not significant (p = 0.18), so the rule is not met. Where the pair is
  active at 6 h the gain is large (+18 points, p = 0.007; reported, not in the rule). Under anti-estrogens in MCF7,
  ARACNe's ESR1 goes from +1.5 to +3.4 at 24 h to about 0 or negative at 6 h. Dose (lowest active vs highest) makes
  no difference to the reading.

## 2026-11-03: early and late-only drug signatures (D89)
- **What:** for 2,854 drug x line pairs with 6 h and 24 h data: early (6 h), late (24 h), late-only (24 h minus what
  it shares with 6 h) and the current aggregated profile, scored with every tumor ARACNe network.
  `scripts/timedose/early_late.R`, `results/lincs/early_late/`.
- **Result:** early and late are different programs (regulator activity correlates 0.37; only 9% of the 24 h signal
  is the 6 h one persisting). The aggregated profile is essentially the late response (0.91 vs 0.58), so the current
  ranking reads mostly secondary effects. Anti-estrogens in MCF7: ESR1, FOXA1, GATA3 and XBP1 read up only in the
  late-only signature; early ESR1 is negative in 6 of 9 profiles.
- **The planned significance test did not work:** it compares each drug's change with other drugs' changes, and since
  every drug changes a lot between 6 and 24 h, almost nothing is unusual (0.01-0.2% of tests). A test against
  replicate noise within the same drug is the right next step.

## 2026-11-03: mechanism profiles in the ranking (D90): keep current
- **What:** the drug ranking (tumor ARACNe only) rebuilt with 6 h profiles at the lowest active dose where a pair is
  active at 6 h, else the lowest active 24 h dose, else the current profile (`scripts/d90/`, `results/lincs/d90/`;
  `PROFILE_FILE` in 79d and `PROFILE_PREFIX` in 79i, defaults unchanged). Literature note: `LITERATURE_TIME_DOSE.md`.
- **Result: rule not met.** Known-use recovery unchanged (Hub AUC 0.60 vs 0.61; the yardstick was already near null).
  Known-drug controls improve (+11 percentile points) and fewer drugs reverse nearly every disease (9.8 vs 12.2%).
  ESR1 in MCF7 is still read as up under most anti-estrogens: most of their 6 h signatures fall below the activity
  cutoff (TAS 0.4), so they used the low-dose 24 h profile, which still shows the late luminal rise. The 2 that used
  6 h profiles read ESR1 down. Dose alone does not fix the reading; time does, where the 6 h profile is used.

## 2026-11-03: exploration: lagged ARACNe across drugs (TimeDelay-ARACNE adapted)
- **What:** TimeDelay-ARACNE (Zoppoli 2010) needs one long time series; LINCS has 2 time points per drug. Adapted:
  samples = the drugs of one line (MCF7 604, PC3 445, A375 325), regulators at 6 h (mrna: 81 measured TFs; act: ~950
  regulator activities from the early signature with the lineage tumor network, incl. ESR1), targets = each landmark
  gene's late change (24 h beyond what its own 6 h value predicts; 86-91% of the 24 h variance). ARACNe-AP as D86.
  `scripts/explore/lagged_aracne_{inputs.py,eval.R}`, `results/lincs/pnet/lagged_*.tsv`.
- **Result: nothing usable.** Lagged links overlap the tumor and DepMap networks at chance (1.0-1.3x) and the
  same-line drug network at 1.7-2.2x. Against knockdown / overexpression signatures every network sits near chance on
  landmark links (lagged 0.46-0.53, drug 0.49-0.51, tumor 0.52-0.56). ESR1 in MCF7: 6 lagged links, all to
  housekeeping / metabolic genes, none to luminal or estrogen-response genes. Across ~600 drugs, early ESR1 activity
  varies for many reasons; the anti-estrogen-specific consequence (8 drugs) cannot dominate a network learned over
  all drugs. Drug-class consequences need per-class contrasts, not one network across drugs.

## 2026-11-03: exploration: how the anti-estrogen misread reaches the ranking
- `scripts/explore/antier_propagation.py`, `results/lincs/early_late/antier_propagation.tsv`: anti-estrogens in MCF7,
  per-regulator contributions to rev (tumor ARACNe, BRCA network). Current (aggregated) profiles: LumA rev -4.4 =
  +6.6 from proliferation-linked regulators (correct direction) and -10.9 from the rest (FOXA1 read +6.5, XBP1 +4.5,
  RARA +5.2, SOX13, CREB3L4: luminal-program regulators read as rising, against their tumor-vs-normal rise); LumB +6.0
  (+11.0 / -5.0); Basal +31.7 (+16.4 / +15.3), where the drug's activity anti-correlates with the basal tumor profile
  (r -0.57). Outcome: anti-estrogens rank at `spec_mech` percentile 96 in Basal, 9 in LumA, 20 in LumB (CUR, `_d85`).
  The late-only signature carries the same pattern (Basal +32.6); the 6 h signature removes it (Basal -0.3) but does
  not reverse LumA either (-8.5, spread over many weak regulators).

## 2026-11-03: gate check (anti-estrogens in basal breast cancer)
- All basal anti-estrogen scores come from MCF7 alone (line match 0.0004, flagged `unmatched_lines`, so tier 3), and
  both gates pass them: `79b` calls ESR1 expressed in Basal (45% of tumors at TPM >= 1, median 0.86 TPM vs 108 in
  LumA), and `79w` requires binding evidence and the receptor in the assay lines (MCF7 has it). Neither asks whether
  the tumor has the target at a level like the assay line. Across diseases, 839 receptor passes (766 cohort x drug)
  rest on a target below 5 TPM and below 5% of its highest cohort (EGFR 218 rows, MET 103, PDGFRA 84, FGFR4 60,
  PPARG 55, FGFR3 52, ESR1 24, ERBB2 22); they rank at a median percentile of 48 (6.5% in the top 10%), so this is a
  hygiene issue for most drugs. A relative-expression rule for the gate would be a separate decision (`79w` is the
  other agent's script).

## 2026-11-03: identity split (D91): keep spec_mech; luminal-only network on the drug side
- **Identity split** (`scripts/d91/`, `results/lincs/d91/`; 79i `IDENTITY_FILE`): regulators whose activity mostly
  separates lineages (DepMap lines) or subtypes (tumors; breast lines) are scored apart, except a drug's own targets.
  Anti-estrogens move up in luminal A / B (positive controls 45 -> 56) and known-use recovery improves slightly, but
  they stay at the 97th percentile in basal and growth-arrest dominance rises (8.2 -> 10.1%), so the rule is not
  met. What keeps basal high is a third component: a quiescence program the drug turns on (MXD4, KLF9, KLF2, NFATC4,
  THRA) that basal tumors have strongly down; it is neither identity- nor proliferation-tagged.
- **Luminal-only network on the drug side** (`scripts/explore/brca_lum_drugside.R`): BRCA_LUM (D65, 769 luminal
  tumors; disease side not materially different, D81) reads anti-estrogens in MCF7 with ESR1 +1.18 (pooled BRCA +2.41,
  positive in 79% vs 100%), FOXA1 / GATA3 / XBP1 about halved, PGR now negative (86%). Subtype networks shrink the
  misread but do not remove it: within luminal tumors the degree of luminal differentiation still co-varies with ESR1.

## 2026-11-03: cell-type breast networks (D92, D92a)
- **What:** networks built per breast cell type before ARACNe, every downstream step rerun under `_d92` / `_d92b` tags
  (`scripts/d92/`; opt-in `NET_MAP`, `MR_OVERRIDE`, `EXTRA_RDS`, `RES_TAG`, `P1_EXTRA` in 50, 79d, 79e, 79h, 79i).
- **D92 (rule as written): keep pooled.** Her2 resembles basal more than luminal, so Basal + Her2 shared a network;
  that put a luminal-vs-basal axis back in, and anti-estrogens stayed at the 94th percentile in basal (needed a
  10-point drop, got 3.5). Luminal subtypes improved a lot (positive controls 21 -> 46).
- **D92a arm (user: > 50 tumors is enough; reported):** Basal (190) and Her2 (81) each on their own network.
  Anti-estrogens fall from the 97th to the 80th percentile in basal and rise in luminal A / B, known-use unchanged;
  dominance +0.2 points (3 drugs), so the strict "does not rise" check is missed by a hair. The anti-estrogen
  inversion is largely corrected: on the basal-only network the luminal regulators barely move under anti-estrogens
  (ESR1 +0.57) and no longer carry the basal disease signature (GATA3 NES +0.18 vs -2.35 pooled).

## 2026-11-03: cell-type networks for all split cohorts (D93): adopted
- 18 networks for the 20 groups of BRCA, CESC, ESCA, LGG, PCPG, SARC, TGCT, UCEC (each group > 50 tumors on its own;
  CESC-AD pooled with UCEC-EEC, PCPG-PGL with PCPG-PHEO). Breast anti-estrogen inversion largely corrected (basal 97
  -> 81, luminal positive controls 21 -> 40), no loss in known-use recovery (41 groups or new cohorts), dominance
  within the 0.5-point tolerance. Splitting LumA / LumB read anti-estrogens a little worse than the shared luminal
  network (positive controls 40 vs 46). Production switch pending.

## 2026-11-03: full-pipeline comparison of network configurations (D94): keep pooled
- Four arms through the whole pipeline (consensus scoring, groups and subgroups, receptor gate, FDA ranking, page,
  known-use, dependency, PRISM / GDSC), each in its own worktree (`../repurposing_d94_<arm>`).
- Cell-type networks help breast clearly (basal anti-estrogens 97 -> 77; luminal controls 26 -> 48 with LumA + LumB
  together, 42 split) and improve known-use specificity and dependency agreement slightly, but drugs that reverse
  nearly every disease rise from 13.2% to 13.9-14.2% (tolerance 0.5 points), mostly through the dedifferentiated
  liposarcoma network (57 tumors, least reproducible). Breast's cell-type networks are less reproducible than its
  very stable pooled network (0.48 vs 0.61), so the uncertainty gate would keep breast pooled anyway.

## 2026-11-03: metaVIPER-style breast (D95) and 85-tumor floor (D96): neither adopted
- D95 helps luminal breast (controls 26 -> 43) but barely moves basal (97 -> 92): combining networks by signal
  strength lets the pooled network's strong anti-estrogen misread dominate. D96 matches the D94 cell-type results
  (basal 97 -> 77, luminal 26 -> 41, ChEMBL known-use 0.500 -> 0.514). Both raise growth-arrest dominance by about
  1 point (151 -> 161-163 drugs). A control rerun on the same network leaves dominance unchanged, so the rise is a
  real effect of the new networks, not noise from rerunning.

## 2026-11-03: production switch to cell-type networks (D97) and queued checks (D98-D100)
- **D97 (user decision):** the 85-tumor cell-type networks are production (`config/group_networks.tsv`, 12 groups).
  Backups `*_pre_d97` (regulon file, drug activity, Step 6 data); last pre-switch results commit `3320b5f`. Rebuilt:
  Step 6 for the 12 groups and 4 subgroups, drug activity, P1 lists, rankings (consensus and TF arm), downstream.
  The 79y2 comparison step first failed (needs `SBX`), resumed. Verified: production equals the D96 test arm on every
  check (basal anti-estrogens 76.5, luminal controls 41.3, known-use 0.615 / 0.514, dominance 14.2%).
- **Written before running, queued after the switch:** D98 selectivity flag for broad reversers at the tier step;
  D99 PROGENy / DoRothEA validation (on-target pathways, ESR1 check, disease contrasts, agreement before vs after
  the switch); D100 early (6 h) and late (24 h) rankings as reported columns, with a test of early profiles as lead.
- Literature notes added: `LITERATURE_NETWORKS.md` (network building), `LITERATURE_SCORING.md` (scoring methods,
  gaps: reversal is largely generic anti-proliferation, Koudijs 2023).

## 2026-11-03: selectivity flag for broad reversers (D98): not adopted
- The 163 drugs that reverse >= 80% of cancers reach tiers 1-2 in 516 cancer x drug rows, but in those rows they
  already rank far above their own typical position (median 36 percentile points above their median), because tiers
  require a top-10% mechanism-specificity score. Requiring 25 points removes only 18 rows; 35 points removes about
  half but loses 4 of 29 known uses (ChEMBL). The existing specificity score already does the selectivity job.

## 2026-11-03: independent pathway check (D99, PROGENy / DoRothEA)
- Drug profiles carry on-target signal: kinase and pathway inhibitors move their pathway as expected in 75% of
  1,446 profiles (PI3K 91%, MAPK 86%, p53 81%, EGFR 68%; JAK-STAT and NFkB at chance). Timing makes no difference
  here (6 h 74% vs 24 h 73%); the lowest active dose reads best (81%).
- PROGENy reads anti-estrogens as blocking estrogen signalling in all 15 MCF7 profiles, while our ARACNe ESR1 reads
  them as activating on every breast network, including the new luminal ones. The ranking improved through the
  disease side, not through a fixed ESR1 reading.
- Disease pathways match known biology in 4 of 5 checks (thyroid MAPK fails). Our lead score agrees only weakly with
  a pathway-based ranking (0.17), slightly more after the switch than before; tier 1-2 drugs reverse pathways a bit
  more than average.

## 2026-11-03: early (6 h) and late (24 h) rankings (D100)
- Both rankings run through production for the 1,123 drugs with both time points. Ranking on early profiles loses
  known-use recovery (Hub AUC 0.61 -> 0.51) and slightly worsens the breast controls, though it cuts broad
  growth-arrest hits (14.2% -> 11.1%); the late ranking is essentially the current one. Aggregated stays the lead.
- Each cancer x drug row now has early and late scores and a timing class (both 4,668; early only 1,871; late only
  2,051; neither 19,408; unknown 22,762) in `results/lincs/d100/timing_scores.tsv`, to be shown on the rankings page.

## 2026-11-03: PROGENy reading for ESR1 / AR (D101), rankings page timing column
- D101: reading ESR1 and AR on the drug side from PROGENy (correct for anti-estrogens) moves anti-estrogens by under
  one percentile point; not adopted. A single receptor carries little weight among ~1,000 regulators; the luminal
  program as a whole drives the anti-estrogen rank.
- Rankings page (`79x`, template): new Timing column (early + late, early only, late only, neither, no 6 h data) with
  early / late percentiles in the tooltip and the score table, a legend entry, and an updated validation note
  (4,483 of 9,823 FDA rows have a timing class). Published as a new page: https://claude.ai/artifact/DkoreTiuHXRkqaEYMG5hfk

## 2026-11-04: relative target-presence gate (D102), breast stage
- Rule (D102): a drug whose ChEMBL-backed targets (79w binding rule) are all absent in a group (median TPM < 1 or
  < 10% of the gene's highest group median) is ranked below every other drug in that group; scores are unchanged.
  Post-ranking step on production (`scripts/d102/d102_gate.py`, `STAGE=1`); 79b / 79w untouched.
- Breast: Basal anti-estrogens 76.5 -> 3.3 (fulvestrant, raloxifene to the bottom; tamoxifen stays at 87 because its
  potent target EBP is expressed), LumA / LumB controls 41.3 -> 43.8, known-use AUC Hub 0.614 -> 0.637, ChEMBL
  0.494 -> 0.488. Stopped by B4: gefitinib (EGFR 1-3 TPM in luminal) and metoprolol (beta receptors < 1 TPM), both
  ChEMBL trial labels for breast cancer, are demoted in LumA / LumB. 53-72 drugs demoted per breast group.
- Issues for an amendment: amplified cohorts set the 10% reference (ERBB2 51 TPM in Basal counts absent vs 1,271
  in Her2); tamoxifen escapes on EBP; anti-estrogens also demoted in Her2 (ESR1 2.1 TPM).
- D102a (amended: 90th-percentile reference, curated mechanism targets first, only approved-use labels stop):
  breast passes (Basal anti-estrogens 76.5 -> 3.0 incl. tamoxifen 3.4; luminal controls 41.3 -> 43.6; Hub AUC
  0.614 -> 0.640) and all 41 groups pass (Hub 0.614 -> 0.617, ChEMBL 0.514 -> 0.530, PRAD controls not lower;
  median 47 drugs demoted per group). Caveat: ALK inhibitors (alectinib, ceritinib) are demoted in lung cancer
  because ALK fusions occur in a few percent of tumors and the median is ~0.15 TPM. Not in production yet.
- D102b (carrier exception: target counts if >= 2% of tumors express it highly) fails both ways: it brings back
  ESR1 in Basal (6 of 190 tumors ER-high, 3.2%) and does not rescue ALK in LUAD (3 of 513 tumors, 0.6%). No
  threshold on expression separates them; a curated biomarker list (fusion targets) would be needed.

## Assumptions register
Moved to [`DECISIONS.md`](DECISIONS.md), section 2 (A1-A43, with current status). The list that used to be here
(numbered 1-16) is in this file's git history.

## Next
Status 2026-10-31. Full list with owners: [`DECISIONS.md`](DECISIONS.md) section 3; step status: [`PIPELINE_V2.md`](PIPELINE_V2.md).
- Answer the Step 7 questions (D18 confirmed rule; TF-only edge confounding), then rerun `86` (variability) for the D79
  subgroups, then Step 9 replication (cohorts to choose; COAD can start with GSE39582 / CPTAC-COAD).
- PI: D26 (recurring module), D60 (kidney reference). User: A37 (LGG, tabled). BRCA_LUM test (D81) running.
- Decide on signaling-only regulons (D82 follow-up); implement lineage-matched scoring (C24).
- Drug side: PRISM / GDSC validation on hold; check the recurring oddities in specificity lists.
- Later: CPTAC-HNSCC as an HNSC replication cohort (GDC CPTAC-3, open STAR counts, GENCODE v36; do not pool with Toil);
  subtype splits for mixed cohorts; recount3 non-cancer diseases.

## 2026-11-02 Receptor evidence gate (79w, C26) and full FDA rankings page (79x)

- Receptor targets of the returned drugs (HGNC receptor families): 251 of 868 expressed targets. Many were annotation
  noise or real binders absent from the cell lines. New gate (thresholds by the user: 1 uM, 100 nM for kinases, receptor
  in >= half of the cohort's matched lines; single-line passes kept and flagged): 22,582 cohort x drug x receptor rows,
  31% pass. Selectivity window of passing pairs median 0.48 log10 (66% within 10x of the drug's best target).
- `results/lincs/disease_drug_rankings_v2.html`: 459 FDA drugs x 41 diseases + 14 subgroups (9,823 rows), all scores,
  flags, receptor gate, literature where 79o searched (225 rows). Runbook steps `receptor_gate`, `rankings_page`.

## 2026-11-02 Closest-compound stand-ins (79y, 79y2, C27) and stand-in page (79z)

- Stand-in = compound tested in a cohort's own lines that shares a target with an FDA drug and changes CollecTRI
  regulator activity alike (r >= 0.4 over >= 2 shared lines; rule fixed before looking). Validation on 503 FDA drugs
  tested in their own lines (similarity measured outside them): Spearman 0.66 vs 0.34 for random target-sharing
  compounds; adopted. Against the drug's own borrowed score (counterfactual rescoring in a sandbox, 485 drugs):
  stand-in 0.65, borrowed 0.56, average 0.68; stand-in better by 0.09 (CI 0.02-0.16).
- Coverage: 443 of 8,419 borrowed / weak FDA rows have a stand-in, in 20 cohorts (none possible in the 27
  borrowed-only cohorts). 1,922 matched compounds have a close FDA drug: `results/lincs/standin_compounds_v2.html`.
- 79o `MODE=fda_all`: PubMed / trial search for every FDA row of the rankings page (7,444 searches, cached).
- 79e: guard for `rev_raw` in consensus mode (raw activity folders for v2 now exist and crashed the rerun).
- 2026-11-03: validation step 1 (protocol pre-registered as C28 before running). (a) Known-use recovery (`79l`,
  `config/known_use.yaml`, Repurposing Hub + ChEMBL 37 indications): lead-score AUC 0.62 on Hub labels but 0.50 on
  ChEMBL, and the "drug scores high in every cohort" baseline matches it (0.62 / 0.60); the specificity test (own
  cancer vs elsewhere) is null (+3.0 points Hub p 0.19; -5.1 ChEMBL p 0.98). Pass rule not met. (b) PRISM / GDSC
  selective killing against the gated v2 rankings (`79m`): receptor gate pass vs fail inconsistent across datasets
  (-0.04 / +0.04 / -0.11), matched evidence no better than borrowed, stand-in average never better than the drug's
  own borrowed score. Pass rule not met -> C16 stands after the gate, FDA filter, evidence split and stand-ins.
  Next: DepMap CRISPR dependency agreement (step 2), binding evidence beyond receptors.
- 2026-11-03: validation step 2 (protocol pre-registered as C29), DepMap 24Q4 CRISPR gene effect added to `05`.
  `79n`: the yardstick holds (AR / PRAD 7.6 SD, ESR1 / LumB 4.9, ERBB2 / Her2 3.6, BRAF / SKCM 2.9; EGFR / LUAD and
  KIT / STAD do not, the line sets lack the mutations). D1 lead score vs target dependency null (rho -0.04);
  D2 the receptor gate carries no dependency information (-0.05 SD, within noise of the 1-8 SD real scale; the split
  tracks expression); D3 positive - drugs approved for a cancer target genes that cancer depends on (+0.24 Hub /
  +0.13 ChEMBL SD, p 0.001), which also shows the C28 labels are sound. 154 rows pair an approved drug with a
  target >= 2 SD depended on and the ranking does not favour them (median percentile 60 vs 53).
