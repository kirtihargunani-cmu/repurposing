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
- **D102c: D102a adopted (user: a subtype's ranking should describe one cohesive group).** New step `79i2` re-ranks
  the production split tables; known-use Hub 0.614 -> 0.617, ChEMBL 0.514 -> 0.531; page republished with the flags
  "Target absent" and "Minority target" (e.g. ALK inhibitors in lung adenocarcinoma).

## 2026-11-04: downstream refresh after the target gate (D102c)
- All outputs that read the lead percentile rerun on the gated ranking (FDA lists, literature, stand-ins and their
  counterfactual, known-use, dependency, PRISM / GDSC sensitivity, the three pages), plus the 79i summaries and the
  D100 timing labels (gate applied to the early / late rankings; D100 decision unchanged).
- Stand-in agreement 0.555 -> 0.593 (random shared-target pairs 0.316 -> 0.422, still below); dependency check D1
  turns slightly positive (matched 0.004 -> 0.056, p 0.02); PRISM / GDSC still no signal.

## 2026-11-04: receptor level in tumor vs assay lines, sizing for D104 (descriptive)
- `scripts/d104/d104_size.py` -> `results/lincs/d104/size_{rows,summary}.tsv`. Per cohort x drug x bound receptor (79w
  mechanism / potent), tumor median TPM vs the weighted mean DepMap TPM of the scoring lines (79e / 79w line rule).
  Provisional classes: line only (line >= 1 TPM, tumor < 10% of line), tumor only (tumor >= 1, line < 1).
- Poorly matched pairs (borrowed or unmatched lines; 4,182 cohort x drug with a bound receptor): 318 have every bound
  receptor line-only (221 already demoted by D102c; 7 of the rest in the top 10%, e.g. cilengitide in PCPG, ITGB3 ~1
  vs 27 TPM; pioglitazone in PCPG, PPARG 3 vs 50); 1,538 have >= 1 receptor the tumor makes but the lines do not
  (median percentile 52, 100 in the top 10%; e.g. antipsychotics in PHEO, DRD2 101 TPM vs < 1 in the lines; genistein
  in luminal breast scored on ER-negative lines). Matched pairs: 15 line-only, 352 tumor-only. Raw TPM scales (Toil vs
  DepMap) are not aligned.

## 2026-11-04: receptor level tumor vs scoring lines (D104), not adopted
- Rule: drugs acting only through receptors are demoted where the scoring lines lack the receptor (unless the drug
  induces it) or the tumor has < 10% of the lines' level. Known-use holds (Hub 0.615 vs 0.617), but prostate
  anti-androgens fall 86.7 -> 5.6: LNCaP makes ~20x more AR than the average prostate tumor. Anti-angiogenic drugs
  (tivozanib, pazopanib in kidney) are demoted because cancer lines lack VEGF receptors, which sit on tumor vessels.
- 928 of 3,587 judged pairs would be demoted; 904 receptor rows rescued because the drug raises the receptor.

## 2026-11-04: D104a (receptor in scoring lines, amended): criteria met
- Dropped the tumor-vs-line abundance test; receptors whose tumor expression rises with stroma / immune content
  (ESTIMATE rho >= 0.3) are treated as tumor-environment targets (flagged, not demoted). 546 of 3,587 judged pairs
  demoted, 286 flagged; known-use Hub 0.622 vs 0.617, prostate controls 86.8 vs 86.7, breast 44.3 vs 43.6.
- Gap: tivozanib in KIRC / KIRP still demoted (ESTIMATE does not capture tumor vessels). Awaiting go-ahead.

## 2026-11-04: D104b (adds a blood-vessel signal): criteria met
- Receptors that track an endothelial score (PECAM1, VWF, CDH5) also count as tumor environment. Tivozanib in kidney
  is rescued (VEGF receptors rho 0.55-0.90); no approved-use drug demoted; known-use Hub 0.621, ChEMBL 0.540.
  The signal is broad (39% of receptor rows reach the cut), so it also rescues some non-vascular cases (NVP-TAE684
  via INSR in lung). 338 pairs demoted, 494 flagged. Awaiting the user's choice between D104a and D104b.

## 2026-11-04: D104c in production (receptor in the scoring lines, rho 0.4)
- `79i3` after `79i2`; 406 demoted, 426 tumor-environment flags (groups). Downstream rerun: known-use Hub 0.623,
  ChEMBL 0.542; dependency D1 on matched lines 0.079 (p 0.002); stand-ins 0.605. Page republished.

## 2026-11-04: D105 (ChEMBL primary targets for every drug): not adopted
- ChEMBL 37 fetched for all 2,698 ranked drugs (1,835 with an id). Known-use Hub 0.642 vs 0.623, toremifene now
  caught in Basal, but criterion 2 fails mechanically (Basal controls 3.0 -> 6.8 as the demoted block grew by 2,694
  rows). Problems: mechanisms stored on salt forms missed (irinotecan, doxorubicin), non-protein mechanisms (DNA
  alkylators) fall back to off-target potency hits, and the relative presence cut demotes proliferation enzymes in
  slow tumors (TOP2A in LGG, TYMS in KICH). Amendment needed.

## 2026-11-04: D105a (FDA drugs, tiered main targets): criteria met
- Curated mechanisms found for 291 of 413 FDA drugs once salt-form records are included; 81 act on DNA or pathogen
  targets and are not gated on human proteins; 37 rely on annotation. Known-use Hub 0.623 -> 0.660, ChEMBL 0.553;
  toremifene caught in Basal; other compounds unchanged. Secondary-target list needs outlier handling (max of single
  measurements). Awaiting go-ahead.

## 2026-11-04: D105b in production (FDA main targets), target-present page
- Production equals the D105a test (Hub 0.660, ChEMBL 0.553; Basal anti-estrogens 2.9). Fixed: 79p / 79o ignored
  the D102c / D104c demotions. Literature search held for the user's review of the pairs (19 new of 611).
- Rankings page: "Target in tumor" column, main-target source, secondary targets. New page with target-present drugs
  only: https://claude.ai/artifact/WXbtdXQzvcxF6fRX35fKSK

## 2026-11-04: target-present page tightened (D105c)
- Presence must rest on measured tumor expression (no rows lost); unverified binding hidden by default with the
  weak-evidence rows (5,173 of 6,864 rows visible); percentiles re-ranked among the page's FDA drugs per cohort.

## 2026-11-07: audit and D110 (fusion targets)
- Full audit of the drug side (5 reviewers: core scoring, gates / flags, validation / outputs, docs vs results,
  methods). Core scores recomputed independently and match production. Main findings: the known-use AUC mostly
  measures drug class (an anticancer yes / no label scores as well; within anticancer drugs the lead score is at
  chance); ALK fusion-partner bug in `79i2`; Cmax lookup, salt list and compound-pick bugs; the runbook misses the
  gate and flag steps; several doc numbers out of date. Fixes are being worked through in order (D110 on).
- D110 in production: a ChEMBL chimeric (fusion) mechanism target counts only through its kinase (Bcr/Abl -> ABL1,
  EML4-ALK / NPM-ALK -> ALK). Alectinib, ceritinib, crizotinib now target-absent where the typical tumor lacks ALK
  (+67 group, +20 subgroup rows); FDA lists 2,426 -> 2,381 rows; two subgroup top-10 changes. Known-use Hub
  0.660 -> 0.631 (reported only).
- D111 in production: the exposure flag (D109b) recovers 24 verified Cmax values lost to a molecular-weight lookup and
  unit spellings (147 of 384 FDA drugs now), and applies Cmax to every LINCS compound of an ingredient. low_exposure
  459 -> 513 FDA rows (52 drugs); bortezomib leaves the PAAD and PCPG-PGL top 10, everolimus COAD-MSI-like.
- D112 in production: Drugs@FDA matching adds 5 salt / ester words, 4 griseofulvin size labels and 4 INN / USAN
  aliases, and drops one wrong match (MRE-269); 542 compounds matched (was 527). New FDA drugs include lapatinib,
  rucaparib, cyclosporine, thioguanine. Missing ChEMBL mechanisms for newly matched drugs fetched (lapatinib now
  EGFR + ERBB2; #6 in BRCA-Her2). Cyclosporine enters 10 top-10 lists without a Cmax reading. Open: 43 FDA drugs
  carry a wrong ChEMBL id (no curated mechanism), proposed as D113.
- D113 in production (11-08): FDA drugs without a ChEMBL mechanism record get their ChEMBL id from the FDA ingredient
  name (44 resolved, 32 applied). FDA drugs on a curated mechanism 295 -> 325 of 422; 'unverified binding' 35 -> 8.
  Tioguanine -> IMPDH1/2 (was annotated as thioguanosine), vincristine -> tubulins, cyclosporine -> PPIA. Meloxicam
  (COX-2 only) now target-absent in 16 groups / 7 subgroups. No top-10 list changed.
- D114 in production (11-08): Cmax read for the FDA drugs added under D112 (cyclosporine 1.5 uM, lapatinib 4.2,
  rucaparib 6.0, griseofulvin 5.7, pentamidine 1.8; 152 drugs now). None reaches the 10x low-exposure cut; no
  ranking change. Thioguanine and clomiphene remain unassessed (no Cmax in label or FDA documents). 79f5 now keeps
  earlier queue entries on a rerun (it would have dropped 50 verified document readings).
- D115 in production (11-08): one eligibility test (`fda_rows.eligible`) for the FDA lists and the pages; pages,
  literature and evaluations now use the same LINCS compound per drug as the FDA lists (110 rows had shown a
  borrowed or flagged sibling). FDA lists unchanged; known-use Hub 0.635, ChEMBL 0.554 (reported only).
- D116 in production (11-08): 'inactive on matched cells' is judged over all compounds of a drug. 132 rows newly
  flagged, 118 cleared (a sibling is active there). Osimertinib enters the LUAD top 10 (temsirolimus out); also PAAD,
  PRAD, SKCM changes. 26 new trials read for the new pairs. Evaluations unchanged.
- D117 in production (11-08): the runbook now runs every gate and flag step in a fixed order (fda_matches ->
  target_presence -> receptor_lines -> evidence_flags -> plate_confined -> exposure -> tiers); a rerun of the chain
  reproduces the FDA lists byte for byte. A flags-only rerun now stops instead of dropping flags. Tiers and evidence
  status are recomputed after all gates (evidence status unchanged everywhere). Fixed a hash-seed dependence in the
  79i2 carrier count.
- D118 (11-08): known-use metric redesigned (drug-level labels and scores, cancer-word regexes, subtype table, primary
  = pooled AUC among eligible anticancer drugs, next to the drug-class baseline). Result: Hub 0.462, ChEMBL 0.430
  (chance 0.5; class baseline 0.66). Within cancer drugs the lead score does not recover each cancer's own drugs; the
  earlier 0.63 came from the gates and from drugs that score high everywhere.
- D119 (11-08): scoring corrections (subgroup regulators in 79d, consistent reference sets, one percentile definition,
  NES-weighted prolif_share, literature counts in the FDA lists, 79k on v2, stop-reason negations). Scores move by
  Spearman >= 0.993; 8 top-10 lists change by one drug. Known use unchanged (0.462 / 0.429); PRISM / GDSC on v2 null.
  Literature lookups are now cached for every mode (a rerun fetches only new pairs).
- D120 (11-08): target-role score, reported beside the ranking (79t). Per cancer x drug: is the drug's main target
  needed by the matched lines (D, CRISPR), altered in the tumors (A, TCGA hotspot mutations / amplifications), or
  turned up (O, regulator activity or expression vs normal), signed by the drug's action; each part shown. Known use:
  Hub 0.674, ChEMBL 0.617 (spec on the same rows 0.489 / 0.474); split halves agree. Roughly half the signal is
  drug-level (0.586 / 0.572 with no cancer information). It rises with the number of targets (best-of-N) and does not
  work for drugs whose target every cell needs (0.46 / 0.44); 194 of 410 top-10 rows rest on minority alterations.
- D121 (11-08): three-layer score (layer 1 target role, layer 2 units / direction / binding, layer 3 mix with
  reversal, weights learned on Hub, tested on ChEMBL minus Hub's pairs). Many-target fix: mean over target units
  (passes the label-free check). Learned weights (alteration-heavy) won on Hub held-out (0.669 vs 0.643) but lost
  on ChEMBL (0.563 vs 0.615) and on the breast controls: baseline (equal weights) stays. Reversal and relative
  binding add nothing; occupancy at Cmax helps where a Cmax exists. Spec doc with flowchart (Claude Docs).
- D121a (11-08, user): lab binding strength stays in the score (gamma 1) even though it does not raise accuracy;
  later fits keep it fixed.
- D122 (11-08): **the layered score is the lead ranking** (groups and the 14 subgroups, subgroup layer 1 on their own
  lines, tumors and normals). Gates: targets expressed, target absent, local route, low exposure; LINCS evidence and
  no reversal are flags; drugs without a protein target listed separately. Known use Hub 0.623, ChEMBL 0.629 (the
  reversal lead on the same rows 0.49 / 0.48); anti-estrogens high in luminal, low in basal. Top 10s almost entirely
  new (median 1 of 10 kept); they cluster by target class (glucocorticoids, statins, glitazones, anthracyclines, EGFR,
  BRAF / MEK). 457 more trials read for the new top 10s. PRISM / GDSC still null.
- D123 (11-08): free fraction read from FDA labels for 196 drugs, 8 new Cmax; free occupancy for 86 drugs (24 below
  50%). It does not beat relative binding on known use; reported, not used.
- Check (11-08, reported, post hoc; `scripts/checks/essential_check.py`): drugs on targets every cell needs
  (chemotherapy). Within a cancer no component separates the approved ones (0.43-0.52). Between cancers, per drug,
  the cancer's proliferation index separates where they are approved (Hub 0.63, ChEMBL 0.57, p <= 0.005);
  alteration points the wrong way (0.37 / 0.38); SLFN11 suggestive for DNA-damaging and nucleotide drugs (few drugs).
  Basis for a possible rule: for these targets drop A, use proliferation for O.
- D124 / D124a / D124b (11-08): targets every cell needs (most genes DepMap common essential) and 13 DNA-acting drugs are
  scored on the cancer's proliferation (PAM50 index), weighted by how closely each target tracks proliferation in
  tumors (rho: TOP2A 0.79, RRM2 0.77, HMGCR 0.29). Known use Hub 0.673, ChEMBL 0.652 (was 0.623 / 0.629); statins and
  cyclosporine leave the top 10s. The tie criterion failed on shared-target ties (glucocorticoids), which the previous
  objective had too; adopted by user decision. In production.

## 2026-11-10: D126 closed — the mechanism-consistency filter fails its own known-failure control
- **B1 ran on all 55 cohorts** (`scripts/d126/b1_consistency.py`; drug activity dumped from 79e itself through a new
  off-by-default `D126_ACT_DIR`, so the vector is production's by construction). 5,914 (drug, cohort) rows, 203 drugs
  with a curated main target and a known action direction.
- **B5.1 fails and D126 is closed, as the spec says it must be.** The flag fires on **0 of 22** scorable
  anti-estrogen rows in the four MCF7 breast cohorts, and on only 7 of 22 even with the coverage precondition
  removed. It does worse than stay silent: in the cohorts where the ESR1 misread lives it calls tamoxifen,
  toremifene and raloxifene **1.00 consistent** in LumA, 0.83 in LumB.
- **Why, and this is the part worth keeping.** ESR1 is reached in **none** of the 22 control rows, because it is
  never among the top 25 regulators by |activity| — those are FOXM1, MXD3, PA2G4, E2F2, MYC, TEAD4, the
  proliferation block. B3's selection rule aims the mechanism check at growth arrest rather than at the drug's own
  axis, and an ESR1 antagonist genuinely does explain lowering FOXM1 and MYC, so the check passes it. The
  proliferation dominance D80 and D85 found in the *score* turns up again inside the diagnostic built to audit it.
- **B5.3 fails by a wide margin:** `mechanism_coverage` median **0.08** (2 of 25 regulators), max 0.48; 37 of 5,914
  rows reach the 0.4 floor and no drug reaches it on its median. The curated signed prior does not connect drug
  targets to tumor-network TFs — this was visible in B0 (46% of targets have an outgoing edge) and is worse in
  practice. **B5.2 passes** (rho 0.165 vs main-target out-degree), so B4 is not triggered, but with coverage this
  low the degree check has little to bite on.
- **Run-order B2 (CARNIVAL, breast): not executable as configured** — no CARNIVAL package, no Gurobi, no CPLEX, no
  highspy (`scipy.optimize.milp` exists but CARNIVAL does not interface with it). Recorded as not executable per
  B3 rather than substituting the heuristic. Moot regardless: its only job was to license a heuristic that has
  failed its own control.
- Nothing in production changed; the flag was reported-only by design and is not emitted. **With Part A closed at
  A1 and Part B closed here, all of `SPEC_graph_metric_scoring.md` is resolved.** A revived D126 would need a
  different regulator-selection rule — the drug's target-proximal regulators rather than its 25 largest shifts —
  which is a new decision with its own control.
## 2026-11-12: precursor — gene-level reversal inside the ESR1 regulon does NOT recover the anti-estrogen axis
- Reported only. `scripts/precursor/{export_regulons.R,esr1_regulon_reversal.py}`, outputs
  `results/lincs/diagnostics/precursor_esr1_regulon/`. No threshold, no pass/fail, no D-number, nothing written to
  production. Premise under test: D126-B1 found anti-estrogens are never examined on their own axis because ESR1 is
  not among a drug's 25 largest regulator shifts, and the dimension diagnostic put ESR1 at median rank #744 on PC1.
  Does keeping the regulon's per-target weights and scoring on its own GENES recover it?
- Quantity, one level below production's `rev = -sum_i w_i a_i / sqrt(sum_i w_i^2)`: `w_ij = |tfmode_ij| *
  likelihood_ij * d_j` over regulator i's targets, `rev_i = -(sum_j w_ij g_kj) / sqrt(sum_j w_ij^2)`. `d` is the
  Step 4 tumor-vs-GTEx z; `g` is the 79c line-specific profile consensused over the cohort's lines exactly as 79e
  builds drug activity. (The signed construction is the same product: `tfmode` cancels in `d_j * g_j`, leaving only
  `tfmode` vs `tfmode^2` as the weight.) Four breast cohorts, 1,511 drugs, ARACNe subtype networks and CollecTRI.
- **Polarity verified first:** the same construction over ALL genes reproduces production's regulator-level `rev` at
  Spearman 0.74-0.80 per cohort, so the restricted numbers are a finding and not an orientation bug.
- **ARACNe arm: inverted.** ESR1-restricted reversal puts the anti-estrogens at the BOTTOM of the cohort's drugs in
  the three ER-relevant cohorts - ESR1 ranks **1082 of 1082** regulons in BRCA-Her2, **1124 of 1134** in LumB, **930
  of 1093** in LumA by how highly it ranks them (median percentile 0.137). It ranks them high only in **BRCA-Basal**
  (ESR1 rank 55 of 1146, mean percentile 0.87) - the ER-negative cohort where 79e's own NEG_CONTROLS say they should
  not be specific. Production meanwhile gives these drugs rev_pct 78-99 in Her2 / LumA / LumB.
- **CollecTRI arm: no inversion, but no lift.** ESR1 is an unremarkable regulon - rank 241 / 484 (Her2), 305 (LumA),
  194 (LumB), median percentile 0.544 - and again its best cohort is BRCA-Basal (rank 55).
- **The regulons that do rank anti-estrogens highest are incoherent** (KDM2B, NKRF, ZNF132, MYC, SSRP1, POU1F1), so
  this is not "a different regulon found it".
- **Reading: this is a reference-choice problem, not a resolution problem.** `d` is tumor vs GTEx normal breast, and
  normal breast is ER-active, so ESR1's targets do not read as up-in-tumor at gene level; an anti-estrogen's push on
  them agrees with the tumor-vs-normal direction instead of opposing it. The aggregate msviper NES is a rank
  enrichment over the whole ranking, not a weighted dot product on the same genes, which is why the sign differs
  when you drop a level. This touches the open item on using TCGA adjacent normals in drug ranking (D48, D57-D60).
- **Not covered:** the raw (non-line-specific) profile arm; estradiol as an agonist anti-control is not among the
  scored drugs; one weighting only. A spec for regulon-level scoring would have to answer the reference question
  first - going finer exposes it rather than fixing it.

## 2026-11-12: drug-side dimension diagnostic — the drug space is ~5 effective dimensions, and its leading axis is potency
- Reported only. `scripts/diagnostics/drug_space_dimension.py`, outputs `results/lincs/diagnostics/`. No threshold, no
  pass/fail, no D-number; nothing in production is read or written. Input is the 79e consensus drug activity dumped
  by D126-B1, i.e. the quantity the score actually uses. Method mirrors `86_subgroup_variability.R` so the drug-side
  and patient-side numbers are comparable (participation ratio of the PCA eigenvalues, regulators centered).
- Two arms, because cohorts differ in drug count (1,144-1,871) and regulator count (1,003-1,223) and both shift a
  spectrum on their own. **1,144 pert_ids and 571 regulators are present in all 55 cohorts**, so the comparable arm
  is an exact common frame — no subsampling and no RNG, which 86 could not manage on the patient side.
- **Effective dimension: median 5.20 (4.09-8.34) in the common frame, 6.04 in each cohort's own space.** PC1 takes a
  median 39% of the variance; 2 PCs reach half, 8 reach 80%, 23 reach 90%. For comparison the patient side is 7.99
  (6.51 at fixed n=30). **The drug space is more collapsed than the tumor space, not less.**
- **The cohorts largely share one leading axis:** pairwise |cos| between cohorts' PC1 loadings has median 0.82, min
  0.64, and 68% of pairs exceed 0.8. It is not cohort-specific biology.
- **And that axis is mostly potency, not direction.** |PC1 score| against the drug's own activity magnitude:
  Spearman **0.685** median in own space, **0.717** in the common frame (range 0.54-0.80). Its top loadings are broad
  chromatin / zinc-finger regulators — PRDM10 (top 25 in 51 of 55 cohorts), ZFP91 (48), ADNP (48), ZBTB11 (46),
  ZBTB1 (42), RLF (40) — **not** the E2F / FOXM1 proliferation block that dominates per-drug |activity| (D126-B1).
- **The number that settles the condensing question: ESR1 ranks median #744 on PC1 by |loading|** (33 cohorts score
  it; best #280). In breast: BRCA-Basal #663/1052, BRCA-Her2 #686/1003, BRCA-LumA #867/1038, BRCA-LumB #937/1054.
  Truncating to the top PCs would delete the target-specific axes outright, and whitening is A1's small-lambda arm,
  which lost monotonically. So the space condenses easily and condensing it for scoring remains the wrong move.
- **What this does open, and it is new:** PC1 looks like a nuisance amplitude axis shared across cohorts rather than
  a biological direction. Removing it is confound removal, not dimension reduction, and is a different proposal from
  D125. It would need its own D-number, the D84 benchmark, and the 22-row anti-estrogen control, because drug potency
  may be real efficacy signal — this diagnostic does not show that removing it helps.

- **Plain-language write-up (2026-11-11):** `docs/REPORT_graph_metric_scoring.html` — a standalone page covering the
  proposed scoring function, the run order, and both failures, with charts drawn from the committed result tables
  (lambda sweep and gain vs the +1.0 bar; the OmniPath filter; path lengths; the B5.1 control scatter with its empty
  flag box; top-25 regulator frequency with ESR1 at zero; the coverage histogram against the 0.4 floor). Also
  published privately at https://claude.ai/artifact/4cXy1rRbMjZE2Y3Kw2qb5K. Report only; no numbers are produced here.

## 2026-11-10: D126 B0 — the signed signaling prior is pulled, filtered, cached and registered
- `scripts/d126/b0_omnipath.py`. **Primary prior: 16,518 signed and directed edges over 5,493 nodes** from 86,189 raw
  OmniPath interactions (literature-curated core, organism 9606, retrieved 2026-11-10); 11,662 activating, 4,856
  inhibiting, 1,199 `sign_disputed`. **SIGNOR sensitivity arm: 10,346 edges over 4,364 nodes.** Cached under
  `data/raw/omnipath/`, checksummed in `data/raw/checksums.sha256`, registered in `config/datasets.yaml`.
- **The open item that blocked B0 is settled (D126-B0):** OmniPath gets its own registry entry, with a `review:` key
  recording that the service publishes no data release label (`/about` returns the server build, 0.1.5), that the
  cached file is the dependency, and that a re-pull is a new decision rather than a refresh.
- **A spec field name was wrong, found by querying the service first as the risks section asked:** there is no
  `consensus_sign` field. The service returns `consensus_stimulation` and `consensus_inhibition`, so B2's condition
  is read as exactly one of the two being true. The intended filter is unchanged; the spec is annotated in place.
- **The support clause, not the sign clause, defines this prior.** 14,939 raw edges fail `consensus_direction` and
  15,073 carry no consensus sign, but **54,598 are dropped for resting on a single resource and a single reference**.
  About two thirds of OmniPath's curated core is single-citation.
- **B2's reason for preferring OmniPath to SIGNOR is weaker than stated.** It chose OmniPath because "coverage is
  the binding constraint"; measured, SIGNOR keeps 63% of the edges and 79% of the nodes. The sensitivity arm is a
  near-equal alternative. The primary arm is left as written.
- **Reported, not acted on:** of the 387 D105b evidence targets on scored FDA rows, 67% are in the prior and only
  **46% have an outgoing signed edge**, which is what a path from the target needs. B5.3's coverage floor is at real
  risk before B1 runs. The B5.1 control is computable: ESR1 and the other control targets are all source nodes.

## 2026-11-10: D125 graph-metric reversal — not met, A2 does not run (A1 early stop)
- The question (D125): `rev` compares the disease and drug regulator vectors as a plain signed dot product, which
  treats regulators as independent when D91 found 7-30% of a network's regulators identity-linked. The proposal
  replaced the identity metric with a graph one, `rev_L = -(d'Kg)/sqrt(d'Kd)` with `K = (L + lambda I)^-1` and `L`
  the normalized Laplacian of the regulon-overlap graph. `K -> I/lambda` as `lambda -> inf`, so the production score
  is nested inside the proposal and only `lambda` had to be chosen.
- **A1 (`scripts/d125/a1_sweep.py`, outputs `results/lincs/d125/a1/`): the sweep says no.** On the D84 perturbation
  benchmark (1,774 signatures, 329 TFs, consensus read) the correct-direction rate is 0.6077 at the `lambda = inf`
  limit and never clears it: 0.573 / 0.596 / 0.604 / 0.605 / **0.609** / 0.608 / 0.607 at lambda 0.1 to 8. Best gain
  **+0.17 points at lambda = 2, against the pre-registered bar of +1.0**, so the spec's early stop fires: D125 is
  recorded as not met and the full-pipeline arm A2 never runs.
- The shape matters more than the verdict. The rate falls monotonically as the graph gains weight, and it falls
  hardest where the benchmark is cleanest: overexpression loses 10.3 points at lambda = 0.1 and is still below the
  limit at lambda = 8. CollecTRI's optimum *is* the limit (+0.00). Per TF, 18 better and 16 worse out of 131.
- The redundancy is real -- 42 graphs, 472-1,266 regulators, median degree 50 after the 0.1 edge cut -- so this is
  not "no overlap to remove". Down-weighting shared-target evidence costs more signal than it removes noise, the
  same way the pnet work (D86) and the latent-factor attempts failed, now on a third kind of input.
- **D125-A1.1, found by the identity check:** the published D84 numbers predate the 2026-11-03 ARACNe v2 TF rebuild
  that added the 12 subtype networks, so ARACNe's benchmark numbers moved (kd 0.597 -> 0.594, oe 0.690 -> 0.694)
  while CollecTRI and DoRothEA reproduce to 1e-15. The sweep reproduces 79n on current inputs at 3.9e-07. D84's own
  `keep_consensus` decision was taken on the old build and is left standing, recorded as open.
- Part B (D126, the mechanism-consistency filter) is untouched and still specified, not run.

## 2026-11-09: validation changed (D127, D125.1) — approval is no longer the primary yardstick
- User: approval status confounds efficacy with patents, market size, rare-cancer economics and sponsor strategy, so
  an unapproved drug-cancer pair is not a negative. The class baseline (anticancer 0/1) already beats the lead score
  on Hub (0.724 vs 0.673), i.e. the metric largely measures drug class.
- **D127 (specified, not run):** biomarker-expectation test. The same drug is compared against itself in two cohorts of
  a biomarker pair, so drug-level confounders cancel. 7 pairs, 4 classes, 49 instances (BRAF/MEK thyroid and melanoma,
  EGFR lung and head/neck, HER2 breast, PARP ovarian). Pass = >= 70% in the expected direction, sign test p < 0.05,
  median gap >= 10 percentile points.
  - Scored on `layered_pct`, which ranks over gate-passing scored FDA compounds and is unaffected by target_absent /
    local_route / low_exposure: the gates otherwise delete the negatives (only 1 of 7 anti-estrogens is eligible in
    Basal), the same trap as the first D124 breast check.
  - Breast anti-estrogens, ALK lung and CDK4/6 are excluded from the primary set and reported separately: their values
    are already known or prior decisions were taken on them. AR PRAD vs BLCA is not computable.
  - Stated limit: it checks that known biology reaches the ranking, not that the objective finds unknown biology.
- **D125.1:** D125's known-use criterion becomes reported only; the controls criterion widens to D127 and leads.
- **D127 run on production (11-09): PASS.** 36 of 49 instances correct (0.735), sign p 0.0007, median gap 15.1.
  Reversal on the same instances 0.612 (p 0.076) and 0.528 on the layer-1 subset: at chance, as everywhere else.
  Layer 1 alone matches the full score's direction (30/36 both) with a much wider gap (39.4 vs 22.8), so layer 2
  compresses separation without reversing it.
- **First finding the new test produced that known use could not:** the thyroid pair fails 3/9 with a median gap of
  exactly 0. BRAF alteration is 86.8% in BRAF-like vs 15.9% in RAS-like thyroid, but layered_A is 100.000 vs 99.996
  and both BRAF inhibitors sit at percentile 100 in both subgroups. The A component is a within-cohort percentile over
  genes, so a top-ranked target cannot discriminate two cohorts however far apart the raw shares are. Open: whether A
  should carry magnitude (raw share or cross-cohort scaling). Needs its own decision; nothing changed.

## 2026-11-09: noise ceiling for the reversal score (D125 section A0)
- Spec `docs/SPEC_graph_metric_scoring.md` (D125 graph-metric reversal, D126 mechanism flag); A0 implemented
  (`scripts/d125/`, outputs `results/lincs/d125/a0/`), Parts A and B not started. Choices before running (D125-A0):
  dose x time cells as A0a strata, active level-5 signatures as replicates (no level-4 profiles on disk), r_c over drugs
  that use a changed pair, production line vectors kept and only changed VIPER columns rerun (spliced).
- The spec's sizing was off: 9,678 drug x line profiles and 40 networks, not 48,274 and 38 (spec annotated).
- Sandbox reproduces production exactly (profiles, VIPER, all 67,174 spec_mech_pct_pregate rows).
- A0.0: 1,019 deep pairs in 32 lines, A0b cohort-level. A0a: r_c median 0.90, but direct (matched-line) drugs 0.61 and
  borrowed drugs 0.905 (inflated: half their lines are unchanged between passes).
- Ceiling: AUC_max 0.97 at r 0.90 (0.89 at r 0.61); reversal known-use AUCs 0.43-0.51 sit near zero attainment
  (-0.15 to +0.03). The reversal score's failure is not the assay noise. A0b dry run: ~3-15 h; waits for a numeric
  stopping rule.
- D125-A0.1 (user): r_c only over drugs whose every scoring line changes between passes (the shared unchanged lines
  inflated borrowed drugs); A0b stops at SE <= 0.01. A0a r_c now 0.81 (old version 0.90 kept as r_c_any); AUC_max 0.95;
  reversal attainment -0.16 to +0.03.
- A0b (complete, 11-09 21:26): stopped at the 20-split floor, SE 0.0045. r_c median 0.881; A0b - A0a = +0.234 in every
  cohort, so dose and time variation is most of what A0a called noise. Only 14 of 55 cohorts have >= 50 drugs under the
  all-lines-changed rule; 41 are r_pooled, 34 keep A0a's r. Final ceiling r 0.873, AUC_max 0.967; reversal attainment
  -0.15 to +0.03. The conclusion holds on the stricter estimator.

## 2026-11-06: trial failure reasons (D108) and evidence flags (D109, D109a, D109b)
- D108: for each cohort's top 10 FDA drugs, trials that list the drug; stop reasons by keyword rule, completed trials
  read from linked abstracts by subagents (631 trials, quotes verified), drug role (single agent / combination /
  comparator / background) and stop attribution. Disease top 10: 147 pairs with a trial, 41 with a failed trial,
  44 with a positive one. Pages show 'Why trials failed'.
- D109: dose-consistency rule not adopted (flagged 45% of judged profiles, incl. doxorubicin); replaced by
  plate_confined (D109a, 781 drugs, incl. metronidazole). inactive_matched_line 2,271 rows; generic_mechanism
  omacetaxine / digoxin / dactinomycin (shown only); low_exposure from FDA label / archived label / review Cmax
  (123 of 384 FDA drugs, D109b), 51 drugs flagged.
- Rebuild: FDA disease rows 2,722 -> 2,426; metronidazole, pentobarbital (PAAD), azacitidine (PAAD), defactinib out
  of the PAAD top 10; fedratinib still #1 (tested 8 uM vs label Cmax 3.4 uM, not flagged).

## 2026-11-06: route flag (D107) and literature rerun
- New flag `local_route` (`scripts/route_flag.py`): no systemic formulation in Drugs@FDA (ChEMBL topical-only for
  non-FDA drugs); bladder exception for intravesical. 49 of 527 FDA ingredients flagged; FDA rankings 2,913 -> 2,722
  disease rows; 8 top-3 placements removed (ciclesonide, halcinonide, ciclopirox, ingenol mebutate).
- Literature (79o pairs, MODE=fda, MODE=fda_all) rerun on the D105b-c rankings: the 11-02 tables were stale (17 of 123
  disease top-3 pairs had changed). Every FDA pair now searched: PubMed >= 1 64.9%, trial >= 1 44.9%, either 69.9%;
  top 3 either 76.4%. A hit means the pair was studied, not that it worked.
- Not covered: systemic drugs with implausible exposure (metronidazole, pentobarbital, cardiac glycosides) and
  class-level trial failures (JAK inhibitors in PDAC: ruxolitinib JANUS 1/2); both still open.

## 2026-11-05: patient-level VIPER on the cell-type networks (D106)
- `61` now reads the network map; rerun for the 12 switched groups (no errors). Patient NES vs old r 0.84-0.97;
  concordance with group MRs up in 6 groups, down in 5 (ESCA-SCC 0.61 -> 0.31); confounding still negligible.
  Subgroups stay fixed; 84 / 86 and the Step 8 inputs are stale on these scores.

## 2026-11-04: disease side left as is (D103)
- User: keep the D79 subgroups fixed and leave the disease side for now. Patient-level VIPER (`61`, corrected patient
  NES) has not been run on the updated pipeline (D97 cell-type networks); patient NES, Step 8 inputs and Step 7-8
  summaries for the 12 switched groups still use the pooled networks. Group-level signatures (drug ranking) are current.

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

## 2026-11-15: D128 implementation plan reviewed against the repo, before A0 (nothing run)
- Notes: `docs/notes/d128_plan_review.md`. Disease effect settled: production limma + voomByGroup `logFC` / `z`
  (the `*_voom` columns are check C12 only), with 79t's contrast rule (GTEx, else adjacent).
- Open, to settle before A0 freezes anything:
  - (1) R as specified repeats the 11-12 precursor calculation, which inverted the anti-estrogens against GTEx
    normal breast; the reference for breast R needs a decision. M is unaffected.
  - (2) No EGFR-sensitizing-mutant LINCS line exists (H1975 only: T790M, 3 profiles), so A4 is expected to be
    not computable in context.
  - (3) Fulvestrant is the only pure ER antagonist (3 pert_ids, 18 MCF7 signatures; corrected from "1 profile"); whether SERMs count as antagonists needs a
    decision; no estradiol profile in MCF7.
  - (4) No primary arm is named; network, gene set, profile arm, hop limit and weights should be frozen together.
  - (5) Generic-response removal may remove ER/proliferation signal in MCF7.
  - (6) `wD` is undefined.
  - (7) The holdouts are new patient groups but reuse D127's drug profiles; PRISM was already read in 79j/79k/79m.

## 2026-11-15: D128-A0 label-blind feasibility audit — run, with one failed run recorded
- `scripts/d128/{config,d128io,a0_audit}.py`, `tests/d128/` (11 pass, pytest added to the env), outputs
  `results/lincs/d128/a0/`, log `logs/d128_a0.log`. Choices fixed first in `DECISIONS.md` (D128-A0). **No M or R value
  computed, no outcome label read, nothing frozen, no `d128-a0-freeze` tag.**
- **The first run failed on performance and is part of the record.** Killed after 43 min without finishing. Diagnostics
  taken while it ran: one core at 100%, CPU time equal to wall time, RSS 1.74 GB, state R — compute-bound and
  single-threaded, not deadlocked or I/O bound. Only `profile_value_audit.tsv` had been written (~3 min in), so
  execution sat in the footprint loop for ~40 min. Cause: `reach()` re-walked the same breadth-first search for target
  genes shared across many drugs and units, and a pandas `.isin()` filter ran once per reached transcription factor
  inside that loop. Progress was invisible because stdout was piped through `tail`, which buffers.
- **Fixed by:** sign bitmasks instead of Python sets; one cached BFS per gene at the largest hop arm, with smaller arms
  filtered from the same result (shortest-path signs do not depend on the cap); regulon gene lists pre-restricted to
  the LINCS gene space so no pandas runs in the inner loop; a reachability cache persisted under
  `results/lincs/d128/a0/cache/` keyed by the OmniPath/SIGNOR checksum, as the plan's A1 section asks; a progress line
  every 25 units with an ETA; and the job run under `nohup` to a real log file. **43+ min unfinished → 141 s**, 202 s
  for the full run including the DepMap and TCGA sections.
- **Feasibility.** 12,283 cohort × compound pairs over 55 cohorts, 113,688 gate-passing (cohort, compound, line, unit)
  rows; 3,931 of 9,288 group pairs and 1,270 of 2,995 subgroup pairs are eligible for at least one of M or R, and every
  other row carries a reason. Of 497 target units, 225 are M-eligible, 226 have no anchored regulator within two hops.
- **The ESR1 control is feasible, which is the point D126 failed on.** Fulvestrant's ESR1 unit resolves 294 genes
  (40 landmark) at zero hops, against a floor of 10. Tamoxifen's ESR1 unit is annotated MODULATOR and is therefore
  ambiguous for M, so only fulvestrant is scorable until D0.2 is decided.
- **A4 (EGFR) confirmed not computable in a matched line** (0 signatures); A2 (ERBB2) is executable but thin, 12
  signatures in SKBR3; A3 (BRAF) has 104 in A375 and HT29; B3 (PARP) has no BRCA1/2-damaging LINCS line.
- **PRISM is not a clean independent label:** it is already read by `79j`, `79k`, `79m` and `d94_decide.py` and cited in
  D94, D97, D102c, D119, D120, D122 and A43. GDSC appears only alongside it.
- **Two rule problems, both found before any score existed.** The ERBB2-amplification rule used a 99th-percentile cut
  and missed SKBR3 at copy number 11.0 (98.7th percentile); it is now an absolute relative-copy-number cut at 2.0,
  where the distribution is bimodal, and the re-run is what made A2 executable. Unsettled: a unit counts as
  TF-anchored when any of its genes is a CollecTRI regulator and is then floored on zero-hop genes alone, so
  palbociclib is below the floor with 0 zero-hop genes despite 1,758 resolved through signaling. That needs a decision
  at the freeze.
- Open and blocking the freeze: D0.1 (reference for R in breast), D0.2 (SERMs), the unit-kind/floor question, and what
  A5 can use for a functional label if not PRISM.

## 2026-11-15: D128-A0.1 — the unit-kind/floor rule corrected (per-gene typing), A0 re-run
- The open item flagged in D128-A0 is fixed and the audit re-run end to end (144 s). Still label-blind: no M or R
  value, nothing frozen, no `d128-a0-freeze` tag, D0.1 and D0.2 still the user's decisions.
- **Was:** a target unit counted as TF-anchored as soon as *any* of its genes was a CollecTRI regulator, and the whole
  unit was then floored on zero-hop genes alone. **Now:** each gene is typed by how it actually anchors — a regulator
  contributes its own regulon at zero hops, a non-regulator contributes a 0-2 hop signed-path footprint — each subset
  is floored on its own genes, and the unit clears if either subset does. Genes are deduplicated across subsets, and a
  gene the two subsets disagree about counts as conflicting. `unit_kind` is now `tf`, `signaling` or `mixed`, with
  per-subset counts and `cleared_by` reported.
- **Why it mattered:** palbociclib (CCND1|CDK4) was below the floor with 0 zero-hop genes although CDK4 resolves 1,781
  genes through signaling, which blocked wave B2 on a typing artefact while 29 in-context signatures exist. B2's unit
  status goes from `below_floor|below_floor` to `eligible|eligible`.
- **Deltas:** units M-eligible 225 → 238; below floor 15 → 2; unit kinds 98 TF → 92 TF + 6 mixed; M-eligible rows
  73,336 → 80,666; rows blocked by the floor 7,849 → 519; M-admitted pairs 4,775 → 4,940. The landmark-only arm gains
  most (185 → 226 units), because per-subset floors no longer import conflicts from TFs a unit reaches only at 1-2
  hops. **No unit lost eligibility**, which is the check the union rule has to pass.
- Of the 13 units admitted, 2 are palbociclib's (cleared by signaling) and 11 are TF-anchored units — imatinib,
  dasatinib, nilotinib, ponatinib, bosutinib, regorafenib, vorinostat — that clear on their own regulon of 10-18 genes
  now that other TFs' conflicting expectations no longer suppress them. The 2 still below the floor fail on data:
  one has a single anchored TF, the other no resolved genes.
- **Correction to a number reported in D128-A0:** fulvestrant's ESR1 unit resolves **478 genes (89 landmark)** on its
  own regulon, not the 294 reported then, which had been intersected with the globally resolved set and lost genes to
  conflicts from TFs reached at 1-2 hops. The conclusion is unchanged and stronger. Cohort × compound pair totals are
  unchanged (3,931 group, 1,270 subgroup).
- `tests/d128/` now 13 tests (2 added for the union rule).

## 2026-11-15: D128 Section 0 frozen, and A1 (signed route index) built
- **Freeze.** All Section 0 rows P1-P14 settled in the `D128-A0 freeze` entry, `scripts/d128/config.py` carries them
  with `FROZEN = True`, tag `d128-a0-freeze` pushed. Two rows changed in application rather than value: P4 (landmark
  genes are a median 15.0% of a resolved footprint after the A0.1 correction; a landmark-only floor admits 226 of 497
  units) and P10 (floors applied per subset, 238 units admitted, 2 below floor on data). P9 was explicitly conditional
  on the audit and is kept: `tanh(z/2)` saturates above 0.9 on only 1.3% of specific-arm values. P14 verified exactly:
  the sandbox reproduces production's direct/borrowed split on 12,283 of 12,283 pairs.
- **A1 built and gated.** `graph.py` (factored from the audited A0 code) and `a1_build_paths.py`; 1.6 s. 14
  hand-calculated path/sign fixtures pass, 27 tests in all. Outputs: 42,432 regulon edges, 31,769 target-to-TF path
  rows over 387 target genes, 367,942 (unit, gene) footprint rows over 257 units.
- Ambiguity is kept rather than resolved: median per-unit TF ambiguity 0.091, gene ambiguity 0.247.
- **The ESR1 footprint is clean**: 478 resolved genes (89 landmark), none conflicted or ambiguous, 405 expected up and
  73 down before the drug's action sign.
- A latent bug surfaced while writing the fixtures and is fixed: the footprint memo was keyed on the network name, so
  two graphs sharing a name could collide. The key now carries the graph checksum.

## 2026-11-15: D128 A2 (M/R engine) and the mandatory ESR1 control
- **A2 built and gated.** `score.py` (pure `score_mechanism` / `score_rescue`), `nulls.py`, `a2_score.py`. 51 tests
  pass, 24 of them A2 fixtures with hand-calculated values. Gate items all met: every defined score in [-1, 1]; a
  score of 0 is distinguishable from NA (0.0 + `scored` vs NaN + `abstained` with a reason); no confident score for a
  unit with no genes; a duplicated gene raises rather than double-counting; nulls seeded, reproducible and
  non-mutating.
- **A null was contaminated and is fixed.** The degree-matched-target null returned a median exactly equal to the
  observation, because the candidate pool held other drugs' units against the *same* target — every ESR1 drug has an
  identical footprint and therefore an identical score. Candidates sharing a target gene are now excluded. Found and
  fixed before the control ran.
- **The ESR1 control is met on the primary arm.** Fulvestrant in MCF7, all three `pert_id`s positive: M = 0.048,
  0.040, 0.025 on the full 478-gene footprint (89 landmark, 0 missing), at the 93rd–100th percentile of both the
  sign-shuffle and gene-label nulls. D126 flagged 0 of 22 anti-estrogen rows; anchoring at the annotated target
  recovers the footprint.
- **The effect is small and we say so:** M ≈ 0.04 on a [-1, 1] scale. It is consistent and reproducible across three
  independently profiled `pert_id`s and clears its nulls, but whether it adds anything over the production layered
  score is an A6 question this control does not answer.
- **The raw-profile arm does not rescue M** (mean 0.0035, one negative), so the plan's "passes only on raw" failure
  mode does not occur: removing each line's generic response helps rather than hides the target signal.
- **SERMs, reported separately (D0.2):** 11 of 12 rows positive on the specific arm; clomifene's single-signature
  profile is the exception at −0.059. No agonist sign-flip control exists (estradiol has 0 MCF7 signatures).
- **Breast R behaved as predicted in advance (D0.1), and the known reference problem reproduced.** R is
  uninformative in the ER-positive cohorts (LumA 0.010, LumB 0.004, Her2 0.001) exactly as recorded, but highest in
  ER-negative Basal (0.049) — the same ordering the 11-12 precursor found. That is the GTEx-normal-breast reference
  problem reappearing, not evidence for fulvestrant in Basal. Nothing was adjusted in response.

## 2026-11-16: D128 A3 (full scoring, aggregation, baselines) and A4 actions 2-5 (controls)
- **The first full A2 run happened.** 113,688 (cohort, drug, source line, target unit) rows, scored twice (two gene
  arms), 227,376 rows in 282 s on 8 threads. The engine was re-indexed first - footprints, units and profile rows
  are looked up instead of re-scanned - and the refactor reproduces the committed A4 ESR1 control **byte for byte**,
  so nothing moved. The two permutation nulls gained a vectorised path used only where it is provably exact, with a
  test asserting element-wise equality against `score_mechanism`.
- **A3's join gate is met.** 10,131 drug-level rows in and out of four production joins with no loss or duplication;
  P14 re-asserted at 12,283 / 12,283 pairs; `layered_pct`'s two independent copies agree to 5e-4 (a 3-dp storage
  difference); the 16 target-role `T` disagreements are the two FDA ingredients that carry two `cmap_name`s each, so
  `target_role_pct` is a drug-level number shared by both `pert_id`s and is now labelled that way. Abstentions
  account for 100% of the 124,168 A0 candidate rows.
- **The scores are small and noisier than they are large.** M median 0.0049 (IQR -0.0089 to 0.0194). Across source
  lines, the weighted sd of M has median 0.0170 against a median |M| of 0.0147: **the disagreement between the lines
  a cohort score is averaged over is larger than the score itself.** 76.7% of unit rows are borrowed.
- **M is nearly independent of the production lead score** (Spearman 0.17 against `layered_pct`); R overlaps VIPER
  reversal (0.43). That means M measures something else, not that the something else is real. **A3 has no outcome
  label and cannot settle whether M/R add value.** That is A5/A6's question and it is open.
- **The landmark sensitivity amendment failed its pre-stated criterion on all three parts, and not marginally.**
  rho(M_all, M_landmark) = 0.677 (needed 0.90), sign flips 25.0% (needed <= 5%), median |delta percentile| 9.87
  (needed <= 5). On the quarter of units with the lowest landmark share, rho 0.423 and flips 37.3%. 567 units clear
  the floor on all BING genes and fall below it on measured genes alone. Fulvestrant's ESR1 M is 2.5x larger on the
  89 measured genes than on all 478, with weaker nulls. **Under the rule written before the result: one unit in
  four changes the SIGN of its mechanism score depending on whether L1000-imputed genes are included, so signal
  cannot be separated from imputation.** Nothing was changed in response - P4 stays frozen, both arms are reported -
  and the conditional population-wide null spend was not made, as the amendment required.
- **The degree-matched-target null is thinly powered nearly everywhere.** Measured structurally, without running a
  null: 95.0% of M-eligible units have fewer than 10 distinct candidate footprints, 27.7% fewer than 3, 9.7% none.
  ESR1's pool collapses to a single distinct footprint (AR), which is why its z came back blank. Only the two
  permutation nulls carry weight.
- **BRAF is the one informative control wave, and its M is MYC-hub-dominated.** Positives (vemurafenib, dabrafenib
  in A375 and HT29) M 0.025-0.041, 3 of 4 confident; the RAF-paradox negatives (BRAF-WT, RAS-active) produce **0**
  confident positives. The negative control passes. But **692 of the BRAF unit's 1,063 footprint genes (65.1%) are
  MYC-regulated and carry 66-69% of the absolute numerator and 69-95% of the signed numerator in the four positive
  rows**; 610 of them are reached by MYC and by nothing else. By the threshold fixed in advance this is hub
  dominance: BRAF's M is largely a MYC-driven proliferation read-out routed through the BRAF anchor.
- **HER2 fails the negative rule** (14.9% confident positives among ERBB2-non-amplified rows against 0% among the
  four SKBR3 positives) on a positive set far too small to carry a verdict - exactly the thinness A0 flagged.
  **EGFR stayed a documented gap:** the only sensitising-mutation line also carries T790M, and no substitute was
  used. Five control compounds have no D121 target unit at all.
- **The D127 replay (descriptive only) does not reproduce the production discrimination.** M points the right way in
  16 of 49 primary instances against `layered_pct`'s 36 of 49 - and the reason is structural: M has no disease term,
  so **15 of 49 instances are exact ties**, including all nine THCA BRAF-like vs RAS-like ones. Among the 34 where M
  can differ at all it is right 47.1% of the time. M as built cannot separate a biomarker pair for a borrowed drug.
- **Failures were analysed by mechanism, not rank.** Missing network path 49.8% of drug rows, line mismatch 42.2%,
  absent target expression 17.3%, sign contradiction 17.0%, disagreeing target units 12.2%. The cleanest internal
  negative control is the one needing no label: where the annotated target is **not expressed in any source line**,
  M is still positive 56.9% of the time (median 0.0006) against 59.8% (median 0.0046) where it is expressed.
  Eight-fold apart in the medians, almost completely overlapping in distribution.
- **Not run, deliberately:** A5, and no CTRP / PRISM / GDSC value was read. **Open for the user:** whether
  `landmark_only` becomes the primary gene space pipeline-wide, and a 55-row regorafenib mismatch where the engine's
  floor is looser than the frozen P10 per-subset floor.

### 2026-11-16 (later the same day) - D128-A2.1: the engine is made to obey the frozen floor

- **The defect A3 flagged is fixed, and it was the engine, not the rule.** P10 as frozen by D128-A0.1 floors each
  subset of a target unit on its own genes (a direct-TF subset on >= 10 zero-hop regulon genes; a signalling subset
  on >= 10 resolved genes **and** >= 2 anchored TFs) and admits the unit if either clears. `score.py` applied only
  the >= 10-gene half, to the merged footprint, and so scored `regorafenib::CHEMBL5122` (DDR2: 173 resolved genes
  behind a single anchored TF). `score.py` now calls the same `graph.floors()` A0 used, `a1_build_paths.py` labels a
  gene reached by both subsets `both` so the subset label is lossless, and `path_ambiguity.tsv` carries the
  signalling subset's TF count. **No frozen value changed and `config.py` is untouched.**
- **The stop rule written before the re-run held.** A2: exactly 754 of 227,376 rows (377 per gene arm, all that one
  unit) moved from scored to abstained with reason `below_floor_M_unit:tf=0,signaling=173/1tf<10/2`; every other
  row's M, R, statuses and coverage counts are bit-identical. A3: the mismatch counter is 0 (was 377 unit rows /
  55 drug rows); 55 drug rows change M, all regorafenib, **none becomes NA** because P13's unit mean simply drops
  DDR2; **R is unchanged on all 10,131 drug rows.**
- **Everything else that moved is mechanical re-ranking:** 138 of 4,940 M percentiles shift by at most 1.18 points
  because those rows left their cohort's ranking. In the D127 replay 4 of 67 rows change - 2 regorafenib (own M),
  lapatinib and olaparib in percentile columns only - and **no instance verdict changes**; the 16/49 and 25/49
  headlines stand. The `absent_target_expression` case report drops 1,748 -> 1,734 rows, which are **exactly** the
  14 regorafenib DDR2 rows.
- **The controls reproduce exactly.** `esr1_control*.tsv`, `esr1_breast_R.tsv`, `control_results.tsv`,
  `control_summary.tsv`, `braf_hub_dominance.tsv` and `braf_hub_genes_all_bing.tsv` are byte-identical, so **the
  BRAF hub numbers do not change and the vemurafenib/dabrafenib positives are untouched** (regorafenib appears in
  that table on its BRAF unit, a different unit that was never below the floor).
- **Landmark sensitivity after the fix:** 10,086 units in both arms (was 10,141), rho 0.674, sign-flip 25.13%,
  median |delta percentile| 9.80 - still **NOT MET** on all three parts. 57 tests pass in `tests/d128/`.
- **Recorded for review:** the engine floors the genes that actually enter the numerator rather than A0's footprint
  count. They coincide today (0 missing responses) but the engine is the stricter of the two if a profile matrix
  ever has gaps.

### 2026-11-16 - D128-A3.1: why M flips sign between the gene arms (diagnostic; nothing changed)

- **The question A3 left open had four candidate answers, written down before the run. Two are supported.**
  Population: 1,852 (unit, `pert_id`, line) rows over 226 units scored in both arms (the cohort does not enter M);
  overall flip rate 27.9%.
- **Control A: flipping is confined to small |M| - but small |M| is where nearly everything lives.** 46.9% flips at
  |M| <= 0.005, 24.2% at 0.02-0.03, 3.3% at 0.05-0.08, 0% above 0.12, so by the pre-stated rule flipping stops at
  |M| >= 0.05. **Only 11.9% of rows, and 12.0% of the 4,940 scored drug rows, get there.** On the conservative
  max(|M_all|, |M_landmark|) measure flipping persists to 15.6% even at 0.08-0.12.
- **Control B: it is NOT imputation.** Random same-size subsets of the all-BING footprint flip 26.0% of the time
  against the landmark split's 28.8%; per unit the gap has median +0.008 (IQR -0.113 to +0.142), with 92 units worse
  and **85 better** on the landmark split. Switching arms would not fix the flipping.
- **H2 is the mechanism, and it is bigger than D128.** Effective independent genes (participation ratio of the
  gene-gene correlation matrix over 2,000 profiles): a median 905-gene footprint carries **n_eff 44**; its inferred
  part 42; its 140 landmark genes 41. The **whole 10,112-gene BING space carries n_eff 41**; all 978 landmarks, 57.
  **A footprint cannot carry more than ~40-55 independent measurements however many genes it lists.** This explains
  the dilution of M, the 2.5x larger landmark-only M on fulvestrant, and - the part that matters beyond D128 - that
  **the gene-label null, which permutes genes as if independent, is too narrow, so its percentiles and z values are
  over-confident wherever LINCS gene counts are treated as sample sizes.**
- **H1 (hub dominance) is not supported once |M| is controlled:** within |M| > 0.03 the flip rate is 9.0% / 9.2% /
  9.6% across hub shares < 0.3 / 0.3-0.5 / > 0.5. The marginal association was a |M| confound.
- **H3:** the median row has 78% of its footprint genes at |z| < 1 and they carry 54% of the absolute numerator, so
  M's sign is decided mostly by near-boundary responses.
- **ESR1 / fulvestrant:** landmark flip rate **0.000** against 14.8% for random 89-gene subsets; n_eff 46 (all) /
  38 (landmark), i.e. the 389 inferred genes add about 8 effective measurements.
- **Nothing was changed:** no arm chosen, no parameter moved, nothing re-ranked, P4 still frozen. Open for the user:
  a minimum-|M| rule before a sign is read, and whether the gene-label null must be replaced project-wide.


## 2026-11-16: D129 — per-gene quality table for the BING gene set (diagnostic only)

`scripts/96_bing_gene_quality.py`, 101 s on 24 threads, log `logs/96_bing_gene_quality.log`, outputs
`results/lincs/gene_quality/bing_gene_quality.tsv` (10,174 rows), `bing_gene_quality_summary.tsv`,
`bing_gene_quality_run.json`. **Nothing was filtered, dropped, re-weighted or curated; no scoring code was touched
and no gene set was chosen.** Every definition was fixed in D129 before anything ran, and every statistic is
label-free: only `level5_bing_hiq.npy` (90,988 signatures x 10,174 genes) and the experiment metadata
(`pert_id`, `cell_iname`, `pert_dose`, `pert_time`) were read. All four statistics came out for all 10,174 genes;
`missing_reason` is empty everywhere.

- **The question D128 raised is answered the same way by every statistic: the inferred genes are not mainly noisy,
  they are mainly redundant and compressed.** Held-out R^2 on the 978 landmarks of the same signature, with the
  split taken by perturbation so replicates cannot leak: **best inferred median 0.875** (p10 0.834, p90 0.917,
  range 0.437-0.964), against the control of each landmark regressed on the other 977, **median 0.630** (p10 0.418,
  p90 0.796). 98.5% of inferred genes are above R^2 0.8 against 9.4% of landmarks; 22.9% above 0.9 against 0.4%.
  No gene reaches 0.95 in the control and only 28 inferred genes do, so the inferred values are not *exactly*
  deterministic functions of the landmarks - they are about seven-eighths of one.
- **The decisive comparison is R^2 against ICC.** For **100.0% of inferred genes the landmarks predict the gene
  better than the gene's own replicate does** (R^2 0.875 vs ICC 0.544 at the medians); for landmarks the same
  comparison is a coin flip (50.9%, 0.630 vs 0.617). Arithmetic on those two pre-registered columns (derived after
  the fact, not a new measurement) says that **at least 73% of an inferred gene's replicate-to-replicate noise is
  itself a function of the same signature's landmark values** (median of the per-gene bound 0.725), against 0.8%
  for landmarks. Inferred genes therefore do not average their noise away against each other the way independent
  measurements would - which is the mechanism that would make permutation nulls too narrow.
- **Reproducibility (ICC(1) over 5,345 replicate groups, 14,906 signatures, k0 = 2.789).** Landmark median 0.617
  (IQR 0.523-0.710), best inferred 0.544 (0.487-0.600); split-half r 0.624 vs 0.551. The inferred distribution is
  shifted down but **narrower**, and it overlaps the landmarks heavily. By the strict reading, **989 of 9,196
  inferred genes (10.8%) fall below the landmarks' 10th percentile of ICC (0.443)**, 273 (3.0%) below the 5th and
  19 (0.2%) below the 1st; by split-half r, 721 (7.8%), 137 (1.5%) and 8 (0.1%). Against the landmark *median*,
  81% of inferred genes are below - but so are 489 landmarks by construction. **There is no large population of
  inferred genes that is simply noise.**
- **Variance compression is present and large.** Median sd over all signatures 0.918 (best inferred) vs 1.569
  (landmark), a factor of **0.585**; p01-p99 range 4.80 vs 9.20 (0.52); mean |z| 0.681 vs 1.039 (0.66); IQR 1.048
  vs 1.331 (0.79, the smallest gap, because level 5 values are z-like near the centre and the compression is in the
  tails). Both parts of the variance shrink by about the same factor - reproducible signal `sd_between` 0.758 vs
  1.383 (0.55) and replicate noise `sd_within` 0.693 vs 1.091 (0.64) - so the *ratio* `sd_between / sd_all` barely
  moves (0.821 vs 0.889) and ICC barely moves with it. **The inference model shrinks the whole response toward the
  mean, it does not add noise.** On signal amplitude the two sets barely overlap: 8,153 of 9,196 inferred genes
  (88.7%) have less reproducible signal than the landmarks' 10th percentile. Level 5 clips at +/-10, which 962 of
  978 landmarks reach but only 818 of 9,196 inferred genes.
- **Effective independent information: adding 9,196 inferred genes does not add effective dimensions, it
  concentrates the ones already there.** Participation ratio of the gene x gene correlation matrix over all 90,988
  signatures: **BING (10,174 genes) 25.1; the 978 landmarks alone 35.8**; best inferred alone 23.6. Entropy
  effective rank 177 vs 215 vs 159; components to 80% of the variance 203 (of 10,174) vs 328 (of 978); to 95%,
  602 vs 699. The sampling-noise reference for the same number of *independent* genes at this sample size
  (Marchenko-Pastur, p / (1 + p/n)) is 9,151 for BING and 968 for the landmarks, so the concentration is real and
  not a small-sample artefact. The covariance version (86's center-only convention) says the same: 23.2 vs 28.5.
  **Read correctly, this does not mean BING holds less information than the landmarks - BING contains them, so its
  span cannot shrink. It means the variance of the 10,174-gene set lies along fewer effective directions than the
  variance of the 978-gene set, because the added genes are copies of directions that already existed.**
- **Per gene, the independent-variance numbers agree.** `uniqueness_all_bing` (1 - R^2 on all 10,173 other genes,
  correlation matrix with a 1e-6 ridge, so an over-estimate of independence): landmark median 0.048, best inferred
  0.0054, nine-fold apart. Summed, the whole BING set is **113.9 independent gene-equivalents, of which the 978
  landmarks hold 50.4 and the 9,196 inferred genes add 63.5.** The marginal number - summing 1 - R^2 on the
  landmarks alone, which ignores that the inferred residuals are correlated with each other - is 1,157, so the
  honest statement is that the 9,196 inferred genes add somewhere between ~60 and ~1,160 gene-equivalents, and the
  conditional measure says the low end.
- **What that means for a footprint, arithmetically.** On the conditional measure a 478-gene inferred footprint is
  worth on the order of 478 x 0.0054 ~ 2.6 gene-equivalents of independent variance, against 478 x 0.048 ~ 23 if
  those genes were landmarks. Both hypotheses D128 raised are therefore supported, but not equally: redundancy is
  the dominant effect by a wide margin, noise is a real but small secondary one confined to ~10% of inferred genes.
- **The surprise, and the thing a filter design has to face: redundancy and reliability are positively correlated
  among inferred genes, not negatively.** Spearman R^2 vs ICC = +0.770; R^2 vs sd 0.852; ICC vs `sd_between` 0.980.
  There is a single amplitude axis - genes with large responses are both well predicted by the landmarks and
  reproducible; genes with small responses are neither. **So the inferred genes that carry the most independent
  information are the least reliable ones, and there is almost no gene that is both.** Of the 9,196 inferred genes
  only 36 have R^2 <= 0.7, and **not one of those has an ICC at or above the landmark median**; 24 of them clear
  even the landmark 10th percentile. The median ICC of the 500 least redundant inferred genes is 0.423, of the 500
  most redundant 0.668. In the other direction, 0 genes have R^2 >= 0.9 together with an ICC below the landmark
  10th percentile. The least redundant inferred genes (RPS4Y1 0.437, HBB 0.511, MAL 0.528, GDF15
  0.543, VIM 0.548) are cell-identity and stress genes with ICC 0.26-0.50, at or below the bottom of the landmark
  distribution. The most redundant are the proliferation and immediate-early block (PPP1R15A 0.964, HSPA6 0.959,
  ASPM 0.959, ATF3 0.957, BUB1/MAD2L1/DLGAP5/CDCA3 ~0.957), all with ICC 0.75-0.81.
- **Assumptions and limits, stated in D129 before the run.** Only linear/second-order dependence is visible, so a
  gene that is a non-linear function of the landmarks still counts here as independent (CMap's inference is a
  linear model, so this is a mild assumption for this particular use). Signatures are treated as exchangeable
  draws, which they are not exactly - compounds repeat across plates - so the effective-dimension numbers describe
  this signature set, not L1000 in general. The ridge on the correlation matrix biases `uniqueness_all_bing`
  upward, so every redundancy statement here is conservative. ICC and split-half r deliberately conflate signal
  with noise: a gene that never responds reproducibly scores low whether the cause is measurement noise or
  genuinely flat biology, and the landmark set is itself enriched for high-variance genes by CMap's own selection,
  which is why the landmarks are used only as an internal benchmark and not as an absolute standard.
- **Not done, deliberately:** no filter proposed, no threshold recommended, no gene set chosen, nothing scored
  twice, no D128 output read. The table is the evidence base; the filter design is the user's separate decision.

### 2026-11-16 - D128-A3.2: the reference-activity flag (computed for every target unit x cohort; one of four checks FAILED)

- **What it is.** `M_disease(unit, cohort, contrast)` is M's own formula fed the disease contrast instead of the
  drug profile - same footprint, same expected signs, same weights, same `tanh(z/2)`, same frozen P10 floor, with
  the P1 columns (`logFC` sign, `z` magnitude). Positive means the tumor ALREADY looks target-inhibited against the
  reference, so R penalises the mechanistically correct drug; negative means the drug reverses the tumor toward
  normal and R is interpretable. Significance from the two well-powered nulls (200 seeded draws each, both must
  clear 97.5 / 2.5); the degenerate degree-matched null was not used. **20,046 (unit, cohort, contrast) rows** over
  257 units, 55 cohorts and 78 disease tables, 296 s. Reporting only: nothing gated, filtered or re-ranked, and M,
  R and every A3/A4 value are unchanged.
- **V2 discrimination: PASSES.** gtex: neutral 70.0%, aligned 14.3%, inverted 4.4%, unavailable 11.3%; eligible-only
  78.9 / 16.2 / 5.0. adjacent: 71.0 / 12.2 / 5.5 / 11.3, eligible-only 80.0 / 13.7 / 6.2. No level over 80%.
- **V4 (expected, not required): PASSES.** ERBB2 x BRCA-Her2 -> **neutral** (M_disease -0.0051, percentiles
  35.5 / 31.0). BRAF x TCGA-SKCM -> **aligned, strongly** (M_disease -0.1424, both percentiles 0.0): the melanoma
  looks BRAF-active against normal skin, so inhibiting BRAF moves toward normal. Both are tumour-acquired events, as
  the pre-statement expected.
- **V3 does the flag predict R? PARTIAL.** On R's level the effect is **weak**: median R 0.0049 (inverted, n = 543)
  vs 0.0100 (aligned, n = 1,691), **Cliff's delta -0.148** (p 2.0e-07); cohort-centred delta -0.140 (p 9.1e-07).
  **The clean evidence is the sign reversal of the M-R relationship within levels: Spearman rho(M, R) = -0.460
  (inverted), +0.393 (aligned), +0.116 (neutral).** Two caveats: that reversal is partly definitional (stated in
  advance - R and M_disease share the disease direction), and a confound control added afterwards shows the flag
  separates **M** more strongly than R (delta -0.351) while the R-level gap **reverses** within units (21 units in
  both kinds of cohort; paired difference +0.0072, 8 of 21 negative). Most of the pooled R gap is composition.
- **V1 ESR1 x breast: FAILED.** BRCA-LumA passes (fulvestrant -> **neutral**, M_disease -0.0322, percentiles
  4.5 / 3.0). **BRCA-LumB FAILS: fulvestrant -> `reference_aligned`** (M_disease -0.0419, sign-shuffle 1.5), which
  asserts R IS readable there - the one verdict the pre-registered criterion forbade. Her2 and Basal: neutral.
  **Nothing was redefined, retuned or reworded.** Two readings are recorded in DECISIONS and neither is chosen; the
  second (that the criterion may have tested the sign of R when the real breast problem was the ordering across
  cohorts) was constructed after the result and is labelled a post-hoc rationalisation.
- **Design limitation the failure exposed, deliberately NOT fixed now.** M_disease carries the drug's annotated
  action sign, so on identical reference data fulvestrant (`lower`) is `aligned` at -0.0419 while estradiol
  (`raise`) is `inverted` at +0.0419. The question "is this pathway active in this reference tissue" belongs to the
  (target, cohort, contrast) triple and must not depend on which drug asks. Proposed fix for the user: compute it on
  a fixed inhibition convention and apply each drug's orientation only when the flag is read.
- **Abstentions (2,262 rows): 1,326 `ambiguous_action`, 858 `below_floor_M_unit`, 78 `no_usable_genes`.** The first
  group is **every SERM row** (tamoxifen 156; raloxifene, bazedoxifene, toremifene, clomifene, estramustine 78
  each), so **the flag abstains on exactly the drug class where the inversion was first noticed.** Coverage gap,
  stated not worked around.
- **Cohort R baselines (companion 1): Basal is high.** Median R per cohort: median over 55 cohorts 0.0068, top
  TCGA-LUSC 0.0245, **BRCA-Basal 0.0167 (4th of 55)**, bottom TCGA-SKCM -0.0082. Breast: Basal 0.0167, Her2 0.0079,
  LumB 0.0068, LumA 0.0050 - Basal's baseline is 3.3x LumA's. Fulvestrant centred: Basal +0.0325, LumA +0.0048, so
  the baseline explains part of the backwards breast ordering, not all of it. `R_centred` now travels on every row.
- **Reference-tissue specificity tau (companion 2):** median +0.105 (IQR -0.057 to +0.235); aligned 0.149, inverted
  0.101, neutral 0.093; breast cohorts 0.347 and the ESR1 footprint itself 0.49, BRAF/SKCM 0.277, ERBB2/Her2 0.224.
  Reported as a column; nothing is gated on it.
- **Property to read the numbers by:** the disease z is much larger than a LINCS z, so the frozen `tanh(z/2)`
  saturates on a median 73.5% of footprint genes and M_disease is close to a signed sign-concordance.
- **Not run:** A5; no CTRP / PRISM / GDSC value was read.


## 2026-11-16: D131 — the reversal score's `1/sqrt(sum w^2)` normalisation, audited (nothing changed)

`scripts/d131/{norm_inflation,impact}.py` -> `results/lincs/d131/`, logs `logs/d131_*.log`, 9 s + 4 s.
Audit only: no production score changed, nothing re-ranked or re-tiered, no fix shipped. Decision rows `D131`
and `D131 result`.

**The error.** `79e_drug_reversal.py:280-285` divides the weighted sum of regulator activities by
`sqrt(sum w_i^2)`, which is the SD of that sum only if the activities are independent unit-variance terms.
Regulons overlap, so they are not. **The overlap is real cascade biology and is not the problem** — if TF A
activates TF B their regulons share targets, and when both are active those shared targets moving is evidence
for both. The fix evaluated keeps all of that in the numerator and corrects only the scaling:
`sqrt(w' Sigma w)` in place of `sqrt(w' w)`.

- **Inflation factor, exact** (it is identically the SD of the shipped `rev` across a cohort's drugs, since the
  score is built to have SD 1 under its own assumption): over 55 cohorts, **median 10.21**, 5–95% 4.21–14.34,
  min 4.10 (PHEO-kinase-signaling), max 15.46 (TGCT-SEM-C1). **The factor varies 3.8-fold between diseases**,
  so a raw `rev` value was never comparable across diseases.
- **Decomposition** (D126-B1 activity dump): correlation part median 7.12 (4.43–13.72), per-regulator spread
  part median 1.25 (0.89–1.72). The independence assumption is the dominant cause.
- `n_eff` of the regulator activities **9.7 of a median 1,052 master regulators**, against an iid sampling
  reference of 555 — real structure, not a small-sample artefact. Read it as *this sum's variance behaves as if
  it were about ten independent terms*, **not** as *there are only ten regulators' worth of information*.
- **The shipped rankings and tiers are not affected.** `norm` is one scalar per cohort, constant across its
  drugs, so every rank and every within-cohort robust z is unchanged. Verified by rescaling each cohort's `rev`
  by an arbitrary positive constant: **55 of 55 cohorts identical** in `rev_pct`, `rev_z` and the
  `(rev > 0) & (rev_pct >= 90)` tier mask. The same argument covers `spec`, `spec_sib`, `spec_pct`, `rank_rev`,
  `rank_spec` and every `tier_*`; the only absolute test on `rev` anywhere is `rev > 0`, a sign test.
- **Affected:** the magnitude of `rev` / `rev_raw`, `rev_sd` and `null_sd_regset` in the summary tables, the
  printed `top_tier1_rev` strings, and any cross-disease or absolute reading of a `rev` value.
- **`prolif_share` moves, but not through the normaliser.** 79i multiplies each part back by
  `sqrt(sum w_part^2)` (D119), so `a` and `b` are raw summed reversals and the share is invariant to `norm`.
  Its exposure is that it compares two raw sums over regulator sets of different size and internal correlation.
  Under a covariance-aware comparison (`k = sd_cov(prolif)/sd_cov(nonprolif)`, median **0.419**, range
  0.041–3.633): median share **0.301 -> 0.503**, **18,050 of 67,174 rows (26.9%) cross 0.5**, and
  `prolif_driven` goes **968 -> 1,539** with **597 rows changing the flag**. Applied to the shipped
  `prolif_share`, so the reconstruction cannot leak in. **Reported as a design question for the user, not a bug
  fix** — which comparison is right is a judgement about what "share of the reversal that is proliferation"
  should mean.
- **VIPER pleiotropy, fact only:** `msviper()` and `viper()` default to `pleiotropy = FALSE` (installed viper
  1.44.0), `shadow()` is a separate post-hoc call, and `50_viper.R:138,148` passes none of them, so it is off in
  production. **Not proposed as a fix for anything here**; it addresses MR inference quality, not variance, and
  enabling it would change production inference. A separate question for the user if wanted at all.
- **Largest remaining exposure, found and not quantified:** `79d_drug_viper.R:96` runs
  `viper(P, reg, method = "none", ...)` with **no `nullmodel`**, so every drug x line regulator activity is an
  analytically normalised aREA NES standardised under the same independence assumption. The disease side (50,
  52, 53, 64, 65, 66, 36, 95, 61) all pass an empirical normal-split null and are clean; the drug side is the
  arm without one. Unlike the `rev` scale this does **not** cancel in a within-cohort rank, because it varies by
  regulator and so changes the weights' relative contributions.
- **Reconstruction honesty:** the reconstructed `rev` matches the shipped values at Spearman >= 0.9989 per
  cohort but not exactly (max deviation 0.19–1.68 on a scale whose SD is ~10), from a small difference in which
  set-tagged regulators enter the consensus mean. The exact inflation does not use the reconstruction at all,
  and the `prolif_share` impact is computed from the shipped values plus a per-cohort SD ratio.
- **Abstentions:** `rev_extra` (79i:163-166, the D91 exception, the one place the normaliser varies by drug) was
  not evaluated. No CTRP / PRISM / GDSC value read. A5 not run. No production file written.

## 2026-11-16: D130 — correlation-preserving nulls; the ESR1 control fails its own corrected null

`scripts/d130/correlated_null.py` -> `results/lincs/d130/`, log `logs/d130_correlated_null.log`, seed 130,
1,000 draws per row (200 also reported and agreeing). Reporting only: M, R, every frozen parameter and every
production ranking are untouched; nothing was re-ranked or gated. Decision rows `D130`, `D130.1`, `D130 result`.

**The method.** Every null in `scripts/d128/nulls.py` shuffles genes as exchangeable units. The corrected
reference replaces the uniform permutation with one drawn from the empirical gene-gene correlation (CRP: draw
`u ~ N(0, R_hat)` as `Z'g/sqrt(m-2)` over all 9,678 LINCS profiles of the same arm, excluding the profile being
scored, and re-assign the observed values by the rank order of `u`). With independent genes it reduces exactly
to the existing null. A self-test then showed a rank permutation is blind to any shift common to every gene —
and that common component is most of why `n_eff` is small — so **D130.1** added CGN, the unconditional
covariance null, rescaled to the observed response magnitude. Both are reported for every row.

- **The ESR1 / fulvestrant control fails.** Specific arm, MCF7, three `pert_id`s. Sign-shuffle z
  **2.81 / 2.92 / 1.52 -> 1.45 / 1.64 / 0.86** (CRP). Gene-label **2.69 / 2.79 / 1.56 -> 2.35 / 2.39 / 1.33**.
  Under CGN, **1.09 / 0.83 / 0.61** (percentiles 86 / 78 / 72). The parent session predicted "roughly
  0.9 / 0.9 / 0.5, i.e. all non-significant" — **confirmed**. Raw arm the same but harder
  (3.53 / 3.64 / 2.68 -> CGN 0.92 / 0.57 / 0.38). **M did not change: 0.048 / 0.040 / 0.025 before and after.**
  What was wrong was the yardstick, not the score.
- **SERM set** (reported separately per D0.2): 7 of 12 rows had old sign-shuffle z >= 1.96; **1 of 12 does under
  CRP**. Only toremifene keeps CGN z above 2.
- **A4 actions 2-5, 164 scored rows.** |z| >= 1.96: sign-shuffle 9/5 -> **CRP 1/0**; gene-label 10/6 -> CRP 5/2;
  CGN 2/1. **All 17 rows D128-A4 marked `confident_positive` lose it.** The BRAF wave D128-A4 called "the only
  informative wave" falls from z 1.82–1.95 to CRP 0.50–0.52 / CGN 0.78–1.28.
- **The two estimation methods disagree, and the cheap one is harsher.** `sqrt(n/n_eff)` is right for an
  *unsigned* mean; M is sign-weighted and the correct analytic factor is `sqrt(e'Re/n)`. Median on the A4 rows:
  VIF **3.47** vs analytic **1.57**, ratio **1.87**. The D130(2a) prediction is confirmed. On fulvestrant it
  happens to be close (3.18 vs 2.65) because that footprint's signs are 85% one-directional (`e_bar` -0.695),
  exactly the case where the unsigned formula is nearly right. **Both shortcuts under-state the surviving z**
  (median |CGN z| 0.62 vs 0.45 analytic and 0.22 VIF) while tracking it in rank (r 0.871 / 0.849). **Trust the
  empirical CGN null; neither shortcut is a substitute, and neither was kinder to ESR1.**
- **The A3 population-wide nulls were never run** — the landmark criterion was not met, so `a3_aggregate.py`
  skipped them by a rule fixed in advance and `a2/nulls.parquet` does not exist. There is nothing there to
  correct. A seeded 1,000-row sample of the population stands in; see the D130 population addendum.
- **Both D128 nulls are affected**, and the sign null worse than the gene-label null (median correction 1.98 vs
  1.21) — the instruction not to assume only the gene-label null was hit was right.
- **Repo-wide inventory** (reported, nothing else re-run): affected gene/module permutation nulls at
  `79e:289-291` (diagnostic only, and already flagged in its own header), `94_replicate_brca_axes.R:61` (a GATE,
  v1 only), `82_breast_axes.R:41` (a GATE, v1 only), `91_replicate_crc.R:81,88` and `:136` (one leg of D33/D34);
  affected parametric-over-genes at **`79h_prolif_regulators.R:44-49`** — a hypergeometric gate (`prolif_linked`)
  **inside the live v2 pipeline** — plus `79d_drug_viper.R:96`, `63:39`, `cameraPR` at `31:149` / `37:59`,
  `62:103,110`, `67:40`, and `|NES| > 1.96` at `50:155`, `61:182`, `66:27` (D81) and `d93_decide.py:62-63`.
  **Not affected:** every msviper / `viperSignature` normal-split null (50, 52, 53, 64, 65, 66, 36, 95, 61),
  which permutes SAMPLES and was built to preserve gene-gene correlation; M3C (80, 91); all ARACNe bootstraps
  (43, 45, d94); all plate / cluster / patient permutations; and **the entire drug-side validation layer** (79k,
  79l, 79m, 79n, 79y, 79y2, 79t2, essential_check, d124, **d127**) — so the D127 biomarker pairs, the primary
  check this project relies on, are **not** implicated.
- **Naming hazard:** `drug_reversal_summary*.tsv` already has an `n_eff` column that is `(sum w^2)^2/sum w^4`
  (~300), a weight-concentration count, not this one. And `n_eff` here means *this sum's variance behaves as if
  it were n_eff independent terms* — **not** that the data hold only n_eff genes' worth of information. The
  correlation is real regulatory biology, which is why the corrected null preserves it.
- **Abstentions / gaps:** A5 not run; no CTRP / PRISM / GDSC value read; no site outside D128 re-run, so the
  inventory magnitudes are structural readings, not measurements; the `degree_matched_target` null stays
  degenerate and CRP does not fix it; CGN assumes Gaussianity given the covariance, which heavy-tailed LINCS z
  violate; `R_hat` is second-order only, so non-linear gene dependence is still under-corrected.

### 2026-11-16 - D130 population addendum: the naive null is ~3.5x too liberal, the corrected ones are near nominal

The A3 population-wide nulls were never run (the landmark criterion was not met), so a seeded sample stands in:
1,000 of the 2,017 distinct (unit, `pert_id`, line) M-scored rows of the frozen `all_bing` arm, 215 units, 1,000
draws per null -> `results/lincs/d130/population_sample_corrected.tsv`. Every M reproduced; 0 rows skipped.

- **Share of rows at |z| >= 1.96** (should sit near 5% two-sided under a calibrated null): naive sign-shuffle
  **18.0%**, naive gene-label **17.1%**; corrected **CRP sign 1.8%, CGN 6.9%, CRP gene-label 7.9%**. CGN and the
  corrected gene-label land near nominal; the corrected sign-shuffle lands below it, exactly as D130(1b) warned
  in advance that rank-mapping a +-1 vector would. **The correction is a calibration fix, not a blanket
  deflation** — that is a stronger claim than "the z's fell".
- Median correction factors: sign-shuffle **2.52x**, gene-label **1.25x**, CGN **1.34x**. The sign null is again
  the worse of the two.
- Footprints: median 1,070 genes, `n_eff` 45.6, `n/n_eff` 22.7.
- **The two estimation methods disagree as predicted:** median VIF factor **4.76** vs analytic **1.68**, ratio
  **1.81** (IQR 1.48–2.96); `e'Re/n` median 2.95. Both under-state the surviving z (median |CGN z| **0.772** vs
  0.467 analytic and 0.248 VIF) while tracking it in rank (r 0.780 / 0.713). **Use the measured corrected null;
  `sqrt(n/n_eff)` is wrong in form for a signed statistic and too harsh in practice, and must not be used to
  rescale a published z.** 2 of 1,000 rows had a non-finite factor and are excluded from the correlations only.
- Nothing was gated, ranked or promoted.

### 2026-11-17 — D133: the activity-dimensionality collapse is not regulon target overlap
Audit only (`scripts/d133/`, `results/lincs/d133/`, `logs/d133_overlap.log`, 22 s). Five cohorts spanning the
full D131 inflation range (PHEO-kinase-signaling 4.10, PCPG-PHEO 4.13, TCGA-LIHC 10.21 = the exact median,
UCEC-EEC 10.26, TGCT-SEM-C1 15.46), four ARACNe networks, group and stratum levels. Decision rows D133
(pre-registration) and D133.1 (the user's circularity amendment) were committed before anything was computed.
- **Gate passed exactly.** The reproduction of `79e`'s per-cohort drug aggregation + `consensus()` matches the
  shipped `d126_act` dump with **max |recon − dump| = 0.000e+00** and element-identical regulator lists, and
  reproduces D131's `n_eff_activity` to four decimals. The consensus-regulon correspondence is confirmed, not
  assumed, and D131's number is independently confirmed.
- **Overlap is strongly predictive where it exists and almost never exists.** Mean activity correlation rises
  0.01 → 0.21 → 0.49 → 0.66 → 0.78 across Jaccard bands 0 → (0.01, 0.05] → (0.05, 0.10] → (0.10, 0.20] →
  (0.20, 0.50]. But **54–63% of all regulator pairs share zero targets** and only ~0.3% exceed Jaccard 0.05.
- **Zero-overlap pairs are not near-uncorrelated:** mean |correlation| **0.214–0.263** against an all-pairs mean
  of 0.247–0.288. About a tenth of them exceed |0.45|.
- **Overlap explains a fifth to a third:** linear `R²` on the weighted cosine 0.160–0.259 (Jaccard 0.094–0.222),
  50-bin nonparametric 0.193–0.355.
- **Headline: removing the overlap-predicted component barely moves the dimensionality.** `n_eff` →
  residual `n_eff`: 9.75 → 12.72, 9.70 → 12.68, 10.36 → 14.48, 9.83 → 12.47, 7.97 → 10.04, against an iid
  reference of 531–582. Pre-registered bands were ≥ 100 (annotation redundancy) and < 30 (real structure):
  **all five land at 10–14.5.**
- **Cleanest form:** a maximal *literally disjoint* regulator set (40–53 regulators) still has `n_eff` 10.8–14.6
  against an iid reference of 38.7–50.7.
- **The circularity arm (D133.1).** On the common regulator subset, `n_eff / n_eff_iid` is ARACNe 0.030–0.037 vs
  curated CollecTRI 0.073–0.098: the learned set collapses **2.1–2.7× more**, so some circularity is real — but
  CollecTRI's own ratio is still a ~12-fold collapse. Inference-circularity is a factor ~2.3 of a 30–60-fold
  collapse. And the cross-check runs the wrong way for the overlap hypothesis: **CollecTRI regulons overlap more
  (23–26% disjoint pairs vs 65–75%) yet their activities correlate less** (mean |r| 0.13–0.17 vs 0.19–0.28).
- **Cascade check.** OmniPath direct-edge pairs are more correlated than pairs in general (+0.12 vs +0.05 raw),
  but **matched on overlap bin they are *less* correlated** (−0.071 to −0.094, permutation p < 0.0001; path ≤ 2
  the same). A documented regulatory relationship buys no correlation beyond shared targets. An OmniPath edge is
  a documented relationship, not causality in these data, and nothing here shows a cascade operating.
- **Reading (ii) of D133(10) holds.** Not annotation redundancy. The one thing the test cannot separate is stated
  in D133 result (11): every NES is an enrichment on the same ranked drug profile, so a dominant component in the
  signature would also do this. That needs a gene-level comparison, not run here.
- Nothing was changed, re-ranked, re-tiered or gated.

### 2026-11-17 — D134-S0: the target-recoverability feasibility audit (no retrieval computed)

Stage S0 of [`SPEC_target_recoverability.md`](SPEC_target_recoverability.md). The question the spec asks is
whether a drug's annotated protein target can be recovered from its transcriptional response. S0 does **not**
answer it and does not try: it converts the spec's annotation-level counts into data-level counts and asks
whether the test can be run at all. No similarity between two different compounds' profiles was computed, no
ranking, no AUROC, no null, no outcome label. Definitions were fixed in `docs/DECISIONS.md` row `D134-S0`,
committed and pushed before the script ran, so the output cannot have been fitted to a result.

`scripts/d134/s0_feasibility.py` (threads 8, 2.8 s) → `results/lincs/d134/`.

**Section 2 reproduces, with two corrections.** On the fixed key (unit = `tid`, 474 of 497 rows non-empty), the
spec's headline numbers come out exactly: 307 compounds, 193 units, 82 / 59 / 30 units with ≥ 2 / 3 / 5
compounds, 111 singletons. Distinct gene targets are **377, not 387**. The coverage line ("median 4 active cell
lines (IQR 1–8); 206 compounds in ≥ 3 lines, 163 in ≥ 5") **does not reproduce** under any of three readings;
the one retrieval can use — distinct lines per `cmap_name` — is median 3 (IQR 1–7), 131 of 247 in ≥ 3 lines,
91 in ≥ 5. Coverage is about a third thinner than stated.

**The attrition.**

| step | testable units | pairs | compounds | lines |
|---|---|---|---|---|
| annotation (spec Section 2) | 82 | 1,061 | 247 | — |
| every compound profiled | 82 | 1,061 | 247 | — |
| **+ pair must share a cell line** | **71** | **651** | **229** | **37** |
| + analog collapse (InChIKey skeleton) | 71 | 647 | 229 | 37 |
| + `nz2` ≥ 200 (illustrative) | 67 | 488 | 202 | 33 |
| + `nz2` ≥ 400 (illustrative) | 49 | 239 | 141 | 26 |
| + `tas_max` ≥ 0.6 (illustrative) | 48 | 206 | 103 | 24 |

**Shared lines are the real constraint.** Requiring that two compounds were actually profiled in the same line
removes 11 units and **39% of the pairs**. The units lost are small and nearly all non-oncology (CACNA1C/D/F/S,
NR3C2, ACE, HTR4, SLC12A3, HTR1D, the GRIN complex, ADRB1, HTR1A, the RAR/RXR complex, F10).

**Analog collapse removes four pairs, and that is a finding about the control.** 307 compounds carry 302
skeletons; only retinoids, betamethasone/dexamethasone, fluocinolone/fluocinonide and levonorgestrel/norgestrel
collapse. The first 14 characters of an InChIKey hash *connectivity*, so two statins or two tinibs are **not**
analogs by this rule — and the compound set is already keyed by `cmap_name`, so near-duplicates were gone
already. Spec control 5.1 as written is a near-duplicate check, not an analog check, and its pre-stated reading
cannot be run: excluding 4 pairs of 651 cannot move an AUROC. S1 must replace it or withdraw it.

**Strength.** TAS is available and was used alongside the gene count. The universe already carries a floor:
`79c` keeps only `tas ≥ 0.4` signatures, which had already cut 35,793 high-quality (compound, line) cells to
9,678 — **73% removed before S0 saw the data**. Over the 1,445 target-annotated profiles, `nz2` (genes at
|z| ≥ 2, raw arm) has median 445 (IQR 238–984); on the `specific` arm it falls to 259, so removing the generic
response removes ~42% of them. `tas_max` median 0.57. Floors to `nz2` ≥ 200 are nearly free; past that the set
halves. **No floor was chosen.** One caution for S1: the best-powered targets are the *weakest* responders —
NR3C1 (22 effective compounds) median `nz2` 246 / `tas_max` 0.46, HTR2A (13) 295 / 0.50, DRD2 (10) 372 / 0.49,
all below the overall medians — so a strict floor deletes the targets the map most depends on.

**Chance.** 563 (line, unit) cells carry ≥ 2 of a unit's compounds. Against the full profiled pool of the line,
chance is median 0.0034 (p95 0.031); **16 cells sit at ≥ 0.05**, including EGFR in three lines whose entire
profiled pool is 2–3 compounds, where chance is 1.0 and retrieval is arithmetically guaranteed. Against the
annotated-only pool, median 0.022 and 115 cells at ≥ 0.05 — the two pools differ about six-fold and the spec
does not say which it means.

**Verdict: S2 is worth running.** 71 units, 229 compounds, 37 lines, median 4 effective compounds and 6 lines
per unit, 53 units at ≥ 3 effective compounds and 29 at ≥ 5. The premise survives its own audit. It still says
nothing about the 111 singletons. Nothing was changed, re-ranked, re-tiered or gated.

## D135-S0 — can a drug's profile be checked against a knockout of its own target? Feasibility and QC (2026-11-17)

Stage S0 of [`SPEC_drug_genetic_concordance.md`](SPEC_drug_genetic_concordance.md), pre-registered as decision
row `D135-S0` and committed before anything ran. **Feasibility and QC only: no concordance, no drug-to-
perturbation comparison, no positive control, no AUROC, no null, no threshold chosen, nothing re-ranked.**
Scripts `scripts/d135/s0a_extract_perturbations.py` (`celligner` env, h5py, 61 s) and
`scripts/d135/s0b_feasibility.py` (`repurposing` env, 20 s), threads 8; outputs in `results/lincs/d135/`.
The 6.5 GB and 11.8 GB gctx matrices were read in 10,000-column blocks and never loaded whole; `trt_sh` and
`trt_sh.cgs` share a file and it is read once.

**Spec Section 2 does not hold as written.** Its three rows are computed under three different filters. The
CRISPR row (208 targets / 837 signatures / 8 lines) reproduces exactly under `is_hiq == 1` on the drug-matrix
lines. The shRNA-consensus row (209 / 1,627 / 10) only reproduces with no `is_hiq` and no line restriction;
under the CRISPR row's filter it is **174 / 515 / 9**. The raw-shRNA row (207 / 6,118 / 12) does not reproduce
under any of sixteen readings — closest 207 / 6,136 / 13 and 207 / 6,109 / 11 — and is reported unreproduced
rather than guessed. And **the drug matrix spans 73 cell lines, not 87**; 87 is close to the 86 lines the
compound set has *before* `79c` applies its `tas ≥ 0.4` floor.

**The paired count is much smaller than the coverage count, as expected.** For CRISPR: 377 target genes → 280
with any `trt_xpr` signature → 279 on a drug-matrix line (8 lines) → 1,936 (gene, line) perturbation cells →
1,049 where an annotated drug was also profiled in that same line → **2,350 paired (drug, gene, line) cells,
201 drugs, 255 genes, 7 lines**. Consensus shRNA gives 2,250 cells / 178 genes / 9 lines, raw shRNA 2,415 /
183 / 10. Five lines carry 389–428 CRISPR cells each; the seventh carries 22.

**The genetic side is much weaker than the drug side, which nobody had written down.** Per-signature `nz2`
(genes at |z| ≥ 2) has median 51 for CRISPR, **14** for consensus shRNA and 88 for raw shRNA; `tas` medians are
0.14 / 0.10 / 0.24. The drug side already requires `tas ≥ 0.4` per signature, and **only 10.5% of CRISPR,
4.2% of consensus-shRNA and 43.9% of raw-shRNA paired cells have a perturbation that would clear that same
bar.** Candidate floors were costed, none chosen: at `nz2 ≥ 100`, CRISPR keeps 1,258 of 2,350 cells but
consensus shRNA keeps 175 of 2,250. The arm the spec designates as the seed-effect mitigation is also the
weakest-responding arm by a factor of six.

**The target-mRNA QC finding, which is the point of the stage.** Per-signature pass rates at z ≤ −1 are
CRISPR **22.5%**, consensus shRNA **19.0%**, raw shRNA **28.0%**. Those numbers are nearly meaningless on
their own, because the determining variable is the L1000 feature space and not the technology. On **landmark**
(directly measured) genes CRISPR passes **70.4%** (target z median −2.65), consensus shRNA 54.8% and raw shRNA
54.8%. On **best-inferred** genes the perturbed gene's own z sits at zero (median +0.01) and the pass rate
collapses to **7.6% / 2.8% / 8.8%** — at or below a shrunken null, i.e. no information at all, because that
value is imputed from the landmarks rather than measured.

The two arms were kept separate, as pre-registered. For shRNA, reduced mRNA is the direct expected readout, so
a landmark failure is real evidence the knockdown did not take. For CRISPR, the pre-stated caution was that a
frameshift transcript escaping nonsense-mediated decay leaves mRNA near normal while the protein is gone, so a
CRISPR mRNA failure must never be counted as a failed perturbation. **That caution stands, but it is not what
limits the measurement here: on measured genes CRISPR passes *more* often than shRNA, the opposite of the
feared ordering.** NMD escape shows up in the tail (29.6% of landmark CRISPR signatures miss z ≤ −1), not in
the bulk. Usable conclusion: mRNA QC is a legitimate perturbation filter for **28 CRISPR genes and 48 shRNA
genes** and for nobody else; applied across the board it would delete ~90% of inferred-gene perturbations for
an artefact of imputation. No cell was excluded on QC here.

**Where the testable set sits in the exculpatory strata (CRISPR, of 2,350 cells):** common-essential target
20.1%, target not expressed in the line 37.2%, drug with more than one annotated target gene **90.7%**, weak
perturbation (illustrative `nz2 < 100`) 46.5%, weak drug profile 5.7%, only one drug timepoint available
69.5%, target mRNA not on the gene space 20.1%. Limitation 2 (non-catalytic scaffolding roles) is **not
computable** — no curated list exists in the repo — so it is reported missing and every clean count below is
an upper bound.

**The bottom line.** Cumulatively for CRISPR: 2,350 paired → 1,878 after dropping common-essential targets →
1,003 after dropping targets not expressed in the line → **80** after dropping multi-target drugs → 70 after
dropping weak drug profiles → **46** after an illustrative perturbation floor. **46 cells, 24 drugs, 14 genes,
6 lines**; 121 cells / 25 genes if drugs with two annotated genes are allowed; 80 cells / 18 genes with no
strength floors; 18–23 cells over 3 genes if the target must also be landmark and pass mRNA QC. Consensus
shRNA ends at **4** cells, raw shRNA at 122. The single dominant cut is polypharmacology (1,003 → 80), which
is as much definitional as biological. The surviving CRISPR set is recognisable: ABL1 (nilotinib, ponatinib),
BRAF (vemurafenib, dabrafenib), EGFR (erlotinib, gefitinib, AZD-9291), TYMS (fluorouracil, floxuridine), ESR1
(fulvestrant, clomifene), VDR, PPARG, AGTR1, DPP4, FKBP1A, PDE5A, HRH2, PTGFR, ACE.

**A problem in the spec's own positive controls, found as a coverage count and not computed.** All three
families are present — HDAC1/2/3/6 (104 CRISPR cells, 5 drugs), PSMB (112 cells, 4 drugs), TOP2A (39 cells,
8 drugs) — but **TOP2A is 100% common-essential and PSMB 57–61%**, so two of the three controls sit inside the
very stratum Section 6 marks as exculpatory. If they pass, the method may be validated by shared cytotoxicity
rather than by mechanism, which would not license reading concordance for non-essential targets. HDAC (25%
essential) is the only control family mostly outside it. S1 must say in advance what such a pass licenses.

**Verdict: S1 and a positive-control S2 are worth running, with a narrowed claim.** CRISPR carries the scope;
46–121 clean-positive-capable cells over 14–25 genes and 6–7 lines is enough for a per-target map of a dozen
well-known pairs and nowhere near enough for a population statistic over 2,350 cells of which 91% are excused
in advance. Consensus shRNA should be reporting-only unless S1 deliberately lowers its floor. Nothing was
scored, re-ranked, re-tiered or gated.

## 2026-11-17: Mechanism-trace procedure, first application — two traces, two different break-points

`docs/PROCEDURE_mechanism_trace.md` Section 10, applied end to end to the two pairs it names:
**fulvestrant → ESR1 → 0 hops → regulon → MCF7 (ACH-000019) → BRCA-LumA** and
**vemurafenib → BRAF → ≤2 signalling hops → 14 regulators → A375 (ACH-000219) → TCGA-SKCM**.
Scripts `scripts/trace/t1_predict.py` … `t4_check_a4.py`; records in `results/lincs/trace/`
(`links_L1_L4.json`, `links_L5_L7.json`, `verdicts.json`, `trace_records.tsv`, rendered prose in
`trace_narrative.txt`). **No DECISIONS row: nothing was decided.** Nothing scored, ranked, nulled or
re-tiered; no production file touched. The seven-field records are the artefact and the prose is
rendered from them (Section 0).

**The freeze was real.** Commit `TRACE FREEZE` carries both predicted footprints and every expected
per-gene direction, and was pushed before any response profile was opened; `t2_observe.py` was written
afterwards and its `expected` strings were written from L1–L4, not from the data. The L6 concentration
was computed from the footprint structure alone and committed **before** `results/lincs/d128/a4/` was
opened at all.

**The two traces break at different links — the Section 10 test passes.**

| | fulvestrant → BRCA-LumA | vemurafenib → TCGA-SKCM |
|---|---|---|
| chain | L1 holds, L2 holds, L3 holds, **L4–L7 all weak** | L1 holds, L2 holds, L3 holds, L4 weak, L5 weak, **L6 breaks**, L7 weak |
| break-point | **L4 — attribution** | **L6 — concentration** |
| effective independent links (Section 3) | **3 of 7** | **4 of 7** |

**L6 found the hub without being told, and reproduces D128-A4 exactly.** From the A1 route index alone:
MYC claims **65.1%** of vemurafenib's 1,063 signed genes and **68.1%** of the observed |contribution|
(58.2% through genes no other regulator of the unit claims); the effective number of independent
regulators is **2.73 of 14**. Removing MYC's 692 genes leaves concordance at **48.8%**, i.e. below
chance, and the mean signed f(z) falls 0.0340 → 0.0137. A4's own columns are 0.651 / 0.6814 / 0.0137 —
identical (`results/lincs/trace/l6_agreement_with_d128_a4.tsv`). MYC is `prolif_linked=True` in the
CollecTRI set, so the proliferation/cytotoxicity alternative survives into the verdict.

**Fulvestrant breaks earlier and for a different reason.** Its footprint is one regulon (ESR1, zero
hops), so L3 contributes no signalling evidence and L6 cannot even run a leave-one-out. Section 5 is
what kills it: only **8.6%** of the 478 genes are claimed by ESR1 alone, a median of **12** CollecTRI
regulators claim each gene — and the pre-stated sub-prediction that uniquely-attributable genes should
move at least as well as shared ones **failed badly**: 39.0% (n=41) vs 52.9% (n=437). The genes that are
specifically ESR1's are the ones moving the wrong way.

**L5 is weak on both.** Sign concordance 51.7% (fulvestrant/MCF7) and 52.6% (vemurafenib/A375); the same
frozen footprint under every compound profiled in the line puts them at the **64.0th** and **82.3th**
percentile, so the generic-response alternative is not ruled out for either. This agrees with D130 from
a different reference (corrected-null z 0.86/0.61 and 0.52/0.95 respectively) — read only after L5 was
committed.

**Two findings about the data, reported rather than worked around.** (1) `d128/a2/gene_contributions.parquet`
is a **seeded 2,000-row sample** (`--gene-contrib sample`, `all_bing` only) covering 186 of 257 scored
units; fulvestrant × MCF7 is not in it. The observation was rebuilt from the frozen footprint and the same
profile matrix with `score.py`'s own formula and reproduces D128-A2's M to 0.00e+00 for both pairs.
(2) **L2's O component has a different source per pair**, which changes the independence accounting:
ESR1's O is a VIPER NES over the *same* CollecTRI regulon L3/L4 use (so L2, L3 and L4 are one source),
while BRAF's O is `logFC_vs_gtex` — the same contrast L7 uses, and **negative** (−0.62), so the
D128-A3.2 reference artefact is live for the melanoma trace and is carried into its verdict.

**Not done, and it matters:** the five negative controls of Section 7 were not run, so the procedure's
false-positive rate is unmeasured and **neither trace counts as a positive**. Section 10's own
instruction is that the controls come before any positive trace is reported.

### The negative controls were run, on ARACNe (D137, 2026-11-17)

They now are. The answer is that **the procedure does not work.**

Pre-stated before any number: a control "closes" when no link L1–L7 is `breaks` or
`not evaluable`, and the fraction of the five that closes is the false-positive rate.
**Three of five closed (0.60) — and both real pairs broke.** The random drug–cohort
pairing scored at the **99.9th percentile** of every compound profiled in its line,
the strongest apparent mechanism of anything tested, real pairs included. A
procedure whose negative controls outscore its positives is not a weak instrument;
it is reading something other than mechanism.

**Why switching to ARACNe mattered so much.** The same target in the same cohort
gets a nearly disjoint footprint from the two networks. Fulvestrant/ESR1 shares
**6 genes** between a 478-gene CollecTRI footprint and a 131-gene ARACNe one;
vemurafenib/BRAF shares 159 of 1063 vs 1223. The whole-footprint response flips
sign for both: fulvestrant +0.025 (64th percentile) on CollecTRI, −0.018 (40th) on
ARACNe; vemurafenib +0.034 (82nd) vs −0.023 (27th). The same code reproduces the
published CollecTRI numbers exactly, so this is the regulon, not a bug. **The
footprint is a property of the network, not of the drug or the target.**

**Two earlier claims are withdrawn.** First, "BRAF is ~65% MYC dependent" was
substantially an artefact of CollecTRI's regulon-size skew — its MYC claims 881
genes against a median regulon of 9, while ARACNe's median is ~100. On SKCM the
BRAF hub is MEF2D at 20.4% over 11 regulators (8.0 effective). Second, and
consequently, L6 did not "surface the MYC dominance unprompted"; it surfaced a
property of the network it was given.

**A domain limit found by the controls:** TOP1, TOP2A, TYMS, DHFR, GART, RRM1 and
TUBB have **zero** OmniPath route rows, so a classical chemotherapeutic acting
through its own effector cannot be traced at all.

**One control is weaker than it should have been.** The random pairing shuffled the
*cohort* and then picked a celligner-matched line for it, and MEK inhibition in
BRAF-V600E HT-29 is real pharmacology. A stronger version randomises the
drug–target–line triple together. The other two closures — a drug whose target is
not expressed in the line, and an alkylator annotated to a nuclear receptor — carry
no such excuse.

**Checked and not supported:** that ARACNe's ESR1 regulon is noise for want of
variance within a subtype cohort. ESR1's sd in BRCA_LUMA is 1.282, the 95.2nd
percentile of all genes.

---

## 2026-11-20 — the proliferation adjustment was removing the wrong half (D157)

The user asked: *"shouldn't we subtract the averaged out proliferation signal so that we don't
subtract the full proliferation that is specific"*. Checking it found the adjustment doing the
opposite of what its own header describes.

`79i` removes growth arrest as `A_adj = A − β·π`. In production `π` is read off the **specific**
profiles — the ones `79c2` builds by projecting out the generic drug response per
`(ModelID, pert_time)`. So the generic part is already gone at the gene level before `π` exists:
`pi_spec` has mean **−0.018 early / −0.014 late** and correlates **0.026** with the group mean. Of
`pi_raw`'s variance, **22.4% is generic and 77.6% drug-specific** — and the drug-specific part is
what was being subtracted.

Splitting that by the drug's cross-line mean within each window:

| window | drug-universal (cytostatic everywhere — caveat C4) | line-selective (candidate mechanism) |
|---|---|---|
| early | 51.6% | 48.4% |
| late | 40.5% | 59.5% |

About half of the discarded signal is a drug arresting growth *selectively* — which under the
discovery goal is a finding, not a nuisance. `PROLIF_MODE` now selects which part to remove;
`full` stays the default so nothing has moved. Full detail in
`docs/NOTE_proliferation_decomposition.md`.

### The D156 recalibration landed, and the window finding held

Re-ran the window scores on the 6 h floor of 0.25:

| | D155 (TAS 0.4) | D156 (6 h floor 0.25) |
|---|---|---|
| (cohort, drug) rows | 12,055 | **13,907** |
| early-only | 1,946 | **3,798** |
| both windows | 2,132 | **3,601** |
| Spearman(early, late) | +0.388 | **+0.385** |
| sign-discordant | 37.0% | **37.4%** |
| hits ≥ 90th percentile | 1,351 | **1,660** |

69% more both-window rows, same answer. Anti-estrogen rows went **20 → 53**, and fulvestrant in
BRCA-LumA now reads **96.3rd percentile early vs 60.2nd late**.

**Still not fixed:** across 33 BRCA anti-estrogen rows the median percentile in **basal** disease
(69.1 early / 89.3 late) is *higher* than in **luminal** (61.3 / 53.0). These drugs still score best
in the subtype they do not treat. The TAS floor bought coverage, not specificity — that is the D151
target-presence gate, still specced and unbuilt.

### β degraded after D156, and it turned out to be dilution

The pooled per-window β check looked worse than D154 (+0.594 / 28.7% sign flips vs +0.687 / 25.2%).
Stratifying the 2,300 both-window cells on whether their early condition would have cleared the old
0.4 floor separates the two explanations:

| stratum | cells | median early TAS | Spearman(β₆, β₂₄) | sign flip |
|---|---|---|---|---|
| **HIGH** | 1,336 | 0.500 | **+0.681** | **25.7%** |
| NEW | 964 | 0.323 | +0.316 | 40.5% |
| *D154 reference* | *1,373* | *—* | *+0.687* | *25.2%* |

`HIGH` reproduces D154 almost exactly on a re-derived sample, so the original finding was not a
quirk of its cells. The drift is the weaker signatures, and it does **not** reach the drug scores
(37.0 → 37.4%) because β is a per-regulator slope while the drug statistic aggregates hundreds of
them. `BETA_TAS_MIN` now allows fitting β on the clean subset and applying it everywhere.

## 2026-11-20 — Route 1 validates Route 2: primary targets yes, secondary targets no (D158)

The full-network per-condition VIPER finished on the 7 cohorts holding 97% of the validation pairs
(2,298–3,253 regulators scored each, 15.5–21.6 min per network). Route 1 then asked the question
directly: when Route 2 ranks every regulator by how much it opposes the disease, where does the
drug's *known* target land?

**2,328 (cohort, drug, window, target) pairs, 169 drugs × 147 targets, 8 cohorts.** Chance@10 = 0.35%.

| arm | recovery@10 | lift | primary | secondary |
|---|---|---|---|---|
| `none` | 1.50% | 4.3× | 2.21% | 0.13% |
| `universal` (D157) | 1.33% | 3.8× | 1.95% | 0.13% |
| `full` | 1.59% | 4.6× | 2.28% | 0.25% |
| `full_filter` | 0.39% | 1.1× | 0.40% | 0.38% |

### Primary targets — a real positive, the first in this project

| quantile | observed | expected | ratio |
|---|---|---|---|
| top 1% | **59** | 15.3 | **3.84×** |
| top 5% | 134 | 76.8 | 1.75× |
| top 10% | 213 | 153.5 | 1.39× |

It holds in both windows (early 1.69×, late 1.79×) and — the part that makes it credible — it
**scales with the evidence behind each score**: top-5% ratio 1.10 at one condition, 2.23 at 2–4,
**2.45 at ≥5**. Artifacts don't do that.

The recoveries are textbook pharmacology and are *not* proliferation-linked: **erlotinib → EGFR at
rank 3 of 2,806**, romidepsin → HDAC10 at 6, panobinostat → HDAC7 at 6, fludrocortisone → NR3C2 at 8
(early), etoposide → TOP2B at 11, belinostat → HDAC1, acitretin → RARA, fingolimod → S1PR1.

### Secondary targets — flat, and the hypothesis is not supported

| quantile | observed | expected | ratio |
|---|---|---|---|
| top 1% | 7 | 7.9 | 0.88 |
| top 5% | 31 | 39.7 | 0.78 |
| top 10% | 75 | 79.3 | 0.95 |

Median rank 1,548 of 2,853 — 53.2%, where chance is 50%. Relative binding strength does **not**
predict rank (Spearman +0.071, if anything the wrong sign). Entries that look plausible individually
— dasatinib → EPHA2 at rank 7, midostaurin → PRKCD at 9, regorafenib → TEK at 13 — are no more
numerous than chance draws from 793 pairs, so they are not claimed as signal.

### "Remove proliferation before or after" — measured, and neither matters here

Subtracting is a near no-op on rank: paired over the same 2,328 pairs, median rank change from
`none` to `full` is **−1 out of ~2,850**, improving 50.9% and worsening 46.4%.

Filtering is actively harmful at the very top, because **81–86% of the rank ≤ 10 hits are
proliferation-linked against a 12.9% base rate** — `full_filter` deletes them and collapses to
chance. But the top-1% primary signal survives filtering, because those recoveries are not
proliferation-linked. So the shallow hits are proliferation and the deep ones are mechanism.

**What this licenses:** the top ~1% of a Route 2 list is trustworthy for primary mechanism.
**What it does not:** any claim of secondary-target discovery from Route 2 as presently formulated.

## 2026-11-21 — The trial-outcome test: a 6.9× volume gradient and a flat null on outcome (D176)

D175 found that `trial_outcomes_v2.tsv` only ever covered `fda_rank` 1–10, so the test everyone
wants — *do top-ranked drugs do better in trials than bottom-ranked ones?* — had no control group.
D176 fetched the missing bands. `79o2` gained a `BANDS` env; bands are taken **per unit**, so
cohorts with different numbers of eligible drugs stay comparable.

**550 top / 572 mid / 550 bottom candidate pairs, 55 units, 9,481 ClinicalTrials.gov records.**

### Volume — a huge gradient, and it is the confound, not the result

| band | mean trials per candidate pair | share with ≥1 trial |
|---|---|---|
| top | **11.73** | 50.4% |
| mid | 3.65 | 41.6% |
| bottom | 1.71 | 32.9% |

Kruskal–Wallis **p = 5.1e-12**; top > bottom p = 4.9e-13. A 6.9× gradient — and it is **fame**.
Curated targets, binding data and alteration frequency are all richer for well-studied drugs, so
"highly ranked" and "heavily trialled" share an upstream cause. This number is not validation and
must not be quoted as such.

### Agreement — the gate that decides whether the contrast means anything

The top band was classified in an earlier session and mid/bottom now, so drift between passes would
map exactly onto the contrast being measured. **60 already-classified top-band trials were shuffled
unlabelled into the new reading queue**, and the veto thresholds were fixed *before* reading:
κ < 0.6, or a positive-rate difference > 0.10.

| | |
|---|---|
| exact outcome agreement | 0.883 (53/60) |
| binary positive-vs-not | 0.917 |
| **Cohen's κ** | **0.834** |
| positive rate, prior vs new pass | 0.517 vs 0.467 |

**Both thresholds met.** Disagreements are boundary calls — 3 positive↔not_reported, 2
negative→mixed, **0 positive↔negative reversals**. All 256 new readings carried a quote verified
character-for-character against that trial's own abstracts (**256 pass, 0 discarded**), and a chunk
integrity check found **0 foreign NCTs, 0 missing** across all four reader chunks.

### Outcome — the actual test, and it is null

Positive rate given ≥1 classified completed trial, with `drug_role` background/comparator excluded
because those trials do not test the drug:

| band | pairs | trials | positive rate |
|---|---|---|---|
| top | 136 | 803 | **0.484** |
| mid | 90 | 302 | **0.513** |
| bottom | 86 | 158 | **0.475** |

Kruskal–Wallis **p = 0.787**; top > bottom p = 0.436; Fisher on trial counts OR = 1.21, p = 0.153.
Not merely non-significant — **not even ordered**.

**What this licenses:** nothing new. **What it refutes:** the claim that the lead ranking predicts
clinical success. Conditional on reaching a completed trial, `layered_pct` does not.

Two limits, stated as limits and not as a rescue: 136 vs 86 pairs needs roughly a **15-point** rate
difference to detect, so a real 5-point edge is invisible; and conditioning on reaching a trial is a
**hard filter** — those pairs already passed human judgement, so little variance survives.

This does not conflict with D127, which asks whether the score tracks biomarker status *within* a
drug (7/7 pairs). Together they bound the claim: the pipeline matches a drug to the cohort where its
target is active; it is **not shown to pick winners**.

## 2026-11-21 — Corticosteroids in LUSC: two real defects, one unfixable assumption (D177)

The user: *"i notice that corticosteriods are highly ranked for lung squamous cell carcinoma which is
just incorrect."* They were right. LUSC's top 10 was **two genes** — 7 corticosteroids on NR3C1,
3 glitazones on PPARG — all with byte-identical components, because `layered_S` is a property of the
target unit, not of the drug.

### Defect 1 — the dependency inversion (removed)

For an agonist (`orient == "raise"`) the score used `100 − v` on the CRISPR dependency percentile.

| target | true dependency percentile | what the score used |
|---|---|---|
| NR3C1 (LUSC) | **7.7** — raw gene effect **+0.052**, knockout slightly *helps* growth | **92.3** |
| PPARG (LUSC) | **1.0** | **99.0** |

Verified exactly: of 5,098 scored rows, **751 matched `100 − dep`, 4,383 matched `dep`, zero matched
neither.** Knockout data cannot be inverted into activation evidence. The one dataset that could
supply it — LINCS overexpression — covers 222 genes and 8 cell lines and contains **none** of our 41
agonist targets, so the axis cannot be built. D is now dropped for `raise` units, not flipped.

### Defect 2 — the logFC fallback (removed)

`comp_overactivity` falls back to `logFC_vs_gtex` for genes absent from the master-regulator table.
For a tissue-foreign receptor that is automatic: GABRA1 and HTR1B are "lower in lung tumour than in
normal lung" because they are neuronal. **182 of 783 `raise` rows rested on it** — the same root
cause as the antidepressant/CNS pattern flagged earlier, not a separate problem.

### What no rescoring fixes — and a correction to my own first diagnosis

I framed this as "absence of evidence inverted into evidence." That is true for the dependency axis
and for the logFC targets. It is **false for NR3C1's overactivity: its NES in LUSC is −2.076, a
genuine measurement.** Within the reversal framework, a suppressed regulator plus a drug that raises
it is a legitimate hit. The framework assumes *reversing a regulator's activity reverses the
disease*, and for a regulator suppressed as a **consequence** of dedifferentiation (GR in squamous
epithelium) that is false. Nothing in the activity data separates the two cases.

### Shipped: flag, don't demote

`layered_n_axes` (1/2/3) and `layered_orient` (`lower`/`raise`) now ride on the split tables, through
79p into the FDA rankings, and are surfaced by 79q in the decode string and as filterable flags.
Census: **3,042 rows on 3 axes, 383 on 2, 231 on 1; 219 `raise`, 4,415 `lower`.**

**LUSC top 4 is now correct** — floxuridine/fluorouracil (TYMS), alpelisib (PIK3CA), palbociclib
(CCND1|CDK6), all 3-axis, all genuine LUSC lesions. **The corticosteroid block remains at slots
5–10**, and across all 41 cohorts the agonist class holds **92 of 410 top-10 slots (22.4%), up from
89 (21.7%)**. The fix did not reduce their share. That residue is expected under this choice: the
user chose flagging over a tuned demotion, because the penalty size would have been fitted to LUSC.

### The cost: D127 no longer meets its criterion

**36/49 (0.7347, p = 0.0007) → 33/49 (0.6735, p = 0.0066)**, median gap 15.1 → 14.2. **VERDICT NOT
MET** (share 0.6735 < 0.70). No threshold was retuned to recover it.

The entire loss is **one pair** — BRAF/MEK, THCA-BRAF-like-classical vs THCA-RAS-like-follicular,
3/9 → 0/9. Six of seven primary pairs are unchanged instance-for-instance. **D177 does not touch any
drug in that pair**: all 18 rows are `layered_n_axes = 3`, `layered_orient = lower`, shrink factor
exactly 1.0. The movement comes from percentile re-ranking as *other* drugs in those strata were
rescored.

And the pair was already failing. 3/9 is below chance (4.5/9), baseline median gap **0.0**, and the
margins are noise — vemurafenib and dabrafenib score **100.0 vs 100.0**, exact ties, because the
percentile saturates in both arms. So D127's headline fell without the lead score getting worse at
what D127 measures. The criterion was pre-set and is not met, and that stands — but it also exposes
a weakness in D127 itself: **a pair with no resolution contributes 9 of 49 instances (18%) of pure
noise.**

**Open, deliberately left to data:** whether raising a suppressed regulator ever helps is empirical,
not a scoring question. D178's phase-progression test covers **211 agonist pairs**.
