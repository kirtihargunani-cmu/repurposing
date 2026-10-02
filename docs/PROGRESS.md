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
  - BRCA GTEx normals female-only (`gtex_sex`; 79 donors). LumA/LumB still contain a few male tumors, so a sex term remains in those models.
  - 40 groups now have a protein-coding signature (58 contrasts).
  - Extended signature (protein-coding + lincRNA + antisense; `GENESET=extended`) in `results/step4_extended/` for mechanism work only; drug matching uses protein-coding. Protein-coding z identical between runs (Spearman >= 0.997); ~5.9k lincRNA + 4.5k antisense genes pass filterByExpr.
- Scope decision (PI, 2026-10-01): Steps 5-7 for all groups with a signature; Steps 8-10 for BRCA and COAD/READ. De-prioritize hematologic tumors (AML/ALL/DLBC): limited and more variable data. Steps 5-10 data requirements in `docs/PROJECT_PLAN_steps5-10.md`.
- Environment moved to R 4.5.3 / Bioconductor 3.22 when the Steps 5-10 packages were installed; Steps 3-4 rerun (ESTIMATE identical; z Pearson >= 0.9998 vs R 4.4).
- Hematologic normal references: recount3 / GDC search recorded in `docs/normal_reference_search.md`; parked per PI.
- Step 5 (regulons):
  - `scripts/40_id_map.py` -> `results/step5/id_map.tsv`: GENCODE v23 Ensembl -> current HGNC (Ensembl match, then symbol, unique previous symbol, unique alias). 97.4% of protein-coding genes mapped (522 unmapped, mostly retired clone names); 5,113 v23 symbols have since been renamed, so symbol joins would have lost them. `primary_for_hgnc` picks one v23 gene per HGNC record (176 duplicates dropped).
  - `tf_universe.tsv`: Lambert 2018 'Is TF? = Yes' (1,639; 1,634 in v23).
  - `scripts/41_regulons.R`: CollecTRI and DoRothEA A-C pulled from the OmniPath REST API (decoupleR/OmnipathR 4.0 wrappers fail against the current server; raw pulls saved in `data/raw/regulons/`), NF-kB / AP-1 complexes collapsed as in decoupleR; aracne.networks (25 TCGA cohorts, Entrez) mapped to Ensembl.
  - Coverage (`regulon_coverage.tsv`): CollecTRI 1,178 regulators / 42,601 edges (99.6% mapped; 614 with < 10 targets); DoRothEA A-C 429 / 32,286 (99.6%; 158 < 10); ARACNe ~6,050 regulators per cohort, ~1,500 of them Lambert TFs, median 95-97% targets mapped.
  - `regulons_core.rds` (CollecTRI + DoRothEA) in git; ARACNe tables and `.rds` in `data/processed/step5/` (1.2 GB). sha256 of both in `regulon_checksums.sha256`.

## Next
- Step 6: VIPER per group (expressed-TF filter, minimum mapped regulon size; ARACNe arm restricted to Lambert TFs).
- CPTAC-HNSCC (GDC CPTAC-3, open STAR counts, GENCODE v36; ~173 tumor / ~59 adjacent-normal RNA-seq files) as a separate HNSC replication cohort. Do not pool with Toil (different pipeline; normals are adjacent, not healthy).
- Decide: subtype splits for mixed cohorts; whether to pull DepMap; recount3 non-cancer diseases.
