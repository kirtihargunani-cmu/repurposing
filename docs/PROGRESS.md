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

## Next
- CPTAC-HNSCC (GDC CPTAC-3, open STAR counts, GENCODE v36; ~173 tumor / ~59 adjacent-normal RNA-seq files) as a separate HNSC replication cohort. Do not pool with Toil (different pipeline; normals are adjacent, not healthy).
- Decide: subtype splits for mixed cohorts; whether to pull DepMap; recount3 non-cancer diseases.
- Step 3 covariate table (purity/ESTIMATE, batch, sex, age).
