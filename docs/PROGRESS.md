# Progress log

## 2026-10-01
- Project scaffolded at `projects/repurposing`.
- Conda env `repurposing` built (R 4.4.3; limma 3.62.1, edgeR 4.4.0, DESeq2 1.46.0, recount3 1.16.0, biomaRt 2.62.0, org.Hs.eg.db 3.20.0, ESTIMATE 1.0.13). Exact versions in `environment.lock.yml`.
- Raw data pulled (~2.9 GB) and pinned in `config/datasets.yaml`: Toil TCGA-TARGET-GTEx (GENCODE v23), BeatAML 2.0 @ `6e403e1`, MONDO v2026-09-01, UBERON 2026-10-01, OncoTree `oncotree_2025_10_03`, GENCODE v23 GTF.
- Note: Toil expected-count file is `TcgaTargetGtex_gene_expected_count.gz` (not `_rsem_`).

- Step 1 (Toil cohorts): `config/toil_disease_curation.tsv` (hand-curated) -> `scripts/10_ground_toil.py` -> `results/step1/disease_grounding_toil.tsv` + `mondo_descendants.tsv`. 40 cohorts (33 TCGA, 7 TARGET), all IDs validated against pinned MONDO/UBERON/OncoTree/GTEx labels; 10,530 tumor samples.
  - GTEx normal match: 13 good, 18 partial, 9 none (CHOL, DLBC, MESO, SARC, THYM, UVM, all leukemias).
  - Mixed-histology cohorts to split before Step 4: LGG, BRCA, CESC, ESCA, PCPG, SARC, TGCT, UCEC, TARGET-ALL.
  - TARGET-RT has no tumor samples in Toil (5 normals only); TARGET-CCSK has 13 tumors.
  - DepMap lineage column empty until DepMap is pulled.

## Next
- Review curation calls (esp. GTEx tissue choices and `partial` matches).
- Decide: subtype splits for mixed cohorts; whether to pull DepMap; recount3 non-cancer diseases.
- Step 3 covariate table (purity/ESTIMATE, batch, sex, age).
