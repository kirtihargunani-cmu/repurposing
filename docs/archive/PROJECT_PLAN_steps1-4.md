> **Status 2026-10-31: original plan** received 2026-10-01, kept for reference. How the pipeline actually runs: [`PIPELINE_V2.md`](../PIPELINE_V2.md).

Step 1 — disease grounding
- MONDO (Monarch): mondo.obo/.owl or OLS4 REST. Per term: ID, label, synonyms, is_a parents/descendants, xrefs (EFO, DOID, NCIt, UMLS, ICD-O-3). Open, CC-BY.
- OncoTree (MSK): flat code list + tissue/lineage mapping (DepMap and Celligner key on OncoTree). Open REST API.
- UBERON: anatomy term for the tissue + parents. Open.
- GTEx / DepMap lineage vocabularies: their own tissue labels, to map UBERON -> each dataset's label.
- (optional) Open Targets Platform: EFO-keyed disease-target associations.
Build one table: mondo_id | label | descendants | efo | oncotree | ncit | umls_cui | uberon | gtex_tissue | depmap_lineage | arm(liquid/solid) | cohort_handles. Pin the MONDO release.
UMLS CUIs come from MONDO xrefs (no UTS license needed unless going to UMLS directly for repoDB mapping).

Step 2 — expression data
Solid: UCSC Xena "TCGA TARGET GTEx" (Toil): TcgaTargetGtex_rsem_gene_tpm (log2(TPM+0.001), Ensembl, ~19,100 samples); TcgaTargetGTEX_phenotype.txt; expected_count for limma-voom. Pull: disease tumors, GTEx normals of matched tissue, TCGA adjacent-normals (sanity check).
Non-cancer: recount3 — available_projects(), create_rse() (raw counts, GENCODE v26), transform_counts()/compute_read_counts(); colData() SRA metadata is free text, curate yourself. GTEx/TCGA also in recount3.
AML: [DECIDED: BeatAML 2.0 GitHub instead of GDC BEATAML1.0 + Tyner 2018 supp.] Normal marrow reference (BloodSpot 3.0 / underlying GEO) [DEFERRED]. No uniform recompute contains both BeatAML and normal marrow: reprocess both, or carry source as covariate. GTEx whole blood is not a marrow substitute. TCGA-LAML is the fallback (~150 tumors, no ex vivo response, can't serve Step 20 validation).

Step 3 — QC and confounders
Phenotype/clinical tables: batch/plate, site, source study, sex, age, sample type. Blast % from TCGA-LAML clinical / BeatAML clinical (AML purity). ESTIMATE (R, computed): stromal/immune/ESTIMATE score per solid sample. (optional) ABSOLUTE / TCGA purity consensus as cross-check. GENCODE GTF matching the matrix release -> restrict to protein-coding. Output: covariate table.

Step 4 — patient signatures
limma (or edgeR/DESeq2 from counts): moderated-t per gene -> signed z; model includes purity, batch/source, sex, age. Ensembl->HGNC map (GENCODE/biomaRt/org.Hs.eg.db) for display only; join on unversioned Ensembl IDs; same GENCODE release as the matrix.

Download queue: Xena TPM + phenotype (solid) · BeatAML (AML) · MONDO release + OncoTree codes · matching GENCODE GTF. Everything else is a package install.
