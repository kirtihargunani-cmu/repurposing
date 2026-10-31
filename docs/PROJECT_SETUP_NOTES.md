# New drug-repurposing project: setup handoff (2026-10-01)

> **Status 2026-10-31: historical** (setup handoff of 2026-10-01). Current state: [`PIPELINE_V2.md`](PIPELINE_V2.md), [`PROGRESS.md`](PROGRESS.md).

Written from a laptop Claude session so the server session can pick up the work. Read this first.

## Decisions already made by the user
- **New project folder, fresh start** (name not chosen yet; ask). It does not depend on the old
  `subset-repurposing` project on the laptop, but reuse its pattern: a `config/datasets.yaml`
  that pins each source's release, URL, access tier and a checksum after the first pull.
- **Data branches to pull now:** solid tumor (UCSC Xena Toil TCGA-TARGET-GTEx), AML (BeatAML 2.0),
  and recount3 (install the package only; pick specific projects later). **No** normal-marrow
  (BloodSpot/GEO) reference for now.
- **BeatAML source: BeatAML 2.0 from GitHub `biodev/beataml2.0_data`**, not GDC BEATAML1.0 plus the
  Tyner 2018 supplementary tables. It is open (CC-BY 4.0; "dbgap" in the filenames is not an access tier) and
  includes ex-vivo probit curve fits:
  - beataml_waves1to4_norm_exp_dbgap.txt  (expression, 707 samples)
  - beataml_probit_curve_fits_v4_dbgap.txt  (drug response, 166 drugs)
  - beataml_wv1to4_clinical.xlsx, beataml_drug_families.xlsx
  - Split on dbgap_subject_id: some subjects have up to 4 RNA-seq samples.

## Known URLs (verified in the old project)
- https://toil.xenahubs.net/download/TcgaTargetGtex_rsem_gene_tpm.gz   (1.3 GB gz, log2(TPM+0.001))
- https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz
- https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap
- expected_count matrix (for limma-voom): TcgaTargetGtex_gene_expected_count.gz on the same hub. The
  user's plan calls it `..._rsem_gene_expected_count`, so check the exact filename.
- https://github.com/biodev/beataml2.0_data/raw/main/<file>

## Pitfalls flagged
- Toil is **GENCODE v23**: the Step-3 GTF and the Step-4 Ensembl->HGNC map must be v23 too. Mixing
  releases silently drops genes. Join on unversioned Ensembl IDs.
- recount3 is GENCODE v26, so it differs from Toil. Carry this if the branches are ever combined.
- TCGA-LAML in Toil: ~173 primary blood samples, with no ex-vivo response data.
- GTEx whole blood is not a marrow reference.

## Software still to install (none present yet)
conda env with R + Bioconductor: limma, edgeR, DESeq2, recount3, biomaRt, org.Hs.eg.db, ESTIMATE
(ESTIMATE is on R-Forge, not CRAN/Bioc, so install it from there), plus Python (pandas, etc.).
Pull MONDO (pin the release; OLS4 also works), OncoTree codes, UBERON.

## Machine
MAIL-10: 64 cores, 503 GB RAM, 3x RTX 6000 Ada, ~497 GB free in home, no SLURM, outbound internet works.

## The user's full plan (Steps 1-4)
Condensed copy (decisions marked inline) in ~/PROJECT_PLAN_steps1-4.md
