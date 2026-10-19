# Spec: tissue composition of tumors vs normal tissue, per disease cohort

Status: **spec written 2026-10-18, not executed.** Companion to the Celligner normal-tissue spec on the `celligner`
branch (`docs/notes/SPEC_celligner_normals.md`). Decision D61 in [`DECISIONS.md`](DECISIONS.md).

## Question

For each disease cohort, how different is the **cell-type make-up** of the tumors from the healthy reference
(GTEx, v2 reference set) and from TCGA adjacent normal tissue? And how much of each tumor-vs-normal signature is a
change in **which cells are present** rather than a change **inside the cells**?

Why it matters: a drug signature measured in a cell line can only reverse changes inside cancer cells. Signature
genes that are high only because normal tissue contains fat, muscle or a different epithelium (e.g. GTEx breast is
fat-rich, GTEx uterus is mostly muscle) cannot be reversed by any drug acting on the tumor cells.

## Inputs

- Samples: one primary tumor per patient (Step 4 rules); the v2 GTEx reference set (`results/step4_v2/gtex_reference.tsv`
  from the v2 spec, Step 4a: one sample per donor per tissue, colon / kidney sub-region filter); TCGA adjacent normals
  (one per patient). The 40 groups with a GTEx signature; adjacent comparisons in the 18 groups that have them.
- Genes: the v2 gene set (mappability >= 0.8, no readthrough loci).
- Expression: TMM log-CPM (prior count 2) computed per group over tumors + its normals, as Step 4.
- Depends on: v2 Step 3 / 4a outputs (gene list, GTEx reference). Can run before the networks are rebuilt.

## Method

### 1. Cell-type scores per sample (two layers)

**1a. Immune and stromal populations: MCP-counter** (Becht et al. 2016, Genome Biology; R package `MCPcounter`,
pinned GitHub commit, installed into the `repurposing` env). Ten populations: T cells, CD8 T cells, cytotoxic
lymphocytes, NK cells, B lineage, monocytic lineage, myeloid dendritic cells, neutrophils, endothelial cells,
fibroblasts. Scores are on log-CPM. ESTIMATE stromal / immune scores (already in Step 3) are reported alongside.

**1b. Tissue-resident cell types: marker panels from the Human Protein Atlas single-cell type data**
(`rna_single_cell_type.tsv`, nTPM per gene per cell type; release pinned and registered in `config/datasets.yaml`).
- Marker definition, fixed before use: a gene is a marker of cell type k if its nTPM in k is >= 4x the highest other
  cell type (HPA "cell type enriched" rule), nTPM in k >= 10, and the gene is in the v2 gene set. Panel = the 20
  markers with the largest fold over the next cell type (fewer if fewer qualify; a panel needs >= 8 markers).
- Cell types scored: those present in the organs of our cohorts (e.g. hepatocytes, proximal tubule cells, collecting
  duct cells, alveolar type 1 / 2 cells, club / ciliated cells, breast glandular / basal cells, adipocytes, smooth
  muscle cells, skeletal myocytes, enterocytes, goblet cells, thyroid follicular cells, prostate
  glandular cells, keratinocytes, melanocytes, neurons, astrocytes, oligodendrocytes, Leydig / Sertoli / germ cells,
  adrenal cortical cells, pancreatic acinar / ductal / islet cells, endometrial stromal cells, ovarian stromal
  cells). The final list = every HPA cell type with a qualifying panel; recorded in `composition_panels.tsv`.
- Score = mean log-CPM of the panel genes (same form as the colon / kidney sub-region scores, D59).
- A proliferation score (Hallmark G2M genes) is reported but is not a cell type.

### 2. Differences per cohort

For each group and cell type: tumor vs GTEx, tumor vs adjacent, and adjacent vs GTEx:
- difference in median score and Cliff's delta (rank-based effect size, robust to skew), with age and sex removed
  within each arm first (linear, as in `scripts/33_normal_variance.R`);
- because marker scores are relative, values are compared **between arms for the same cell type**, never between
  cell types.
- GTEx scores are also reported after the v2 collection correction (Step 4a covariates set to the anchor), to show
  whether the correction changes composition estimates.

### 3. How much of the signature is composition

Per group, using the Step 4 tumor-vs-GTEx z (v1 now; v2 when built):
- Build a gene x cell-type matrix from HPA: each gene's share of its total nTPM in each cell type (cell-type
  specificity).
- Regress the signature z on that matrix across genes (one coefficient per cell type, ridge penalty chosen by
  10-fold CV). R^2 = the share of the signature explained by which cell types genes belong to ("composition share").
- Report the top cell types driving it (largest coefficients) and check that their sign agrees with the score
  differences in part 2.
- Repeat on the purity-unadjusted model (Step 4 without the ESTIMATE term) to show what the ESTIMATE adjustment
  already removes.

## Positive controls (fixed before running)

Expected differences, from known biology and our earlier findings; reported whatever they show:

| Cohort | Expected |
|---|---|
| BRCA (all subtypes) | Adipocytes much lower in tumors than GTEx breast (GTEx breast is fat-rich) |
| UCEC-EEC / SER, UCS | Smooth muscle much higher in GTEx uterus (myometrium) than tumors |
| TCGA-LUAD | Alveolar type 1 cells lower in tumors than lung |
| TCGA-COAD / READ | Epithelial (enterocyte) scores similar in tumors and the filtered GTEx colon; smooth muscle lower than unfiltered GTEx would show |
| TCGA-KIRC / KIRP | Proximal tubule markers lower in tumors than kidney cortex (ccRCC dedifferentiates) |
| TGCT-SEM | T cells / B lineage higher in tumors (lymphocyte-rich by definition) |
| TCGA-SKCM | Melanocytes higher, keratinocytes lower in tumors (mostly metastases) |

## Outputs

`scripts/35_tissue_composition.R` -> `results/step4/composition/`:
- `composition_panels.tsv` (cell type, markers, HPA release)
- `composition_scores.tsv.gz` (sample x cell type; group, role)
- `composition_differences.tsv` (group x cell type x comparison: medians, Cliff's delta, p, FDR)
- `composition_share.tsv` (group: composition share of the signature, top cell types, with / without ESTIMATE term)
- `composition_controls.tsv` (the positive controls)

## Assumptions

- Marker-panel scores track the abundance of a cell type, and a cell type's markers are not strongly regulated
  inside tumor cells (violated for e.g. hepatocyte genes in LIHC, proximal-tubule genes in KIRC: the tumor cells
  themselves derive from these cells, so lower scores mix "fewer cells" and "dedifferentiated cells").
- HPA single-cell reference cell types cover the cell types in our tissues; rare or tumor-specific states are missing.
- MCP-counter and marker scores give relative abundance across samples, not proportions.

## Possible extension (not in this spec)

Reference-based proportions (e.g. BayesPrism with Tabula Sapiens references) for 3-4 organs, to calibrate the marker
scores against estimated cell fractions.

## Runtime

Minutes per group; under 1 h in total once the v2 gene list and GTEx reference exist. New installs: `MCPcounter`;
new download: HPA single-cell type data.
