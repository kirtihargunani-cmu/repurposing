# Data Requirements: Steps 5–10

> **Status 2026-10-31: original plan** received 2026-10-01, kept for reference. How the pipeline actually runs: [`PIPELINE_V2.md`](../PIPELINE_V2.md).

Mechanism-grounded Subset-Level Drug Repurposing — Phase B (regulons, VIPER, master regulators) and Phase C (stratification). Received 2026-10-01; condensed copy. Project-specific notes marked **[NOTE]**.

## Step 5 — Regulons (only step in 5–10 with substantial new downloads)
- **CollecTRI** (open): `decoupleR::get_collectri()` / `OmnipathR::collectri()`; columns source_genesymbol, target_genesymbol, mor (±1), optionally n_resources/references. ~1,200 TFs; count interactions after download.
- **DoRothEA** (academic/non-commercial — record in license table): `decoupleR::get_dorothea()`, confidence A–C.
- **ARACNe, run ourselves** (open): ARACNe-AP or ARACNe3; expression matrix (genes × samples, log) + TF list; ~100 bootstraps, MI p 1e-8, DPI; `viper::aracne2regulon()`.
- **ARACNe precomputed** (open): `aracne.networks` (TCGA cohorts, Entrez IDs). Check whether laml is included; built on ~150 TCGA-LAML tumors (thin-sample problem).
- **Human TF list**, Lambert et al. 2018 (humantfs.ccbr.utoronto.ca), "Is TF?" + DBD family; TF universe and ARACNe input.
- ID-space friction: CollecTRI/DoRothEA = symbols, aracne.networks = Entrez, Step 4 = Ensembl. Build one mapping table; record unmapped targets per regulon (regulon-coverage statistic).
- AML: a real ARACNe arm needs a combined TCGA-LAML + BeatAML compendium, harmonized at Step 2.

## Step 6 — VIPER
Software only: viper (aREA, aracne2regulon) or decoupleR (run_ulm / run_viper / run_aucell). Filters: expressed-TF floor (from own matrix); minimum regulon size after the Step 5 ID join.

## Step 7 — Master regulators
`viper::msviper()` + `ledge()` + `shadow()` on the Step 4 signature. Optional genomic support: TCGA mutations (GDC MAF / Xena), TCGA GISTIC2 gene-level CNV, BeatAML mutations (doc says full calls are dbGaP phs001657; Tyner 2018 supplement has a curated table).
**[NOTE]** BeatAML 2.0 GitHub (pinned commit) also ships `beataml_wes_wv1to4_mutations_dbgap.txt` and raw `beataml_waves1to4_counts_dbgap.txt` openly.

## Step 8 — Consensus clustering
ConsensusClusterPlus; M3C (Monte-Carlo null for k) or diceR; cluster (silhouette). PAC computed from the consensus CDF or taken from M3C.

## Step 9 — Cluster validation
(a) Independent cohort processed the same way (AML: the other of TCGA-LAML/BeatAML; AMLCG GSE37642, Leucegene; solid: ICGC/cBioPortal). (b) Confounders: Step 3 table. (c) Known schemes: AML ELN 2022 / WHO 2022 / ICC 2022 (derived per patient), TCGA-LAML 2013 NEJM subgroups; breast PAM50 (genefu); colorectal CMS (CMSclassifier); pan-cancer PanCanAtlas subtypes. (d) silhouette; ARI (clue/mclust); chi-square / logistic regression for confounders; SigClust or M3C null.

## Step 10 — Cluster signatures
limma cluster-vs-rest on the NES matrix; msviper cluster-vs-normal re-aggregating Step 4 signatures; the same (byte-identical) Step 5 regulon object, reused in Step 16.

## Net new downloads
CollecTRI (+ DoRothEA) · Lambert TF list · GENCODE/Entrez/symbol mapping table · TCGA MAF + GISTIC for the disease · one independent cohort · subtype-label source. Build the ID mapping first; record unmapped fraction per regulon; drop TFs with too few mapped targets before Step 8.
