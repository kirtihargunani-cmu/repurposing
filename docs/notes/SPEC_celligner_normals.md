# Spec: where healthy tissue lands among cell lines (Celligner), per disease cohort

Status: **spec written 2026-10-18; first part run provisionally 2026-10-19** (`scripts/72b_celligner_gtex_full.py`,
commit e091894): all v2 GTEx reference normals projected with the tumors in chunks of 100 ->
`results/celligner/normal_projection_{gtex_v2_samples,group,lineage,drift}.tsv`. Anchor drift check: same nearest line
in 82-98% of samples (median 92%); 17 of 33 chunks pass the pre-set rule, so acceptance (weight level vs chunks of 50)
is still open and the analysis questions below are not answered yet (status 2026-10-31). Companion to the disease-side composition spec
(`docs/SPEC_tissue_composition.md` on `main`). Builds on `scripts/72_celligner_normal_control.py` (normal-tissue
negative controls, 2026-10-08).

## Question

For each disease cohort: which cell lines (and lineages) do the cohort's **tumors** resemble, which do its **healthy
GTEx tissue** (and TCGA adjacent normal tissue) resemble, and how different are the two? In Celligner terms, how much
of a tumor's line match is "the organ" and how much is "the cancer"?

## What exists already (script 72, GTEx arm)

- At most 50 GTEx donors per group (seeded sample), v1 GTEx tissue choices, no colon / kidney filter.
- Per group: nearest-line distance tumors vs normals (median AUC 0.88), share of normals whose nearest line is among
  the tumors' 5 most frequent (LIHC 100%, THCA 74%, PRAD 72%), tissue match of the nearest line.
- Not measured: full line-weight profiles of normals, lineage composition, overlap with the disease line weights
  used for drug matching.

## Inputs

- **GTEx samples: the v2 reference set** (v2 rebuild spec Step 4a: one sample per donor per tissue, the colon and
  kidney sub-region filter, only donors with collection data), read from `main` with
  `git show main:results/step4_v2/gtex_reference.tsv` once that file is committed. All donors, no cap
  (about 3,400 unique samples; a sample used by several groups, e.g. breast for the four BRCA subtypes, is projected
  once).
- TCGA adjacent normals: reuse the existing adjacent arm (`normal_control_adjacent_*`; all 737).
- Celligner model and DepMap reference: unchanged (`data/processed/celligner/model.pkl`, DepMap 24Q4, 1,588 cancer
  lines as candidates).

## Method

### 1. Projection in chunks

Celligner centers each batch on its own mean, so normals must stay a small share of the batch they are projected
with (script 72 design: all main-run tumors + normals as one batch, `transform(compute_cPCs=False)`).
- Split the GTEx samples into chunks of at most 700 (about 6% of each batch, as in the v1 GTEx arm at ~9%), balanced
  by tissue; each chunk is projected together with all 10,211 tumors.
- **Drift check:** a fixed set of 50 GTEx samples (seeded, mixed tissues) is added to every chunk. Their nearest-line
  distance and identity across chunks measure chunk-to-chunk drift (accept if the same nearest line in >= 90% and
  distance rho >= 0.98; otherwise reduce chunk size).
- Tumor reproducibility vs the main run, per chunk (as script 72: d1 rho, same nearest line).

### 2. Line weights for normals

Same as for tumors (`scripts/71`, `74`): softmax weights over negative distance with tau 7.8; hybrid mode (organ-key
lines plus any closer line) as in the disease line-weight matrix. Per group: the normals' line-weight vector = mean
of its GTEx donors' weights (and of its adjacent normals' weights).

### 3. Per-cohort comparisons

For each group with >= 8 GTEx donors (and, separately, >= 8 adjacent normals):
- **Distance:** nearest-line distance tumors vs normals (AUC, medians), as script 72.
- **Shared weight:** sum over lines of min(tumor weight, normal weight); 1 = tumors and healthy tissue resemble the
  same lines in the same proportions, 0 = no overlap.
- **Lineage composition:** weights summed by DepMap `OncotreeLineage` for tumors and normals; difference per lineage
  and total variation distance. This is the "composition" view in cell-line terms: e.g. uterus tumors on
  uterine lines vs GTEx uterus on soft-tissue / sarcoma lines.
- **Organ-driven lines:** lines carrying >= 5% of the tumors' weight whose normal-tissue weight is >= 50% of the
  tumor weight (the line matches the organ, not the cancer; LIHC / Hep G2 is the known case).
- **Cancer-specific lines:** lines with tumor weight >= 5% and normal weight < 1%.
- Link to composition (disease-side spec): per group, correlate the lineage shifts with the cell-type score
  differences (e.g. a high smooth-muscle score in GTEx uterus should go with weight on soft-tissue lines).

## Outputs

`scripts/72b_celligner_gtex_full.py` -> `results/celligner/`:
- `normal_projection_gtex_v2_samples.tsv` (sample: tissue, group(s), nearest lines, distances, chunk)
- `normal_projection_gtex_v2_weights.tsv.gz` (sample x line weights; in `data/processed/celligner/`)
- `normal_projection_group.tsv` (group: distances, AUC, shared weight, total variation distance, organ-driven and
  cancer-specific lines; GTEx and adjacent arms)
- `normal_projection_lineage.tsv` (group x lineage: tumor weight, GTEx weight, adjacent weight)
- `normal_projection_drift.tsv` (anchor samples across chunks; tumor reproducibility per chunk)

## Use for drug matching (proposed, not decided)

- Flag organ-driven lines in the disease line-weight matrix: drug evidence from those lines may reflect organ
  identity rather than cancer (already accepted in principle, 2026-10-08 decision: tissue-level matches are fine;
  the flag makes the share visible per drug).
- The healthy-tissue line weights are the natural input for the safety filter (see open questions in
  `drug_side_decisions_and_caveats.md`): lines resembling the patient's normal tissue stand in for normal-tissue
  drug response.

## Assumptions

- Projecting normals with all tumors in one batch keeps the tumor-defined centering (script 72 design).
- Chunking does not change normal positions beyond the drift check limits.
- Weights with the tumor-derived tau are meaningful for normal tissue, which sits farther from all lines.

## Runtime

About 5 chunks x ~12 min (the main fit took 12 min) plus metrics: about 1-1.5 h. Depends on the v2 GTEx reference
file from `main` (or on rebuilding its selection rules here if needed earlier).
