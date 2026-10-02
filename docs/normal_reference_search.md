# External normal references for hematologic groups (search, 2026-10-01)

Status: **parked**. The PI de-prioritized hematologic tumors (AML / ALL / DLBC) on 2026-10-01. This records what exists if they are revisited.

Method: `recount3::available_projects()` (human); per-project SRA metadata (`sra.sra.<project>.MD.gz`) for all 8,802 human SRA projects keyword-scanned (marrow, CD34, HSPC, tonsil, lymph node, germinal center, naive/memory B × healthy/normal/control), hits curated by hand. GEO series from NCBI. TARGET normal marrow via the GDC API. "Verified" = checked against recount3 metadata.

Caveat for all recount3 SRA projects: a recount3 sample is a run (SRR). Several studies split one biological sample across lanes (6 runs per sample in SRP149374; 8 in SRP098651), so sum by `sample_acc`/GSM before voom.

## A. TCGA / TARGET in recount3
| Item | Result (verified) |
|---|---|
| TCGA LAML | recount3 `LAML` (tcga), 178 samples (Toil: 173) |
| TCGA DLBC | recount3 `DLBC` (tcga), 48 samples (Toil: 47) |
| NCI TARGET | Not in recount3 ("TARGET_ALS" is the ALS consortium) |
| GTEx `BONE_MARROW` | All 204 = K562 CML cell line; no GTEx lymph node. SPLEEN (255) exists |

Pediatric normal marrow outside recount3: GDC TARGET-AML has 88 "Bone Marrow Normal" STAR-count files from 70 normal donors (`TARGET-00-BM####`, `TARGET-00-RO#####`; Bolouri 2018 NBM). There are also 324 under TARGET-20 patient IDs, probably remission marrow (unverified). No donor ages in GDC. GDC STAR / GENCODE v36, so pipeline-matched only with GDC-STAR TARGET tumors. TARGET-ALL has 0 normal marrow RNA-seq files.

## B. Healthy bone marrow in recount3 SRA (verified; n = biological samples)
| SRA | GEO | Healthy n | Cells | Library |
|---|---|---|---|---|
| SRP133442 | GSE111085 | 23 (+44 MDS) | CD34+ BM, adult | SMARTer Ultra Low, PE |
| SRP158628 | GSE118944 | 12 (6F/6M) | CD34+ BM, young adult | PrepX polyA, SE |
| SRP149374 | GSE114922 | 8 CD34+ (+ myeloid/erythroid fractions; 13 donors) | CD34+ BM | SMARTer, PE; 6 runs per sample |
| SRP050146 | GSE63569 | 5 | CD34+ BM | NEBNext polyA stranded, PE; same group as SRP149374 |
| SRP098651 | Leucegene | 3 (ages 33-51) | whole BM (only whole-marrow set) | TruSeq polyA, PE; 8 runs per sample |
| SRP062025 | GSE55689 | 8 donors / 22 samples | sorted HSC/GMP/MEP, older adult | SMARTer Ultra Low, SE |

Smaller or sorted: SRP065216 (Corces 2016 hematopoietic hierarchy), ERP006121 (E-MTAB-2672), SRP105369 (GSE98310; includes 4 normal early pre-B / pre-B samples, the only normal B-progenitors found), SRP058719, ERP003613 (HPA, 4 whole BM). Not suitable: SRP139763 (single-cell), SRP163252 BM (CD138+ plasma cells), SRP058414 (pre-B ALL patients).

Bottom line: there's no healthy whole-marrow cohort of 10 or more donors. CD34+ pooled across SRP133442 + SRP158628 + SRP149374 + SRP050146 gives about 48 adult samples from 4 library batches (study as a covariate). CD34+ vs bulk blasts mixes differentiation stage into the signature. No pediatric healthy marrow in recount3.

## C. DLBCL comparators in recount3 SRA (verified)
Tissue: ERP003613 (HPA: 5 lymph node, 4 spleen), SRP163252 (7 healthy tonsil, rRNA-depleted), ERP006650 / ERP023317 (3 + 3 tonsil), ERP000546 (BodyMap, 1 LN), GTEx spleen (255). Rejected: SRP193844 (trial patients), SRP072528 / SRP123359 / SRP198243.
Sorted B cells: SRP021509 (GSE45982; tonsil naive 4, GC 4), SRP133767 (GSE111310; tonsil GC 4, naive 3, memory 2), SRP132924 (GSE110669; tonsil naive / GC / memory 3 each), SRP077918 (naive 2, GC 2), SRP048820 (centrocytes 2), SRP156583 (blood B subsets, 37 healthy; inferred), SRP125125 (Monaco 2019; inferred). Totals about 13 GC-B, 12 naive and 5 memory tonsil B samples across small studies, so study is confounded with cell type.

## D. Uniform alternatives
BloodSpot / HemaExplorer: microarray or normalized values, not usable for voom. Haemopedia: TPM only (raw = SRP150419). ARCHS4: kallisto counts for most GEO series, voom-compatible but a third pipeline. TARGET NBM via GDC STAR: see A. Leucegene / BeatAML normals: not uniformly reprocessed with TCGA.

Working files from the search (per-project metadata, keyword scan) were left in the session scratchpad and are not part of the repository.
