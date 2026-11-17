# SPEC — Does a drug's profile carry its own target's signature?

**Status:** design only. Nothing implemented, no decision row written.
**Proposed ID:** D135. **Drafted:** 2026-11-17.
**Runs BEFORE** `SPEC_target_recoverability.md` (D134) stage S2.

## 1. The question, and why it comes first

Compare a drug's transcriptional profile against a **genetic perturbation of its own annotated target** in the same cell line. An inhibitor should resemble the knockout; an agonist should resemble overexpression or anti-correlate with the knockout.

This is the most fundamental version of the question and it precedes everything else we have built:

- **No network.** No regulon, no signed graph, no hops. Two measured profiles compared directly. D128's ESR1 failure is ambiguous between "the signal is not in the data" and "our network does not capture it"; this test cannot be confounded that way.
- **No second drug.** Each drug is tested alone, so the 111 singleton targets that D134 cannot touch become testable.
- **The direction is predicted in advance**, from D121's `orient` (lower/raise) and `actions` (INHIBITOR/AGONIST). A result with the wrong sign is interpretable, not merely null.

If a drug's profile does not carry its own target's signature, then no target-anchored score can work, and D134's drug-drug retrieval is testing something derivative.

## 2. Feasibility, measured 2026-11-17

Our drug matrix spans **73** cell lines (63 for the 307 annotated compounds). Restricted to lines shared with it:

| perturbation | targets covered (of 377) | signatures | shared lines |
|---|---|---|---|
| **CRISPR knockout (`trt_xpr`)** | **208** | 837 | **8** |
| shRNA consensus (`trt_sh.cgs`) | 174 | 515 | 9 |
| shRNA raw (`trt_sh`) | 207 | 6,109 | 11 |

**Corrected by S0 (2026-11-17).** The draft said 87 drug-matrix lines; that was the count *before* `79c`'s TAS >= 0.4 floor, and the real figure is 73. The draft's three rows were also computed under three different filters; the table above is now consistent (`is_hiq`, restricted to drug-matrix lines), which moves the cgs row from 209/1,627/10 to 174/515/9. The raw-shRNA row did not reproduce under any of sixteen readings S0 tried and is given as S0 measured it.

**Coverage is not scope.** Paired (drug, gene, line) cells - the real unit - are **2,350** for CRISPR over 255 genes, 201 drugs and **7** lines, with five lines carrying 389-428 cells each. See §6a.

## 3. Primary and secondary arms, and why

**CRISPR is primary.** The CMap L1000 paper reports that shRNAs targeting the same gene are only *slightly* more similar to each other than random, because shared-seed off-target effects outweigh on-target ones; a PLOS Biology analysis of the same data finds RNAi off-target effects far larger than usually appreciated while CRISPR off-target activity is negligible. Using raw shRNA as ground truth would measure seed effects.

- **Primary:** CRISPR knockout (`trt_xpr`).
- **Secondary:** shRNA consensus (`trt_sh.cgs`) - CMap's own consensus across shRNAs for a gene, which is the published mitigation for seed effects.
- **Reported, not relied on:** raw shRNA. **Pre-stated:** disagreement between CRISPR and raw shRNA is expected and is evidence about shRNA, not about the drug.

**S0 correction to this arm ordering.** The consensus-shRNA arm designated here as the seed-effect mitigation is the **weakest-responding arm by about 6x** (median `nz2` 14, against 51 for CRISPR and 88 for raw shRNA), and a strength floor of `nz2 >= 100` leaves it just 175 of 2,250 paired cells. It is reporting-only unless S1 deliberately lowers its floor. Separately, the NMD-escape caution drafted for CRISPR did not materialise in the feared direction: on measured target genes **CRISPR passes mRNA QC more often than shRNA** (81% vs 57%), so NMD escape is a tail effect, not the bulk.

## 4. The test and its null

**Statistic:** signed concordance between the drug profile and the target-perturbation profile in the same line, with the sign predicted in advance by the drug's annotated action.

**Null (the part that gives it meaning):** the same drug against perturbations of **other genes**, matched on perturbation signature strength. The question is not "does the drug correlate with this knockout" but "does the drug resemble **its own** target's knockout more than it resembles other knockouts." Gene-level resampling must preserve gene-gene correlation (D130/D130.1); do not reintroduce the independence error.

**Also test the target's signed neighbourhood.** The literature reports that drug correlations appear for knockdowns of the target *and/or* of genes up- and downstream. Report target-only concordance as primary, neighbourhood concordance separately, and never merge them.

## 5. Validating the test itself — required before any drug result is read

**Positive controls, fixed in advance:** drug-target pairs whose concordance should be unmistakable (HDAC inhibitors vs HDAC knockout; proteasome inhibitors vs PSMB knockout; topoisomerase poisons vs TOP2A).

**A flaw in that control set, found at S0 and not repaired here.** All three families have data (HDAC 104 CRISPR cells / 5 drugs; PSMB 112 / 4; TOP2A 39 / 8), but **TOP2A is 100% common-essential and PSMB 57-61%** - two of the three sit inside limitation 3's exculpatory stratum, where a knockout yields a death signature. A pass there may demonstrate shared cytotoxicity rather than shared mechanism, and would **not** license reading concordance for non-essential targets. **HDAC (25% essential) is the only mostly-clean family.** S1 must state in advance what a pass on an essential-gene control licenses, or replace those two controls.



**Pre-stated rule: if the positive-control set fails, the method is reported as non-functional and no drug-level conclusion is drawn.** A failing control means the comparison is broken, not that the drugs are non-concordant.

**Perturbation QC:** a knockout that did not knock down is not evidence. Where the target's own mRNA is measured, require it to move in the expected direction (the check D84 already used). Perturbations failing QC are excluded and counted, never scored as non-concordant.

## 6. Limitations, annotated — read this before reading any result

**This test is asymmetric. A positive result is strong evidence; a negative result is weak evidence.** Many independent things can break the comparison, and the literature expectation is that drug-to-genetic concordance is *difficult to achieve* even for true targets. So failure is uninformative unless the positive controls pass AND the stratifications below rule out the listed cause. **A null result must never be reported as "the drug does not engage its target."**

| # | Limitation | What it does | How we detect it | How to read a negative result if it is biting |
|---|---|---|---|---|
| 1 | **Kinetics mismatch** | Drug read at 6-24h; CRISPR over days. Different cellular states. | Concordance by drug timepoint. | If concordance rises with drug exposure time, kinetics explains the failure. **Excused.** |
| 2 | **Knockout is not inhibition** | Deleting a protein removes scaffolding functions that an active-site inhibitor leaves intact. | Stratify by targets with known non-catalytic roles. | A genuine inhibitor can legitimately fail. **Excused for those targets.** |
| 3 | **Essential targets** | Knocking out an essential gene yields a death/arrest signature, not a mechanism signature; the drug at sub-lethal dose does not. 1,842 rows already carry `essential_target`. | Stratify by DepMap common-essential. | Essential targets failing is **expected and excused**. |
| 4 | **shRNA seed effects** | Off-target activity can exceed on-target. | CRISPR vs raw-shRNA agreement. | shRNA-only failure is evidence about shRNA. **Excused.** |
| 5 | **Weak perturbation** | A knockout with no transcriptional consequence cannot be matched. | Perturbation signature strength; target-mRNA QC (§5), landmark targets only. | **Excluded, not scored.** |
| 5b | **mRNA QC is only available for 48 of 377 targets** (28 with CRISPR) | The target's own transcript is a *measured* landmark gene for 48 targets; for the rest it is imputed, and imputed target z sits on zero (median -0.39 CRISPR, -0.02 cgs) - no information about whether the perturbation worked. On measured genes the QC is excellent (CRISPR median target z **-3.90**, 81% pass). | Feature space of the target gene. | **QC is usable for 28 CRISPR genes and nobody else.** Applied across the board it would delete ~90% of inferred-gene perturbations for an imputation artefact. Use a strength floor or DepMap gene effect elsewhere. |
| 6 | **Weak drug response** | A drug that moved nothing cannot be matched. Note `79c` already keeps only TAS >= 0.4, removing 73% of cells. | Drug signature strength. | **Conditioned on; excused below the floor.** |
| 7 | **Polypharmacology** | The drug hits many proteins; the knockout hits one gene. | Number of annotated targets per drug. | Multi-target drugs are **excused**. |
| 8 | **Target not active in the line** | No concordance is possible if the target is not expressed or the pathway is not running. Literature: connectivity is expected where the target is mutated or highly expressed. | Target expression and mutation status per line. | **Excused.** |
| 9 | **Thin line overlap (8-12)** | Power is concentrated in a few lines. | Per-line result counts. | A **power** limit, not a biological one. |
| 10 | **Annotation may be wrong** | D121 targets are ChEMBL annotations above a binding threshold; 230 of 307 compounds carry exactly one, which reflects curation. | Not detectable here. | **Irreducible.** A "wrong" result may be right about biology. |

**The corollary that matters for how we act:** because so much of the attribution above is exculpatory, this test can establish that target signal **is** present but cannot establish that it is **absent**. Design every downstream decision to depend on the positive finding, not the negative one.

## 7. Pre-registered expectation

Record before running: **the literature expectation is weak concordance.** Drug-to-knockdown correlation is described as difficult to achieve in practice even for true targets. A broadly negative result is therefore *consistent with prior work* and is not news; a broadly positive result would be the notable outcome. State this in advance so a negative result is not written up as a discovery about this pipeline, and so a positive one is not discounted.

## 8. Stages

Standing routine: a `D135-<stage>` decision row committed and pushed **before** each stage runs; dated result row, `PROGRESS.md` and plain-language update after; push to main, merge into `celligner`.

- **S0 — feasibility and QC.** Perturbation coverage per target and line, perturbation strength, target-mRNA QC pass rates, essentiality and expression strata. **No concordance computed.**
- **S1 — freeze.** Statistic, null, strength floors, positive-control set, verdict thresholds. Tagged.
- **S2 — positive controls only.** If these fail, stop and report the method as non-functional.
- **S3 — the drug set**, CRISPR primary, with all strata from §6.
- **S4 — the per-target concordance map**, with every negative annotated by which limitation applies.
- **S5 — disposition:** what it implies for D128, for D134's S2, and for target-anchored scoring generally.

## 9. Outputs

`results/lincs/d135/`: `feasibility.tsv`, `perturbation_qc.tsv`, `positive_controls.tsv`, `concordance_by_pair.tsv`, `concordance_by_target.tsv` (the map, with limitation annotations), `d135_run.json`.

## 10. Constraints

- **Audit only.** No production score changed, nothing re-ranked.
- **No outcome labels.** No CTRP, PRISM or GDSC. This test needs none.
- Thresholds fixed before results are seen; later changes are dated amendments, never edits of frozen rows.
- Threads < 32, set explicitly. The `trt_sh` and `trt_xpr` matrices are 11.8 GB and 6.5 GB: memory-map, never load whole.

## 11. Sources for the claims in §3, §4 and §7

- Subramanian et al., *Cell* 2017, "A Next Generation Connectivity Map" — shared-seed similarity exceeding same-gene similarity; CRISPR recovering small-molecule connections missed by CGS.
- Evaluation of RNAi and CRISPR technologies by large-scale gene expression profiling in the Connectivity Map, *PLOS Biology* — comparable on-target efficacy, RNAi off-target effects far larger than appreciated, CRISPR off-target negligible.
- Predicting protein targets for drug-like compounds using transcriptomics, *PLOS Comp Biol* — drug-to-target-knockdown correlation difficult to achieve in practice; correlations arise for the target and for up/downstream genes.
