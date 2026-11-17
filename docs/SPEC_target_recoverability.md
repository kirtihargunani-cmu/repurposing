# SPEC — Target recoverability from transcriptional profiles

**Status:** design only. Nothing implemented. No decision row written yet.
**Proposed ID:** D134. **Drafted:** 2026-11-17.

## 1. The question

**Can a drug's annotated protein target be identified from its transcriptional response, and for which targets?**

This is foundational and should have preceded the scoring work. Retrieval accuracy is an **upper bound on any target-anchored score**: if a drug's target cannot be recovered from its profile, no mechanism model built on that profile can do better, however well specified.

It is also the generalisation of a test already run once. The D128 ESR1 control was exactly this question for a single target. It passed under a null that assumed independent genes and failed under the corrected null (D130). Running it across all testable targets is the powerful version, and it is label-free — it uses target annotation only, no efficacy data, no outcome labels.

## 2. Scope, and the hard limits

Measured from `results/lincs/target_role/d121/units.tsv` and `results/lincs/drug_annotations.tsv` on 2026-11-17:

| | count |
|---|---|
| FDA compounds with a D121 target unit | 307 |
| Distinct target units | 193 |
| Distinct gene targets | 387 |
| **Targets with >= 2 compounds (testable)** | **82**, covering 247 compounds |
| Targets with >= 3 compounds | 59, covering 224 compounds |
| Targets with >= 5 compounds (well powered) | 30, covering 156 compounds |
| **Singleton targets (untestable by retrieval)** | **111 (58% of units)** |

Compound coverage: median 4 active cell lines (IQR 1-8); 206 of the testable compounds are active in >= 3 lines, 163 in >= 5.

**State this limit prominently in any result.** Retrieval cannot say anything about 111 of 193 targets. The domain map covers at most 82, and only 30 are well powered. A target absent from the map is untested, not undetectable.

## 3. The core test

**Unit of analysis:** a (compound, cell line) consensus profile. Level 5 signatures are already MODZ-consolidated within (pert_id, line, dose, time); aggregate across dose and time within (compound, line) for the primary arm, and keep them separate as a sensitivity arm.

**Query and candidates:** for each (compound, line) profile, candidates are all *other compounds* profiled in the **same cell line**. Other `pert_id`s of the same `cmap_name` are **excluded** — retrieving a different batch of the same drug is self-retrieval, not target recovery.

**Primary metric: per-target AUROC.** For a target T, rank all candidate pairs by profile similarity and ask whether same-T pairs rank above different-T pairs. AUROC is used as the primary because it is calibrated against class imbalance by construction, which matters here: NR3C1 carries 25 compounds while most targets carry one or two, so raw top-1 accuracy would be dominated by target frequency.

**Secondary metrics:** top-1 and top-5 retrieval, mean reciprocal rank — reported with their frequency-aware chance baselines, never bare.

**Similarity:** Spearman correlation over the gene profile as primary; cosine as sensitivity. Both arms (`specific`, generic-response-removed; and `raw`) and both gene sets (all-BING; landmark-only) are run.

## 4. The chance baseline, and the null

Two separate things, both required.

1. **Analytic per-query chance:** for a query in line L with target T, chance = (compounds in L sharing T, excluding self) / (candidate compounds in L - 1). Report observed against this expectation per query, never pooled without it.
2. **Label-permutation null:** shuffle target assignments across compounds **preserving the target-size distribution**, recompute, repeat. This is the reference the verdicts are read against.

**Use the corrected nulls from D130.** The permutation operates on compound labels, not genes, so it does not inherit the gene-independence defect — but any gene-level resampling introduced anywhere in this spec must preserve gene-gene correlation (D130/D130.1). Do not reintroduce the error this project has just spent three stages correcting.

## 5. Confound controls — all four are required

A positive result is uninterpretable without these.

1. **Chemical analogy.** Drugs sharing a target are frequently structural analogs, so similarity may reflect chemistry rather than mechanism. Flag close analogs by the **first block of the InChIKey (14 characters = skeleton hash)**, which is already in `drug_annotations.tsv` and needs no new dependency. Report AUROC with analog pairs **excluded**. *Pre-stated reading: if analog-excluded AUROC drops substantially, retrieval is measuring chemistry, not mechanism.*
2. **Signature strength.** A compound that produced no transcriptional response cannot be retrieved, and must not be scored as "target undetectable". Condition on a strength measure (count of genes at |z| >= 2, or TAS where available) and report a strong-signature-only arm. *This separates "the assay cannot resolve this target" from "this drug did nothing at this dose".*
3. **Cell line.** Comparisons are within-line by construction. Additionally report whether recoverability varies by line, and whether it tracks line-level response amplitude.
4. **Generic response.** Both arms are run. *Pre-stated reading: if retrieval works only on `raw`, it is the shared stress-and-arrest program being recovered, not target identity.*

## 6. The main deliverable — the domain map

Not a global pass/fail. A **per-target table**: target unit, genes, compounds tested, lines, AUROC, null percentile, analog-excluded AUROC, strength-conditioned AUROC, and a verdict.

**Pre-registered verdict rule, to be fixed before any retrieval is computed:**

- **recoverable** — AUROC >= 0.65, above the 95th percentile of the label-permutation null, **and** >= 3 compounds
- **not recoverable** — AUROC not above the null, with >= 3 compounds and adequate signature strength
- **underpowered** — 2 compounds, or insufficient strength, or too few shared lines; reported, never counted as either verdict
- **untested** — the 111 singletons

This map is the point of the exercise. It defines the **domain where target-anchored scoring can work**, which is far more useful than a verdict on the method, and it would let any future method be scoped to targets where its premise holds rather than applied everywhere.

## 7. Link to the dimensionality findings

D130 measured the gene-level LINCS space at ~41 effective dimensions and D131/D133 put regulator activities at ~10. Against 193 distinct targets, that gap is the open question this spec addresses.

Compute the effective dimensionality of **between-target** variation against **within-target** variation. If targets are distinguishable, between-target variance should occupy materially more dimensions than within-target. Report both against the ~41 ceiling.

Note in advance: 193 targets does **not** imply 193 dimensions. Inhibitors of KIT, FLT3, ABL1 and PDGFRA converge on much the same growth-arrest programme, and that convergence is real biology rather than measurement failure. The test distinguishes convergence from non-detection; it does not assume either.

## 8. What this cannot establish

State these in any result; they are not hedges, they bound the claim.

- **Not causality.** Retrieval could succeed through off-target effects shared among analogs, which control 5.1 bounds but does not eliminate.
- **Not clinical relevance.** A detectable transcriptional fingerprint says nothing about therapeutic benefit.
- **Not generality.** Findings apply to these cell lines, these doses and ~24h readouts. The ESR1 case already showed 6h and 24h profiles pointing in opposite directions.
- **Not the 111 singletons**, which remain untested.
- **Annotation is not truth.** D121 target units are ChEMBL annotations above a binding threshold; 230 of 307 compounds carry exactly one, which reflects curation rather than pharmacology. Real polypharmacology is broader, so "wrong" retrievals may be right about biology and wrong about the annotation.

## 9. Stages

Follows the standing per-stage routine: a `D134-<stage>` decision row committed and pushed **before** the stage runs; a dated result row, `PROGRESS.md` entry and plain-language update after; push to main, merge into `celligner`.

- **S0 — feasibility audit.** Testable targets, shared-line counts, signature strength distributions, chance baselines. **No retrieval computed.** Confirms or revises the Section 2 numbers.
- **S1 — freeze.** Similarity metric, aggregation, verdict thresholds, null design. Tagged.
- **S2 — core retrieval** plus both baselines.
- **S3 — confound controls** (Section 5), all four.
- **S4 — the domain map** (Section 6).
- **S5 — dimensionality link** (Section 7).
- **S6 — disposition.** What the domain implies for D128 and for any future target-anchored scoring.

## 10. Outputs

`results/lincs/d134/`: `feasibility.tsv`, `retrieval_by_query.tsv`, `retrieval_by_target.tsv` (**the domain map**), `confound_controls.tsv`, `dimensionality.tsv`, `d134_run.json`.

## 11. Constraints

- **Audit only.** No production score changed, nothing re-ranked, no fix shipped.
- **No outcome labels.** No CTRP, PRISM or GDSC value read. This test needs none, which is what makes it cheap to run now.
- Thresholds fixed before results are seen; a change after is a dated amendment with its reason, never an edit of a frozen row.
- Threads < 32, set explicitly.
