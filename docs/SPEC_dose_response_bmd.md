# SPEC — Intrinsic mechanism evidence from dose-response

**Status:** design only. Nothing implemented, no decision row written.
**Proposed ID:** D136. **Drafted:** 2026-11-17.
**Read first:** `docs/notes/settled_measurement_facts.md`.

## 1. Why this is different from everything else we have run

Every measurement in this project so far is **comparative**: `n_eff` is variance across drugs, `rev_pct` is a rank within a cohort, D134 retrieval is drug-versus-drug, the nulls compare a unit against other units. None asks whether **one drug, on its own terms**, shows coherent mechanism.

Dose-response is the classical intrinsic criterion. A drug either shows graded, reproducible engagement of its target across its own dose series, or it does not - and the answer does not depend on what else is in the panel, which cohort it is scored against, or which cell lines were borrowed.

It is also immune to the artefacts that have been biting us. No cohort composition, no borrowed lines, no degree-matched null, no comparison population. And it sidesteps the amplitude axis found at D131/D133 (|PC1| correlates 0.827 with a drug's overall response magnitude), because **potency is not magnitude**: where a curve crosses a threshold is a different quantity from how high it goes.

## 2. Grounding: this is established method, applied to a new unit

Benchmark-dose modelling is the standard for transcriptomic dose-response (BMDExpress, toxicology, since 2007): fit a curve per gene, derive a gene-level BMD, then **summarise gene-level BMDs into gene-set-level BMDs**. Set-level BMDs are validated - they recapitulate points of departure from traditional apical endpoints.

What is new here is the **unit of aggregation**: regulons and signed target footprints rather than pathways, and a **coherence** term that the standard pathway summary discards. The literature gap is explicit - Novartis's mvAC50 paper states there are no generally accepted methods to estimate EC50/IC50 from multi-parametric readouts.

**Why BMD rather than EC50.** L1000 dose ranges often do not reach a plateau, so EC50 would be extrapolated from unsaturated curves and would drift with how far the dose series happens to extend. BMD needs no plateau. More importantly, it lets us set the threshold from our own measurement error (Section 3).

## 3. The BMR is defined from replicate noise (user, 2026-11-17)

Rather than an arbitrary 10% or 1 SD, the benchmark response is **per gene**, from the replicate structure D129 already measured:

> **BMR_g = k x sd_within_g**

where `sd_within_g` is that gene's replicate standard deviation across signatures of the same (pert_id, cell line, dose, time) - the one-way ANOVA within-group term in `results/lincs/gene_quality/bing_gene_quality.tsv`, present for all 10,174 genes (landmark median 1.091, best-inferred 0.693).

A gene-level BMD is then interpretable in a unit that means something: **the dose at which this gene's response becomes distinguishable from the drug's own batch-to-batch variation.**

**This also solves the compression problem for free.** Inferred genes have lower replicate noise (0.693 vs 1.091) because their whole scale is compressed (D129: sd_all 0.918 vs 1.569). A per-gene noise-scaled threshold rescales each gene by its own variability, so no separate correction for feature space is needed and none should be added.

`k` is **not chosen here** - S1 fixes it, having costed candidates (1.0, 1.5, 2.0) at S0.

## 4. The nested structure

Three levels, each retaining the level below rather than collapsing it.

**Level 0 - gene.** Per (drug, line, gene) across the dose series: fit a dose-response curve, record direction, and derive the BMD as the lowest dose at which the fitted curve crosses BMR_g.

**Level 1 - regulon.** Aggregate member genes using the regulon's weights and expected signs. Report two things, not one:
- a weighted regulon-level BMD (potency)
- **coherence**: the dispersion of member-gene BMDs

Coherence is the point of the nesting. A regulon whose genes share a common BMD is evidence of **one regulatory event** driving them. A regulon whose genes have scattered BMDs is noise being averaged into a number that looks respectable. Collapsing to a single activity per dose - what the pipeline does now - destroys exactly this distinction.

Genes moving **against** their expected sign contribute to a discordance measure, never to the BMD.

**Level 2 - target footprint.** Aggregate regulators along the signed path, as D128 does, carrying both potency and coherence upward.

## 5. Right-censoring - the subtlety that would otherwise bias everything

Most genes will **never cross BMR at any tested dose**. Dropping them and averaging the rest biases the regulon BMD downward, because only the responsive tail contributes, and it makes a weakly-responding regulon look potent.

Non-crossing genes are **right-censored at the highest tested dose**, and the regulon aggregate must be censoring-aware (a Kaplan-Meier style median, or an explicit censored-regression fit). **The censored fraction is reported beside every BMD, always.** A regulon BMD computed on 8% of its genes is not comparable to one computed on 70%, and no table may present them side by side without that column.

## 6. The intrinsic null

Random gene sets drawn from **the same drug's own profiles**, matched on set size **and** on the `sd_within` distribution of the real regulon - otherwise the null is easier simply by containing quieter genes.

**The discriminating prediction, fixed in advance:** on-target engagement should appear at **lower dose** than the generic stress-and-arrest programme, which arrives near the top of the dose range. So the target's regulon should show a **lower and tighter** BMD than matched random sets from the same drug.

That is an intrinsic specificity test: one drug, its own doses, its own genes, no other compound involved.

## 7. Feasibility, measured 2026-11-17

Across the 307 annotated compounds: 104,292 signatures, 13,993 (compound, line) cells.

| | cells | compounds |
|---|---|---|
| >= 3 distinct doses | **8,778** | **284** |
| >= 5 distinct doses | 4,172 | 265 |
| >= 2 timepoints | 2,736 | 272 |

Timepoints present: 2, 3, 4, 6, 12, 24, 48, 72, 96, 120 h. **This is far better coverage than any comparative test** - 284 compounds against D134's 71 testable targets or D135's 46 clean cells.

**Note:** `79c` keeps one consensus profile per (drug, line) at TAS >= 0.4, so the production pipeline has discarded this entire dose structure. D136 reads level 5 directly.

## 8. Timepoint, and a hypothesis worth testing early

The anti-estrogen diagnosis was that estrogen-response genes fall by 6h and **rise again by 24h** as cells shift to a differentiated state, and that a rule preferring early timepoints fell back to 24h for most anti-estrogens anyway.

So there is a live possibility that **the pipeline has been reading at a timepoint where the direct mechanism is already overwritten by secondary consequences.** If regulon coherence and BMD are systematically better at 6h than 24h, that reframes the D128 ESR1 failure as a timepoint choice rather than a method failure. **Test this early** - it changes how everything downstream should be built.

## 9. Curve fitting

Hook effects are real in L1000: DOSE-L1000 fits GAMs rather than sigmoids above four concentrations precisely because sigmoidal fits silently misfit non-monotone shapes. **Non-monotone genes are flagged and reported, never forced through a monotone curve.** S1 fixes the model (monotone fit with a non-monotone flag, or GAM where doses allow) and the minimum dose count.

## 10. What this cannot establish

- **Not benefit.** Coherent dose-response is evidence of target engagement, not of therapeutic effect.
- **Not causality.** A coherent response proves the regulon moves with dose, not that the annotated target mediates it; a strongly correlated off-target would look the same.
- **Not the target annotation.** D121 units are ChEMBL annotations above a binding threshold (see the settled-facts note).
- **Not generality** beyond these lines, doses and timepoints.
- **Direction only where annotated.** Units flagged `direction_assumed` (every SERM) have no predicted sign, so the signed test abstains for them.

## 11. Stages

Standing routine: a `D136-<stage>` row committed and pushed **before** each stage runs; dated result row, `PROGRESS.md` and plain-language update after; push to main, merge into `celligner`.

- **S0 - feasibility and noise floor.** Dose-series coverage per (drug, line); `sd_within` distributions; what candidate `k` values cost in crossing rates and censoring; dose-count and timepoint availability. **No curve fitted, no BMD computed.**
- **S1 - freeze.** `k`, curve model, minimum dose count, censoring-aware aggregate, coherence statistic, null matching. Tagged.
- **S2 - gene-level curves and BMDs**, with censored fractions.
- **S3 - regulon aggregation**, potency and coherence.
- **S4 - target footprint and the intrinsic null.**
- **S5 - the timepoint question** (Section 8), if not already answered at S0.
- **S6 - disposition:** what it implies for D128, for D134/D135, and for target-anchored scoring.

## 12. Outputs

`results/lincs/d136/`: `feasibility.tsv`, `gene_bmd.parquet`, `regulon_bmd.tsv`, `target_bmd.tsv`, `null_comparison.tsv`, `timepoint.tsv`, `d136_run.json`.

## 13. Constraints

- **Audit only.** No production score changed, nothing re-ranked.
- **No outcome labels.** No CTRP, PRISM or GDSC. This test needs none.
- Thresholds fixed before results are seen; changes afterward are dated amendments.
- Threads < 32, set explicitly. Level 5 `trt_cp` is 35.5 GB - memory-map and read in blocks, never whole.
- Do **not** add a feature-space correction; Section 3 handles it.

## 14. Sources

- BMDExpress / transcriptomic benchmark dose: gene-level fit then gene-set-level summary; set-level BMDs recapitulate apical points of departure.
- DOSE-L1000 (*Bioinformatics* 2023, doi:10.1093/bioinformatics/btad683) - GAMs over sigmoids above four concentrations; hook effects in L1000 dose-response.
- Renner et al., mvAC50 (*Sci Rep* 2020) - no generally accepted method for EC50/IC50 from multi-parametric readouts.
