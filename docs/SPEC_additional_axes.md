# Additional scoring axes: exposure, normal-tissue selectivity, bypass state, and R's mechanism class

**Status:** design, 2026-11-15. Nothing implemented or run. Written after the review of what T (target
relevance), M (mechanism concordance), R (target-local rescue) and V (functional benefit) do **not** cover
(see `notes/d128_plan_review.md` and the D128 plan). Decision rows (D129-D132) are added to `DECISIONS.md`
once the D128-A0 audit settles, so this file does not collide with that run.

Each axis below is **reported separately** and never folded into the lead score without its own approval, as
with M and R. Where data are missing the output is an abstention with a reason, never a pass and never a zero.

---

## E · Achievable exposure (D129) — a filter

**Why it is not in T/M/R/V.** All four can be satisfied by a drug that never reaches the target at an active
free concentration. Exposure is currently split across eligibility rules (D107-D114) and a reported-only
diagnostic (D123); it is not a stated axis.

**What exists on disk**

| Input | File | Coverage |
|:--|:--|:--|
| Route / formulation | `results/lincs/drug_annotations.tsv` (`route_topical_only`, `atc`, `oral`, `parenteral`) | all FDA rows |
| Free fraction `fu` | `results/lincs/d123/occupancy.tsv` | 165 of 232 rows |
| Plasma `cmax_uM` | same | 101 of 232 rows |
| Free occupancy, `low_free_occupancy` | same | 24 True / 207 False |
| Structures (SMILES) for descriptors | `data/raw/chembl/ChEMBL_37/molecules.json` | 1,851 of 1,851 |

**Proposed filter, three tiers.** The gate is *conservative by construction*: it may exclude only where
evidence exists, and abstains otherwise.

1. **E1 route gate (hard, already in production).** Local-only formulations are excluded from systemic
   rankings. No change beyond recording it as part of this axis.
2. **E2 free-occupancy gate (new, narrow).** Where both `fu` and `cmax_uM` exist, compute free occupancy
   against the unit's potency as D123 does. `low_free_occupancy` becomes a **demotion flag with a stated
   reason**, not a silent exclusion. **Coverage is the binding constraint: only 101 drugs have `cmax`,** so
   this can never be a universal gate. Drugs without exposure data are `exposure_unknown`, ranked as now and
   flagged.
3. **E3 compartment gate (new, cohort-conditional).** For cohorts whose tumors sit behind a barrier — the
   LGG subtypes and TCGA-GBM — a systemically dosed drug that does not cross is not a candidate however good
   its T/M/R. Proposed implementation: CNS-MPO (Wager et al.) computed from ChEMBL SMILES, which needs
   **RDKit installed into the `repurposing` env** (not currently present). CNS-MPO is a crude, physicochemical
   prior, so it is a **flag and a reported score, never a hard exclusion**, and it applies only to the brain
   cohorts. An alternative, hand-curating CNS penetration for the FDA list, is more accurate and does not
   scale; decide in D129.

**Output:** `exposure_status` ∈ {ok, low_free_occupancy, route_excluded, barrier_risk, unknown} plus the
numeric `occ_free` and, where computed, `cns_mpo`. One row per (drug, cohort) because E3 is cohort-conditional.

**What this axis cannot do.** Tumor penetration is not plasma exposure; we have no tissue PK. E2 uses plasma
free concentration as a stand-in and should be described that way.

---

## S · Selectivity against normal tissue (D130) — a new pipeline

**Why it is not in T/M/R/V.** Nothing in the four asks what the drug does to cells that are not the tumor. A
target can be depended on, engaged and rescuing, and still be required in normal tissue. This is the one gap
with **no representation anywhere in the current pipeline**, and toxicity is a leading reason repurposing
candidates fail.

**The honest constraint.** DepMap screens cancer lines only. There is no normal-cell dependency measurement,
so "needed in normal cells" must be approached by proxy and labelled as such.

**Three components, each reported separately.**

- **S1 · Essentiality burden.** Is the target a DepMap common essential
  (`data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv`), and what is its Chronos effect distribution
  across all screened lines rather than only the cohort's? A target that every line needs is a toxicity
  signal. Reuses `79n_dependency.py`'s matrix; **must be computed on the full line panel, not the cohort
  subset**, or it duplicates T's dependency component.
- **S2 · Normal-tissue expression breadth.** From GTEx (`data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz`,
  already used in Step 4): the number of normal tissues expressing the target above a fixed threshold, its
  maximum normal TPM, and the tumor-to-normal ratio for the cohort in question. Report breadth and maximum
  separately — a target confined to one dispensable tissue differs from one expressed everywhere.
  **Do not weight tissues by "criticality";** that is a clinical judgement, not a measurement, and it would
  bury an assumption inside a number.
- **S3 · Window proxy.** The cohort's tumor dependency (T's D component) minus the pan-line essentiality of
  S1, reported as a difference with both parts visible. Positive means the cancer needs the target more than
  cells in general do.

**Validation before use.** Check S against drugs whose toxicity is established: a broad cytotoxic such as a
topoisomerase or tubulin agent should score poorly on S, and a selective agent against a lineage-restricted
target (ESR1 in breast, AR in prostate) should score well. If S cannot separate those, it is not measuring
what it claims.

**Circularity to avoid.** T already contains cohort dependency. S must use the **pan-line** distribution and
**normal-tissue** expression, so the two are not the same number twice. Report `corr(T, S)` across rows; a
high correlation means S is adding nothing.

---

## B · Bypass and resistance state (D131) — can it be learned?

**Yes, with one design that keeps the discovery and the application on different data.** The failure mode to
avoid: dependency already feeds T, so learning bypass from dependency and then adjusting T with it uses the
same evidence twice.

> **Learn the rule in cell lines (DepMap). Apply it in patients (TCGA).**

**B1 · Discovery, in lines.** For a target gene `g`, take the lines that express `g` and regress Chronos
dependency on `g` against the genomic state of candidate modifier genes (mutation from
`OmicsSomaticMutationsMatrix*`, copy number from `OmicsCNGene.csv`, expression from the TPM matrix). A
modifier whose altered state predicts **low dependency despite target expression** is a candidate bypass.

- Restrict candidates to the target's signed pathway neighbourhood in the cached OmniPath graph (already
  built for D126/D128) rather than testing every gene. This controls multiple testing and keeps every hit
  mechanistically interpretable.
- Seeded, with permutation nulls, and effect sizes reported with confidence intervals; `n` is about 1,100
  lines overall but far fewer per lineage, so most tests will be underpowered. Report power, not just p.

**B2 · Known-answer validation, fixed before running.** The method must recover **RB1 loss weakening CDK4/6
dependency**. Secondary checks: NF1 or RAS state and RAF dependency; PTEN state and PI3K-pathway dependency.
If B1 cannot recover the RB1/CDK4/6 relationship, which is well established and present in these data, it is
not ready to propose anything new. This is the same discipline as the anti-estrogen gate for D128.

**B3 · Application, in patients.** For a cohort and a drug whose target has a validated bypass modifier,
compute the **carrier fraction of the bypass state** in that cohort's tumors from MC3 and GISTIC
(`data/raw/tcga_pancan/`). Output `bypass_fraction` and a flag; high carrier fraction demotes the drug in
that cohort and states why. This uses patient genomics, not the cell-line dependency that trained the rule.

**B4 · Graph redundancy (secondary, cheap).** From the signed OmniPath graph, count alternative paths from
the target's downstream regulator set that do not pass through the target. Mechanistically motivated but
speculative; reported as a diagnostic only, never a gate.

**Not used for discovery:** PRISM and GDSC drug response. Those screens were already read by `79j`, `79k` and
`79m`, so using them to find bypass rules and then to validate them would be circular. They may serve as a
**held-out check** on B1 hits, subject to the prior-use audit in D128-A0.

**Expected yield, stated in advance.** Few targets will have enough lines and enough altered modifiers for a
powered test. The realistic output is a short, well-evidenced list (CDK4/6-RB1 and perhaps two or three
others), not a general resistance model. That is still worth having, because those few cases are exactly the
ones where a confident ranking is most misleading.

---

## R's mechanism class (D132) — an amendment to D128

**The problem.** R assumes "moving the tumor toward normal is good." That is false for whole mechanism
classes, so a low R is sometimes a category error rather than evidence against a drug.

**Proposed label,** attached to every R value and derived from inputs already present:

| `r_class` | Meaning | Derived from | How R is read |
|:--|:--|:--|:--|
| `interpretable` | Targeted agent against an over-active driver; normalisation is the expected effect | ChEMBL action type (inhibitor/antagonist/degrader) on a non-essential, non-DNA-acting unit | R is meaningful and may be gated |
| `not_expected` | Benefit does not take the form of a return to normal expression: synthetic lethality (PARP/HRD), cytotoxic and DNA-acting agents, differentiation agents | D124 `dna_acting` and `essential_unit` flags, plus a short curated synthetic-lethality list | R reported, **never** read as evidence against the drug |
| `ambiguous` | Mixed, unknown or modulator action | 79t `orient()` flags `direction_assumed` / `mixed` | R reported with the flag; no gate |
| `out_of_scope` | Mechanism not visible in a cancer-cell transcriptional assay: immune, angiogenic, microenvironment, antibody | Target compartment (D104 exception set) and modality | R and M both abstain with this reason |

**Consequences.** (1) Any pooled statistic over R is computed **within `r_class`**, never across all drugs,
because the classes are not comparable. (2) The D128 promotion criteria apply to `interpretable` rows only.
(3) `out_of_scope` makes explicit something the pipeline currently handles implicitly — that a cancer-cell
assay cannot see an immune mechanism, so silence there is not negative evidence.

This is the only one of the four that belongs **inside** D128, since it changes how an existing output is
read rather than adding a measurement. It is an amendment to Section 0, added as a D128 amendment row once
A0 settles, in the manner of D125-A0.1.

---

## Sequencing

1. **D132 (R's mechanism class)** first and with D128: small, uses existing flags, and prevents
   misinterpretation of the first R results.
2. **D129 (exposure)** next: mostly promoting and recording existing D123 and route evidence. Decide the E3
   question (CNS-MPO versus curation) before any implementation.
3. **D130 (selectivity)** next: new but straightforward, and it fills the only axis with no current
   representation.
4. **D131 (bypass)** last: the most research-like, gated on recovering the RB1/CDK4/6 relationship first.

Do not let any of these widen D128 while its controls are unresolved. D128's value is that it is narrow and
falsifiable; these axes are reported beside it, and each carries its own decision entry and its own gate.
