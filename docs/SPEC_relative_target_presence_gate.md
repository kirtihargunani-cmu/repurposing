# SPEC — relative target-presence gate (D151, proposed)

**Status:** pre-registration. Nothing changed in production by this document.
**Date:** 2026-11-19 · **Raised by:** the user — *"if we gate by target prevalence between cell lines
and tumour tissue, anti-estrogens should not even be scored for basal when that target is not present."*

## 1. The defect

The gate (`config/drug_gate.yaml`, `scripts/79b_target_gate.py`) calls a gene present in a disease if
**TPM ≥ 1 in ≥ 20% of that cohort's tumours**. It is an absolute floor, so a gene can clear it while
being effectively absent relative to where it actually matters:

| cohort | ESR1 median TPM | % tumours ≥ 1 TPM | current gate | share of ESR1's own maximum |
|---|---|---|---|---|
| BRCA-Basal | **0.86** | 44.7% | **passes** | **0.6%** |
| BRCA-Her2 | 2.08 | 67.9% | passes | 1.5% |
| BRCA-LumA | 108.41 | 99.8% | passes | 78.1% |
| BRCA-LumB | 138.76 | 100% | passes | 100% |

A **126-fold** difference between LumA and Basal, and all four pass. D91 recorded this on 2026-11-03
("ESR1 counts as present in Basal, median 0.86 vs 108 TPM in LumA") and it was never acted on.

Consequence, measured in D150: anti-estrogens score at the **76th–96th percentile in BRCA-Basal** on
the target-regulon cosine — a false positive the gate should have removed before scoring.

## 2. The proposed rule

Add a **relative** criterion alongside the absolute one. A gene is present in a cohort only if it
clears the existing floor **and**

> `median_tpm(gene, cohort) ≥ rel_min × max over cohorts of median_tpm(gene, ·)`

i.e. the gene must reach a stated fraction of its own highest cohort-level expression. This is
scale-free per gene, so it needs no per-gene tuning, and it directly encodes the user's point: a
target is "present" where it is expressed *like it is where it matters*, not merely above a floor.

## 3. Threshold, chosen on coverage and frozen before scoring

Sweep over the 30,215 (gene, cohort) pairs that currently pass:

| `rel_min` | pairs blocked | share |
|---|---|---|
| 0.05 | 3,817 | 12.6% |
| **0.10** | **6,486** | **21.5%** |
| 0.20 | 10,667 | 35.3% |
| 0.33 | 15,841 | 52.4% |

**Proposed `rel_min` = 0.10.** It removes about a fifth of currently-passing pairs — a real tightening
without gutting coverage — and blocks ESR1 in Basal (0.6%) and Her2 (1.5%) while keeping LumA (78%)
and LumB (100%). Every value from 0.05 to 0.20 does that, so the ESR1 case does **not** discriminate
between them and the choice is made on coverage alone, not on making a control pass.

## 4. Pre-stated acceptance tests

1. **Anti-estrogens become ineligible in BRCA-Basal and BRCA-Her2**, and remain eligible in LumA/LumB.
2. **Known-use recovery does not fall**: median known-use AUC across cohorts (Hub and ChEMBL labels,
   C28) drops by no more than 0.01 on either label set.
3. **Coverage floor**: no cohort loses more than 35% of its gate-passing drugs; any that does is
   reported individually rather than absorbed into a median.
4. **The biomarker pairs of D127 do not degrade** — they are the project's primary check per
   `approval-is-weak-validation`.

Fail any of these and the rule is reported and not adopted.

## 5. What this does and does not fix

**Does:** removes a false positive — anti-estrogens scored in a disease that does not express their
target. That is a win on its own terms and it is what the user asked for.

**Does not:** create the corresponding true positive. D150 also found anti-estrogens at the
**0.9th–12.7th percentile in BRCA-LumA/LumB**, the disease they actually treat, and no prevalence
filter touches that. That failure is now traced to timing, not presence: in MCF7 — an ER-positive
line that matches luminal disease well — the ESR1 programme responds **−0.056 at 6 h (correct)** and
**+0.094 at 24 h (inverted)**, consistently across all five anti-estrogens, and the cosine against the
LumA disease vector is **−0.132 at 6 h (opposes, ranks high)** versus **+0.093 at 24 h (agrees, ranks
low)**. The pipeline reads 24 h. That is a separate pre-registration (D152).

## 6. Implementation

One rule in `config/drug_gate.yaml` (`expression.rel_min`) and one condition in
`scripts/79b_target_gate.py` alongside the existing `min_fraction` test, plus a new
`gate_reason = "target_low_relative_to_max"` so blocked rows state why. Re-runs `79b` → `79w` → the
gate-dependent downstream steps; no scoring formula changes.
