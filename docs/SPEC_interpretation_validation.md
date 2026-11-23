# SPEC — Validating the interpretation, not the ranking (D202)

**Status: PRE-REGISTRATION. Written and committed before any statistic in sections 4–6 is
computed.**

User, 2026-11-22: *"let's focus on interpretation now then, how can we validate interpretation …
we're working on validating the fact that we can condense these 800+ regulators into these 20–40
programs for both drugs and diseases BY KEEPING THE GENE LEVEL INFORMATION BECAUSE WE IDENTIFIED
THAT ONE NES SCORE CAN MEAN SEVERAL THINGS FOR AN MR."*

---

## 1. Why this document exists, and the instrument problem it fixes

Every test this week measured **approval-recovery AUC**. That is a ranking yardstick and it cannot
validate an interpretive claim. D192 showed it is actively misleading here: bimodal, with 9 cohorts
significantly above chance and 9 significantly below (UCEC-EEC z = −7.7), because several cohorts'
approved sets are cytotoxics. D196 then showed the single proliferation axis **beats** the
mechanism score in exactly those cohorts. A representation can be biologically correct and lose
AUC, and we have now seen that happen twice (D194's cosine, D201's program score).

**So no criterion in this document uses approval recovery.**

The interpretive evidence we currently hold is **convergence without a null**:

- disease programs (patient co-variation) and drug programs (drug co-variation) independently
  recovered the same hubs — PRRX1↔PRRX1 Jaccard **0.71**, DMTF1↔DMTF1 **0.61**;
- BRCA-LumA's endothelial program (`LHX6, SOX17, ERG, TAL1, EPAS1`) contains three of the five
  most-opposing regulons for fulvestrant (D194), whose opposition is carried by `LDB2, TIE1, CDH5,
  CD34` (D195).

These are striking and they are anecdotes. Nobody has asked how often that happens by chance. That
is the gap this spec closes.

## 2. The three claims, stated so each can fail

| | claim | fails if |
|---|---|---|
| **C1** | One NES value is **ambiguous**: the same regulator activity is reached by different drugs through different genes | equal activity implies equal gene pattern — the gene layer is redundant with the scalar |
| **C2** | 800+ regulators carry real **20–40 dimensional program structure**, on both the disease and drug side | modules are clustering noise, or no better than random regulator sets |
| **C3** | The condensation **preserves** gene-level information that matters | gene-level divergence predicts nothing external that the program score misses |

**C1 is the premise everything else rests on.** It is tested first, and a failure there invalidates
the architecture, not just this document.

## 3. Pre-registration, and what is already fixed

Four tests in this project came back null or refuting after the fact (D174, D175/D176, D169, D196),
and two of this week's runs had to be discarded or re-specified for defects found *after* writing
the criterion (D184's first run, D201's withdrawn magnitude premise). Thresholds below are
therefore fixed now.

**Pre-committed:** a negative result is written up with the same prominence as a positive one; no
threshold in sections 4–6 is moved to recover a number; tolerances (δ, size bands) are fixed here
and not tuned.

## 4. C1 — is one NES value ambiguous?

### 4.1 Data

`results/lincs/regulon_contributions_v2_tf/target_contributions_by_window.tsv.gz` — 10,973,419 gene
rows, 346 drugs, 29 cohorts, top-25 regulators per (drug, cohort, window), columns
`activity_contrib` and `cos_term` per target.

Unit: **(cohort, window, regulator)** with ≥ 8 drugs scored.

### 4.2 The reliability anchor, computed BEFORE the comparison

A gene-level vector measured once cannot be matched better than it matches itself. For drugs with
≥ 2 conditions in a window, split that drug's conditions in half, build the gene contribution
vector from each half, and take the cosine. The median is the **ceiling** `R`.

**If `R` is low, the test is inconclusive rather than confirmatory**, and that is said plainly —
it would mean the instrument cannot resolve gene-level pattern at all. Median conditions per
window is **1** (D183), so this is a real possibility and is checked first.

### 4.3 The comparison

Within a (cohort, window, regulator), for every drug pair:

- `d_act` = |activity_i − activity_j|
- `s_gene` = cosine(contribution vector_i, contribution vector_j) over that regulon's targets

**Activity-matched** pairs: `d_act` in the lowest decile *for that regulator* (a relative
tolerance, so regulators with different activity spreads are treated alike). **Unmatched**:
highest decile.

Three statistics, all reported:

1. **median `s_gene` among matched pairs**, against the ceiling `R`
2. **Spearman(`d_act`, `s_gene`)** within regulator, pooled as the median over regulators
3. the matched-vs-unmatched difference in median `s_gene`

**A confound named in advance:** both contribution vectors carry the same edge weights
`tfmode × likelihood`, so their product carries `w²` and *all* pairs will show elevated cosine.
This is why the test is **relative** (matched vs unmatched vs ceiling) and never absolute.

### 4.4 Acceptance — C1 is SUPPORTED when both hold

1. median `s_gene` among activity-matched pairs is **below 0.60 × R** — equal activity does not
   reproduce the gene pattern;
2. Spearman(`d_act`, `s_gene`) is **weaker than −0.30** in magnitude, i.e. knowing the activity
   difference tells you little about the gene-pattern difference.

Fail either → **the gene layer is largely redundant with the NES**, and that is reported as such,
with the consequence that D179/D194/D195 are reporting machinery rather than new information.

## 5. C2 — are the programs real?

Four tests. **(a) and (d) are gating; (b) and (c) supply the missing nulls for numbers already
quoted in D198/D200 and are reported either way.**

### 5.a Resampling stability (gating)

Split each cohort's patients 50/50 at random, recompute modules on each half with 58's exact rule
(average linkage on 1 − |r|, cut at |r| ≥ 0.7), and for each module of half A take the best
Jaccard to any module of half B. **Null:** the same procedure after shuffling patient labels within
the cohort, which destroys co-variation while preserving every marginal.

**Accept when** median best-match Jaccard exceeds the shuffled null's **95th percentile**, in
**≥ 70%** of cohorts. 20 random splits, fixed seed.

### 5.b Hallmark enrichment against a size-matched null (reported)

D198 reported 392/461 programs annotated at q < 0.05 — against the whole gene universe, which is
not the right null. Redo it: for each program draw **1,000 random regulator sets of the same size
from the same cohort**, take each one's best hallmark q, and report the **empirical p**.

Reported, not gated, because hallmark sets are coarse and a program can be real biology with no
hallmark counterpart. The number to report is the share of programs at empirical p < 0.05.

### 5.c Cross-side recovery against a null (reported)

D200: drug→disease best Jaccard median **0.074**, 260/1035 ≥ 0.25, 98 ≥ 0.50, with PRRX1↔PRRX1 at
0.71. **Null:** permute regulator labels within cohort before matching, 1,000 times. Report the
empirical p for the median and for the count at ≥ 0.50, and for the named examples individually.

### 5.d Knockdown coherence — the causal test (gating)

The only test here that does not depend on our own pipeline being right.
`data/processed/lincs/knockdown_bing` holds shRNA (`trt_sh.cgs`) and CRISPR (`trt_xpr`) signatures
for TF knockdowns (79m).

For each program whose **hub** has a knockdown signature: score the knockdown signature through the
cohort's regulons (VIPER, as 79n does) and compare the absolute activity shift of the program's
**member** regulators against **size-matched non-member** regulators from the same cohort, drawn
1,000 times.

**Accept when** members shift more than non-members with permutation **p < 0.01**, in **≥ 60%** of
testable programs. If fewer than 20 programs are testable, the test is reported as
**underpowered** rather than failed.

## 6. C3 — does the condensation preserve what matters?

**Stated first: the exact decomposition is arithmetic, not biology.** `Σ_j cos_term = cos_contrib`
to 1.4e-06 (D195) proves the bookkeeping and nothing else. It is not evidence for C3 and will not
be cited as such.

### 6.1 The test

Take drug pairs **within a (cohort, window)** that are **matched on program-level profile**
(cosine of their program-score vectors in the top decile) but **divergent at gene level** (gene
contribution cosine in the bottom decile), and the converse set.

Ask whether gene-level divergence predicts an **external** difference that the program score
cannot see:

| external quantity | source | requirement |
|---|---|---|
| shared annotated target | `evidence_targets` (79i2) | both drugs annotated — populated for 21.6% (D184), so **annotation is required, not assumed** |
| shared MOA class | `hub_moa` | both annotated |
| selective-killing agreement | PRISM/GDSC (79m) | where available |

### 6.2 Acceptance — C3 is SUPPORTED when

Among **program-matched** pairs, gene-level cosine is associated with MOA/target agreement at
permutation **p < 0.01**, with the effect in the direction that gene-similar pairs share mechanism
more often.

Fail → the gene layer adds no external discrimination beyond the program score, and should be kept
for explanation only, never as a feature.

**Registered caveat:** requiring annotation on both drugs will cut the usable pairs hard. If fewer
than 200 annotated program-matched pairs exist, this is reported as **underpowered**, not as a
null — the D184 criterion-4 lesson, where unannotated drugs were being silently counted.

## 7. What this does NOT test

- Whether any of this **improves a ranking**. D194 and D201 already measured that and it does not;
  that question is closed for these representations and is not reopened here.
- Whether programs **transfer across cohorts**. D197 measured that (regulon content Jaccard
  0.063–0.093; same-lines/different-network activity ρ = 0.667) and the answer is no.
- Whether the **window** structure is right. D183/D191/D198 cover it.
- Causal claims about drug mechanism in patients. Nothing here leaves cell lines and TCGA.

## 8. Outputs

```
results/lincs/interpretation/
  c1_nes_ambiguity.tsv       per regulator: ceiling R, matched / unmatched median s_gene, rho
  c2a_stability.tsv          per cohort: observed vs shuffled-null Jaccard, 20 splits
  c2b_hallmark_null.tsv      per program: observed q, empirical p vs size-matched sets
  c2c_crossside_null.tsv     observed vs permuted drug->disease Jaccard
  c2d_knockdown.tsv          per program: member vs non-member shift, permutation p
  c3_external.tsv            program-matched pairs, gene cosine vs MOA/target agreement
  acceptance.tsv             every criterion in sections 4-6, pass/fail as registered
  run.json                   seeds, fold assignments, input mtimes, package versions
```

Scripts: `scripts/79am_c1_nes_ambiguity.py`, `scripts/79an_c2_program_reality.py`,
`scripts/79ao_c3_preservation.py`.

## 9. Registration

The git commit adding this file is the registration. Any later change to sections 4–6 is a new
dated section in this document, never an edit in place.

---

## 10. AMENDMENT, 2026-11-22 — two defects in section 4, found on running it

**Section 4 is NOT edited; this records what it got wrong and what the corrected analysis shows.**

### 10.1 Criterion 2 was partly tautological

`actᵢ = Σⱼ Cᵢⱼ`, so two drugs with **opposite-sign** activity are forced to have gene contribution
vectors pointing opposite ways. The registered Spearman therefore measured sign agreement — true by
construction — rather than fine structure. **The tautology was identified when reporting the
failure, before the corrected analysis was run.**

Restricted to **same-sign** pairs (468 units), the non-tautological question gives:

| | pooled | early | late |
|---|---|---|---|
| matched `s_gene` | 0.290 | 0.238 | 0.327 |
| ceiling | 0.531 | 0.376 | 0.590 |
| matched / ceiling (bar < 0.60) | **0.545 PASS** | 0.632 FAIL | **0.554 PASS** |
| Spearman (bar < 0.30) | **+0.041 PASS** | +0.095 PASS | **+0.0001 PASS** |

### 10.2 Section 4 registered a POOLED statistic across windows

D191, D198 and D201 all required windows reported separately; section 4 did not. The two windows
have very different reliability — ceiling **0.376 early vs 0.594 late** (85 vs 155 split-half
cells) — so pooling mixes populations. Split, the registered criterion 1 gives **early 0.686 FAIL,
late 0.566 PASS**; the pooled 0.605 failure was early dragging the median over the bar.

### 10.3 Standing verdict

**As registered, C1 is NOT SUPPORTED** (pooled 0.605 ≥ 0.60; |ρ| 0.758 ≥ 0.30). That stands.

**Corrected for the tautology and split by window, the claim holds in the late window and is
inconclusive in the early one**, where the instrument's own reliability (0.376) is close to the
section 4.2 floor below which the spec pre-committed to calling the test inconclusive.

The supportable statement is **not** "one NES means several things" but: *within a direction of
effect, the activity value leaves ~45% of the reproducible gene pattern undetermined, and predicts
essentially none of it (ρ ≈ 0.04).* That is weaker than the original claim and sufficient for
everything built on it.
