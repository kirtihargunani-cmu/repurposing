# SPEC — D127 with one vote per target unit (D209)

**Status: PRE-REGISTRATION.** Committed before the deduplicated statistic is computed.
User, 2026-11-23: approved "register a target-deduplicated D127 and rerun".

---

## 1. Why

D127's sign test treats every compound as an independent instance. But the lead score
(`layered_S`, hence `layered_pct`) is a property of the drug's **target unit**, not of the drug
(D177). Compounds that share a target unit in both cohorts get identical percentiles, so they are
one observation counted several times. D207's negative control showed it plainly: 4 of 5 HDAC
inhibitors carry byte-identical scores (48.3 vs 19.7), so "5/5, p = 0.031" is really about 2
observations.

The sign test's p-value assumes independence. Counting copies inflates n and shrinks p.

## 2. What is already known — stated, not hidden

The **instance-level** results on the production arm are known (D208): primary 36/49, PASS;
negative control 5/5. The per-instance scores are in `results/lincs/d127/instances.tsv`. What has
**not** been computed is any statistic over deduplicated units. The rule below is fixed on
independence grounds alone, before that computation.

## 3. The rule

- **Unit of independence** (within one pair): the tuple
  `(layered_genes in the high cohort, layered_genes in the low cohort)`. Compounds that share it
  share the scored target unit on both sides.
- **One vote per unit:** the unit's gap = the **median** of its compounds' gaps
  (`layered_pct(high) − layered_pct(low)`); correct = gap > 0; a gap of exactly 0 is a tie and
  counts as not-correct (as in D127).
- A compound with `layered_genes` missing on either side is its own unit (it cannot be matched).
- **Same thresholds as D127, not re-opened:** share ≥ 0.70, one-sided sign test p < 0.05, median
  gap ≥ 10.0 percentile points, over the 7 primary pairs pooled.
- Same pair sets: primary (decision), seen (reported), D207 expansion sets and the HDAC negative
  control (reported, with the same s6.3 reading: control share ≥ 0.70 = the instrument is not
  trustworthy at this resolution).

## 4. Acceptance and what each outcome means

1. **Primary, deduplicated**, against the thresholds above. **This becomes D127's headline**,
   with the instance-level figure reported beside it, because it is the version whose p-value
   means what it says.
2. **Negative control, deduplicated:** reported with its unit count. If it collapses to 1–2 units,
   it is reported as **underpowered**, not as "passing" or "failing".
3. A fail at the unit level is reported as prominently as a pass, and the instance-level PASS is
   then described as resting on repeated observations.

**Pre-committed:** no threshold, unit definition or pair changes after seeing the result.

## 5. Registration

The git commit adding this file is the registration. Later changes go in a new dated section.
