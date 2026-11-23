# SPEC — Does gene-level detail separate drugs that `rev` scores identically? (D203)

**Status: PRE-REGISTRATION. Written and committed before `rev_gene` is computed for any pair.**

User, 2026-11-22: *"a regulator reaches its activity through genes unrelated to the disease pattern
gets down-weighted while one hitting the disease-relevant genes counts fully — that is exactly the
idea that I was going for"*, then: *"let's focus on separating drugs that rev scores identically."*

---

## 1. The claim, and why it is NOT a ranking claim

Two scoring attempts have already returned null (D194 contribution cosine, D201 program score), and
both were judged on approval-recovery AUC — an instrument D192 showed to be bimodal and D196 showed
the bare proliferation axis beats in half the cohorts.

This document registers **two** claims, both fixed before any number is computed.

**Claim A (user's prediction, 2026-11-22):** *"it should raise reversal AUC significantly because
of the downweighting."* Down-weighting regulators that reach their activity through
disease-irrelevant genes should remove false reversal and improve approval recovery.

**Claim B:** the narrower, instrument-independent one:

> Among drug pairs that `rev` scores **the same**, does the gene-level projection separate them,
> and does that separation correspond to a real mechanistic difference?

A score that merely reshuffles tied drugs at random would pass a "separation" test trivially. So
separation alone is **not** sufficient and §5 criterion 2 is the one that carries Claim B.

**The two claims are registered together but judged separately, and either can fail alone.** They
are not redundant: AUC answers "does it rank better on known uses", separation answers "does it
carry information `rev` cannot". A score can do the second without the first — D192 measured that
approval recovery is bimodal (9 cohorts significantly inverted) and D196 that the bare
proliferation axis beats the mechanism score in COAD/PRAD/LUSC, whose approved sets are cytotoxic.
**So §5 criterion 4 is paired with an attribution guard: an AUC gain that comes from increased
proliferation alignment is not a win, it is the confound returning.**

## 2. The score

```
act_geneᵢ = ⟨dᵢ , Cᵢ⟩ / ‖dᵢ‖          dᵢ = disease_contrib_nes over regulon i's genes
                                      Cᵢ = activity_contrib, the drug's gene contributions
rev_gene  = −Σᵢ nesᵢ · act_geneᵢ / ‖s‖
```

`rev` is identical except that `actᵢ = Σⱼ Cᵢⱼ` sums the gene contributions **blind to which genes
they are**. `act_geneᵢ` instead projects them onto the disease gene pattern.

**The normalisation is the whole design, and it is fixed here.** The divisor is **‖dᵢ‖ only** — the
disease vector, which is identical for every drug in a cohort. It removes the regulon-size effect
without touching drug magnitude. It is **not** divided by ‖Cᵢ‖: that is what `cos_contrib` does,
and the D202 decomposition attributed the entire loss of the cosine score to exactly that
per-regulon, per-drug normalisation (`rev / ‖D‖` cost +0.004; adding the per-regulon normalisation
cost −0.035). Repeating it would repeat a measured mistake.

A regulator whose activity comes from genes orthogonal to the disease pattern gets
`⟨dᵢ, Cᵢ⟩ ≈ 0` and drops out **by geometry, with no threshold and no tuned parameter**.

## 3. Data

`target_contributions_by_window.tsv.gz` — 10,973,419 gene rows, 346 drugs, 29 cohorts, top-25
regulators per (cohort, drug, window); `disease_contrib_nes` and `activity_contrib` are columns of
the same table, aligned by construction.

Unit of pairing: **(cohort, window)**. Pairs are formed only within a unit; nothing is compared
across cohorts (D197).

## 4. The noise floor, computed BEFORE any pair is examined

Separation is meaningless without knowing how much `rev_gene` moves for a drug compared **with
itself**. Using the D202/C1 split-half machinery (`79e3 SPLITHALF=1`), for every drug with ≥ 2
conditions in a window compute `rev_gene` from each half and record

```
Δ_self = |rev_gene(half A) − rev_gene(half B)|       expressed in within-cell percentile points
```

`Δ_self` median is the **floor F**. Any separation at or below F is noise.

**If F is so large that it exceeds the spread of `rev_gene` across drugs, the test is INCONCLUSIVE
rather than negative** — that would mean the instrument cannot resolve the quantity at all, and it
is said plainly.

## 5. Design and acceptance — fixed here

**`rev`-tied pairs:** within a (cohort, window), drug pairs whose `|Δ rev_pct|` falls in the
**lowest decile**. Reported alongside: pairs in the highest decile, as a contrast.

**Criterion 1 — separation exceeds noise.**
Median `|Δ rev_gene_pct|` among `rev`-tied pairs ≥ **2 × F**.
*Fails → the gene projection does not distinguish tied drugs beyond measurement error.*

**Criterion 2 — the separation is mechanistically real. THIS IS THE RESULT.**
Among `rev`-tied pairs where **both drugs are annotated** (`evidence_targets` / `hub_moa`),
pairs that share **no** annotated target and **no** MOA class must show larger `|Δ rev_gene_pct|`
than pairs that share either, with **permutation p < 0.01** (10,000 shuffles of the
shares/does-not-share label within cell).
*Fails → `rev_gene` separates tied drugs, but the separation tracks nothing mechanistic, and the
honest conclusion is that it is noise with structure.*

**Criterion 3 — the separation is not just `rev` leaking back.**
Spearman(`|Δ rev_pct|`, `|Δ rev_gene_pct|`) among tied pairs must be **< 0.30**; otherwise
`rev_gene` is re-expressing `rev` rather than adding gene information.

**Criterion 4 — Claim A: ranking.**
Approval-recovery AUC of `rev_gene` vs `rev`, per (cohort, window), **identical rows**:
median delta **≥ +0.03** and better in **≥ 60%** of cells — the same bar D201 was held to.
Reported per cohort, so the D192 bimodality is visible rather than averaged away.

**Criterion 5 — the attribution guard on criterion 4.**
If criterion 4 passes, Spearman(`rev_gene`, proliferation-program score) must **not** exceed
Spearman(`rev`, proliferation-program score) by more than **0.10**.
*A gain bought by aligning harder with proliferation is reported as the confound returning, not as
a win. Criterion 4 passing while 5 fails is written up as "recovers cytotoxics better", explicitly
not as "sharper mechanism".*

**Pre-committed:** annotation is **required, not assumed** — `evidence_targets` is populated for
21.6% of rows (D184), and D184's criterion 4 initially passed vacuously by counting unannotated
drugs as disjoint. If fewer than **200** annotated tied pairs exist, criterion 2 is reported as
**UNDERPOWERED**, never as a null. No threshold here moves after a number is seen.

## 6. What this does NOT claim

- Adoption. Even if criterion 4 passes, `rev` remains production until a promotion decision is
  taken on its own terms (the D186 pattern); passing earns the right to propose, nothing more.
- Not a causal claim about mechanism in patients.
- Not transferable across cohorts (D197).
- Does not revisit the window structure (D183/D191/D198).

## 7. Outputs

```
results/lincs/interpretation/
  d203_rev_gene.tsv.gz      per (cohort, drug, window): rev, rev_gene, both as percentiles
  d203_noise_floor.tsv      per cohort/window: F from the split halves
  d203_pairs.tsv.gz         tied pairs with |Δrev|, |Δrev_gene|, annotation agreement
  d203_acceptance.tsv       the three criteria, pass/fail as registered
  d203_run.json             seeds, counts, verdict
```

Scripts: `scripts/79an_rev_gene.py` (score + floor), `scripts/79ao_rev_gene_separation.py` (pairs).

## 8. Registration

The git commit adding this file is the registration. Later changes to §2–§5 go in a new dated
section, never an edit in place.

---

## 9. AMENDMENT 1, 2026-11-22 — sign error in the section 2 formula, found on the first smoke test

**Section 2 is NOT edited.** As registered it read

```
act_geneᵢ = ⟨dᵢ , Cᵢ⟩ / ‖dᵢ‖ ;   rev_gene = −Σᵢ nesᵢ · act_geneᵢ / ‖s‖
```

which **double-counts the disease direction**: `dᵢ` already carries it, since `Σⱼ dᵢⱼ = nesᵢ`
(D179). Multiplying by `−nesᵢ` again flips the sign. Caught on the first smoke test — TAL1 in
BRCA-LumA read `rev_contrib` **+0.29** (correct reversal: disease NES −3.27, drug activity +4.47)
against `rev_gene_contrib` **−0.032**.

**Corrected, and simpler than the registered version:**

```
rev       = −Σᵢ (Σⱼ dᵢⱼ)(Σⱼ Cᵢⱼ) / ‖s‖      a product of sums
rev_gene  = −Σᵢ  Σⱼ dᵢⱼ · Cᵢⱼ    / ‖s‖      a sum of products
```

`rev` multiplies the totals; `rev_gene` requires the genes to match one by one. A regulator whose
activity comes from off-pattern genes has `Σⱼ dᵢⱼCᵢⱼ ≈ 0` even when both totals are large. The
`‖dᵢ‖` divisor is dropped from the score (it was there to make regulons comparable, which the dot
product already does) and `act_gene = ⟨dᵢ,Cᵢ⟩/‖dᵢ‖` is **reported only**. `‖s‖` is a per-cohort
constant and does not affect within-cohort ranking.

**Nothing in sections 3–6 changes**; the criteria, bars, nulls and the attribution guard all stand
as registered. Post-fix sanity: Spearman(`rev_contrib`, `rev_gene_contrib`) = **+0.70** per drug,
and the score demotes near-orthogonal regulators (TGIF1, `cos_contrib` 0.016) while promoting
gene-concordant ones the totals miss (KLF4, `cos_contrib` −0.264).
