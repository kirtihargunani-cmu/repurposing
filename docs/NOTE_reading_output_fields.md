# Reading the output fields: which columns must be read together

**Why this page exists.** On 2026-11-20 the gate reported `gate_reason = target_expressed` for
tamoxifen in BRCA-Basal while the same row carried `target_absent = True` and a `spec_mech_pct` of
**4.1**. Both are correct. Read alone, the first says "this drug's target is present in basal breast
cancer", which is false. The fields answer different questions at different stages, and a reader who
takes one as a verdict will be wrong in a way no single column reveals (D168).

This is not a defect to be fixed by adding rules — D168 measured three versions of that fix and
rejected all of them. It is a **reading rule**, and it applies to humans, to downstream scripts and
to the LLM coordination layer.

## The rule

> **No gate or flag column is a verdict on its own. A drug's standing in a disease is the pair
> (stage-1 gate, stage-2 evidence), and the lead ranking is a third thing again.**

## The stages, and what each column actually means

| column | written by | question it answers | what it does NOT mean |
|---|---|---|---|
| `passes_gate` / `gate_reason` | `79b` | was this drug excluded at the **coarse** stage? | **not** "this drug is a candidate". It passes if **any** annotated non-ADME target is expressed. |
| `target_absent`, `evidence_targets`, `target_source` | `79i2` (D102a/D102c) | is a **curated mechanism** target of this drug expressed in this cohort? | only computed where mechanism knowledge exists; silence is not a pass |
| `tier_spec_mech` | `79i2` | 3 = demoted below every other drug of the group | — |
| `spec_mech_pct`, `rev_mech_pct` | `79i`/`79i2` | within-cohort percentiles of the mechanism-split reversal | a high `rev_mech_pct` with `target_absent = True` is the **D150 false positive being suppressed**, not a lead |
| `layered_pct` / `layered_S` | `79t3` | the **lead** ranking (D122) | it never reads the reversal; it is target-evidence driven (D165) |
| `rev`, `rev_raw` | `79e` | reversal on adjusted / unadjusted profiles | `rev_raw` is a diagnostic, not a score |

## The worked example, in full

Tamoxifen, BRCA-Basal, current build:

```
gate_reason      target_expressed          <- PKC targets are expressed (PLD1|PRKCD|PRKCE|PRKCQ|PRKCZ)
evidence_targets ESR1:0.86/11.1; ESR2; ESRRG
target_absent    True                      <- every MECHANISM target is absent
tier_spec_mech   3                          <- demoted below every other drug
rev_mech_pct     91.4                       <- strong reversal signal ...
spec_mech_pct    4.1                        <- ... correctly suppressed
```

Read `gate_reason` alone and tamoxifen looks viable in basal breast cancer. Read `rev_mech_pct`
alone and it looks excellent. Read the row, and the pipeline's actual answer is "no" — which is the
right answer: PKC is **not** a master regulator in BRCA-Basal (absent from the TF-only production
regulator set entirely, and non-significant in the full set, best FDR 0.103 with incoherent signs).

## Scale: this is normal, not rare

**1,192 rows across 171 drugs and all 41 cohorts** pass `79b` while their mechanism target is absent
— fluoxetine (40 cohorts), loperamide (40), paroxetine (38), tamoxifen (37), fluphenazine,
perphenazine, ceritinib, triflupromazine. All 1,192 are `tier_spec_mech = 3`, and **0 of the 558 FDA
top-10 rows are among them**. The layered design works; the single-column reading does not.

## Consequences for anyone consuming these tables

1. **Never filter on `passes_gate` alone** to build a candidate list. Use the lead ranking
   (`layered_pct`, or `fda_rank` from `79p`), which already incorporates the later stages.
2. **Never report `gate_reason` as evidence of target presence.** Quote `evidence_targets` and
   `target_absent` with it, or neither.
3. **A high `rev*_pct` beside `target_absent = True` is a suppressed false positive**, and saying
   so is the correct report — not "promising but gated".
4. **`passes_gate` counts are a coverage statistic, not a candidate count.**

## Related

- D168 (why no rule was added), D150 (the measured false positive), D102a/D102c (the `79i2` gate),
  D122 (`layered_pct` is the lead), D165 (the lead score never reads the reversal).
- `docs/LLM_COORDINATION_NOTES.md` R6 carries this as a reporting requirement.

---

## `layered_n_axes` and `layered_orient` — how much evidence, and in which direction (D177, 2026-11-21)

**Why these exist.** The user flagged that corticosteroids ranked at the top of LUSC, which is wrong.
Decomposing it exposed two real defects and one assumption that no column can fix.

**Defect 1 — the dependency inversion (removed).** For an agonist (`layered_orient == "raise"`) the
score used to invert the CRISPR dependency percentile, `100 - v`. NR3C1 sits at the **7.7th**
dependency percentile in LUSC — raw gene effect **+0.052**, i.e. knocking it out very slightly helps
growth — and was scored **92.3**. PPARG, true percentile **1.0**, was scored **99.0**. Knockout data
cannot be inverted into activation evidence; that needs CRISPRa/ORF, which we do not have (the LINCS
overexpression set covers 222 genes and **none** of our 41 agonist targets). D is now dropped, not
flipped, for `raise` units. Alteration frequency was already dropped.

**Defect 2 — the logFC fallback (removed).** `comp_overactivity` falls back to `logFC_vs_gtex` for
genes absent from the master-regulator table. For a tissue-foreign receptor that fallback is
automatic: GABRA1 and HTR1B are "lower in lung tumour than in normal lung" because they are
neuronal, not because anything regulatory was measured. **182 of 783** `raise` rows rested on it.
Overactivity now counts for a `raise` unit only when it is a real VIPER NES.

**What no column fixes.** NR3C1's NES in LUSC is **−2.076** — a genuine measurement. Within the
reversal framework, "suppressed regulator + drug that raises it" is a legitimate hit. The framework
assumes *reversing a regulator's activity reverses the disease*, and for a regulator suppressed as a
**consequence** of dedifferentiation (GR in squamous epithelium) that assumption is false. Nothing in
the activity data distinguishes driver-suppression from consequence-suppression. So corticosteroids
are **flagged, not demoted** — the user's call, over demoting them with a tuned penalty, because the
penalty size would have been fitted to LUSC and that is how we would fool ourselves.

### The reading rule

> **`layered_pct` is a rank, not a weight of evidence. Always read it beside `layered_n_axes`.**
> A score built on one axis is not a score built on three, and the percentile does not say which.

| column | values | read it as |
|---|---|---|
| `layered_n_axes` | 1, 2, 3 (NaN on the essential-unit / DNA-acting path, which scores from potency instead) | how many of dependency / alteration / overactivity actually contributed |
| `layered_orient` | `lower` (inhibitor) / `raise` (agonist) | an **agonist can only ever reach 1 axis**, because D and A are structurally unavailable to it |

Census after D177: **4,465 rows on 3 axes, 648 on 2, 617 on 1.** `T` is shrunk toward the null by
`sqrt(n_axes / 3)`, so a one-axis score is pulled 42% of the way to 50 — but that shrink does **not**
make a one-axis hit equivalent to a three-axis one, and it was never meant to. It reordered LUSC's
top 4 (floxuridine/fluorouracil on TYMS, alpelisib on PIK3CA, palbociclib on CCND1|CDK6 — all
3-axis, all genuine LUSC lesions) and left the corticosteroid block at slots 5–10. **That residue is
expected and is what the flags are for.**

On the rankings page both appear in the decode string (`[1 of 3 axes; AGONIST: direction inferred,
no dependency/alteration evidence possible]`) and as filterable flags (`agonist (direction
inferred)`, `1-axis evidence`).

**Open question, deliberately left to data:** whether raising a suppressed regulator ever helps is
an empirical question, not a scoring one. D178's phase-progression test covers **211 agonist pairs**
with ChEMBL indication records and can answer it directly.
