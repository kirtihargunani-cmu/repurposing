# D128 implementation plan: issues found before A0 (2026-11-15)

Review of `docs/IMPLEMENTATION_target_mechanism_scoring.md` and its parent spec
(`docs/SPEC_target_mechanism_scoring_benchmark.md`, draft PR #2, branch
`spec/target-mechanism-scoring-20261114`) against the repository's data and earlier results.
**Nothing has been implemented or run.** The issues below should be settled, or recorded as A0 gaps, before A0
freezes anything.

## Settled: which disease gene effect to use

Production's primary Step 4 signature is limma with voomByGroup weights (D69): the columns `logFC`, `t`, `z` and
`adj.P.Val` in `data/processed/step4_v2/<group>__<contrast>.tsv.gz`. The `*_voom` columns are only a cross-check
(check C12, `scripts/30_signatures.R:34-36`). Step 6 VIPER reads the `z` matrix (`scripts/50_viper.R:80`).
`79t_target_role.py:115` reads `logFC`, using the GTEx contrast when its status is `ok` and the adjacent contrast
otherwise. **D128 does the same:** the disease direction comes from the sign of `logFC`, any disease weight is built
from `z`, and the contrast is chosen with 79t's rule.

## Open issues

### 1. R is likely to fail the mandatory ESR1 control: the reference problem (needs a decision)

Neither the spec nor the plan cites the 2026-11-12 precursor result (commit `55c80d7`, `docs/PROGRESS.md`,
`results/lincs/diagnostics/precursor_esr1_regulon/`). That test scored gene-level reversal on the genes in ESR1's
regulon:

- **ARACNe arm: inverted.** The anti-estrogens ranked near the bottom in BRCA-LumA, LumB and Her2, and ranked high
  only in BRCA-Basal, the ER-negative cohort.
- **CollecTRI arm:** no inversion, but no lift either.
- **Cause: the reference.** GTEx normal breast is ER-active, so ESR1's targets do not read as "up in tumor".

R as specified, `-sum w * sign(logFC) * f(z)` over target-linked genes, is essentially the same calculation. So A1's
same-drug comparison (LumA/LumB high vs Basal low) is predicted to come out backwards for R. M has no disease term
and is not affected, which makes **M the real ESR1 test**.

Options, to be fixed in A0 before any score is computed:
- Use adjacent normals instead. BRCA-LumA has 62 and LumB 22; Her2 has 9, so its adjacent contrast was skipped.
- Use a within-tissue contrast (LumA vs Basal) as the disease direction for breast.
- Report R under several references and mark them all exploratory.

This is the same open question as D48 / D57-D60 (adjacent normals in drug ranking).

### 2. A4 (EGFR) cannot be tested in a genotype-matched line

No EGFR-sensitizing-mutant line is profiled in LINCS (no HCC827, PC9 or HCC4006). H1975 carries L858R plus T790M,
which makes it resistant to gefitinib and erlotinib. It has 3 profiles, and only afatinib and neratinib were run
there. Every EGFR positive would come from a mismatched line, so under E16 those rows abstain or are flagged. A4 is
expected to be **not computable in context**, and A0 should record that.

What the other waves have (DepMap 24Q4 hotspot calls):

| Wave | In-context lines | Compounds there | Negatives |
|:--|:--|:--|:--|
| A1 ESR1 | MCF7 (PIK3CA) | fulvestrant, tamoxifen, raloxifene, bazedoxifene, toremifene, clomifene | anti-estrogens in non-ER lines (A375, HT29, A549) |
| A2 ERBB2 | SKBR3 (amplification still to confirm in A0) | lapatinib, neratinib, afatinib, canertinib, erlotinib | the same drugs in non-amplified lines |
| A3 BRAF | A375 (BRAF), SKMEL5 (BRAF) | vemurafenib, dabrafenib, PLX-4720, SB-590885; MEK inhibitors kept separate | RAF inhibitors in KRAS-mutant A549 / HCC515 (paradox) |
| B2 CDK4/6 | MCF7 | palbociclib | an RB1-null line is still to find |
| B3 PARP | MCF7 | olaparib, rucaparib, talazoparib, veliparib | |

### 3. The ESR1 positive control is thin (needs a decision)

- Fulvestrant is the only pure antagonist (`ANTAGONIST|DEGRADER` on ESR1). In MCF7 it has 3 `pert_id`s with 18 signatures in total (corrected 2026-11-15; an earlier count of one profile was wrong).
- The other anti-estrogens are SERMs, i.e. partial agonists, annotated `MODULATOR` on ESR1 (production orients them "lower" with a `direction_assumed` flag; tamoxifen's one `ANTAGONIST` row is on ESRRG). Under E07 (no default inhibitor sign) they abstain
  unless we decide to treat SERMs as antagonists in breast tissue. That should be a written decision.
- Estradiol is not profiled in MCF7, so there is no agonist control that should flip the sign.

### 4. Too many free choices; no primary arm named

The plan lists many sensitivity arms but does not say which one is primary. The following should be frozen
together in A0, before any score is read, so that the arm that passes is not the one that gets reported:

- primary network (ARACNe vs CollecTRI; they disagreed sharply for ESR1 in the precursor);
- landmark-only vs all genes;
- generic-removed vs raw profiles;
- hop limit;
- the weights `wM` and `wD`.

### 5. Removing the generic response may remove the ER signal

79c projects each line's generic response out of its profiles, and that response is largely proliferation. In MCF7,
ER drives proliferation, so part of the anti-estrogen footprint may be removed with it. The precursor did not test
raw profiles. Run both arms for the ESR1 control, with the primary named in advance.

### 6. `sign(logFC)` gives every gene a full vote

Every gene counts ±1 in R, however small its change. R therefore depends entirely on `wD`, which neither document
defines (e.g. `|z|`, a capped `|z|`, or an FDR cutoff). Fix it in A0.

### 7. The holdouts are less independent and more costly than they look

- A5 means rebuilding the Step 4/6 signatures and regulator scores for each new genotype group. The data needed are
  on disk: MC3 MAF, GISTIC thresholded copy number (`data/raw/tcga_pancan/`) and TCGA subtypes.
- But the drug-side evidence is not new. For example, BRAF inhibitors in A375 are the same drug-line profiles that
  D127 already used; only the patient grouping is new. Label these holdouts as new patient groups only.
- PRISM 20Q2 has already been read in `79j`, `79k` and `79m`, and possibly in D94. A0 must work out whether it is
  still clean enough to serve as independent validation. GDSC 8.5 looks less used; check it the same way.

### 8. Minor mismatches between the spec and the plan

The two documents number the stages differently (in the spec `a1_score` is scoring; in the plan A1 is the graph),
place the tests differently (`tests/test_d128_*` vs `tests/d128/`) and use different output formats (tsv vs
parquet). The implementation follows the plan. pyarrow 25.0.0 is in the `repurposing` env.

## Suggested order

1. Merge the spec branch into main (standing workflow: work on main, push after each step).
2. A0: loaders, the audit, and the gaps above recorded as data.
3. ESR1 M check early. It needs no disease term, so it is the cheapest test of whether D128 can work at all.

## Resolution (2026-11-15)

The plan now has a Section 0 that names a primary arm for every free choice (P1-P14) and sets out the two user
decisions (D0.1 reference for breast R; D0.2 SERMs) with recommended options. The workflow has been aligned with
the project's practice: commits on main, a D128-<stage> entry before each run, and a freeze tag. Issues 1, 3, 4, 5
and 6 are addressed by Section 0; issue 2 by the Section 5 status table; issue 7 by A0 action 7 and A5 action 5;
issue 8 by Section 1. The spec branch is merged into main. DECISIONS: row D128 (proposed; nothing run).
