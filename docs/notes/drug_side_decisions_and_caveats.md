# Drug side: decisions and caveats

An ongoing list for the drug side (branch `celligner`, worktree `~/projects/repurposing-celligner`). Add to it as
decisions are made and caveats come up. Detailed numbers and run history are in `docs/PROGRESS.md` and the commit
messages; this file is the short list.

- **Part 1:** your decisions.
- **Part 2:** method choices made while building, not yet reviewed by you.
- **Part 3:** caveats in the data and results.
- **Part 4:** open questions.

## 1. Decisions (yours)

| Date | Decision | Where it lives |
|---|---|---|
| 2026-10-08 | Cell-line matching is for drug applicability: a patient is matched to lines of the same organ (tissue key), not to the closest line anywhere | `config/celligner_tissue_keys.tsv`, `scripts/73` |
| 2026-10-08 | Substitute keys where the organ has no usable lines: thymoma -> T-cell lines; PCPG -> neuroblastoma lines | `config/celligner_tissue_keys.tsv` |
| 2026-10-08 | Hybrid matching for weights: the organ's lines, plus any other line closer than a patient's nearest organ line | `scripts/73` (MODE=hybrid), `scripts/74` |
| 2026-10-09 | A drug signature counts as active if it is high quality (`is_hiq`) and TAS >= 0.4 | `scripts/77`, `79b`, `79c` |
| 2026-10-09 | A drug not tested in the cancer's own tissue is borrowed if it was tested in >= 3 other cancer lines | `scripts/79b` |
| 2026-10-09 | Target check: a drug stays on a cancer's list if one of its own annotated targets is expressed (TPM >= 1) in >= 20% of the cancer's tumors | `config/drug_gate.yaml`, `scripts/79b` |
| 2026-10-09 | Drug-handling genes (metabolism enzymes, transporters) never count as targets | `config/drug_gate.yaml` |
| 2026-10-09 | No family rescue: a related gene (paralog, isoform) cannot stand in for the annotated target; drugs should be specific to their target | `config/drug_gate.yaml` |
| 2026-10-09 | Locally acting drugs (topical only, gut-lumen, eye, etc.; 50 drugs) are flagged, not removed | `scripts/79a` |
| 2026-10-12 | Blood cancers match myeloid or lymphoid lines separately, not one blood key | `config/celligner_tissue_keys.tsv` |
| 2026-10-12 | HA1E (non-cancer kidney line) stays excluded: it is never among any tumor's 10 closest lines | `scripts/76`, `76b` |
| 2026-10-12 | No tissue-tiered borrowing; the >= 3-line rule stays as it is | `scripts/79b` |
| 2026-10-12 | Drugs are scored with the same project ARACNe regulons as the disease side, full regulator set (TFs, kinases, signaling proteins; no Lambert-TF restriction), on measured + best-inferred LINCS genes | `scripts/76b`, `79c`, `79d` |
| 2026-10-12 | Approved drugs are the top-tier output; the list should not be narrowed by much | `scripts/79e` (tiers) |
| 2026-10-12 | Multi-target drugs are flagged (> 5 annotated targets), not removed | `scripts/79e` (`multi_target`) |
| 2026-10-12 | Disease-specificity score (option A) added next to the reversal score, both kept for comparison | `scripts/79e` (`spec`, `tier_spec`) |
| 2026-10-12 | Strata: Celligner line weights from each stratum's own patients | `scripts/74b` |
| 2026-10-12 | Strata: the target check stays as computed on the parent group's tumors | `scripts/79b` |
| 2026-10-12 | Strata ranked by reversal and specificity, including specificity against sibling strata (`spec_sib`, to build) | `main:docs/notes/subgrouping_for_drug_side.md` |
| 2026-10-12 | LumA PC3 cell-line issue (C1 below): noted, revisit only if it distorts rankings | this file |

## 2. Method choices made while building (for your review)

These were chosen to make the scoring work. None was a decision of yours; each is a setting that can be changed.

| Choice | Why | Setting |
|---|---|---|
| Each cell line's average drug response is removed from every drug profile (gene level) before VIPER | Almost every active drug slows proliferation; without this, 87% of drugs "reversed" bladder cancer and VIPER misread cytostatic drugs (see C3) | `scripts/79c` (specific profiles) |
| Drugs are ranked by percentile within each cancer; no p-values | The permutation null was 3-5x too narrow (~38% of drugs at p < 0.05) | `scripts/79e` |
| "Reversing" = top 10% of the cancer's drugs (and score > 0) | Simple and transparent | `TOP_PCT = 10` in `scripts/79e` |
| Master regulators used for scoring: FDR < 0.10 in the disease run | Same threshold as the disease-side summaries | `MR_FDR` in `scripts/79e` |
| Tested in the cancer's own tissue: average over those lines, weighted by Celligner; borrowed: plain average over all lines | The borrowing test found distance weighting adds nothing | `scripts/79e` |
| Specificity compares against cancers of other tissues only, each with >= 20 master regulators | Sibling subtypes would otherwise count against each other | `MIN_MRS_REF` in `scripts/79e` |
| Tier 1 = approved + passes target check + reversing; tier 2 = not approved + passes + reversing; tier 3 = rest | Follows "approved drugs first" | `scripts/79e` |

## 3. Caveats in the data and results

### C1. LumA PC3 low end gets its drug data mostly from HER2-type cell lines (2026-10-12)

**Status:** noted, no change. Revisit if it distorts the breast stratum rankings (user, 2026-10-12).

**The axis.** LumA PC3 splits luminal A patients by proliferation:
- High end: fast-growing tumors (AURKA, MYBL2 high).
- Low end: slow-growing tumors with strong estrogen-receptor activity (ESR1, GATA3 high), the most classically
  ER-driven luminal A tumors.

**The issue.** Each stratum's drug data comes from the cell lines its patients most resemble (Celligner weights,
`scripts/74b_stratum_line_weights.py`). Only six breast lines have LINCS drug data: MCF7, BT-474, SK-BR-3, BT-20,
MDA-MB-231 and Hs 578T. Among those six:

| Stratum | Leans mainly on |
|---|---|
| `BRCA-LumA-PC3-high` (fast-growing) | MCF7 (46%), ER-positive |
| `BRCA-LumA-PC3-low` (ER-driven) | BT-474 and SK-BR-3 (80%), both HER2-amplified; SK-BR-3 is ER-negative |

So the most ER-driven patients get their drug data mostly from HER2-type cells, one of them without the estrogen
receptor. That is the opposite of what the biology suggests.

**Likely cause.** The low end's closest lines overall are mostly lines without drug data (MFM-223, MDA-MB-453,
CAMA-1). Only 8% of its in-tissue weight falls on the six drug-tested lines (15% for the high end), so the split
among them rests on a thin slice of similarity. MCF7 is a fast-growing line, which likely ties it to the high end.

**How much it matters.** A drug's score only uses the lines where that drug was tested. Many drugs were tested only
in MCF7, so both ends use MCF7 for them and nothing changes. It matters for drugs tested in several breast lines:
for those, the low end's score leans on HER2-type cells.

**Possible fix if needed.** Let ER-driven strata use only ER-positive drug-tested lines (MCF7, BT-474).

Numbers: `results/celligner/stratum_line_weights_long.tsv` (`in_key` lines), `stratum_line_summary.tsv`.

### C2. The reversal scores do not yet validate on known drugs (2026-10-12)

Drugs with a known mechanism rank only moderately: on average at the 73rd percentile by reversal and the 53rd by
specificity (14 drugs; 50 = chance). Examples: erlotinib 96th in lung, 5-FU 93rd in colon, flutamide 88th in prostate,
but afatinib 18th in lung and irinotecan 41st in colon. Anti-estrogens rank at the very top in basal (ER-negative)
breast cancer and in the middle or lower in luminal A, also with a plain gene-level comparison. Likely cause: the
tumor-vs-normal comparison is dominated by what all tumors share (mainly proliferation), so strong growth-arresting
drugs "reverse" almost any cancer. Treat both current rankings as exploratory. Strata (subgroup contrasts) are the
planned fix. Numbers: `results/lincs/drug_reversal_controls.tsv`, `results/lincs/reversal_baseline/`.

### C3. Tumor-derived regulons misread growth-arresting drugs (2026-10-12)

The ARACNe networks are built from tumors, where lineage regulators are anti-correlated with proliferation (e.g.
luminal breast tumors grow slowly). Their target lists therefore include proliferation genes with a negative sign.
When a drug stops proliferation, VIPER reads it as switching those lineage regulators on: fulvestrant in MCF7 scores
FOXA1 +6 and GATA3 +5, although at the gene level it lowers the classic estrogen-response genes (GREB1, MYB, XBP1).
Removing each line's average drug response (part 2) reduces this but does not remove it; dropping Hallmark
proliferation genes did not help either.

### C4. General growth inhibitors dominate the reversal score (2026-10-12)

44 drugs (mainly mTOR, PI3K, MEK, IGF1R and CDK4/6 inhibitors) are in the top 10% by reversal in >= 20 of the 41
cancers. The specificity score removes this (0 drugs), but it also lowers the known-drug controls (C2).

### C5. Three cancers have too few master regulators to score (2026-10-12)

UCS (2), ESCA-SCC (3) and CESC-AD (10) master regulators on the drug genes. Their scores rest on very few
regulators; they are flagged `low_mrs` and left out of the specificity reference set.

### C6. Drug data is thin and uneven across tissues

- Only 73 cancer lines have active drug data. Many cancers have almost no drugs tested in their own tissue and rely
  on borrowed drugs: e.g. kidney, thyroid, testis, head and neck and soft-tissue cancers 0, glioblastoma 1.
  14 cancers have >= 100 drugs tested in their own tissue.
- Breast drug data is mostly MCF7; all anti-estrogen data is MCF7.
- Only 8-18% of a stratum's in-tissue Celligner weight falls on lines with drug data (see C1).

### C7. Borrowing drug effects across lines is approximate

A drug's effect in an untested line, predicted from other lines, correlates at r 0.46 with the real effect when the
drug was tested in >= 3 lines (best possible ~0.55). Matching lines by expression similarity did not improve this.
The test used the 978 measured genes, pooled doses and times, and mostly lines from different tissues.

### C8. Most drug genes are inferred, not measured

LINCS measures 978 genes; the other 9,196 used for regulator scoring are computationally inferred from them
("best inferred"). They are needed (without them only 11-20% of regulators can be scored, with them 89-92%), but they
carry less independent information than measured genes.

### C9. Target annotations are uneven

- Targets come from LINCS and the Drug Repurposing Hub and are direct binding targets, not downstream effects.
- Some drugs list every subunit of one target (e.g. tubulins for paclitaxel, proteasome subunits for bortezomib,
  GABA-A subunits for benzodiazepines), so `n_targets` overstates promiscuity for these; others are genuinely
  multi-target (multi-kinase inhibitors, antipsychotics) or noisy (caffeine 36 targets).
- 152 passing drugs have > 5 targets; 90 of them are approved (20% of the 460 approved passing drugs).
- 466-945 drugs per cancer have no annotated target and pass the target check by default.

### C10. The target check is a low bar

"Expressed in >= 20% of tumors at TPM >= 1" lets low-level expression pass: ESR1 passes in basal (ER-negative)
breast cancer, so anti-estrogens stay on the basal list (and rank at the top there, C2).

### C11. Blood cancers

- The four blood groups (B-ALL, AML x2, DLBCL) have no drug lists: they lack normal tissue for a disease signature.
- Before the myeloid / lymphoid split, ~30% of B-ALL and ~60% of DLBCL patients were closest to a line of the other
  lineage; they now match only their own lineage.

### C12. Cell-line matching measures similarity, not composition

A Celligner weight means "share of the cancer's patients that resemble this line", not the share of tumor cells
like it. Substitute keys (thymoma -> T-cell lines, PCPG -> neuroblastoma) are approximations.

### C13. Strata are close to their parents

Stratum signatures correlate 0.97-0.99 with their parent's (disease side). In the 41-group baseline, COAD and READ
already get near-identical drug rankings (rho 0.96), so reversal alone will mostly reproduce the parent's list;
sibling specificity (`spec_sib`) is what can separate strata. READ clusters are tentative and some breast axis ends
may be technical (see the disease-side note).

## 4. Open questions

- Should the top-10% "reversing" cutoff change (e.g. 20-25%) to keep more approved drugs in tier 1?
- If C1 distorts breast results: restrict ER-driven strata to ER-positive drug-tested lines?
- Should the target check be stricter (C10), e.g. a higher TPM or fraction for hormone receptors?
- Should subunit-type targets (C9) be collapsed before counting `n_targets` for the multi-target flag?

## 5. Likely questions and current answers

Questions a PI, reviewer or audience is likely to ask, with the answer as of 2026-10-12. Update the answers as the
work moves on.

### Cell-line matching
- **Why cell lines, and why Celligner?** Drug response data (LINCS) exists only for cell lines. Lines differ
  systematically from tumors (no immune cells or stroma, culture effects), so Celligner aligns the two first.
- **Why restrict matches to the same organ?** The question is whether a drug can act in that tissue, not which line
  is closest anywhere (decision, part 1). Hybrid matching also admits a closer non-organ line.
- **How well do the lines represent the tumors?** Partly. A weight means "share of patients that resemble this line"
  (C12). Only 8-18% of a stratum's matching weight falls on lines with drug data (C6).
- **Cancers without matching lines?** Substitutes: thymoma -> T-cell lines, PCPG -> neuroblastoma. Blood cancers
  match myeloid or lymphoid lines separately.

### Drug data (LINCS)
- **Which signatures count?** High quality with a clear effect (TAS >= 0.4): 29,812 of 90,988, 3,058 drugs in
  73 lines.
- **Doses and time points?** Averaged into one profile per drug x line; dose-response is lost.
- **Most genes are inferred.** Less independent information, but needed: with the 978 measured genes only 11-20% of
  regulators can be scored, 89-92% with inferred genes (C8).
- **Drugs never tested in the tissue?** Borrowed if tested in >= 3 other lines; predicted vs real effect r 0.46,
  best possible ~0.55 (C7). Weighting by line similarity added nothing.
- **Do dish effects reflect patients?** Not fully: no immune system, stroma or drug metabolism, and doses may exceed
  achievable levels. Standard limitation.

### Target check
- **Why 20% of tumors at TPM >= 1?** A low bar so a target present in a subset of patients keeps the drug; trade-off:
  low-level expression passes, e.g. ESR1 in basal breast cancer (C10).
- **Why no family rescue?** Target-specific drugs, to limit off-target binding (decision).
- **Drugs without an annotated target?** 466-945 per cancer pass by default; a missing annotation is not evidence of
  no target (C9).
- **Multi-target drugs?** Flagged, not removed; 90 of the 152 are approved; some just list every subunit of one
  target (C9).

### Reversal scoring
- **Why regulators, not genes?** Same mechanism-level view as the disease side; regulators summarize many genes.
- **Tumor networks applied to cell lines?** A known risk; it causes the growth-arrest misreading (C3).
- **Why remove each line's average drug response?** Almost every active drug slows proliferation; without the
  correction 87% of drugs "reverse" bladder cancer (part 2).
- **Where are the p-values?** The permutation null was 3-5x too narrow, so drugs are ranked by percentile instead.
  A calibrated null is unsolved.
- **Why the top 10%?** A transparent cutoff and one setting, not a statistical threshold.
- **Does it recover known drugs?** Partly: 73rd percentile on average; erlotinib, 5-FU, flutamide high; afatinib,
  irinotecan not; anti-estrogens wrongly high in basal breast cancer (C2). Expect this to be the main question.
- **Isn't it just finding anti-proliferative drugs?** Largely yes for the reversal score: 44 growth inhibitors are
  top 10% in >= 20 cancers. Specificity removes them but does not recover known drugs better (C4).

### Strata
- **Do subgroups get different drug lists?** Still to be shown: stratum signatures correlate 0.97-0.99 with the
  parent's, so sibling specificity is the test (C13).
- **Why line weights from the stratum's own patients?** Subgroups may resemble different lines, e.g. CMS1/3-like vs
  CMS2/4-like colon tumors lean on different colon lines.

### Overall
- **Better than CMap / LINCS connectivity?** Regulator level, tumor-matched cell lines, a tumor-expressed target
  check, subgroup focus, and a disease-specificity score.
- **How would you validate it?** Known-drug controls (current); drug sensitivity in the matched lines (DepMap PRISM,
  GDSC); clinical trial / approval concordance; eventually experiments. PRISM / GDSC is the obvious next check.
- **Why prioritize approved drugs?** Faster translation: an approved drug has known safety.
