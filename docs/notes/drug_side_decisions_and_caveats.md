# Drug side: decisions and caveats

An ongoing list for the drug side (branch `celligner`, worktree `~/projects/repurposing-celligner`). Add to it as
decisions are made and caveats come up. Detailed numbers and run history are in `docs/PROGRESS.md` and the commit
messages; this file is the short list.

How the steps run (v2): [`../PIPELINE_DRUG_SIDE_V2.md`](../PIPELINE_DRUG_SIDE_V2.md) and `scripts/run_v2_drug.sh`.

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
| 2026-10-12 | Strata: the target check is judged on each stratum's own tumors (same rule: TPM >= 1 in >= 20%); tissue key and direct / borrowed split are the parent's | `scripts/79b` (STRATA=1) |
| 2026-10-12 | Strata ranked by reversal and specificity, including specificity against sibling strata (`spec_sib`) | `scripts/79g` |
| 2026-10-12 | Her2 PC4 kept but flagged weak; flag carried into the drug outputs | `scripts/79g` (`weak`) |
| 2026-10-12 | LumA PC3 cell-line issue (C1 below): noted, revisit only if it distorts rankings | this file |
| 2026-10-20 | Borrowing stays; no same-organ or same-body-area restriction. Instead, a confidence for each borrowed drug (predicted transfer accuracy), first version diagnostic only (C17) | `scripts/79l` |
| 2026-10-28 | v2: the drug side scores the TF-only regulator set (build `v2_tf_xmap`), so drug and disease are compared on the same TFs. Replaces "full regulator set" of 2026-10-12 for v2 | `scripts/79d` REGSET=tf |
| 2026-10-28 | Every regulator counts, weighted by its disease NES; no master-regulator FDR cutoff (main D74). Why: the FDR lists changed a lot with small changes to the null (colon 425 / 0 / 1 master regulators) while the ranking of regulators did not | `scripts/79e`, `79i` (SCORE=weighted) |
| 2026-10-28 | Drugs are scored with a consensus of two regulon sets, ARACNe TF and CollecTRI (main D75; C18, C19). Why: the tumor network misreads anti-estrogens; CollecTRI reads them correctly; the knockdown benchmark was too thin to pick one set per TF, so both are averaged | `scripts/79e`, `79i` (REGULON=consensus) |
| 2026-10-28 | Flag drugs tested only in cell lines that do not resemble the disease (`unmatched_lines`, line match < 0.1); they stay ranked but not in tiers 1-2 (main D76; C20) | `scripts/79e` |
| 2026-10-28 | Borrowed drugs: transfer confidence refit with a distance term; flag `low_transfer` when confidence < 0.5; same treatment as `unmatched_lines` (main D77; C21) | `scripts/79l` VERSION=v2, `79e` |
| 2026-10-28 | Every disease and subgroup reports its evidence status: matched / borrowed only / no adequate compound, with a sentence for reports (main D78; C22). Requirement R1 for the LLM layer (`main:docs/LLM_COORDINATION_NOTES.md`) | `scripts/79e` summary |
| 2026-10-29 | Reports lead with disease specificity (`spec_mech`, then `spec`); overall reversal is context (main D80; R2) | reports / LLM layer |
| 2026-10-29 | Subgroups (strata) for v2 are the 14 of main D79; ranked by `79e` STRATA=1 (replaces `79g` for v2) with sibling specificity | `scripts/79e`, `74b`, `79b` |
| 2026-10-29 | PRISM / GDSC validation of the v2 scores on hold | - |

## 2. Method choices made while building (for your review)

These were chosen to make the scoring work. None was a decision of yours; each is a setting that can be changed.

| Choice | Why | Setting |
|---|---|---|
| Each cell line's average drug response is removed from every drug profile (gene level) before VIPER | Almost every active drug slows proliferation; without this, 87% of drugs "reversed" bladder cancer and VIPER misread cytostatic drugs (see C3) | `scripts/79c` (specific profiles) |
| Drugs are ranked by percentile within each cancer; no p-values | The permutation null was 3-5x too narrow (~38% of drugs at p < 0.05) | `scripts/79e` |
| "Reversing" = top 10% of the cancer's drugs (and score > 0) | Simple and transparent | `TOP_PCT = 10` in `scripts/79e` |
| Master regulators used for scoring: FDR < 0.10 in the disease run | Same threshold as the disease-side summaries (v1; v2 weights every regulator by NES instead, D74) | `MR_FDR` in `scripts/79e` |
| Consensus: a regulator is averaged over the sets where it is scored on both the disease and the drug side (so drug activity and disease NES always come from the same regulons) | Avoids comparing a drug's CollecTRI reading with a disease's ARACNe reading | `consensus()` in `scripts/79e`, `79i` |
| Consensus split: a regulator is proliferation-linked if it is linked in either set; the proliferation adjustment is done per set before averaging | Conservative: anything that may carry growth arrest is kept out of the mechanism score | `scripts/79i` |
| Line-match cutoff 0.1 (drug's best tested line has < 10% of the Celligner weight of the disease's best line) | Chosen from the distribution before any ranking was looked at: it is bimodal, a quarter of direct pairs below 0.01, half above 0.19 | `LINE_MATCH_MIN` in `scripts/79e` |
| Transfer-confidence cutoff 0.5 ("transfers less than half as well as a replicate in the same line") | Chosen from the distribution before rankings (no natural break; 11.5% of borrowed pairs below) | `CONF_MIN` in `scripts/79l` |
| Subgroup borrowed-drug confidence = the parent group's | Same drug data and lines; the distance term carries no information (C21) | `scripts/79e` STRATA=1 |
| Tested in the cancer's own tissue: average over those lines, weighted by Celligner; borrowed: plain average over all lines | The borrowing test found distance weighting adds nothing | `scripts/79e` |
| Specificity compares against cancers of other tissues only, each with >= 20 master regulators | Sibling subtypes would otherwise count against each other | `MIN_MRS_REF` in `scripts/79e` |
| Tier 1 = approved + passes target check + reversing; tier 2 = not approved + passes + reversing; tier 3 = rest | Follows "approved drugs first" | `scripts/79e` |
| Proliferation-linked regulators: targets enriched for Hallmark E2F + G2M (hypergeometric FDR < 0.05, fold >= 2) or drug-induced activity tracking the drug's proliferation effect (abs r >= 0.5) | Fixed before looking at drug results; diagnostic split (C15) | `scripts/79h`, `R_MIN` in `scripts/79i` |

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

**v2 update (2026-10-29):** with consensus regulons, NES weighting and the flags: erlotinib 96th in lung (osimertinib
#1), flutamide 85th in prostate, irinotecan 98.5th on specificity in colon. Basal-breast anti-estrogens still score
high but are flagged (C20) and kept out of the tiers. Luminal A anti-estrogens stay mid-ranked once proliferation is
removed (fulvestrant 40th on the mechanism part). `results/lincs/split_controls_v2_consensus.tsv`.

### C3. Tumor-derived regulons misread growth-arresting drugs (2026-10-12; cause corrected 2026-10-12)

The ARACNe networks are built from tumors, so a regulator's target list is the genes that co-vary with it across
tumors, not the genes it switches on and off in a few hours. Fulvestrant in MCF7 scores FOXA1 +6 and GATA3 +5,
although at the gene level it lowers the classic estrogen-response genes (GREB1, MYB, XBP1).
- First explanation (proliferation genes as negative targets) is not the main cause: FOXA1, GATA3, ESR1 and XBP1
  in the breast network have almost no Hallmark proliferation genes among their targets (6, 3, 1, 5 of 483, 354,
  201, 290), and their drug-induced activity correlates only -0.2 to -0.4 with the drug's proliferation effect.
- Actual cause: the breast ESR1 regulon is a luminal-identity module (positive targets CA12, FOXA1, XBP1, SCUBE2,
  SLC39A6; negative targets basal / stress genes S100A8, KRT7, PHGDH, PSAT1). GREB1, PGR, TFF1 and MYB are largely
  not in it. In MCF7, fulvestrant leaves the luminal genes about flat (GATA3 even rises) and lowers the negative
  targets, so VIPER reads "ER program up".
- Consequence: anti-estrogens look anti-reversing in luminal tumors (ER program up in the tumor) and reversing in
  basal tumors (ER program down). Removing proliferation does not change this (C15).
Removing each line's average drug response (part 2) reduces this but does not remove it.
**v2 update (2026-10-28):** see C18 and C19. The knockdown benchmark shows the ESR1 regulon reads ESR1 *knockdown*
correctly, so the regulon is not wrong in general: anti-estrogen drugs raise GATA3 / FOXA1 mRNA, which a knockdown
does not. The consensus with CollecTRI (direct ER targets) corrects the direction in luminal breast.

### C4. General growth inhibitors dominate the reversal score (2026-10-12)

44 drugs (mainly mTOR, PI3K, MEK, IGF1R and CDK4/6 inhibitors) are in the top 10% by reversal in >= 20 of the 41
cancers. The specificity score removes this (0 drugs), but it also lowers the known-drug controls (C2).

### C5. Three cancers have too few master regulators to score (2026-10-12)

UCS (2), ESCA-SCC (3) and CESC-AD (10) master regulators on the drug genes. Their scores rest on very few
regulators; they are flagged `low_mrs` and left out of the specificity reference set.
**Resolved in v2 (2026-10-28):** with NES weighting (D74) every tested regulator counts, so no cancer lacks regulators
(317-436 regulators effectively carry each score).

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

### C14. Sibling specificity separates only some strata, and the ER check fails (2026-10-12)

- Most stratum rankings are near their parent's: median rho 0.97 by reversal (top-10% overlap 0.75); siblings 0.92
  (0.55). Basal strata are essentially the parent (rho 0.99-1.00; `spec_sib` spread only 0.2-0.5).
- The most distinct strata: READ-C2 (rho vs parent 0.76), LumA PC2 / PC3 (0.84-0.89), COAD CMS1/3-like (0.91),
  Her2 PC3, LumB PC1.
- Check on the LumA proliferation axis (PC3): palbociclib is most specific to the proliferation end (100th
  percentile), as expected. But anti-estrogens are also more specific to the proliferation end than to the ER end
  (`spec_sib` 1st-5th percentile at the ER end), the opposite of the biology. Same pattern for the Her2 PC3 ESR1 end.
  Likely the same cause as C2 / C3: in MCF7 the anti-estrogen effect is mostly growth arrest, which the scoring
  reads as reversing proliferation, not ER signaling.
- Tier overlap between two siblings is ~0 by construction: with two siblings, `spec_sib` of one is roughly minus
  that of the other.
Numbers: `results/lincs/stratum_reversal_summary.tsv`, `stratum_reversal_controls.tsv`.

### C15. Splitting off proliferation does not fix the anti-estrogen result (2026-10-12)

Test (`scripts/79h`, `79i`): each cancer's master regulators split into proliferation-linked (targets enriched for
Hallmark E2F / G2M genes, or drug-induced activity tracking the drug's proliferation effect, |r| >= 0.5; median 17%
of a cancer's MRs) and the rest. New columns, all kept: `rev_prolif` (reversal of the proliferation MRs), `rev_mech`
(the other MRs, after removing the proliferation-explained part of every regulator's activity), `rev_nonprolif`
(the other MRs, no adjustment), `prolif_effect` (how much the drug lowers proliferation genes), `spec_mech`,
strata `spec_sib_mech`, `prolif_share`, and a `prolif_driven` flag.
- Generic drugs: top 10% in >= 20 cancers: 36 by `rev` -> 29 by `rev_mech`; 75 by `rev_prolif`; 0 by either
  specificity score. Only a modest reduction.
- Known drugs (mean percentile, expected-high controls): `rev` 73, `rev_prolif` 85, `rev_mech` 55, `spec` 53,
  `spec_mech` 41. Erlotinib (96), gefitinib (87), flutamide (87) keep their ranks on the mechanism part;
  irinotecan rises (41 -> 72); anti-estrogens and palbociclib drop in luminal A (fulvestrant 58 -> 22).
- Anti-estrogens in basal stay at the 91st-98th percentile by `rev_mech`; at the LumA PC3 ER end `spec_sib_mech` is
  3rd-11th percentile (was 1st-5th). Not fixed.
- Reason: the ER misreading is a regulon-content problem, not a proliferation one (C3, corrected).
- What it does show: proliferation reversal carries most of the known-drug signal (`rev_prolif` 85) and of the
  anti-estrogen and palbociclib signal in luminal tumors, in line with your point that effective drugs often work by
  stopping proliferation. Lineage-specific drugs whose signal is not proliferation (erlotinib, flutamide) are
  carried by the mechanism part.
Numbers: `results/lincs/split_overview.tsv`, `split_controls.tsv`, `drug_reversal_split.tsv`,
`stratum_reversal_split.tsv`, `prolif_regulators.tsv`, `prolif_regulators_p2.tsv`.

### C16. The scores do not predict measured drug sensitivity (PRISM / GDSC, 2026-10-12)

Test (`scripts/09b`, `79j`, `79k`; protocol fixed before results): for each cancer and stratum, the drug's selective
killing of the cancer's own cell lines (Celligner-weighted, z-scored per drug across all lines so general toxicity
is removed) in PRISM Repurposing 20Q2 (secondary 8-dose, 681 scored drugs; primary 1-dose, 1,135) and GDSC1 / GDSC2
release 8.5 (271 / 184 scored drugs).
- The data works as a yardstick: known drugs are selective where expected (BRAF / MEK inhibitors in melanoma
  97th-100th percentile, HER2 drugs in Her2 breast 86-98, fulvestrant in luminal breast 96-99 and in basal 4-84,
  EGFR inhibitors in lung adenocarcinoma mostly 80-99).
- Our scores do not predict it. Median within-cancer Spearman with selective killing: -0.04 to +0.04 for every score
  (rev, rev_prolif, rev_nonprolif, rev_mech, prolif_effect, spec, spec_mech) in every dataset; the range across
  cancers is about -0.3 to +0.3. Top-10% enrichment is about 1x (chance). The few nominal p < 0.05 (PRISM primary
  spec 0.03, rev_nonprolif 0.02) are tiny effects and do not repeat across datasets.
- Across cancers per drug (does a drug score highest where it kills most): median rho 0.00-0.06; same picture.
- Strata: `spec_sib` does not predict which sibling's lines are more sensitive (median rho 0.00-0.03; sibling line
  weights overlap ~0.72, so the contrast is weak anyway).
- Target expression alone (fraction of targets expressed in the tumors) does no better within cancers.
- Anti-estrogens: PRISM / GDSC show the biology our scores invert (C3): selective in luminal, not in basal.
Caveats of the test: viability in 3-5 day assays is not reversal; LINCS signatures are 6-24 h; most drugs have
LINCS data in few lines (borrowed), while PRISM / GDSC measure many lines.
Numbers: `results/lincs/sensitivity_validation/` (`t1_summary.tsv`, `t1_per_unit.tsv`, `t2_across_units.tsv`,
`t3_siblings.tsv`, `controls.tsv`, `overlap.tsv`).

### C17. Confidence for borrowed drugs (2026-10-20; first version, diagnostic)

Borrowed drugs are 73% of scored disease x drug rows (37,017 of 50,760), 67% of tier 1 and a median 9 of each
cancer's top 10 by `rev`. 12 cancers have no drug tested in their own organ (kidney x3, Wilms, thyroid, head and neck,
adrenal, testis x2, sarcoma x3); the organ's tumor lines in LINCS have no high-quality signature, and thyroid, testis,
adrenal and soft tissue have no LINCS line at all. Restricting to the same body area was considered and not taken
(user 2026-10-20): drugs cross organs through shared targets / pathways (guselkumab: psoriasis -> ulcerative colitis
and Crohn's), and the target check already asks whether the target is there.
Test (`scripts/79l`): predict a borrowed effect's accuracy from the leave-one-line-out data (`scripts/77`, 7,661
held-out pairs). Predictors: the drug's consistency across its other held-out lines, how many lines that rests on,
and how many lines are averaged. OLS, cross-validated R² 0.60 (Spearman 0.76) by drug and by held-out line alike,
0.63 where the nearest tested line is in another lineage; lines only 0.17; distance to the nearest tested line adds
0.006 (not used: the committed Celligner tables keep only 10 neighbors per line). Calibrated within 0.03 in every
decile. `conf` = predicted r / replicate level (0.553), capped at 1; direct drugs set to 1 (assumed). `rev_conf` =
`rev` x `conf`, with percentile and tier as in 79e.
- Borrowed `conf`: 11.5% below 0.5 (4,256 pairs), 9,810 at 1.
- Rankings move modestly: 752 of 877 tier-1 rows stay, 333 of 410 top-10 places; most change in the borrowed-heavy
  cancers (PCPG, pancreas, adrenal, head and neck, sarcoma UPS, bladder 5-6 of 10 kept). 45 tier-1 borrowed rows had
  `conf` < 0.5.
- Caveat: consistency favors drugs that act the same everywhere (top: dactinomycin, panobinostat, epirubicin,
  emetine, vorinostat) and penalizes target-dependent drugs (tamoxifen 0.11, bazedoxifene 0.15, idelalisib 0.19,
  gefitinib 0.33), whose effect varies with the target, not only with noise. Using `rev_conf` as is would add to the
  generic-drug problem (C4). Proposed v2: fit transfer accuracy on lines that express the drug's targets (needs
  DepMap expression per line).
- Run with system Python (pandas 1.3.5); the `celligner` env is not on MAIL-15.
Numbers: `results/lincs/transfer_confidence_model.tsv`, `transfer_confidence_calibration.tsv`,
`drug_transfer_confidence.tsv`, `transfer_confidence_summary.tsv`.
**v2 (2026-10-28):** now used as a flag, see C21.

### C18. Anti-estrogens: the tumor-network ESR1 regulon misreads the drugs, not ESR1 itself (2026-10-28)

**What we saw.** In v2 (TF-only network) anti-estrogens in MCF7 still raise ESR1, FOXA1, GATA3 and XBP1 activity
(ESR1 z +1.19 against 150 random breast-line drug profiles). At the gene level the drugs do what they should (GREB1,
CCND1, MYB down) but GATA3 (+3.9) and FOXA1 also go up, and those are ESR1's strongest positive targets in the tumor
network.
**What we tried** (scratch tests, 2026-10-28):

| Scoring of ESR1 | z | Right direction? |
|---|---|---|
| ARACNe pooled BRCA network | +1.19 | no |
| + shadow / pleiotropy correction | +1.07 | no |
| dropping targets that are themselves regulators | +1.19 | no |
| ARACNe luminal-only network (BRCA_LUM) | +0.76 | no, half the error |
| ARACNe targets also in CollecTRI | -0.32 | weakly |
| CollecTRI ESR1 regulon | -1.69 | yes |

**Why.** A tumor network learns which genes rise and fall together across patients, and in breast tumors the lineage
genes (GATA3, FOXA1) move with ER. A drug acting for hours changes the direct ER targets, not lineage identity.
Shadow correction cannot help: it reassigns targets shared between regulators, while here the targets respond in
opposite directions.
**What we did:** consensus scoring with CollecTRI (C19).

### C19. Knockdown benchmark: no regulon set is clearly better; scored with a consensus (2026-10-28)

**Test** (`scripts/79m`, `79n`): LINCS shRNA (consensus of several hairpins) and CRISPR knockdowns of TFs in cancer lines.
Rules fixed before the run: the TF must be expressed in the line and its own mRNA down where measured. 1,803 of 7,608
signatures usable, 336 TFs, 18 lines. A regulon set is "right" for a TF if the knocked-down TF scores among the 25% most
decreased regulators (median) and negative in >= 60% of signatures (>= 3 signatures).
**Result.** All sets detect knockdowns only weakly and about equally (AUC 0.55 shRNA / 0.63 CRISPR; TFs passing 22%
ARACNe, 24% CollecTRI, 30% DoRothEA; ARACNe vs CollecTRI better in 51% of shared TFs, p = 0.6). ESR1 and FOXA1
knockdowns are read correctly by ARACNe (C18). With 1-5 signatures per TF, choosing a set per TF would rest on noise
(the pre-set rule would have dropped 143 TFs), so the per-TF choice was not applied.
**Decision (user, D75):** consensus of ARACNe TF and CollecTRI for every regulator. Result on breast controls: luminal A
fulvestrant 57th -> 72nd percentile, luminal B anti-estrogens 66-93 -> 75-95.
Numbers: `results/lincs/knockdown_benchmark_{summary,tf}.tsv`, `knockdown_regulon_choice.tsv` (not used).

### C20. Drugs tested only in cell lines unlike the disease (line-match flag, 2026-10-28)

**What we saw.** Basal-breast anti-estrogens ranked 97-100th. LINCS tested them only in MCF7 (luminal); basal breast
gives MCF7 a Celligner weight of 0.0001. When no tested line carries weight, the score falls back to equal weights over
the tissue's lines, so basal breast was scored on a luminal line's response.
**Measure.** `line_match` = Celligner weight of the drug's best tested line in the tissue / the disease's top line in
the tissue. Absolute coverage did not work (Celligner spreads weight over ~1,000 lines, so half of all pairs cover
< 1%); the ratio separates the cases (bimodal).
**Rule (user, D76):** `unmatched_lines` = direct drug with line_match < 0.1: kept in the ranking, kept out of tiers 1-2.
Flagged most: basal and Her2 breast (MCF7), brain (LN-229), ESCA-SCC (AGS, gastric), LUSC (A549, adenocarcinoma),
PCPG (SK-N-SH, neuroblastoma), prostate (986 of 1,587 drugs, mostly PC-3, AR-negative).

### C21. Borrowed drugs: transfer confidence as a flag; distance adds nothing (2026-10-28)

**Why not flag all borrowed drugs:** 73% of rows are borrowed and many cancers (and 12 of the 14 subgroups) have no
tested line of their own; their tiers would be empty although many drugs act the same in every line.
**v2 model** (`scripts/79l` VERSION=v2, runs before `79e`): the C17 consistency model plus an excess-distance term (how
much further the drug's tested lines are from the disease than its nearest DepMap line, Celligner 70-PC space, the same
definition for training and for disease patients). **Distance adds nothing** once the drug's consistency is known
(cross-validated R^2 0.602 -> 0.603 by drug, 0.596 -> 0.586 by held-out line): how well a borrowed drug carries over
is a property of the drug.
**Rule (user, D77):** `low_transfer` = borrowed with confidence < 0.5 (11.5% of borrowed pairs, 132 drugs): kept in
the ranking, kept out of tiers 1-2. The C17 caveat stays: consistency favors drugs that act the same everywhere.

### C22. Evidence status per disease (2026-10-28)

**Why (user):** gives a number-based way to say that no compound was tested in a line resembling the disease, so an
adequate compound cannot be reported, instead of showing the top of a ranking built on dissimilar lines.
**Columns** (`drug_reversal_summary*.tsv`, `stratum_reversal_summary*.tsv`): `best_line`, `best_line_match` (the best
any tested drug could reach), `n_t12_matched`, `n_t12_borrowed`, `evidence_status` (matched / borrowed only / no
adequate compound), `evidence_note` (sentence for reports).
**v2 result:** 14 of 41 groups matched (breast LumA / LumB / Her2, colon, rectum, lung adenocarcinoma, liver, prostate,
melanoma, pancreas, stomach, esophageal adenocarcinoma, cervical adenocarcinoma, neuroblastoma); 27 borrowed only.
Subgroups: only the two COAD subgroups matched. Uterine groups have a well-matched line (MFE-319) but none of its drugs
reaches tiers 1-2.

### C23. v2 results: what the rankings say (2026-10-29)

- By overall reversal (`rev`) the same drugs top most cancers (32 drugs in the top 10% of >= 20 groups: statins,
  palbociclib, MEK inhibitors, anthracyclines). By specificity with proliferation removed (`spec_mech`) none does.
  Hence reports lead with specificity (D80).
- Subgroups: COAD MSI-like is led by BRAF / MEK inhibitors, and HSP90 inhibitors are the most MSI-specific (MSI colon
  cancer is enriched for BRAF V600E, an HSP90 client); THCA BRAF-like by HDAC inhibitors (studied for re-differentiating
  BRAF thyroid cancer); HNSC atypical-classical (HPV-enriched) by the MDM2 inhibitor AMG-232; BLCA basal-squamous by
  Src / MEK inhibitors. LGG and TGCT-SEM subgroups barely differ in drug response (`spec_sib` <= 0.6 and <= 0.2), so
  they should be reported as one drug list.
- Not yet validated against measured drug sensitivity (PRISM / GDSC on hold).
- Subgroups with proliferation removed (2026-10-29, user; `79i` STRATA=1 -> `stratum_reversal_split_v2_consensus.tsv`:
  `spec_mech`, `spec_sib_mech`): sibling rankings barely change in most subgroups (Spearman 0.98-0.99 with `spec_sib`),
  so their sibling differences are not about proliferation. They change in HNSC (0.72-0.80): AMG-232 drops out of the
  atypical-classical top 5 (now vinca alkaloids and rigosertib lead), and basal-mesenchymal is led by topoisomerase I
  (topotecan, SN-38) and HSP90 inhibitors. Subgroup reports lead with `spec_sib_mech` / `spec_mech` (D80).

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
- The sibling-specificity formula (`spec_sib`, note section 4) was "to confirm with the user"; it is implemented as proposed.
- **Disease-side v2 rebuild (2026-10-18; drug side rerun on v2 2026-10-28 / 29, PRISM / GDSC validation on hold):** the Step 4 signatures, ARACNe networks and master regulators this branch scores against are being rebuilt with a TCGA-vs-GTEx lab correction, GTEx post-mortem covariates, the mappability/readthrough gene filter and ESTIMATE purity everywhere (spec `docs/SPEC_rebuild_v2.md` and decisions D45-D56 on `main`). All drug-side scores and the PRISM/GDSC validation (C16) are v1 and are to be rerun on v2 once its checks are reviewed. The GC tilt of v1 signatures is enough to reproduce group MR rankings (median rho 0.20, KIRC 0.66), which may contribute to the poor validation.
- **TCGA adjacent normal tissue in drug ranking (user 2026-10-18: the two most interesting uses; not designed yet).** On the disease side GTEx stays the reference (D57-D58 on `main`); adjacent normals can still serve the drug side:
  - *Safety filter.* The tumor's own surrounding tissue is what a drug also hits. Candidate components: (a) target window: a drug's targets higher in the tumor than in the same patients' adjacent normal tissue (paired where TCGA has pairs, ~10% of patients) and than in GTEx; (b) normal-tissue impact: penalize drugs whose LINCS effect suppresses the organ's normal-function genes (genes high in both GTEx and adjacent normal tissue of that organ); (c) lines resembling healthy tissue (`SPEC_celligner_normals.md`) as stand-ins for normal-tissue drug response. Open: thresholds, and whether it filters or down-weights.
  - *Inflammation around the tumor as a target.* Tissue next to tumors carries a shared pro-inflammatory state (Aran et al. 2017, Nat Commun 8:1077). Candidate: score drugs on reversing the inflammation part of adjacent-vs-GTEx (Hallmark inflammatory response, TNFA/NF-kB, IL6/JAK/STAT3, interferon; Aran's NAT genes), using only the most normal-looking adjacent samples (low same-lab tumor-likeness score, `scripts/33_normal_variance.R` on `main`). Main caveat: surgery itself induces stress / immediate-early genes (FOS, JUN, EGR1, ATF3) and the adjacent-vs-GTEx gap also holds lab, surgery-vs-death and donor-population differences; those genes must be separated (e.g. drop genes that respond to ischemia within GTEx, and immediate-early genes) before the inflammation part is a target.
- ~~Borrowed-drug confidence (C17): shrink, cutoff or column?~~ Decided 2026-10-28: flag at conf < 0.5 (C21). Still open: condition it on target expression.
- ~~How to fix the ER misreading~~ Addressed 2026-10-28: consensus regulons (C18, C19) and the line-match flag (C20). Original note: Splitting off proliferation did not (C15). Options: score drug-side TF activity with consensus / perturbation regulons (CollecTRI or DoRothEA are already in `data/raw/regulons`) as a sensitivity check next to the ARACNe regulons; or let drug-sensitivity data decide. PRISM / GDSC (C16) support none of the current scores, so the next step is choosing what to change (regulons, gene-level scoring, or using sensitivity data directly).

## 5. Likely questions and current answers

Questions a PI, reviewer or audience is likely to ask, with the answer as of 2026-10-12 (v2 changes: C18-C23). Update the answers as the
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
  best possible ~0.55 (C7). Weighting by line similarity added nothing. Per-drug confidence: `conf` in
  `drug_transfer_confidence.tsv` (C17).
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
- **Is it validated against drug-sensitivity data?** Tested on PRISM and GDSC: none of the scores predicts which
  drugs selectively kill a cancer's cell lines (median rho about 0), although that data clearly recovers known
  drugs (C16). This is the main open problem.
- **Does it recover known drugs?** Partly: 73rd percentile on average; erlotinib, 5-FU, flutamide high; afatinib,
  irinotecan not; anti-estrogens wrongly high in basal breast cancer (C2). Expect this to be the main question.
- **Isn't it just finding anti-proliferative drugs?** Largely yes for the reversal score: 44 growth inhibitors are
  top 10% in >= 20 cancers. Specificity removes them but does not recover known drugs better (C4). The score is
  now also split into a proliferation part and a mechanism part (C15): proliferation carries most of the
  known-drug signal; EGFR and androgen-receptor drugs are carried by the mechanism part.
- **Why do anti-estrogens rank high in basal breast cancer?** The tumor-derived ESR1 regulon is a luminal-identity
  module that lacks the acute estrogen-response genes, so fulvestrant reads as raising the ER program (C3). It is
  not a proliferation effect (C15).

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
