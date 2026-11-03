# Drug side in plain language, with results

Written 2026-11-02 from the v2 results. This is the non-technical version of
[`PIPELINE_DRUG_SIDE_V2.md`](PIPELINE_DRUG_SIDE_V2.md), which has the commands, files and formulas. It starts where
[`DISEASE_SIDE_PLAIN_LANGUAGE.md`](DISEASE_SIDE_PLAIN_LANGUAGE.md) ends. Decision and caveat numbers (D.., C..) point
to [`DECISIONS.md`](DECISIONS.md) and [`notes/drug_side_decisions_and_caveats.md`](notes/drug_side_decisions_and_caveats.md).

## What we are trying to do

The disease side measured which control switches (transcription factors) are stuck on or off in each of 41 cancers
and 14 subgroups.

Drugs also flip these switches. A public database (LINCS) measured what about 3,000 drugs do to gene activity in
cancer cell lines grown in the lab. From that we can read which switches each drug flips.

A drug is a candidate for a cancer if it flips that cancer's switches **back toward healthy**. We then ask three
follow-up questions:

1. Does it do this **for this cancer in particular**, or for every cancer? (A drug that "fixes" everything is usually
   just a general poison.)
2. Does it work through something other than **just stopping cells from dividing**? (Almost any toxic drug slows
   division, so that alone is not interesting.)
3. **How trustworthy is the lab data** behind it? Was the drug tested on cells that actually resemble this cancer?

## How it works, step by step

1. **Pick lab cell lines that resemble each cancer.** A tool called Celligner compares patient tumors with cell lines and
   says which lines look most like them. We only use lines from the same tissue, because the question is whether a
   drug's effect will carry over to this disease.
2. **Clean up the drug data.** We keep good-quality, clearly active drug measurements. Each cell line's "average drug
   response" is removed first. Otherwise nearly every drug looked like it reversed every cancer, because nearly every
   drug slows growth.
3. **Borrow carefully.** If a drug was never tested on this cancer's tissue, we use its average effect in at least 3
   other cancer lines. We also estimate how well that borrowed effect will carry over. This mostly depends on how
   consistent the drug is across cell lines, not on how far away the other tissues are (C21). Borrowed drugs with
   confidence below 0.5 are flagged (`low_transfer`; about 11.5% of borrowed pairs).
4. **Check that the drug's target is there.** A drug stays on a cancer's list only if the protein it is designed to hit
   is actually made in that cancer's tumors (in at least 20% of them). A related protein from the same family does
   not count.
5. **Read the drug's effect on each switch.** We use two "maps" of which genes each switch controls. The first is
   learned from the tumors themselves (ARACNe). The second is a curated map of direct targets (CollecTRI). We average
   the two. The tumor-learned map mixes in "cell identity" genes. For example, it reads anti-estrogen breast cancer
   drugs as *switching estrogen on* (C18, C25). The curated map reads them correctly, so we use both.
6. **Score each drug against each cancer.**
   - **Reversal (`rev`)**: how strongly the drug pushes the cancer's switches back toward healthy. Strongly changed
     switches count more.
   - **Specificity (`spec`)**: reversal for this cancer minus the drug's average reversal for cancers from other
     tissues.
   - **Mechanism specificity (`spec_mech`)**: specificity after removing the part that comes only from stopping cell
     division. **This is the score we lead with (D80).**
   - For subgroups, `spec_sib_mech` compares the subgroup with its sibling (for example MSI-like vs CIN-like colon
     cancer).
7. **Flag weak evidence.** If a drug's best test line looks nothing like the cancer, it gets the flag
   `unmatched_lines` (C20). Flagged drugs stay in the table but are never put in the top tiers.
8. **Label each cancer's evidence.** Each cancer is labeled `matched` (good candidates were tested in similar lines),
   `borrowed only`, or `no adequate compound`. Reports must state this label next to any list.
9. **Check which receptors the drug really acts through.** A drug's annotated receptor counts for a cancer only if a
   drug database (ChEMBL) shows it binds that receptor strongly enough (1 µM; 100 nM for kinase drugs), and the
   receptor is actually made in the lab cell lines the cancer's score came from. Otherwise the score cannot be
   coming through that receptor. We also report how much more strongly the drug binds its main target (C26).
10. **Find stand-ins for drugs never tested on this cancer's cells.** Most FDA drugs were only tested on other cancers'
    cells. For these we look for another compound that was tested on this cancer's cells, hits the same target and
    changes the cell's switches the same way. Its score is shown next to the drug's own as a second opinion. Checked
    first on drugs that were tested here: a stand-in's rank follows the drug's real rank fairly well (C27). Only
    cancers with their own lab data (about 20 of 55 groups) can have stand-ins.
11. **Link to the outside world.** We match drugs to FDA approvals, look up PubMed papers and clinical trials for the
    top drugs, and build two browsable pages: the FDA top lists (`results/lincs/fda_rankings_v2.html`), every FDA
    drug with all its scores, flags, receptors, stand-ins and literature (`results/lincs/disease_drug_rankings_v2.html`),
    and the compounds tested on each cancer's own cells with their closest FDA drugs (`results/lincs/standin_compounds_v2.html`).

## Results

### Does it find drugs we already know work? Mostly yes.

These are percentiles: 99 means better than 99% of drugs for that cancer.

| Cancer | Known drug | Reversal | Specificity | Mechanism specificity | Read |
|---|---|---|---|---|---|
| Lung adenocarcinoma | erlotinib (EGFR) | 96 | 96 | 94 | found; osimertinib is #1 overall |
| Lung adenocarcinoma | gefitinib (EGFR) | 87 | 88 | 83 | found |
| Prostate | flutamide (anti-androgen) | 85 | 84 | 87 | found |
| Colon | irinotecan | 63 | 98.5 | 97 | found by specificity |
| Colon | fluorouracil | 81 | 57 | 59 | middling |
| Breast luminal B | tamoxifen | 95 | 73 | 62 | found by reversal, weaker on specificity |
| Breast luminal A | tamoxifen / fulvestrant | 85 / 72 | 46 / 18 | 49 / 15 | **missed** (see below) |

**A worked example.** Erlotinib in lung adenocarcinoma pushes 647 of 1,061 switches the right way. Its score is not
carried by one or two genes. Palbociclib in luminal B breast cancer works mostly through cell-division switches (E2F7,
E2F8, MYBL1). That matches how the drug is known to act, and it is why its mechanism score is lower than its reversal
score. Table: `results/lincs/worked_example_regulators.tsv`.

### What it gets wrong, honestly

- **Anti-estrogens in luminal A breast cancer rank low.** The tumor-learned map misreads them (C18). Removing
  "identity" links from that map helped only a little (luminal A 9.7 -> 12.3, luminal B 23.4 -> 26.0 percentile; D83,
  C25).
- **Anti-estrogens rank at the top in basal (triple-negative) breast cancer, where they should not work.** They were
  only tested in MCF7, a luminal line that looks nothing like basal tumors. The `unmatched_lines` flag catches this and
  keeps them out of the top tiers (C20).
- **Some drug classes still show up everywhere.** Anthracyclines (doxorubicin, daunorubicin, epirubicin) and
  topotecan top the FDA lists of many cancers. Odd hits such as metronidazole, pentobarbital and cardiac glycosides
  need checking.
- **Before the specificity step, 32 drugs** (statins, palbociclib, MEK inhibitors, anthracyclines) were in the top 10%
  for 20 or more cancers. With mechanism specificity, none are.

### How good is the data behind each cancer?

- **14 of 41 cancers are `matched`**: breast HER2 / luminal A / luminal B, cervical adenocarcinoma, esophageal and
  stomach adenocarcinoma, neuroblastoma, colon, rectum, liver, lung adenocarcinoma, pancreas, prostate and melanoma.
- **27 are `borrowed only`.** Every candidate for these rests on data from other tissues.
- **Stricter check (same cell type, at least 10 drugs; C24): 14 of 41 cancers and 2 of 14 subgroups pass.** The ones
  that fail:
  - **No lab line of that cell type exists (4):** pheochromocytoma / paraganglioma, chromophobe kidney, seminoma.
  - **Lines exist but have no good-quality LINCS drug data (12):** for example sarcomas, kidney clear-cell and papillary,
    thyroid, adrenal, Wilms, IDH-mutant glioma.
  - **Data only in lines that don't resemble the tumors (8):** basal breast, all squamous cancers (only one lung
    squamous line, 3 drugs), glioblastoma, serous endometrial.
  - **Too few drugs tested (3):** bladder (1 drug), ovarian (3), endometrioid endometrial (7).

For these cancers the honest report is "no good-quality data at this step", not a ranked list.

### Top FDA-approved candidates (by mechanism specificity)

The top 3 for some `matched` cancers:

| Cancer | Top 3 |
|---|---|
| Breast HER2 | temsirolimus, osimertinib, gedatolisib |
| Lung adenocarcinoma | tofacitinib, temsirolimus, dasatinib |
| Prostate | panobinostat, azacitidine, belinostat |
| Liver | digitoxin, niclosamide, omacetaxine |
| Neuroblastoma | floxuridine, trifluridine, vinblastine |
| Melanoma | clofarabine, artesunate, ciclopirox |
| Breast luminal A / B | epirubicin, daunorubicin, idarubicin (the recurring anthracyclines) |

Subgroups, against their sibling:

| Subgroup | Top 3 | Read |
|---|---|---|
| Colon MSI-like | trametinib, binimetinib, mirdametinib (all MEK inhibitors) | fits the BRAF / MAPK biology of MSI colon cancer |
| Colon CIN-like | cytarabine, gemcitabine, cladribine | DNA-synthesis drugs |
| Thyroid BRAF-like | panobinostat, romidepsin (HDAC inhibitors), teniposide | HDAC inhibitors are in thyroid trials |
| Thyroid RAS-like | ciclopirox, trifluridine, pitavastatin | |
| HNSC basal-mesenchymal | topotecan, azathioprine, regorafenib | |

The LGG and seminoma subgroups barely differ in drug response, so report one list per parent cancer. All rankings,
with paper and trial links, are on the FDA rankings page.

## What is not done yet

- **No test yet against measured drug killing of cell lines** (PRISM / GDSC). This is on hold. Until then these are
  ranked hypotheses, not validated predictions.
- Scoring on same-cell-type lines only (instead of same-organ) is not implemented. The data-gap check reports it but
  does not change rankings.
- The odd recurring hits above still need to be reviewed.
