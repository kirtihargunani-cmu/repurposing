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
   `unmatched_lines` (C20). If a drug only exists as a cream, eye drop, inhaler or other local form, it cannot reach a
   tumor inside the body at a useful dose, so it gets the flag `local_route` ("Local use only"; D107). The one
   exception is a bladder rinse (intravesical) for bladder cancer. Flagged drugs stay in the table but are never put
   in the top tiers.

   Three more checks (D109) send a drug to the weak-evidence pile: its effect was seen on only one lab plate and a
   second plate testing the same doses saw nothing (a sign of a plate mix-up, as with metronidazole); its score was
   borrowed from other cancers although it was tested on this cancer's own cells and did nothing there; or the lab
   dose that worked is more than 10 times the highest level the drug reaches in patients' blood (taken from its FDA
   label or FDA review). Drugs that simply shut down the whole cell (protein synthesis, the sodium pump) are marked
   'generic mechanism' but kept. (Fix D111, 2026-11-07: the blood-level check had skipped 24 drugs whose FDA values
   were found but not converted to the lab's units, among them bortezomib and everolimus; they are now checked, and
   the check covers every lab version of a drug, not only one. Bortezomib leaves two top-10 lists.)
   (Fix D112, 2026-11-07: some FDA drugs were missed because their FDA name carries an unusual salt word, for
   example "lapatinib ditosylate", or a different spelling, for example "tioguanine" vs "thioguanine". Fifteen lab
   compounds now match, among them lapatinib, a standard HER2 breast cancer drug, which now ranks 6th of 20 FDA
   drugs in HER2 breast cancer. Cyclosporine, an immune-suppressing transplant drug, also appears and lands in ten
   top-10 lists; whether lab doses are reachable in patients has not been checked for it yet.)
   (Check D114, 2026-11-08: blood levels were then read from the FDA labels of cyclosporine, lapatinib, rucaparib,
   griseofulvin and pentamidine. For all five, the lab doses that worked are within reach of patients' blood
   levels, so none is marked. Thioguanine and clomiphene have no blood level in their FDA documents and stay
   unchecked.)
   (Fix D115, 2026-11-08: some drugs were tested as more than one lab sample. The FDA list used the best sample
   that passed all checks, but the web pages sometimes showed a different sample of the same drug that had failed
   one, so a drug could look weaker on the page than in the list. Both now use the same sample; the lists did
   not change.)
   (Fix D116, 2026-11-08: a drug borrowed from other cancers is now marked weak when any lab sample of it was tested
   on this cancer's own cells and did nothing, not only the same sample. In lung adenocarcinoma this moves
   temsirolimus out of the top 10 and osimertinib, a standard EGFR lung cancer drug, in.)
   (Housekeeping D117, 2026-11-08: every step that applies a check now runs from one script in a fixed order, and
   rerunning it gives exactly the same lists. One display count, how many tumors carry a missing target, could
   come out differently from run to run; that is fixed. The lists did not change.)
   (Fixes D119, 2026-11-08: seven small scoring and bookkeeping errors from the audit were corrected. The scores
   barely moved, eight top-10 lists swapped one drug, and the checks against known uses and against measured cell
   killing came out the same as before: the fixes made the numbers correct, not better.)
   (New score D120, 2026-11-08, shown beside the list, not used to rank: instead of asking whether the drug's
   whole effect on the cell undoes the cancer's pattern, it asks whether the drug's own target matters in this
   cancer, in three ways that are shown separately: do this cancer's cells need it, is it mutated or amplified in
   the patients' tumors, and is it turned up compared with normal tissue. This simpler question does pick out each
   cancer's approved drugs (about 62-67 times in 100 against 50 by chance, where the main list is at about 47-49),
   and the breast cancer checks come out right: anti-estrogens rank high in hormone-driven breast cancer and low in
   basal. Limits: about half of it comes from some drugs scoring well in every cancer; drugs with many targets get
   an unfair boost; it misses chemotherapy aimed at things every cell needs; and much of its top 10 rests on
   mutations that only a minority of patients carry.)
   (Follow-up D121, 2026-11-08: the target score was rebuilt in three layers (does the target matter, how does
   the drug act on it, how is that combined with the main list's score) and we let the data choose how much each
   part counts. The weights it learned from one list of approved uses did not hold up on the other list, so all
   three parts keep equal weight. Adding the main list's score did not help, and neither did how tightly the drug
   binds its targets in lab tests; how much of the target the drug covers at patients' blood levels helped a little
   where that is known.)
   (Decision D121a: binding strength stays in the score anyway. A drug that only weakly grabs a target should not
   get full credit for that target, even if the check against approved uses cannot see the difference.)
   (Switch D122, 2026-11-08: **the lists are now ordered by the target score**, not by how well a drug's lab
   signature undoes the cancer's pattern. The old score is still shown beside it, and the lab-data warnings (cells
   not well matched, effect not reproduced, drug did nothing in this cancer's own cells) are shown as flags instead
   of removing drugs. Drugs that do not act on a protein, such as many classic chemotherapies, are listed in a separate
   table. On the check against approved uses the new order does clearly better (about 62-63 times in 100 against
   about 48-49 for the old order on the same drugs), and anti-hormone drugs come out high in hormone-driven breast
   cancer and low in basal breast cancer. Almost every top 10 changed. Drugs that share a target now sit together,
   for example steroid hormones, cholesterol drugs (statins) and diabetes drugs (glitazones), which deserves a close
   look before anything is acted on.)
   (D123: how much of each drug circulates unbound in the blood was read from the FDA labels for 196 drugs. Using it to
   estimate how much of the target the drug occupies in patients did not improve the check, so it is reported only.)
   (D124: drugs aimed at things every dividing cell needs, like classic chemotherapy, are now scored by how fast the
   cancer divides, weighted by how closely the drug's target follows cell division in patients' tumors. Chemotherapy now
   rises in fast-growing cancers, cholesterol drugs (statins) and cyclosporine no longer reach any top 10, and the check
   against approved uses improved to about 65-67 in 100.)
12. **Why did trials fail?** For each cancer's top 10 FDA drugs we look up the trials of that drug in that cancer,
    read why stopped trials stopped and what the published results of finished trials said (with a quote), and
    count a failure only when the trial actually tested the drug (D108).
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

### What changed on 2026-11-03: separate networks for different cell types

**The problem we traced.** The breast "map" of which genes each switch controls was learned from all breast tumors
mixed together. In that mix, ER-positive (luminal) and ER-negative (basal) tumors are different kinds of cells, so the
map learned "ESR1 goes with being luminal" rather than "ESR1 switches these genes on". Anti-estrogens make MCF7 cells
look a little more luminal after a day, which that map read as ESR1 switching on. The ranking then did the opposite
of what doctors do: anti-estrogens near the top for basal breast cancer and near the bottom for luminal.

**What we tried, in order.** Removing "identity" links from the map, learning maps from drug experiments, relearning
link directions, using early (6 h) drug profiles, splitting the score into identity and mechanism parts, and blending
several maps together (metaVIPER). Each helped a little or only on one side.

**What fixed most of it.** Building a separate map for each cell type that has enough tumors (85 or more): luminal A,
luminal B and basal breast cancer each get their own, and so do the clearly different subtypes of cervical,
esophageal, glioma, pheochromocytoma, sarcoma and endometrial cancer. Smaller subtypes keep the shared map.

| Anti-estrogens (rank percentile) | Before | After |
|---|---|---|
| Luminal breast cancer (should be high) | about 26th | about 41st |
| Basal breast cancer (should be low) | about 97th | about 77th |

Checks against known drug uses and cell-line gene dependencies stayed the same or improved slightly. The cost: about
12 more drugs (of 1,144) now look like they help almost every cancer, mostly broad growth blockers such as mTOR
inhibitors. Our main score already pushes those to the middle of each list, and a check (D98) confirmed that where they reach the
top tiers they are already clearly specific to that cancer (on average 36 percentile points above their usual rank), so
no extra filter was added. This change was adopted by the
user's decision, knowing that trade-off (D97).

**Fixing the ESR1 reading alone is not enough (D101).** An independent pathway map (PROGENy) reads anti-estrogens
correctly, so we tried using it for ESR1 and AR on the drug side. The rankings barely moved (under one percentile
point): ESR1 is one of about 1,000 switches in each cancer's score, and what holds anti-estrogens back is the whole
luminal program being read as rising, not ESR1 itself.

The rankings page now shows a **Timing** label for each drug where 6 h data exist: works early and late, early only
(a direct effect the cell adapts to), late only (a slower or secondary effect), or neither.

**Still open.** Basal breast cancer still ranks anti-estrogens fairly high, because the drugs switch on a "slow down"
program that basal tumors have switched off. We also ranked drugs separately on their early (6 h) and late (24 h)
responses (D100), because the two are different programs: the early one reflects what the drug does directly, the late
one also includes how the cell adapts. Ranking on the early response finds fewer broad growth blockers but recovers
fewer drugs already used for each cancer (0.61 -> 0.51), so the main ranking keeps using all time points. Each drug
now carries an early score, a late score and a label (works early and late, early only, late only, neither, or no
early data), which shows whether its effect is direct or develops later.

### How much did the cell-type maps gain? The evidence in numbers

**In short:** the cell-type maps fix most of the breast anti-estrogen inversion, which is solidly supported. Outside
breast we can only show they do no harm. Before = pooled maps (results commit `3320b5f`); after = production since D97.

**Evidence 1: anti-estrogens now rank more like clinical practice.** Median rank of the 6 anti-estrogens
(percentile within each subtype's list; higher = recommended more strongly):

| Breast subtype | Should be | Before | After |
|---|---|---|---|
| Luminal A (ER+) | high | 15th | **38th** |
| Luminal B (ER+) | high | 30th | **43rd** |
| Basal / triple-negative (ER-) | low | 95th | **70th** |
| Her2 (map unchanged, too few tumors) | middle | 56th | 56th |

Before, anti-estrogens ranked 80 points higher in basal (where they do not work) than in luminal A (where they are
standard care). Now the gap is 32 points: still the wrong way round, but less than half as large. Examples:
fulvestrant luminal A 14.5 -> 37.2 and basal 98.7 -> 76.5; tamoxifen luminal A 49 -> 74; toremifene basal 93 -> 68.
Not every drug improved: tamoxifen in luminal B dipped slightly (62 -> 57).

**Why: the basal disease profile is no longer defined by "not luminal".** The drug side did not get better at reading
anti-estrogens; our maps still read them as switching ESR1 on (D99). The gain came from the disease side. On the old
mixed map, basal tumors were described largely as "the luminal regulators are switched off" (ESR1 -2.5, GATA3 -2.4,
XBP1 -2.6). So any drug that nudged luminal genes up looked like it was fixing basal cancer. With a basal-only map,
that luminal axis largely disappears (ESR1 -0.7, GATA3 +0.2, XBP1 -0.2), and the raw reversal score of anti-estrogens
in basal falls by a third (30.3 -> 20.1).

**Evidence 2: approved drugs stand out a little more in their own cancer.** Comparing each approved drug's rank in the
cancer it is approved for with its rank elsewhere, the median advantage rose from +2.9 to +6.6 percentile points (Hub
labels); on ChEMBL labels it moved from -5.0 to -4.0. Overall known-use recovery stayed about the same (AUC Hub 0.620
-> 0.615, ChEMBL 0.500 -> 0.514; 0.5 is chance). Agreement with cell-line gene dependencies moved from -0.038 to
-0.029: slightly better, still essentially zero.

**What it cost, and what is uncertain.**

- Twelve more broad growth blockers look good for almost every cancer (151 -> 163 of 1,144 drugs). They stay
  mid-ranked on the main score (D98), so they do not flood the top of the lists.
- The new maps are noisier. Built twice from independent halves of their runs, they share 0.40-0.57 of their links,
  against 0.61 for the old combined breast map, because each has fewer tumors.
- Outside breast there is no ground truth for a benefit. In the other 7 cancers the rankings moved (only 10 of the top
  20 stayed for leiomyosarcoma, 12-13 for esophageal cancer) and known-use recovery did not drop: no harm, but no shown
  gain either.
- Basal breast cancer still ranks anti-estrogens too high (70th). Part of that is a "slow down" program the drugs switch
  on, which basal tumors have switched off. Fixing the ESR1 reading alone did not help (D101), because the whole
  luminal program is still misread on the drug side.

**Bottom line.** The cell-type maps fixed the clearest known error, the breast inversion, by about 60%, and made
approved drugs stand out somewhat more in their own cancer. The cost is slightly noisier maps and about 12 more broad
growth blockers. The evidence for a gain is concrete only in breast; elsewhere the change is neutral as far as we can
measure.

### Is the old reversal score limited by noisy lab data? No (2026-11-09)

Before trying a new version of the reversal score, we measured how repeatable the lab data behind it are. We scored
every drug twice, each time from a different experiment on the same drug and cell line (another dose or time), and
compared the two rankings. For drugs measured in a matching cell line the two rankings agree with a correlation of
about 0.6. We then ran a stricter version that splits each drug's experiments into two halves and compares the two
rankings, repeated twenty times: across cancers the typical agreement is about 0.87. The gap between the two methods
shows that most of what the first one counted as noise is really the drug behaving differently at different doses and
times. At 0.87 the data would allow a score to put a cancer's approved drug above another cancer drug about 97 times
in 100, and even at the pessimistic 0.6 about 89 times in 100. The reversal score does it about 43-51 times in 100,
which is chance. So the problem is the score, not the noise in the data. The current lead ranking (target role) does not use these lab data,
so this limit does not apply to it.

### Why "is it already approved?" is a weak way to check the rankings (2026-11-09)

Until now the main check was whether approved drugs rank highly in the cancers they are approved for. That test has a
built-in flaw: a drug can go unapproved for a cancer for reasons that have nothing to do with whether it works, such as
patents, the cost of a trial, too few patients to make one worthwhile, or simply that nobody tried. Those drugs count
as failures in the test even when they might work, which is exactly the kind of drug this project exists to find. The
test also rewards being a cancer drug at all: simply labelling every cancer drug as a hit scores higher than our
ranking does.

So we are changing the main check. The new one compares a drug against itself in two groups of patients where biology
says the answer: BRAF drugs should rank higher in BRAF-driven thyroid cancer than in the RAS-driven kind, EGFR drugs
higher in lung adenocarcinoma than in squamous lung cancer, HER2 drugs higher in HER2-positive breast cancer than in
the other subtypes. Because it is the same drug on both sides, patents, company decisions and market size cancel out.
Approval data is still reported, but it no longer decides anything. The honest limit of the new check is that it
confirms the rankings carry biology we already know; it cannot show they find biology nobody knows yet.

**First result (2026-11-09): the rankings pass, and the new check immediately found something.** Across 49 drug
comparisons the ranking put the drug higher in the group where biology says it should be 36 times, about 73%, which is
unlikely to happen by chance. The old reversal score managed 61%, no better than a coin toss. One comparison failed
badly: BRAF drugs in thyroid cancer, where the ranking gives almost exactly the same score to both patient groups even
though BRAF is altered in 87% of one group and only 16% of the other. The reason is that the score ranks a drug's
target against all other genes in that cancer rather than measuring how much the target differs between the two groups
of patients, so a target that comes top in both groups looks identical in both. The approval-based check could never
have shown this. How to fix it is not yet decided.

### What it gets wrong, honestly

- **Anti-estrogens in luminal A breast cancer ranked low** because the tumor-learned map misread them (C18). Removing
  "identity" links helped only a little (D83, C25); separate cell-type maps (D97, above) lift them to about the 41st
  percentile, which is better but still not near the top.
- **Anti-estrogens rank at the top in basal (triple-negative) breast cancer, where they should not work.** They were
  only tested in MCF7, a luminal line that looks nothing like basal tumors. The `unmatched_lines` flag catches this and
  keeps them out of the top tiers (C20). Tiers alone do not move the rank, though: with the cell-type maps they still
  sit at about the 70th percentile. **Test on 2026-11-04 (D102):** a drug gets no credit in a cancer whose tumors
  barely make its target (here ESR1: 0.86 units in basal tumors vs 108 in luminal A). In breast this sends
  fulvestrant and raloxifene from the 76th to the 3rd percentile in basal and nudges them up in luminal, without
  hurting the known-use check. It is not adopted yet: it also pushed down two drugs that have been tried in breast
  cancer trials (gefitinib, metoprolol, whose targets are nearly absent in luminal tumors), which the rule set in
  advance did not allow, and tamoxifen escapes because it also binds a second protein that basal tumors do make.
  A revised version (D102a) judges each drug on its known mechanism (so tamoxifen is judged on the estrogen
  receptor), compares against a typical high level instead of the single highest cancer, and only blocks on
  approved uses. It passes: in basal all four anti-estrogens drop to about the 3rd percentile, luminal controls rise
  a little, and the known-use checks hold across all 41 cancers. One catch before switching it on: drugs for targets
  that only a small share of patients carry, such as ALK inhibitors in lung cancer (ALK fusions in a few percent of
  tumors), are pushed down too, because the typical tumor has almost none of the target.
  We tried letting a target count if even a small share of tumors make a lot of it (D102b). It backfired: 6 of 190
  basal tumors make a lot of the estrogen receptor (3%), more than the share of lung tumors in our data with high ALK
  (3 of 513), so the rule brought the anti-estrogens back in basal and still did not rescue the ALK drugs. Expression
  alone cannot tell a rare but targetable mutation from leftover expression.
  **Decision (D102c, now in use):** each subtype's list describes its typical tumor, since that is the point of
  splitting patients into subtypes. So the revised rule is switched on: drugs whose target the typical tumor lacks
  go to the bottom of that list and are marked "Target absent". Where a handful of tumors do make the target (like
  ALK in lung cancer), the page adds "Minority target" so a biomarker-selected use is not lost from view.
  **Fix (D110, 2026-11-07):** an audit found the ALK drugs were not actually being pushed down. ChEMBL lists their
  target as the fusion protein (for example EML4-ALK), and the code counted both halves, so EML4, a gene almost every
  cell makes, made the target look present everywhere. The drug only binds the ALK half, so now only ALK counts.
  Alectinib, ceritinib and crizotinib now drop to the bottom wherever the typical tumor lacks ALK (67 disease and 20
  subgroup rows). The same applies to imatinib-type drugs (only ABL1 counts), which changes nothing because ABL1 is
  made almost everywhere.
- **The cell lines must have the receptor too (D104c, 2026-11-04).** A drug's score comes from cell lines. If the
  drug works through a receptor and none of those lines makes it (nor does the drug switch it on), the experiment
  could not show the drug's real effect, so the drug moves near the bottom and is marked "Receptor absent in tested
  lines". Exception: some receptors sit on blood vessels or immune cells inside the tumor, not on the cancer cells,
  and cancer lines never have them (VEGF receptors, targets of kidney-cancer drugs like tivozanib). Those are kept
  and marked "Acts on tumor environment". About 400 drug-cancer pairs were moved down; the known-use check rose
  slightly (0.617 -> 0.623). A first version compared tumor and line levels directly and failed: one prostate cell
  line makes 20 times more androgen receptor than a typical tumor, which pushed the prostate drugs to the bottom.
- **Each FDA drug is judged on its real main target (D105b, 2026-11-04).** We now take the target from ChEMBL's
  hand-curated record of how the drug works (found for 291 of 413 FDA drugs), or else from the drug's listed target,
  checked against lab binding data where possible. Lab binding data never picks the target on its own: when we tried
  that, irinotecan's "target" came out as an enzyme it only hits in a side assay. Drugs that act on DNA or on germs
  are not checked against human proteins. This gave the biggest improvement so far in finding drugs already used for
  each cancer (0.623 -> 0.660). A separate page now shows only drugs whose main target is actually made by the
  cancer's tumors.
  (Fix D113, 2026-11-08: 32 FDA drugs, among them vincristine, prednisone and thioguanine, had been linked to a
  duplicate or wrong ChEMBL record that has no curated mechanism, so they were judged on their listed targets
  instead. They are now looked up by their FDA name; 325 of 422 FDA drugs now use the curated mechanism. The
  top-10 lists did not change; meloxicam, whose curated target is COX-2 only, now drops where tumors lack COX-2.)
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

- **Tests against measured drug killing of cell lines (PRISM / GDSC), known drug uses and gene dependencies** have
  been run (C28, C29) and do not support the scores. **Known uses (redone carefully as D118, 2026-11-08):** if we
  only compare cancer drugs with each other, the drugs already approved for a cancer do not rank higher for that
  cancer than other cancer drugs do; they rank slightly lower (0.46 and 0.43 on the two label sources, where 0.5 is a
  coin toss). The earlier "better than chance" (about 0.63) mostly reflected that cancer drugs rank above non-cancer
  drugs, and that the target checks push some drugs down; simply knowing whether a drug is a cancer drug scores 0.66.
  The scores also do not predict which cell lines a drug kills. A published check of this whole family of methods (Koudijs 2023) found the same: reversal
  scores mostly pick up general anti-growth effects. These remain ranked hypotheses, not validated predictions.
- **A signalling map was downloaded to check drug mechanisms (D126, step 1 of 4, 2026-11-10).** The plan: when a
  drug moves a set of regulators, ask whether a signal could actually travel from the drug's own target to those
  regulators, in the right direction. That needs a map of which protein switches which other protein on or off. We
  pulled one (OmniPath), kept only the links where the curators agree on both the direction and the on/off sign and
  where at least two sources or two papers back it, and saved our own copy with a fingerprint, so the result cannot
  change under us: **16,518 links between 5,493 proteins.**
  Two things worth knowing before the check is built. First, what shrank the map was not disagreement about
  direction or sign but thin evidence: about two thirds of the starting links rest on a single paper from a single
  database. Second, the check may not reach far enough -- only 46% of the drug targets we score have any outgoing
  link in the map, and a drug whose target is a dead end cannot be checked at all. The rule written in advance
  requires the check to cover at least 40% of the regulators for at least 60% of drugs, and that is now in doubt.
  Nothing is flagged yet; this step only built and pinned the map.
- **A smarter way of comparing the two maps was tried and did not work (D125, 2026-11-10).** When we compare a
  tumor's regulator pattern with a drug's, we currently treat every regulator as a separate piece of evidence. That
  is not true: many regulators control overlapping sets of genes, so a block of twenty that largely share targets is
  counted twenty times. The fix was to build a map of which regulators overlap, and down-weight evidence that is
  shared. It is a one-knob change: turn the knob all the way down and you get exactly today's score back.
  We tested the knob on a dataset where the answer is known - thousands of experiments where a single regulator was
  switched off or forced on, so the right reading is not in doubt. **Every setting read the experiments either no
  better or worse than today's score**, and the more weight we gave the overlap map, the worse it got: the cleanest
  experiments (forcing a regulator on) lost up to 10 points. The rule written before the test required a gain of at
  least 1 point to go on; the best was 0.17. So the idea was stopped there and nothing in the ranking changed.
  Worth knowing: the overlap is genuinely there (the typical regulator shares targets with about fifty others), so
  the answer is not "nothing to fix" - it is that the overlap carries real signal, not just repetition. This is the
  third different attempt to strip out shared structure (the others used tumor expression patterns) and the third to
  remove signal along with it.
- **An independent check with a different kind of map (PROGENy pathways, D99).** Good news: the drug data are sound.
  Drugs that block a pathway (PI3K, MAPK, EGFR, p53 activators) move that pathway the right way in 75% of profiles.
  Revealing news: PROGENy reads every anti-estrogen profile as switching estrogen signalling off, while all our
  ARACNe breast maps, even the new luminal ones, still read ESR1 as on. The anti-estrogen ranking improved because the
  disease side changed, not because the ESR1 reading was fixed. Our ranking and a pathway-based ranking agree only
  weakly (0.17), so they capture different things.
- Scoring on same-cell-type lines only (instead of same-organ) is not implemented. The data-gap check reports it but
  does not change rankings.
- The odd recurring hits above still need to be reviewed.
