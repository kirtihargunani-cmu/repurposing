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

## How the decisions developed, D81 to D128 (added 2026-11-15)

One story in four phases: the project kept trying to fix one score, gradually realised the problem was the question
that score asks, and changed the question.

**In one sentence:** the pipeline started by ranking drugs on whether they push a tumor's switch activity back toward
normal ("reversal"). The anti-estrogen test case kept showing that reversal misreads the biology. Patching the
networks, timing and profiles did not fix it, so the project (1) tightened which cancers a drug is even eligible for,
(2) switched the lead ranking to "does this cancer need the drug's target?", and (3) is now testing whether reversal
can be rebuilt in a narrower, more honest form (D128).

### Phase 1: can reversal be made to read the right biology? (D81-D101)

**The problem everything revolved around.** Anti-estrogens shut down estrogen signalling in MCF7 breast cancer cells.
That is unambiguous, and it is the easiest test a drug-ranking method could face. Yet the tumor-trained network
(ARACNe) read ESR1, the estrogen receptor, as going *up* under these drugs. If the method cannot get this case right,
it cannot be trusted on harder ones.

| What was tried | What happened, plainly |
|:--|:--|
| Narrower tumor maps; drop links that just mark the subtype (D81-D83) | Barely helped. Within luminal tumors, "being luminal" and "ESR1 active" travel together, so deleting obvious markers cannot untangle them. |
| Pick the better network (D84) | Neither ARACNe nor CollecTRI won clearly. The blend of the two was kept. |
| Learn network links from drug experiments instead (D86-D87) | No better than the existing networks, sometimes worse. Fixing ESR1 tended to break other things. |
| Raw vs "generic response removed" profiles (D85) | Removing each cell line's generic response is essential. Without it, 82% of drugs "reverse" 80% of cancers, because general growth arrest looks like reversal everywhere. Raw profiles helped ESR1 a little but not overall. |
| Early timepoints (D88-D90, D100) | **The key insight of this phase:** at 6 hours, estrogen genes have already dropped. By 24 hours the cells shift toward a luminal identity, and *that* is what the network misread as ESR1 rising. Timing explained the misread better than dose. But most anti-estrogen 6-hour profiles were too weak to use, so early and late became reporting columns only. |
| Separate "identity" from "reversal" (D91, D98) | Helped luminal controls, but anti-estrogens still looked great in Basal breast cancer, which has no estrogen receptor. A "cells going quiet" program made them look good there. |
| Separate tumor maps per cancer type (D92-D97) | The biggest win of the phase. The Basal anti-estrogen false positive fell from the 97th to the 76th percentile, and real breast positive controls rose from 25 to 41. Adopted by explicit decision despite slightly missing one tolerance. It improved the *ranking* more than it fixed ESR1 itself. |
| Independent pathway readout, PROGENy (D99, D101) | PROGENy got it right: estrogen signalling down in 15 of 15 profiles. Swapping it in barely moved the ranks, because the misleading luminal signal is spread across many switches. |

**What Phase 1 taught:** the misread was not one bad link or one bad setting. It is built into how a tumor-trained
network sees a cell changing identity. Many patches each helped a little; none fixed it.

### Phase 2: stop the obvious mistakes, whatever the score says (D102-D119)

Since the score could not be made to read the biology, the project put guardrails around it: rules for when a drug
is even eligible in a cancer.

- **The target must be present in the typical tumor (D102).** This alone dropped the Basal anti-estrogen false
  positive from about 76 to about 3. The trade-off: rankings describe the *typical* tumor, and rare carriers get a
  flag, not a rescue.
- **The target must be present in the test cells (D104).** Done carefully: a too-strict comparison unfairly punished
  flutamide, and anti-angiogenic drugs whose targets sit in blood vessels the cell assay does not have. Exceptions
  were added for the tumor environment.
- **Use each drug's real main target (D105, D110).** The old "strongest binder" rule assigned wrong targets, and ALK
  inhibitors got false credit through fusion partners.
- **Clean up the evidence (D107-D114).** Creams and other local-only drugs were removed from systemic rankings; a
  dose rule that flagged 45% of profiles was dropped; profiles that do not reproduce are caught; exposure, name and
  salt-form problems were fixed.
- **Make everything consistent and reproducible (D103, D106, D115-D119).** The FDA lists now reproduce byte for
  byte. The D119 corrections barely changed the rankings, and that mattered: it showed **the reversal failure was not
  a bug**.

### Phase 3: change what is measured, and how it is judged (D118-D127)

This is the turning point.

1. **The yardstick was flattering reversal (D118).** Once the "known use" benchmark was cleaned up, reversal scored
   0.46 and 0.43, where 0.50 is a coin flip. The earlier ~0.63 was an artefact.
2. **Ask a different question: does this cancer need the drug's target? (D120)** Is the target essential, altered, or
   overactive? This scored 0.67 and 0.62 on the same rows and recovered the right patterns, such as ESR1 in luminal
   breast and BRAF in melanoma.
3. **Make it fair across drugs (D121).** Families and complexes count as one target, so drugs with long target lists
   do not get more chances. Equal weights were kept, and **reversal got zero weight**.
4. **Promote it (D122).** The layered score became the lead ranking for all 41 cancer groups and 14 subgroups.
   Weaknesses in the LINCS evidence became flags rather than vetoes.
5. **Handle chemo-like drugs (D124).** Weighting proliferation by how closely a target tracks cell division lifted
   TOP2A-type drugs and dropped statins from the top 10s. This was confirmed on the data it was developed on, so it is
   a seen-data gain.
6. **Change the judge (D127).** Approval history mixes biology with commercial history, so the test became "does the
   same drug rank higher in the cohort that carries its biomarker?" Production passed 36 of 49 comparisons
   (p = 0.0007). The clear miss was **thyroid**: BRAF alterations were far more common in one group (87% vs 16%), but
   converting to a within-cohort percentile erased that difference.

**What Phase 3 taught:** "is the target relevant?" works. "Does the drug reverse the tumor?" does not, at least not
as currently measured.

### Phase 4: is reversal dead, or just measured badly? (D125-D128)

- **Ceiling test (D125 A0).** The measurements are repeatable enough that a good score *could* reach an accuracy
  (AUC) of about 0.97. Reversal sits near chance anyway, so **the problem is the scoring, not noisy data**. This is
  the most important finding of this phase.
- **Discount switches that overlap (D125 A1).** Did not work: +0.17 points against the +1.0 threshold. The shared
  targets were real signal.
- **Trace the path from a drug's target (D126).** Failed its own control: 0 of 22 anti-estrogens flagged. It looked
  only at each drug's 25 biggest switch changes, and ESR1 is never among them; growth-arrest switches are.
- **Score only ESR1's own target genes (2026-11-12).** The anti-estrogens came out *backwards*: worst in ER-positive
  cancers, best in Basal. The likely reason is that **normal breast tissue already has estrogen activity**, so
  comparing a tumor with normal breast points the wrong way for estrogen blocking. Not yet tested with another
  reference.
- **D128 (2026-11-15, planned, nothing run).** Built directly on those lessons:
  - start from the drug's annotated target, not its biggest changes (fixes D126);
  - keep **M** ("is the drug doing what its target predicts?") separate from **R** ("does it push disease genes
    back?");
  - write down in advance that R is expected to fail in breast because of the reference problem.

  Plan: [`IMPLEMENTATION_target_mechanism_scoring.md`](IMPLEMENTATION_target_mechanism_scoring.md); issues found
  before it starts: [`notes/d128_plan_review.md`](notes/d128_plan_review.md).

  **First step done (2026-11-15): the stock-take, before any scoring.** Nothing was scored and no answer was looked
  at. The question was only: for which cancers and drugs could this be worked out at all? Of about 12,300
  cancer-and-drug pairs, roughly 5,200 can be scored on at least one of the two new measures; every one of the rest
  now carries a stated reason, most often that no route could be traced from the drug's target to any switch.

  The encouraging part is the anti-estrogen test case. Starting from ESR1 itself, we reach **294 of its own genes**,
  where ten would have been enough, so the footprint that the earlier path check (D126) could never find is clearly
  there. The catch is that only fulvestrant counts for now: the other anti-estrogens are recorded as "modulators",
  which does not say whether they raise or lower the receptor's activity, so they are set aside until that is decided.

  Three limits were confirmed rather than assumed. There is **no lung cell line carrying the EGFR mutation** that the
  matching drugs were made for, so that comparison cannot be done honestly and is marked as such. The HER2 comparison
  survives on one cell line and twelve experiments. And **PRISM**, the drug-sensitivity data we had hoped to validate
  against, turns out to have been used repeatedly in earlier decisions, so it can no longer serve as an independent
  check.

  Two problems with our own rules surfaced while there were still no results to bias them, and both are now fixed.
  One was a cut-off that mislabelled a strongly HER2-amplified cell line as ordinary; correcting it is what made the
  HER2 comparison possible at all. The other was that a drug such as palbociclib was being judged on too narrow a
  slice of its evidence: because one of its two targets happens to be a switch, the whole drug was judged on that
  switch alone and failed, even though its other target reaches nearly 1,800 genes. Each target is now judged on its
  own terms and the drug counts if either one suffices, which admitted 13 more targets and unblocked the CDK4/6
  comparison without taking anything away. A first attempt at the run also had to be abandoned after 43 minutes for being too slow; caching the
  repeated graph searches brought it down to about two minutes, and both the failure and the fix are in the record.

  **The first real test is in (2026-11-15).** The scoring machinery now exists, and the anti-estrogen case that the
  whole project keeps returning to was the first thing run through it. Starting from ESR1 itself, fulvestrant's
  measured response does line up with the direction its target predicts — in all three independent experiments, and
  more than it does for 200 scrambled versions of the same footprint. The earlier path-tracing attempt found nothing
  at all here, so this is real progress on the question that mattered.

  Two honest qualifications. The agreement is **weak**: about 4% of the strongest possible match, averaged over 478
  genes. It is consistent and clearly better than chance, but small, and whether it adds anything to the ranking we
  already have is a later question we have not answered. And the drug-versus-disease half of the new scoring behaved
  in breast exactly as we predicted in writing beforehand: uninformative in the oestrogen-driven cancers, and
  highest in the one breast cancer that lacks the receptor entirely. That is the normal-tissue comparison problem
  resurfacing, not a reason to give the drug to those patients, and we changed nothing to make it look better.

  **The whole thing was then run, on every drug and cancer it could be run on (2026-11-16), and the honest summary
  is that it does not yet work well enough to be trusted.** Four results matter, and none of them was tuned.

  *First, half the measurement may be an artefact of how the data are made.* The machine that measures drug
  responses only measures about a thousand genes directly; the other nine thousand are **guessed by a model** from
  those thousand. We asked, for the first time, whether the new score survives if we keep only the genes that were
  really measured. It does not. **One drug-target pair in four flips sign** — the drug goes from "doing what its
  target predicts" to "doing the opposite" — depending only on whether the guessed genes are included. The question
  and the pass mark were written down before we looked. The pass mark was missed on all three parts, and badly. We
  changed nothing: both versions are now reported side by side, and the choice of which to trust is a decision for
  the whole project, not just for this score, because every score in the pipeline is built on the same guessed
  genes.

  *Second, the one control that works is measuring the wrong thing.* Melanoma drugs that target BRAF were the
  richest test case, and they behave as they should: the score is positive in the two cell lines that carry the BRAF
  mutation, and the same drugs in cells **without** it never reach the confidence bar. That is a real pass. But when
  we opened up where the score comes from, **65% of the genes behind it are controlled by MYC**, a master switch for
  cell growth that almost every effective drug disturbs, and those genes carry 69-95% of the actual number. So the
  score is mostly reading "this drug slowed the cells down", dressed up as "this drug hit BRAF". We are reporting
  that as a finding rather than patching it.

  *Third, the score cannot tell the two halves of a biomarker pair apart.* The project's main fairness test takes
  one drug and asks whether it ranks higher in the cancer that should respond than in the cancer that should not.
  The existing ranking gets that right about three times in four. The new score gets it right once in three — and
  the reason is built into it: it contains **nothing about the disease at all**, so when a drug's cell-line evidence
  is the same for both cancers, the two answers are **identical by construction**. That happened in 15 of the 49
  comparisons, including every thyroid one. On the rest it was right 47% of the time, which is a coin toss.

  *Fourth, a check that needs no outside data, and the score only just survives it.* If a drug's target is **not
  switched on at all** in the cells the measurement came from, the drug cannot have hit it there, so the score
  should be flat. It is still positive 57% of the time, against 60% when the target is present. The typical value is
  eight times smaller when the target is absent, which is the right direction, but the two cases overlap so heavily
  that no single result can be told apart.

  Two further things worth knowing. The score is **smaller than the disagreement between the cell lines it is
  averaged over** — the spread across lines is larger than the typical score — so a cohort-level number is thinner
  than it looks. And the "is this specific to this target?" check we built turns out to have almost nothing to
  compare against: **95% of targets have fewer than ten genuinely different alternatives** in our list, and a tenth
  have none, so that particular reassurance is largely unavailable. For the oestrogen-receptor case there was
  exactly one alternative (the androgen receptor), which is why that check came back blank.

  **A rule the program was not obeying, fixed (2026-11-16).** Before any of this is taken further, one bookkeeping
  error was corrected. There is a written minimum for how much evidence a target needs before we are willing to
  score it at all: ten genes reached, and - for targets reached indirectly, through a signalling chain - at least
  two switches at the end of that chain. The program was only checking the ten-gene half. One target (the DDR2
  protein, reached by the drug regorafenib) had 173 genes but only **one** switch behind all of them, so it was
  being scored when the rule says it should not be. The rule was right and the program was wrong, so the program
  was changed, not the rule. Re-running everything moved exactly what was predicted in writing beforehand and
  nothing else: those rows now abstain with a stated reason, regorafenib's own score changes in 55 cancer-drug
  rows, a hundred-odd other drugs shift by at most one percentile place because the ranking renumbers, and every
  control - the anti-estrogen test, the BRAF melanoma test and its MYC finding - comes back byte for byte identical.

  **Why one pair in four flipped, worked out (2026-11-16).** The sign-flip finding above needed an explanation
  before anyone decides which version of the gene list to trust, so we tested four possible causes, all written
  down in advance. Two survive. First, flipping happens almost entirely where the score is close to zero: below
  0.005 it flips half the time, above 0.05 it essentially stops - but **only about one row in eight is ever above
  0.05**, so for most of the list the sign is not something to read. Second, and more important: we counted how
  many genuinely *independent* measurements a list of genes actually contains, and the answer is about **45,
  whatever the length of the list**. A 478-gene footprint carries 45; a 1,900-gene one still carries 45; the entire
  ten-thousand-gene space carries 41. The genes move together, so adding more of them adds almost no new
  information - it only dilutes the average and, worse, makes our chance comparisons look more impressive than
  they are, because those comparisons shuffle genes as though each were a separate experiment. **That is a caution
  about every score in this project that counts genes, not just this one.** The two causes we did NOT find: it is
  not the guessed genes (a random slice of the same size is just as unstable as the measured slice - so switching
  to measured-only would not fix it), and it is not the MYC hub problem (once you compare like with like, the hub
  makes no difference). Nothing was changed in response; the two things the evidence asks for - a minimum size
  below which a sign is not reported, and a fairer chance comparison - are decisions for the user.

  **Can the "pushes it back toward normal" score even be read? Now computed, and it fails one of its own tests
  (2026-11-16).** The second half of the new scoring compares the tumor with healthy tissue and asks whether the
  drug moves it back. That question only makes sense if the healthy tissue is not already running the pathway the
  drug blocks - the breast/oestrogen problem this project keeps hitting. Until now that caution was a handwritten
  note about one drug. It is now a number computed for **every target and every cancer**: we feed the *disease*
  comparison through exactly the same machinery that scores the *drug*, which tells us whether the tumor already
  looks as though the target were switched off. Three verdicts come out - "do not read this one", "safe to read",
  and "no signal either way" - plus an honest "cannot tell" where the data is missing. **Three checks of four came
  out well.** The verdicts are spread sensibly (no single answer swallows more than 80% of cases). Melanoma versus
  normal skin reads "safe to read" for BRAF drugs, and HER2 in HER2 breast cancer reads "no signal", both of which
  are what you would expect for faults the tumor acquires rather than features of healthy tissue. And where the flag
  says "do not read", the relationship between the two halves of the score genuinely turns upside down (-0.46
  against +0.39 where it says "safe"). **The fourth check failed.** We had written down beforehand that the
  oestrogen case in luminal-B breast cancer must *not* come back "safe to read", because we know from experience
  that it is not. It came back "safe to read". We have not touched the definition or the threshold to make that go
  away. The failure also exposed a real design flaw: the answer currently depends on which drug is asking - the
  blocker and the hormone itself get opposite verdicts from the same data - when it should be a property of the
  tissue alone. That fix is written down for the user to approve rather than applied now, because changing a
  definition after seeing it fail is how you fool yourself. Two further honest notes: the flag **abstains on every
  SERM** (tamoxifen and its relatives), which is precisely the class where the problem was first seen; and some
  cancers simply sit further from healthy tissue than others - basal breast cancer's typical drug scores three
  times higher than luminal-A's whatever the drug does - so every score now carries its cohort's baseline beside it.

  What this does not say: the new score is not wrong, and it is not finished. It is near-independent of the ranking
  we already use, which is what we wanted — but "different" is not "better", and nothing here can tell the two apart
  because no outcome data was looked at. That test exists, and it is deliberately **not** run yet: it needs its own
  written-in-advance design, and until that is written no drug-sensitivity data is opened.

### Threads that run through the whole history

1. **The anti-estrogen case is the project's canary.** Almost every phase was triggered by it, and D128's first gate
   is still that test.
2. **Patches that help the ranking do not necessarily fix the mechanism.** The cell-type maps (D92-D97) and the
   PROGENy swap (D99, D101) both showed this.
3. **When the benchmark cannot separate cause from coincidence, change the benchmark:** known use, then cleaned known
   use, then biomarker pairs.
4. **Decisions are rule-first, and failures get recorded.** Several approaches (early-only timepoints, D125 A1, D126)
   were closed by their own pre-set rules rather than rescued.

## What is not done yet

- **Tests against measured drug killing of cell lines (PRISM / GDSC), known drug uses and gene dependencies** have
  been run (C28, C29) and do not support the scores. **Known uses (redone carefully as D118, 2026-11-08):** if we
  only compare cancer drugs with each other, the drugs already approved for a cancer do not rank higher for that
  cancer than other cancer drugs do; they rank slightly lower (0.46 and 0.43 on the two label sources, where 0.5 is a
  coin toss). The earlier "better than chance" (about 0.63) mostly reflected that cancer drugs rank above non-cancer
  drugs, and that the target checks push some drugs down; simply knowing whether a drug is a cancer drug scores 0.66.
  The scores also do not predict which cell lines a drug kills. A published check of this whole family of methods (Koudijs 2023) found the same: reversal
  scores mostly pick up general anti-growth effects. These remain ranked hypotheses, not validated predictions.
- **The mechanism check was built, tested and dropped (D126, 2026-11-10).** Using the signalling map below, we
  asked of every drug: could a signal actually travel from the drug's own target to the regulators it moved, in the
  right direction? The test we had written in advance was whether it catches a failure we already know about —
  anti-estrogens in luminal breast cancer, where our maps read the drug as switching ESR1 *on* when it switches it
  off. **It did not catch it. It did the opposite**, calling tamoxifen and toremifene fully consistent there. The
  reason is worth keeping: the check looks at the 25 regulators each drug moves most, and those are almost always
  the growth regulators, not the drug's own target. ESR1 never appeared in the list at all, so the question being
  answered was "can this drug explain slowing growth?" — to which the honest answer is yes, and which tells us
  nothing about mechanism. The same growth-signal dominance we have been fighting in the ranking turned up inside
  the tool built to audit it.
  It also could not reach far enough: on average the map connected a drug's target to only 2 of the 25 regulators,
  against the 10 the rule required. Nothing was added to the rankings, and the idea is closed. Reviving it would
  mean looking at regulators near the drug's target instead of the ones it moves hardest — a different test, which
  would need its own advance rule.
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

### Two measuring errors found on 2026-11-16, and what they do and do not break

Both come from the same mistake, made in two different places: treating things that move together as if they
moved independently. Genes do not move independently — switch on one master regulator and hundreds of its
targets move together. Regulators do not move independently either, because one regulator often controls
another, so their target lists overlap. **That shared movement is real biology, not noise, and neither fix
below removes any of it.** What was wrong was only the yardstick used to say how surprising a number is.

**Error one: our "is this better than chance?" tests were too easy to pass (D130).**
To ask whether a drug's effect on its target's genes is real, we compare it against shuffled versions of itself.
The shuffling treated each gene as its own independent piece of evidence. A 478-gene list behaves more like
about 46 independent pieces, so the shuffled comparisons were far too similar to each other and almost anything
looked unusual beside them. We rebuilt the comparison so the shuffles keep the real pattern of which genes move
together. The headline casualty is our one positive control: **fulvestrant, an anti-oestrogen, in a breast
cancer cell line. It cleared the old test comfortably and does not clear the corrected one.** The effect itself
did not change at all — it was always small — but we can no longer say it is more than chance would give. The
same correction removes almost all the other "confident" signals in that family of tests. Nothing in the drug
rankings was built on those tests, so no ranking moved; what changed is that a result we had been treating as
evidence is now not evidence.

We also checked the corrected test on a thousand randomly chosen drug-target-cell-line rows. Under a test that
is working properly, about 5% of them should look "unusual" by chance. Under the old test **18%** did. Under the
corrected tests, 7% and 8% — close to where they should be. So this was not a case of making everything look
less impressive: the old yardstick was simply about three and a half times too easy to pass, and the new one is
about right.

**Error two: the reversal score is on a scale about ten times too big (D131).**
The reversal score adds up how much a drug pushes back on each of a disease's master regulators, then divides by
a figure that would be correct only if the regulators were independent. They are not, so the division is too
small and the score comes out roughly **ten times** larger than the scale it claims to be on — and the amount of
over-statement is different for every disease, ranging from about four times to about fifteen. **The good news,
and it is the important part: this does not move any ranking.** The dividing figure is one number per disease,
the same for every drug in that disease, so it cancels out of every ranking, percentile and tier we publish. We
checked all 55 diseases and subgroups and the answer was identical in every one. What it does break is comparing
a raw reversal number between two diseases, or reading one as though it meant something on its own. We have not
changed the score; this is written down for the user to decide on.

One related number does move. When we split a drug's effect into "general slowing of growth" and "specific
mechanism", the split compares two raw totals without allowing for the two sets being different sizes and
differently interlinked. Comparing them fairly shifts about **27%** of drug-disease rows across the halfway line.
That is a genuine judgement call about what the split should mean, not a bug, and nothing was changed.

**Was the "ten mechanisms" finding just an artefact of our regulator lists? No (D133).**
The previous audit found that when we add up a thousand regulator readings for a drug, the total behaves as
though it were built from only about **ten** independent pieces rather than a thousand. The obvious objection,
and the right one to raise, is that this may say nothing about biology: we read each regulator by looking at the
genes it controls, and if two regulators control mostly the same genes then of course they give the same answer.
That would mean our instrument cannot tell a thousand regulators apart — a worse problem than the one we thought
we had. We tested it on five diseases chosen to span the whole range of the effect.

The answer is no, and it is not close. Regulators that share genes *do* give near-identical readings — the
effect is real where it happens — but **over half of all regulator pairs share no gene at all**, and those pairs
are just as correlated as everything else (typical strength 0.21–0.26, against 0.25–0.29 for all pairs). When we
subtract out everything that gene-sharing can predict, the "about ten pieces" figure moves only to **about
twelve or fourteen** — nowhere near the ~550 we would expect if the readings were independent. If we go further
and pick out a set of regulators that share *literally no gene* with each other, 40 to 53 of them, their
readings still behave like about 11 to 15 independent pieces instead of 40 to 50.

We also checked a second worry the user raised: our main regulator lists are *learned from which genes move
together in tumours*, so they could have inherited the very correlation we are measuring. We re-ran everything
on a hand-curated list built from published literature instead, which cannot have inherited it. It collapses too
— about two to three times less than the learned list, so some of the worry is justified, but it still collapses
roughly twelve-fold. And the telling detail runs the opposite way to the gene-sharing story: the curated lists
**share more genes** with each other, yet their readings are **less** correlated.

Finally, regulators that are documented as regulating each other are more correlated than regulators in general
— but no more correlated than other pairs that share the same number of genes. So a documented link buys us
nothing beyond gene-sharing, and we are not claiming to have seen a signalling cascade at work.

**What we still cannot rule out.** Every regulator reading is taken from the same ranked list of genes for a
given drug. If a drug's effect on the genome has one very dominant direction, every regulator will pick that up,
whatever genes it controls. So the honest statement is not "there are only ten mechanisms" and not "our
instrument is broken", but: **drug responses share a small number of dominant directions, and our regulator
readings all reflect them.** Nothing was changed, re-ranked or re-scored.

### Can we even tell what a drug targets, just from what it does to the genes? (first a head-count, 2026-11-17)

There is a question we should have asked at the very beginning: **if you only saw what a drug did to a cell's
genes, could you work out which protein it was built to block?** If the answer is no, then no amount of clever
scoring on top of those gene readings can be right about mechanism, because the information simply is not
there. We have written down how to test this properly. Before running it, we did a head-count — just counting
what data exists, with no test performed and no result produced. The point of doing it this way round is that
the counting cannot be bent to flatter a result, because there is no result yet.

**The idea of the test.** Two different drugs that block the same protein should leave similar fingerprints on
the genes. So, for each protein target, take the drugs that hit it and ask whether their fingerprints really do
look alike. That needs at least two different drugs per target. On paper we have 193 targets, of which 82 have
two or more drugs — and 111 have only one, which can never be tested at all.

**What the head-count found.**

- **The fingerprints have to come from the same cell line to be comparable.** Comparing a drug tested in a
  breast cell to one tested in a lung cell mostly measures the cell, not the drug. Insisting both drugs were
  tested in the *same* cell takes us from **82 testable targets down to 71**, and removes **39% of the drug
  pairs**. The targets we lose are small and mostly non-cancer ones — blood-pressure channels, a few serotonin
  receptors, the clotting target.
- **The planned check for "these are really the same molecule in disguise" barely does anything, and that is
  worth knowing.** Three drugs that are near-copies of each other are not three independent pieces of evidence.
  Our plan was to spot near-copies using a standard chemical fingerprint. It turns out that fingerprint only
  catches *the same molecule* in a different salt or mirror-image form — it does not catch, say, two different
  statins or two different kinase drugs, which is what we actually meant by near-copies. It merged only four
  small groups (the vitamin-A drugs, two steroids, two steroid creams, two contraceptives) and changed nothing.
  **So that safeguard, as written, is not the safeguard we thought it was**, and we have flagged that it needs
  a real chemical-similarity measure before the test runs rather than after.
- **A drug that did nothing cannot be recognised.** We measured how strong each drug's gene fingerprint is. Our
  data already throws away weak responses — three-quarters of the raw measurements were dropped long ago by an
  existing quality rule — and among what remains, a typical fingerprint moves about 450 genes noticeably. We
  deliberately did **not** pick a strength cut-off here; we only costed several. A gentle one is nearly free
  (71 targets become 67). A strict one halves the set (down to 35–49). **One uncomfortable detail:** the
  targets with the *most* drugs — the steroid receptor, two brain-signalling receptors — have some of the
  *weakest* fingerprints, so a strict cut-off would delete exactly the targets the whole exercise leans on.
- **How often would you be right by luck?** Usually almost never: in a typical cell line you are choosing among
  dozens to hundreds of drugs, so luck gives you about a 0.3% chance. But **sixteen cases are hopeless**: cell
  lines where only two or three drugs were ever tested, so "finding the right one" is arithmetically guaranteed
  and means nothing. Those are now flagged.

**The verdict: the test is worth running.** After all of the above we are left with **71 targets, 229 drugs and
37 cell lines**, with a typical target carrying 4 independent drugs across 6 cell lines, and 29 targets carrying
5 or more drugs. That is a genuinely useful map of where mechanism-based scoring could work. It is smaller than
the plan assumed, and the plan's own coverage figures were about a third too optimistic, which we have
corrected. And it still says nothing whatsoever about the 111 targets with only one drug — those stay untested,
which is not the same as untestable. Nothing was scored, re-ranked or changed.

### Does a drug's fingerprint look like switching off its own target? A second head-count (2026-11-17)

The previous head-count asked whether two drugs that hit the same protein look alike. This one asks something
more direct, and arguably more basic: **if a drug is built to block protein X, does the mark it leaves on a
cell's genes resemble what happens when you simply delete the gene for X?** That comparison needs no network,
no regulon, no second drug — just two measured fingerprints side by side. CMap profiled thousands of gene
knockouts (CRISPR) and knockdowns (shRNA) in the same way they profiled drugs, so in principle we can look.

Again we counted first and tested nothing. **No drug was compared with any knockout in this step.** What
follows is a stock-take of what could be compared, and of whether the knockouts themselves are trustworthy.

**The plan's own coverage table was wrong in three ways.** It claimed our drug data spans 87 cell lines; it
spans **73**. Its three rows (CRISPR, consensus shRNA, raw shRNA) were each counted with a different quality
filter, so they are not comparable with one another, and one of the three cannot be reproduced at all under
any of sixteen readings we tried. We report that as unreproduced rather than quietly substituting our own
number. This is the second plan in a row whose feasibility table did not survive being recomputed.

**The real scope is much smaller than "208 targets".** Coverage counts how many targets were knocked out
somewhere. What the test actually needs is a *cell*: a target, a cell line, a knockout of that target in that
line, **and** a drug for that target profiled in the same line. Counting those properly: 377 target genes →
280 have a CRISPR knockout → 279 of those in a line we also have drugs for (only **8** such lines) → **2,350
usable (drug, target, line) combinations**, covering 201 drugs and 255 targets. Five cell lines carry almost
all of it.

**The knockouts are far weaker than the drug profiles.** Our drug side keeps only strong signatures (a quality
score of 0.4 or better). The knockouts mostly cannot clear that bar: only **one CRISPR knockout in ten** does,
and only one consensus-shRNA knockdown in twenty-five. A typical CRISPR knockout moves about 51 genes
appreciably; a typical consensus shRNA knockdown moves 14. That matters, because the arm the plan nominated as
its *safest* source of knockdowns — the consensus shRNA, designed to cancel out shRNA's notorious off-target
effects — turns out to be the weakest-responding arm by a factor of six.

**The most important measurement: did the knockouts actually work?** The obvious check is to look at the
target's own gene: after knocking it out, its mRNA should go down. Headline pass rates are poor — 22% for
CRISPR, 19% for consensus shRNA, 28% for raw shRNA. **But those headline numbers are an illusion created by
how L1000 works.** The platform physically measures about 1,000 genes and *predicts* the rest from them. On
the genes it really measures, the check passes **70% of the time for CRISPR** and 55% for shRNA. On the
predicted genes, the target's own value sits at zero no matter what happened, and the pass rate falls to 3–9%
— which is no information at all, not a failed experiment.

So the honest reading has two halves, and they are not the same. For **shRNA**, which works by destroying the
transcript, reduced mRNA is the direct expected readout, so a failure on a measured gene is genuine evidence
that the knockdown did not take. For **CRISPR**, which cuts the gene rather than the message, a broken gene
can still produce a near-normal amount of mRNA while the protein is gone — so a CRISPR mRNA failure must never
be counted as a failed knockout. We wrote that caution down before looking, and it still stands; what we did
not expect is that it is not the limiting problem here, because CRISPR actually passes *more* often than
shRNA on measured genes. The limiting problem is that **only 48 of our 377 targets are genes L1000 measures at
all.** For the rest, this quality check simply does not exist and must not be applied.

**How much of what remains could give an interpretable answer?** The plan lists ten reasons a negative result
would be excusable — the target is essential to survival, the target is not switched on in that cell line, the
drug hits many proteins, the knockout is too weak, and so on. If nearly everything falls into one of those
buckets, the test can only ever confirm, never refute. Counting them: of 2,350 CRISPR combinations, 20% are
essential targets, 37% have a target that is not expressed in that line, and **91% involve a drug with more
than one annotated target**. Stacking the filters leaves **46 combinations** that are clean on every count we
can compute (121 if we allow drugs with two targets, 80 if we drop the strength requirements). One of the ten
reasons — that deleting a protein also removes its scaffolding role, which a drug would leave intact — cannot
be checked at all, because no list of such proteins exists in our data. We report it missing rather than
invent one, which means even the 46 is an optimistic count.

**What those 46 actually are is the encouraging part.** They are the textbook pairs: ABL1 with nilotinib and
ponatinib, BRAF with vemurafenib and dabrafenib, EGFR with erlotinib, gefitinib and osimertinib, TYMS with
fluorouracil, ESR1 with fulvestrant, plus PPARG, VDR, DPP4 and a handful more. If the comparison works
anywhere, it should work there.

**One flaw in the plan's own safety check.** It proposed three "this must obviously work" controls: HDAC
inhibitors, proteasome inhibitors and topoisomerase poisons. All three exist in our data — but TOP2A is
*always* essential for survival and the proteasome subunits nearly always are, which is exactly the
circumstance the plan itself marks as unreadable, because deleting an essential gene mainly tells you the cell
died. If those controls pass, it may be because dying cells look like dying cells, not because the method
recognises mechanism. HDAC is the only one of the three that is mostly clean.

**The verdict: worth continuing, with a smaller claim.** The CRISPR arm is the one with real scope, and
46–121 clean combinations across 14–25 targets is enough for a careful target-by-target map of a dozen famous
drug–target pairs. It is nowhere near enough to make a sweeping statement about all 2,350 combinations, nine
in ten of which are excusable in advance. And the direction of the logic is fixed from the start: this test
can show that target signal **is** there; it can never show that it is absent. Nothing was scored, changed or
re-ranked.

---

## Why we stopped subtracting all of "proliferation" (2026-11-20)

Nearly every cancer drug slows growth. That is a problem for ranking, because if we just ask "which
drug most reverses the tumour's gene activity", the answer is always *whichever drug is most toxic* —
in every disease, for every patient group. So the pipeline estimates how much of each drug's effect
is simply growth arrest and subtracts it out.

The question the user asked was whether we were subtracting too much. We were — and in a way worth
spelling out, because the fix changes what the pipeline can find.

Think of a drug's effect on growth as three separate things stacked on top of each other:

1. **The background.** Cells in a dish slow down for all sorts of reasons that have nothing to do
   with the drug — the handling, the medium, the time elapsed. Every drug in the same plate shares it.
2. **The drug's general toxicity.** Some compounds slow *every* cell line they touch. These are the
   ones that look brilliant against every disease and treat none of them specifically.
3. **Selective growth arrest.** This drug stops *this* kind of cancer from growing, and leaves others
   alone.

Layer 1 is noise. Layer 2 is the trap. **Layer 3 is the whole point** — a drug that selectively
arrests one tumour type is exactly the kind of discovery this project exists to make.

It turns out an earlier step already removes layer 1, at the level of individual genes, before the
growth-arrest correction ever runs. We measured this: the growth signal going into the correction has
an average of about **−0.015**, essentially zero, and almost no relationship to the background
(correlation 0.03). The background was gone. So the correction, which was written to remove layer 1,
was in practice eating layers 2 **and** 3.

How much of layer 3 was being thrown away? About half — **48% early and 59% late** of the growth
signal is a drug arresting some cell lines and not others.

So there is now a switch with three settings: remove everything (what we did before, still the
default so no result has silently changed), remove only the drug's general toxicity and keep the
selective part, or remove nothing. Which one is right is not something to argue about — the plan is
to run the upcoming validation under all three and see which one best recovers drug targets we
already know about.

### A smaller thing that went right

Two weeks of results rested on a finding that growth arrest behaves differently at 6 hours than at
24 hours. When we loosened a data-quality filter and let in roughly twice as much early data, that
finding looked shakier. The worry was that we had been fooled by a small sample.

We checked by splitting the new data into the part that would have passed the old filter and the part
that only just scraped in. The old-filter part reproduced the original result almost exactly — a
correlation of **0.681 against the original 0.687**, and a disagreement rate of **25.7% against
25.2%**. The finding was real. The wobble came entirely from the weaker new measurements, and it
never reached the drug rankings themselves, which barely moved (37.0% to 37.4%).

More data, same answer, on a sample 69% larger. That is the outcome you want from this kind of check.

---

## Can the pipeline find a drug's target without being told it? (2026-11-20)

This was the big test, and it came back split — one clear yes and one clear no.

The pipeline has a step that ignores drug labels entirely. It takes a drug, looks at everything the
drug changes in a cancer cell, and asks: which of these changes most opposes what the tumour is
doing? Whatever comes out on top is the pipeline's guess at how the drug works. The guess is made
without ever being told what the drug binds.

To check whether the guess is any good, we ran it on drugs whose targets we already know, and asked
where the known target landed in a list of roughly 2,850 candidates.

**For a drug's main target, it works.** Among the top 1% of each list we found **59 known targets
where chance would give about 15** — nearly four times better than guessing. Some individual results
are striking: for erlotinib, a lung cancer drug, its actual target EGFR came back **3rd out of
2,806**. Romidepsin, panobinostat and belinostat all surfaced the HDAC enzymes they are designed to
block. Etoposide surfaced TOP2B. Fingolimod surfaced S1PR1.

The detail that makes this convincing is not any single hit — it's that **the method does better when
it has more data**. Drugs measured in one experiment land 1.1× better than chance; drugs measured in
five or more land **2.45×** better. A fluke does not improve when you feed it more evidence.

**For a drug's secondary targets, it does not work.** Many drugs hit more than one protein, and the
second or third one is often where a new use would come from. That is exactly what we hoped to
surface. Here the result is flat: **7 known secondary targets in the top 1%, where chance gives about
8**. Strength of binding did not help either — the secondary targets a drug grips most tightly ranked
no better than the ones it barely touches.

A few entries do look right — dasatinib finding EPHA2, midostaurin finding PRKCD — but with 793
pairs being scored, that is about how many good-looking results you would get by luck alone. So we
are not counting them.

**And the proliferation question was settled along the way.** The worry was that growth arrest, which
nearly every cancer drug causes, was drowning out the real secondary mechanisms. We tested that
directly by running the whole thing four ways: no removal, partial removal, full removal, and full
removal plus discarding every growth-related regulator. The secondary result stayed flat in all four.
Growth arrest was not what was hiding them.

One thing we did learn: the shallowest hits are mostly growth-related, while the deepest and best
ones are not. Of the results that made the top 10, **more than 80% were growth regulators**, against
a background rate of 13%. But the genuinely good recoveries — EGFR, the HDACs, the hormone receptors
— sit in the top 1% and have nothing to do with growth. The lesson is to trust the very top of the
list and be sceptical of the merely-near-the-top.

## Does our ranking predict whether a drug actually works in patients? (D176, 2026-11-21)

Short answer: **no — not on the evidence we have, and we tested it properly.**

Here is the test. For each cancer cohort we took the drugs we rank near the **top**, drugs near the
**middle**, and drugs near the **bottom** of our list. We pulled every clinical trial on record for
each of those drug-plus-disease combinations, read the published results, and asked two questions.

**Question 1: are top-ranked drugs trialled more often?** Yes, enormously. Top-ranked drug-disease
pairs have an average of **11.7 trials** on record; middle-ranked have 3.7; bottom-ranked have 1.7.
That is a 6.9× difference and it is not a fluke.

But this is **not** good news, and it is important to understand why. Our score is built from
curated drug-target lists, laboratory binding measurements, and mutation frequencies. All of those
sources are *much* richer for drugs that scientists have studied heavily. A famous drug gets a high
score partly *because* it is famous — and famous drugs are also the ones that get trialled. So the
gradient mostly measures attention, not biology. We report it, and we do not count it as a win.

**Question 2: among drugs that did get trialled, do our top-ranked ones succeed more often?**
**No.** The share of trials reporting a positive result is essentially identical in all three bands:

| where we ranked the drug | share of trials that worked |
|---|---|
| top | 48.4% |
| middle | 51.3% |
| bottom | 47.5% |

Not only is the difference not significant (p = 0.79) — it is not even in the right order. The
middle band is nominally highest.

**How do we know the reading was fair?** The top band was read in an earlier session and the other
two now, so if our standard for "positive" had drifted between sessions, that drift alone could have
produced — or erased — a difference. So we secretly mixed **60 already-read trials** back into the
new queue, unlabelled, and checked whether we gave them the same answer twice. We did: **88% exact
agreement, κ = 0.83**, with zero cases of a trial flipping between positive and negative. We also
required every single reading to quote the actual abstract word-for-word, and all 256 passed. The
thresholds for calling the check a failure were written down *before* we looked.

**What this means in plain terms.** We can say the pipeline is good at pointing to the disease where
a drug's target is actually switched on — that is what the biomarker test (D127) shows, 7 out of 7.
We **cannot** say it picks the drugs that will survive a clinical trial. Those are different claims,
and only the first one is supported.

Two honest caveats. The test is **small**: with 136 top-band and 86 bottom-band pairs, we could only
have detected a fairly large gap, so a modest real advantage would be invisible here. And it is a
**hard test by design**: a drug only reaches a trial after human experts already judged it promising,
so almost all the obviously-bad candidates have been filtered out before we start counting. Neither
caveat changes the conclusion; they describe how strong a statement the null supports.

## Why corticosteroids topped the lung-cancer list, and what we did (D177, 2026-11-21)

Steroids like dexamethasone were ranking at the very top for lung squamous cell carcinoma. They do
not treat that cancer — they are supportive care. So something was wrong, and two things were.

**Problem 1: we were rewarding drugs for targeting genes the tumour doesn't need.** Our score asks,
among other things, "does the tumour depend on this drug's target?" We measure that by deleting the
gene in cancer cells and seeing whether they die. For most drugs — which *block* their target — more
dependence is better. But steroids *activate* their target rather than blocking it, and for those
the code flipped the number upside down: it treated "the tumour does **not** need this gene" as
strong evidence. NR3C1, the steroid receptor, sits in the **bottom 8%** of dependency in lung
squamous cancer — deleting it slightly *helps* the cells grow — and the score read that as **92 out
of 100**. We checked this exactly: 751 of 5,098 rows were being flipped this way. Deletion
experiments simply cannot tell you what happens when you *turn a gene up*, so we stopped using them
for activating drugs.

**Problem 2: "lower in tumour than in normal tissue" was being treated as a finding.** For genes we
hadn't properly measured, the code substituted a cruder comparison — is this gene lower in the
tumour than in healthy tissue? For a brain receptor in a lung tumour the answer is automatically yes,
because brain receptors aren't in lungs. That is not a discovery. **182 rows** were resting on it,
and it is the same thing that made antidepressants show up across unrelated cancers.

Both are now fixed.

**The part we cannot fix by changing the score.** The steroid receptor really *is* switched down in
these tumours — that measurement is solid. Our whole method says: find what the tumour has pushed out
of balance, and find a drug that pushes it back. By that logic a steroid is a reasonable suggestion.
The flaw is in the logic itself: some things are switched down *because* a tumour has lost its normal
identity, not because switching them down is what makes it a tumour. Pushing those back doesn't cure
anything. Our data cannot tell those two situations apart.

**So we label instead of hiding.** Every drug now carries two new pieces of information: how many
independent lines of evidence support it (out of three), and whether it works by blocking or
activating its target. An activating drug can only ever have **one** line of evidence, because the
other two are impossible to collect for it. The rankings page says so in plain words next to the
drug. The rule is: **a high rank is a position in a list, not a measure of how much evidence there
is — always read the two together.**

**What improved and what didn't.** The top four drugs for lung squamous cancer are now
fluorouracil/floxuridine, alpelisib and palbociclib — all backed by three lines of evidence, and
alpelisib and palbociclib target genes that are genuinely broken in this cancer. **The steroids are
still there, at positions 5 to 10.** Across all 41 cancers, this family of drugs actually holds
*slightly more* top-10 spots than before (22.4%, up from 21.7%). We are reporting that rather than
claiming a win.

**And it cost us something.** Our main accuracy check — does the score correctly rank a drug higher
in the cancer that carries its biomarker? — fell from 36 out of 49 to **33 out of 49**, which is
below the 70% bar we set in advance. We did not adjust the bar. Looking closely, the whole drop comes
from a single comparison in thyroid cancer that was **already failing before this change** and that
our method cannot actually resolve: two of the drugs score 100.0 in *both* groups, an exact tie. So
the fall says more about that test than about the change. But the bar was set beforehand and it is
not met, and that is how we are recording it.


## Seeing what a drug's score is made of (2026-11-22, D179)

A drug's score for a cancer is a sum over hundreds of switches (regulators), and each switch's reading is itself a sum
over the genes it controls. We used to keep only the final number. Now we keep every layer: for a drug in a cancer, which
switches it pushes the right way and which the wrong way, and under each switch, which genes moved. The layers add up
exactly to the score, so this is a true breakdown, not an estimate.

Example: palbociclib in lung squamous cancer scores well because it shuts down the E2F cell-division switches, through
genes such as FOXM1, MYBL2 and MKI67 that are strongly up in the tumor. That is exactly how palbociclib is meant to
work, and it is also the general "stop dividing" signal that many cancer drugs share. The breakdown makes that
visible for any drug. Across drugs, scores are usually spread thin: the ten biggest switches typically carry under a
tenth of the total.

## We now use one set of gene networks, and that fixed three hidden problems (D208, 2026-11-23)

The pipeline used to be able to run on a blend of two kinds of gene networks: ones learned from
each cancer's own tumors (ARACNe) and a general, literature-curated set (CollecTRI). We had already
decided to use only the tumor-learned ones. Removing the blended version completely showed that
some steps were still quietly using it:

- **Our main validation check was reading the old version.** Re-run on the current pipeline, the
  biomarker test (do drugs score higher in the cancer that has their target biomarker?) now
  **passes**: 36 of 49 comparisons go the right way (73%), against a 70% bar.
- **Part of the main score still mixed in the curated networks.** It no longer does. Most rankings
  barely move, but about 900 drug-disease scores shift by more than 10 percentile points.
- **About a fifth of drug-disease pairs were missing their warning flags** (for example, "the dose
  needed in the lab is far above what patients can reach"). Every pair now has them.

We also tried adding new biomarker comparisons. The test included a deliberate "should fail" check,
and it passed anyway. Looking closer, several drugs in that check hit the same target, so they get
the same score and count as one piece of evidence rather than five. So the new comparisons are not
counted. The same issue applies to the original test, so we will re-run it counting each target once.

**Update, same day (D209).** We re-ran the biomarker test counting each drug target once instead of
each drug. The 49 comparisons come from 20 distinct targets; 16 of the 20 go the right way (80%), so
the test still passes. The "should fail" check turned out to rest on only 2 distinct targets, which is
too few to tell us anything, so it neither passes nor fails.
