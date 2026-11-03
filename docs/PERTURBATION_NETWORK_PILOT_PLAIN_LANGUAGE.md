# Perturbation-network pilot in plain language

Written 2026-11-02 before the pilot ran; the result is added at the end. The rule is decision D86 in [`DECISIONS.md`](DECISIONS.md); the
findings that led to it are in [`PROGRESS.md`](PROGRESS.md) (2026-11-02 entries). Background on the drug side:
[`DRUG_SIDE_PLAIN_LANGUAGE.md`](DRUG_SIDE_PLAIN_LANGUAGE.md).

## The problem this tries to fix

To read which control switches (transcription factors) a drug flips, we need to know which genes each switch
controls. That list is the switch's "regulon". We get it from ARACNe, which looks at hundreds of patient tumors and
asks which genes rise and fall together with each switch.

The drawback is that tumors mostly differ in **what kind of cell they are**: luminal or basal breast cancer, more or
less mature, more or less immune infiltration. Those differences built up over years. So ARACNe learns which genes go
with which cell state. That is the right thing for describing a disease, and the disease side works well with it.
It is not the same thing as what a drug changes in a few hours.

The clearest example is ESR1, the estrogen receptor. Its tumor-learned target list is mostly the luminal program
(GATA3, FOXA1 and similar genes). Anti-estrogen drugs block ESR1, but in MCF7 breast cancer cells they raise many of
those luminal genes. So the tumor map reads anti-estrogens as **switching estrogen on**, the opposite of what they do
(caveats C18, C25).

## What we already tried

- **Removing cell-identity links** (D82, D83): dropping links that exist only between subtypes helped the disease
  side, but drug readings barely changed.
- **Checking links against experiments** (D84, exploration 2026-11-02): we tested every ARACNe link against LINCS
  experiments where a switch was knocked down or overexpressed. About 55% of links move the predicted way. Links
  found in many tumor types are more reliable (up to 68%) than links found in one (52%).
- **Down-weighting the unreliable links, or removing hidden factors** (the main patterns of variation across
  tumors): this made individual links more accurate but did **not** make the switch readings more accurate. The
  biggest patterns in tumors are the main switches' own programs. In breast cancer, 92% of ESR1's links run through
  them, so removing the patterns removes the switches too.
- **The curated list (CollecTRI)** is not better link by link either (52% of its links move the predicted way). It
  does read anti-estrogens correctly, which is why the drug side currently averages the two lists.

The conclusion was that tumor data cannot, even in principle, separate "genes ESR1 switches on" from "genes that are
part of the luminal state". A better list has to come from a different kind of data.

## The idea

Use the same algorithm (ARACNe), but feed it **drug experiments in one cell line** instead of patient tumors.

LINCS measured what about 1,100 to 1,200 different drugs do to each of several cell lines. In one line the cell type
never changes. The only thing that differs between experiments is the drug's short-term effect. Genes that move
together across many drug treatments are linked by the wiring that responds within hours, which is what we need to
read drugs.

An analogy: the tumor map compares many different cities and notices that cities with more hospitals have more
doctors. The perturbation map takes one city, changes things one at a time, and watches what moves.

The new lists would be used **only on the drug side**. The disease side keeps the tumor maps.

## What the pilot does

1. **Three cell lines:** MCF7 (breast), PC3 (prostate) and A375 (melanoma). They have the most drug experiments and
   the most knockdown / overexpression experiments to check against. MCF7 alone has too few checks to decide.
2. **Measured genes only.** LINCS measures 978 genes directly and calculates the other ~9,000 from those. We checked
   that the calculated genes are about 88% a fixed mix of the measured ones, so a network built on them would mostly
   rediscover the calculation formula. The pilot therefore uses the 978 measured genes, and the 81 of them that are
   switches.
3. **Two versions of the drug data:** as measured ("raw"), and with each line's average drug response removed
   ("line-specific", what the drug side uses now). This also feeds the open question D85 about which version to use.
4. **No peeking.** Drugs are split into five groups by how they work. A drug is only ever scored with a map built
   without its group, so the map cannot simply memorise that drug's own effect.
5. **36 maps in total** (3 lines x 2 versions x 6 builds: one with all drugs and five leaving one group out), each
   with the same ARACNe settings as the tumor maps.

## How the pilot decides

The check uses experiments that were never used to build any map: LINCS knocked down or overexpressed many switches
directly in these three lines. For each experiment we know the right answer:

- knockdown: the switch's activity should read **lower**
- overexpression: it should read **higher**

Three lists read every experiment, using the same measured genes and only the switches all three can score:

| List | Learned from |
|---|---|
| Perturbation map (new) | drug experiments in the same cell line |
| Tumor map (current) | tumors of the matching cancer (breast, prostate, melanoma) |
| Curated list (CollecTRI) | published literature |

We count how often each list gets the direction right, comparing them experiment by experiment.

**The pilot passes** if, in either data version, the perturbation map

1. gets the direction right at least **5 points more often** than the tumor map,
2. and that lead is **not chance**: among experiments where the two give different answers, the perturbation map
   wins clearly more often (p < 0.05),
3. and it is **no more than 2 points worse** than the curated list.

Pass: build maps for the other cell lines with enough data. Fail: stop and report. Either way, using the maps in the
actual drug ranking needs its own rule later. All of this was written down before running anything.

**Also reported, but not part of the decision:** how many switches each list can score; whether the maps stay the
same when a fifth of the drugs is left out; and the few drugs that act directly on a measured switch (for example
steroid drugs on the glucocorticoid receptor NR3C1).

## What the pilot cannot tell us

- **ESR1 is not tested.** ESR1, AR, FOXA1, E2F1, SP1 and several other important switches are among the calculated
  genes, so the pilot cannot build maps for them. If the approach works for the measured switches, a follow-up could
  keep the curated target list for these and learn each target's direction from the drug experiments.
- **It is a verdict on the approach, not on each switch.** There are too few experiments per switch to choose a list
  switch by switch.
- **The test is small.** About 200 experiments, probably fewer after the "all three lists can score it" filter.
  A real gain of about 5 points could fail to reach significance; the test reliably detects gains of about 10 points.
  A fail would mean "no clear advantage at this size", not proof that the idea is wrong. The size of the difference
  is reported either way.
- **A known leak in the drug groups.** Grouping by how drugs work uses the wording of the mechanism labels, which
  differ for drugs that work the same way: tamoxifen and fulvestrant ended up in different groups. This only affects
  the secondary drug check, not the main knockdown / overexpression test, and will be fixed before any scale-up.

## Where things are

| What | Where |
|---|---|
| Rule | `docs/DECISIONS.md`, D86 |
| Inputs (profiles, regulators, drug groups) | `scripts/pnet/pnet_inputs.py`; groups in `results/lincs/pnet/folds.tsv` |
| Building the maps | `scripts/pnet/pnet_run_aracne.sh` (24 threads); networks in `data/processed/lincs/pnet/` (not in git) |
| Scoring and the decision | `scripts/pnet/pnet_eval.R`; results in `results/lincs/pnet/` |
| Earlier explorations | `scripts/explore/` |

## Result (2026-11-02): no advantage, so we stopped

- **The main test came out much smaller than planned.** Only 45 experiments on 13 switches could be scored by all
  three lists, because on the 978 measured genes the curated list has enough targets for few switches.
- **The perturbation map did not do better.** It got the direction right 62% of the time, the tumor map 69% and the
  curated list 62%. The gap is not significant either way. Both data versions (raw and line-specific) gave the same
  totals.
- **A wider, informal comparison agrees.** Dropping the curated-list requirement gives 133 experiments on 37 switches:
  57% for the perturbation map vs 59% for the tumor map. It did better in MCF7 and worse in A375.
- **The maps are only moderately stable.** About two thirds of the links reappear when a fifth of the drugs is left
  out.
- **The drug check was too small to read** (30 drug profiles).

What this means: for the switches LINCS measures, maps learned from drug experiments read switches no better than the
tumor maps, at least at this size. The pilot could not test ESR1 and the other unmeasured switches, so the
anti-estrogen problem is still open; for now the drug side keeps averaging the tumor map with the curated list, which
reads anti-estrogens correctly.
