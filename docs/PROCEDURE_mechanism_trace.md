# Procedure — building a dependable mechanism trace

**Status:** procedure, not yet applied. **Drafted:** 2026-11-17.
**Read first:** `docs/notes/settled_measurement_facts.md`.

## 0. The problem this procedure exists to solve

A trace is a causal narrative: this drug binds this target, which sits above these
regulators, whose genes move this way, and that programme is dysregulated in this
tumour. Narratives are **infinitely flexible**. For any drug and any tumour a
plausible chain can be assembled, and unlike a score - which fails loudly when
miscalibrated - a well-told mechanistic story simply sounds right.

So almost all of this procedure is about making traces **capable of failing**.

**The governing rule: the trace is a DATA STRUCTURE; the narrative is a rendering of
it.** Prose is generated from the structured record, never the other way round.
Anything asserted in prose that has no field behind it does not belong in the trace.

## 1. The chain, declared before any data is read

Seven links, evaluated in **forward order only** - drug to target to disease to
regulators to genes to observation to closure. Never constructed backwards from a
conclusion ("this drug works, what explains it?"). The forward constraint is what
stops the chain being fitted to its endpoint.

| # | Link | Question |
|---|---|---|
| L1 | drug → target | Does the drug bind the target, and on what evidence? |
| L2 | target → disease | Is the target altered, over-expressed or depended on in this cohort, in what fraction of patients? |
| L3 | target → regulators | What signed path reaches which regulators, in how many hops? |
| L4 | regulators → genes | Which genes, with what weight and expected direction? |
| L5 | predicted → observed | Did those genes move as predicted, in which cell line, matched how? |
| L6 | concentration | Is the signal carried by many independent pieces of evidence, or one? |
| L7 | disease closure | Are those same genes dysregulated in the tumour in the direction the drug reverses? |

L6 is not a step in the biology. It is the audit that the BRAF case showed we need:
an apparently healthy L4-L5 that rests entirely on one hub is not five hundred
pieces of evidence, it is one. **L6 counts; it does not subtract** - see §5. Report
the effective number of regulators carrying the footprint, never a score recomputed
with a regulator's genes removed.

## 2. Every link carries the same seven fields

No link may be recorded as prose. Each is a record:

| field | what it holds |
|---|---|
| `expected` | what we predict **before** observing, from the link above |
| `observed` | what the data shows |
| `verdict` | holds / weak / **breaks** / not evaluable |
| `provenance` | which file and which decision row the evidence comes from |
| `independence` | which other links share that source (see §3) |
| `concentration` | the largest single contributor's share of the link's evidence |
| `alternatives` | competing explanations that survive this observation |

`expected` must be written from the link above, before `observed` is read. A link
whose expectation is written after its observation is a description, not evidence,
and is marked `not evaluable`.

## 3. Independence is declared, not assumed

A chain of seven links is only seven pieces of evidence if the links have
independent sources. Where two links draw on the same underlying data they are
**one link reported twice**, and the record must say so.

Known entanglements in this pipeline, to be declared wherever they occur:

- L3 and L4 both derive from the regulon. A regulon error moves both.
- L2 (dependency) and any essentiality-based reasoning in L6 share DepMap.
- L5 and L7 share the cell line only when the line is a direct tissue match; with a
  **borrowed** line (76.7% of unit rows) they do not, and the trace must say which.
- ARACNe regulons are learned from tumour co-expression, so L3/L4 and L7 are not
  fully independent for ARACNe-derived regulators (D133.1).

**A trace reports its effective number of independent links, not its total.**

## 4. The weakest link is the verdict

No averaging, no scoring, no compensating. A chain with six strong links and one
broken link is a **broken chain**, and the report leads with the break.

This is the whole reason for tracing. A score of 0.03 averages a dead link into a
live one and returns a number that looks mildly encouraging.

## 5. Attribution is DESCRIBED, never used to subset the response

**Withdrawn 2026-11-17 (user).** An earlier version of this section told the trace
to compare "uniquely attributable" genes against shared ones, and to recompute the
score with a dominant regulator's genes removed. **Both operations are withdrawn,
and results produced by them are not to be cited.**

**Why the unique-vs-shared test was confounded.** A gene is "uniquely attributable"
because no other regulator claims it. Where the median gene is claimed by ~12
regulators, being claimed by one is unusual - and what makes a gene unusual that way
is a small regulon with thin evidence, not biological specificity. The first
application found fulvestrant's 41 unique genes 39.0% concordant against 52.9% for
437 shared ones, and read it as "the ESR1-specific genes move the wrong way". It is
at least as likely to be **well-annotated genes versus poorly-annotated genes** -
a property of how CollecTRI was built, not of ESR1.

**Why leave-one-out subtraction is worse.** BRAF reaches transcription *through*
MEK-ERK-MYC. Removing MYC's genes and recomputing does not give a purer BRAF signal;
it gives BRAF with its principal effector deleted. The same holds for any cascade:
shared structure is real biology (D133 established the overlap pairs are obligate
complexes and close paralogs - NFKB1/RELA, RFX5/RFXANK/RFXAP, BMAL1/CLOCK), so
subtracting it removes mechanism, not confound.

**What is kept.** The *structural* description, which counts rather than tests:
how many regulators claim each gene, and how many **effective** regulators the
footprint represents. "1,063 genes arrive through ~2.7 effective regulators" is a
statement about how much independent evidence exists. It requires no subsetting,
no removal and no comparison between gene subsets.

**The standing rule: the transcriptional response is holistic.** Do not carve a
subset out of it and judge the subset. Where the question is whether the response
is consistent with the mechanism, compare **whole response to whole response** -
against the target's genetic perturbation (`SPEC_drug_genetic_concordance.md`) or
across the drug's own dose series (`SPEC_dose_response_bmd.md`).

## 6. Alternatives are enumerated at each link, not at the end

At minimum, at each link, and carried forward if they survive:

- **generic response** - would these genes move under most drugs in this line?
- **cytotoxicity** - is the cell dying? Report alongside the proliferation and
  `dna_acting` indicators the pipeline already computes.
- **off-target** - does the drug hit something else that explains this better?
- **amplitude** - is this just a large overall response rather than a specific one?
  (D131/D133: |PC1| correlates 0.827 with overall response magnitude.)
- **reference artefact** at L7 - is the target's pathway active in the normal
  reference tissue, which inverts the comparison? (D128-A3.2.)

An alternative that is not ruled out is **carried into the verdict**, not dropped.

## 7. Negative controls — the procedure's own validation

**Before any real trace is believed, the procedure must break on pairs where it
should.** If a chain closes on these, it is generating stories and no positive
trace from it means anything.

Required negative controls:

1. **Target not expressed** in the source line - L5 must break.
2. **RAF inhibitor in BRAF-wild-type, RAS-active** context - already available from
   D128-A4, where it behaved correctly.
3. **Random drug-cohort pairing** - shuffle drugs across cohorts; L2 or L7 must break.
4. **A cytotoxic** where the "mechanism" is cell death - L6 and the cytotoxicity
   alternative must fire.
5. **A drug whose target has no anchored TF** - L3 must return `not evaluable`
   rather than finding a path anyway.

**Pre-stated: the fraction of negative controls on which the chain closes is the
procedure's false-positive rate, and it is reported with every set of traces.**

## 8. Order of operations

1. Declare the pair and the chain structure. Nothing read yet.
2. L1, L2 - drug and disease side, independent of any response data.
3. L3, L4 - the predicted footprint, with expected directions written down. **Still
   no observation read.**
4. Freeze the prediction (file written, committed).
5. L5 - read the observed response and compare.
6. L6, L7.
7. Verdict from the weakest link; render prose from the record.

Step 4 is the discipline that makes the whole thing work: the predicted footprint
and its directions exist as a committed artefact before any drug profile is read.

## 9. What a trace is not

- **Not a score.** It produces a verdict and a break-point, not a number, and
  traces are not ranked against each other.
- **Not evidence of benefit.** A closed chain means the mechanistic account is
  consistent, not that the drug helps a patient.
- **Not independent of the annotation.** L1 rests on ChEMBL binding calls.
- **Not a substitute for a holdout.** Nothing here reads an outcome label, which is
  what makes it cheap - and also what limits it.

## 10. First application

Two contrasting pairs, built end to end, to test whether the procedure
discriminates where scores did not:

- **TF-anchored:** fulvestrant → ESR1 → zero hops → regulon → MCF7 → BRCA-LumA.
- **Kinase:** vemurafenib → BRAF → signalling hops → regulators → A375 → SKCM.

We already know the second is ~65% MYC and the first fails its corrected null. **The
test of this procedure is whether L6 surfaces the MYC dominance without being told,
and whether the two traces break at different links.** If both produce the same
verdict, the procedure adds nothing over the score and should be abandoned.

Then the five negative controls of §7, before any positive trace is reported.

## 11. Result of §7 — the procedure's false-positive rate (D137, 2026-11-17)

The five negative controls were run for the first time, on cohort-matched ARACNe
regulons. Pre-stated before any number was computed: a control **closes** when no
link L1–L7 has verdict `breaks` or `not evaluable`, and the fraction that closes
is the procedure's false-positive rate.

| | pair | network | break | closes |
|---|---|---|---|---|
| real | fulvestrant / ESR1 / BRCA-LumA | BRCA_LUMA | L6 | no |
| real | vemurafenib / BRAF / TCGA-SKCM | SKCM | L5 | no |
| NC1 | regorafenib / RET / TCGA-LUAD — target not expressed | LUAD | L4 | **YES** |
| NC2 | vemurafenib / BRAF / BRCA-Basal — RAF paradox | BRCA_BAS | L5 | no |
| NC3 | MEK-162 / MAP2K2 / TCGA-READ — random pairing | COADREAD | L2 | **YES** |
| NC4 | estramustine / ESR2 / BRCA-LumB — cytotoxic | BRCA_LUMB | L1 | **YES** |
| NC5 | midostaurin / KDR / TCGA-LIHC — no anchored TF | LIHC | L3 | no |

**False-positive rate 3 / 5 = 0.60, and both real pairs break.** The procedure is
not merely noisy on ARACNe; it is anti-correlated with the thing it is meant to
detect. The random pairing scores at the **99.9th percentile** of all compounds in
its line — the strongest apparent mechanism of anything tested.

Per §7 this is binding: **no trace produced by this procedure, on either network,
counts as a positive.** That includes the two CollecTRI traces of §10, whose
surviving findings are hereby reported as uncalibrated.

### What the controls established about the procedure itself

1. **The footprint is a property of the network, not of the drug or the target.**
   The same target in the same cohort gives nearly disjoint gene sets under the two
   regulons — fulvestrant/ESR1 shares 6 genes between a 478-gene CollecTRI and a
   131-gene ARACNe footprint — and the mean signed response flips sign for both real
   pairs. L4 has no network-agnostic definition, so §4's "weakest link" verdict is
   reading the regulon choice.
2. **L6's concentration test was partly measuring CollecTRI's regulon-size skew.**
   CollecTRI's MYC claims 881 genes against a median regulon of 9; ARACNe's median
   is ~100. On SKCM the BRAF footprint's hub is MEF2D at 20.4% over 11 regulators
   (8.0 effective), not MYC at 65.1%. The §10 claim that L6 "surfaced the MYC
   dominance unprompted" is withdrawn: it surfaced a property of the network.
3. **L6 is undefined for a zero-hop TF anchor.** Hub share is 1.00 by construction
   when L3 reaches one regulator, so fulvestrant breaks at L6 on both networks for a
   structural reason rather than an observed one.
4. **Classical chemotherapy is outside the procedure's domain.** TOP1, TOP2A, TYMS,
   DHFR, GART, RRM1 and TUBB have zero OmniPath route rows, so a cytotoxic acting
   through its own effector is `not evaluable` at L3 — not a verdict about the drug.
5. **NC3 is a weak control and its closure is partly our fault.** The shuffle
   randomised the *cohort*, then chose a celligner-matched line for it; MEK
   inhibition in BRAF-V600E HT-29 is real pharmacology. A stronger version must
   randomise the drug–target–line triple together. The other two closures (NC1, NC4)
   carry no such excuse.

**Not supported on checking:** that ARACNe's ESR1 regulon is noise for want of
variance in a subtype cohort. ESR1's sd in BRCA_LUMA is 1.282, the 95.2nd
percentile of all genes.

## 12. ARACNe-only, and the stability test it forced (D138, 2026-11-17)

CollecTRI is retired by decision; only ARACNe regulon sets are used. That removes
the "two-network agreement" repair §11 proposed, so the equivalent test has to run
*inside* ARACNe: rebuild the same footprint once per cohort network, change nothing
else, replay the identical L5.

| | fulvestrant / ESR1 / MCF7 | vemurafenib / BRAF / A375 |
|---|---|---|
| networks with a footprint | 41 | 41 |
| footprint size | 3 – 2,295 genes | 476 – 1,451 |
| pairwise Jaccard (median) | **0.011** | 0.112 |
| mean signed f(z) | −0.2524 … +0.0492 | −0.0374 … +0.0986 |
| **positive in** | **5 of 41** | 34 of 41 |
| generic percentile | 1.3 – 72.0 | 16.1 – 97.5 |
| **its own cohort network** | BRCA_LUMA: −0.018, 39.6th | **SKCM: −0.023, 26.6th — worst of the 41** |

**The instability is internal to ARACNe.** Two ARACNe networks share ~1% of their
ESR1 footprint. Changing which cohort network is used moves vemurafenib's answer
from the 16th to the 97.5th percentile — and the cohort-matched network, the one we
would defend as biologically correct, is the *worst* of the 41.

**Rejected as a rescue:** that low-variance cohorts give noise regulons. Restricting
ESR1 to the 17 networks where ESR1 is above the 75th percentile of expression
variance leaves 3 of 17 positive; all-BRCA, the best-powered ESR1 regulon in the
project (sd at the 99.8th percentile, 251 genes), gives −0.077 at the 18.3rd
percentile.

**This is not new.** The ARACNe ESR1 regulon reading anti-estrogens in MCF7 as ESR1
*activation* is recorded in D83 (+2.2 to +2.7 across all 10 profiles, against
CollecTRI's −2.4), D85, D88, D89, D90 and D92. The trace reproduced a known
property of the regulon rather than discovering a new artefact — which is the one
encouraging thing here, because it means the measurement is behaving consistently.

**Consequence for this procedure.** Under ARACNe-only, L4 cannot be defined from a
single cohort network: for a TF-anchored target every well-powered network reads the
drug backwards, and for a signalling target the cohort-matched network is an
outlier. §11's verdict stands and is now unconditional — **no trace from this
procedure counts as a positive.**
