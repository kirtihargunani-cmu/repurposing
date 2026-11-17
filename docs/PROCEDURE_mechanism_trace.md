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
