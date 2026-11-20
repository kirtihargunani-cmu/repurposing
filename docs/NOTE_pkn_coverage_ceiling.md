# Would a bigger knowledge base help? Measuring the PKN ceiling before pricing IPA

**D159, 2026-11-20.** Prompted by the user's question about using Ingenuity Pathway Analysis (IPA)
to trace transcriptional routes from drug targets to dysregulated tumour activity.

The case for a commercial knowledge base rests on one claim: that our coverage gap is caused by the
strictness of *our* filtering, which a larger curated base would fix. This tests that claim using
sources we already hold, before anything is priced.

## First, a correction to how D143 has been quoted

D143 recorded "only ~33% of tested regulators exist in the prior knowledge network". That figure is
for the **TF-only** regulator set. On the **full** set actually scored in D158 — 4,233 regulators
across 7 networks — coverage is **54.1%**. The gap is real but materially smaller than the number
we have been repeating.

Also: of the 1,942 regulators OmniPath misses, **84.6% are named, studied genes** (ABCG1, the
adhesion GPCR family, AGO1, AFF1/3/4, ADNP) and only 15.2% are zinc-fingers. These are not genes
without literature.

## What each gate of our filter actually costs

Our registered filter (`config/datasets.yaml`) is three gates, and an edge must clear all three:

1. `consensus_direction` — the aggregated resources agree the edge points A→B
2. exactly one of `consensus_stimulation` / `consensus_inhibition` — they agree on the sign
3. ≥ 2 distinct resources **or** ≥ 2 distinct PubMed ids

Measured over the 33,343 raw gene-level OmniPath interactions (self-loops and complex nodes removed):

| gate | edges kept | share |
|---|---|---|
| 1 — direction consensus | 20,167 | **60.5%** |
| 2 — + one consensus sign | 20,042 | 60.1% |
| 3 — + evidence floor | 15,860 | 47.6% |

**Gate 1 is the dominant filter**, costing 39.5% of all interactions. Gate 2 is nearly free (125
edges). Gate 3 — the one proposed for relaxation — deletes 4,182 edges, 20.9% of what survives
gates 1–2, and only **213 of them are regulators we score** and **31 are drug-target genes**.

## What relaxing buys

Five arms, scored on node coverage and on directed reachability within ≤ 3 hops from the 387
annotated drug-target genes to the 4,233 scored regulators (1,638,171 possible pairs):

| arm | edges | nodes | % activating | regulators | drug targets | reachable pairs |
|---|---|---|---|---|---|---|
| curated (production) | 15,646 | 5,264 | 70.2% | 2,291 (54.1%) | 258 (67%) | 66,340 |
| relaxed (drop gate 3) | 19,634 | 6,158 | 68.5% | 2,504 (59.2%) | 289 (75%) | 85,862 |
| **permissive** (drop consensus, keep directed) | 19,353 | 6,115 | 70.1% | 2,486 (58.7%) | 288 (74%) | **100,173** |
| curated + SIGNOR | 19,438 | 6,130 | 68.6% | 2,500 (59.1%) | 289 (75%) | 84,847 |
| relaxed + SIGNOR | 19,634 | 6,158 | 68.5% | 2,504 (59.2%) | 289 (75%) | 85,862 |

Three things follow.

**1. Relaxing is free money, and the feared tradeoff does not materialise.** `permissive` gives
**1.5× the reachable pairs** at essentially unchanged activation balance (70.1% vs 70.2%). D145's
failure mode — an empty CARNIVAL solution at 96–97% activating edges — is not triggered by any arm
here. The builder's docstring warning about a "90.5% activating, 70,389-edge" network describes a
different construction and does not apply to these.

**2. SIGNOR adds nothing.** It is already inside OmniPath; the union is identical to `relaxed` to
within 0 regulators.

**3. There is a hard ceiling at ~59%, and it is OmniPath's, not ours.** Every arm lands between
58.7% and 59.2% regulator coverage. **1,729 of 4,233 regulators (40.8%) are absent from OmniPath
under any filtering we can apply.** No amount of loosening reaches them.

## What this means for IPA

The honest reading is split, and it does not settle the question on its own.

**For IPA:** the 40.8% residual is a genuine resource ceiling, not an artefact of our strictness, and
84.6% of the missing regulators are studied genes that a 24M-relationship curated base would
plausibly cover. That is a real argument that a larger knowledge base buys something OmniPath cannot.

**Against buying yet:** D144–D147 found F-3 inverse CARNIVAL recovered **0 of 4** drivers *on genes
it already covered*, and the ARACNe causal map returned an empty solution. If the inference step is
what fails, more edges will not fix it. The decisive experiment is therefore free: re-run F-3 on the
`permissive` arm, which supplies 1.5× the reachable target→regulator pairs at no cost. **If recovery
is still zero with 50% more routes available, coverage was never the binding constraint** and IPA
would be solving the wrong problem.

That test is the pre-registered next step, and it must run before any licence is priced.

## Caveats

1. Reachability here is **structural** — a directed path of ≤ 3 hops exists. It does not mean
   CARNIVAL can find a consistent signed solution along it; sign coherence is a strictly harder
   condition and is what D145 failed.
2. The `permissive` arm admits edges whose sign is **not** consensus across resources. That is a
   deliberate loosening of the builder's "nothing is guessed" rule and must be reported as a
   sensitivity arm, never silently promoted to production.
3. IPA's own knowledge base is also literature-curated and carries the same study bias in kind, if
   not in degree. The 40.8% residual bounds what is missing from OmniPath; it does **not** promise
   IPA holds it. Only a coverage check against IPA's gene list could establish that, and that needs
   access we do not have.

## Guard against a conflation (user, 2026-11-20)

**The D158 3.8× result does not use the PKN, and more edges cannot improve it.** Verified directly:
neither `route1_validation.py` nor `candidate_mechanisms.py` contains any reference to the PKN,
OmniPath or CARNIVAL. Their inputs are the condition profiles, ARACNe VIPER activity, `group_mr.tsv`,
celligner line weights and the ChEMBL annotations. That result would read identically against a PKN
of zero edges.

Two separate machines, and only one consumes edges:

| | **A — identification** | **B — tracing** |
|---|---|---|
| question | which regulator does this drug move opposite to the tumour? | by what route does binding this protein cause that change? |
| inputs | ARACNe regulons, VIPER activity, disease NES | drug target → signed directed paths → CARNIVAL ILP |
| reads the PKN | **no** | yes |
| status | **works: 3.84× in the top 1% for primary targets (D158)** | **0 of 4 recoveries (D144–D147)** |

So a larger knowledge base is an investment in **B**, which does not currently work, not an
improvement to **A**, which does. Any case for IPA must be argued on B's terms.

**The strongest form of that case** is not "more edges, better targets". It is this: D158's secondary
-target failure (0.88×, chance) was measured by **A**, an activity test — does the drug move that
regulator in the opposing direction? **B** would score a secondary target by an entirely different
criterion — is there a coherent signed *path* from it to the dysregulated regulators? That has never
been tested on secondary targets, because tracing does not yet work at all. If B can be made to work
it is a new line of evidence on the secondary-target problem, not a restatement of A's.
