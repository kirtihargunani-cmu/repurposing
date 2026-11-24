# Experiment C — can a drug's targets be recovered from its regulator response?

Protocol `5cc2227` + C settings `cbd2fb5` (frozen before any C computation). Run on the 2026-11-23 inputs: 55,892
(drug, target) pairs over cancers and subgroups, both windows, 84 s. Null: each cell's target lists permuted among its
drugs 1,000 times. Uncertainty: cluster bootstrap over target genes. Resolved = interval excludes the permutation
null, in either direction (A1.7).

## What resolved

**1. Primary targets stand out in the drug's own regulator response at 24 h.** Where the primary target is itself a
scored regulator, it sits at the 58th percentile of |activity| against 52 for shuffled targets (lift 1.11; 95% CI
0.54-0.62; p = 0.001; 357 pairs, 14 targets). In subgroups at 6 h the lift is larger (74th vs 48th percentile; 7
targets). Both rest on few distinct targets (feasibility-limited: most drug targets are kinases or enzymes, not
regulators), which is why Experiment A's version of this check did not resolve. The effect is there; the target set is
small.

**2. The primary target's wider downstream neighbourhood moves.** Over 2-hop OmniPath neighbourhoods of primary
targets (83 targets, 1,874 pairs), regulators sit at the 52.2nd percentile against 51.2 for shuffled targets at 6 h
(CI 51.6-52.8, p = 0.001); 24 h is in the same direction (p = 0.01) and just reaches the null. Small, but measured over
many targets, so it is a broad signal rather than a few strong cases.

**3. At 24 h, the downstream regulators move against the naive expected direction.** Signing each 1-hop downstream
regulator by drug orientation x edge sign (an inhibitor of a stimulating target should lower it), the late window is
consistently *below* the null: 45.4 vs 48.5 in cancers (CI 42.7-48.0; 35 targets) and 42.1 vs 48.9 in subgroups (CI
38.8-46.2; 25 targets); the early window is at the null (50.4 vs 49.5). The 6 h response follows the network's signs
no better than chance and the 24 h response inverts them. This matches Experiment B (delayed, program-wide changes) and
the anti-estrogen inversion first seen in D152, and it is the pattern adaptive feedback would produce. It is a
hypothesis for follow-up, not evidence of feedback signalling (that needs a phosphoproteomic or time-course readout).

## Not resolved at current data depth

| Question | Result | Why it stays open |
|---|---|---|
| Primary target, exact, 6 h | 54th vs 51st percentile, p = 0.058 (signed: 57 vs 50, p = 0.026) | 14 targets; intervals include the null |
| Primary target, 1-hop downstream (unsigned) | 51.3 vs 51.4 early; 52.1 vs 51.4 late (p = 0.07) | small effect at this depth |
| Secondary targets with measured binding, exact | 52 vs 52 early, 51 vs 52 late | 15-16 targets are regulators: feasibility-limited |
| Secondary targets, 1-hop downstream | 51.0 vs 50.0 late (p = 0.036), interval reaches the null | 39 targets |
| Annotation-only (Repurposing Hub) targets | at the null throughout | annotation alone carries no recoverable signal |
| Exposure-to-potency (dose / pChEMBL) vs recovery, secondary | Spearman -0.08 to +0.13, all intervals span 0 | 60-62 drugs per window |

The primary-target yardstick resolves at the exact level (24 h) and at 2 hops (6 h) but not at 1 hop unsigned, so the
secondary-target nulls are not yet informative: per the frozen rule they are read as unresolved, not negative.

## What this supports and what it does not

It supports that the regulator response carries target information for drugs whose primary target is
binding-confirmed (exact at 24 h, neighbourhood at 6 h) and that the 24 h response has a structured direction relative
to the signalling network. It does not show target engagement, and it does not show that secondary targets are
engaged. The SIGNOR-only network gives the same picture with fewer targets (descriptive). Pinned networks only
(OmniPath and SIGNOR retrieved 2026-11-10); no curated transcription regulons used (D208).

**Implementation note:** the summary's `resolved` column was computed one-sided (interval above the null); the frozen
rule is two-sided, so the two below-null results above were identified by checking `ci_hi < null`. The 2-hop cap now
sorts its set before truncating, so the selection is reproducible.

**Files:** `C_summary.tsv`, `C_coverage.tsv`, `C_potency.tsv`, `C_pairs.tsv.gz`, `run.json`.
