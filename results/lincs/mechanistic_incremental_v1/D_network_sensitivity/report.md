# Experiment D — why one regulator reads differently through different networks

Protocol `5cc2227` + D settings `d3763b7` (frozen before any D computation). Run `20261124T175241-3873623`, 223 s,
inputs = the 2026-11-23 run. Identical drug conditions read through two networks of one tumour family, so the network
is the only thing that changes. Decomposition into shared-target and unique-target parts is exact (max error 8e-7).

## What resolved

**1. Target overlap explains most of the cross-network disagreement, in every evaluable pair.** Across regulators,
the more targets two networks share for a regulator, the more alike they read it across drugs:

| Pair | Window | Regulators | Median overlap (Jaccard) | Median agreement | Overlap vs agreement (95% CI) | Null |
|---|---|---|---|---|---|---|
| BRCA luminal A vs basal | early / late | 782 | 0.10 | 0.71 / 0.73 | **0.82** (0.80-0.84) / **0.81** (0.78-0.83) | 0.06 |
| BRCA luminal B vs basal | early / late | 806 | 0.09 | 0.69 / 0.70 | **0.84** (0.82-0.86) / **0.83** (0.80-0.85) | 0.11 |
| BRCA luminal A vs luminal B | early / late | 801 | 0.15 | 0.79 / 0.81 | **0.82** / **0.82** | 0.10 |
| BRCA subtypes vs pooled BRCA (3 pairs) | both | 760-776 | 0.15-0.30 | 0.77-0.88 | **0.79-0.83** | 0.07-0.08 |
| Cervical vs cervical squamous | early / late | 839 | 0.38 | 0.88 / 0.90 | **0.77** / **0.72** | 0.07-0.11 |
| Esophageal adeno vs squamous | early | 764 | 0.05 | 0.57 | **0.73** (0.69-0.76) | 0.15 |
| PCPG vs pheochromocytoma | late | 802 | 0.39 | 0.91 | **0.68** (0.64-0.72) | 0.10 |
| Endometrial vs endometrioid | early | 745 | 0.21 | 0.84 | **0.75** (0.71-0.78) | 0.30 |

All 16 evaluable pair-windows resolve (null = overlap shuffled among regulators of similar regulon size). The result
holds on held-out drugs: computed separately on two halves of the drugs, the overlap-agreement correlation is 0.63-0.84
in each half, and per-regulator agreement reproduces across halves at 0.75-0.99.

**2. The disagreement sits on the targets one network has and the other lacks.** 83-95% of the absolute difference
in activity comes from unique targets, above their share of edges wherever overlap is substantial (cervical 84% of the
difference from 45% of edges; PCPG 84% from 44%). Shared targets almost always carry the same sign (median sign
agreement 1.0), and sign agreement adds little beyond overlap (0.09-0.29).

**3. Most regulators read the same way in related networks, and the ones that do not are explained by rewiring.**

| Pair (both windows) | Context-robust (agreement >= 0.7 in both drug halves) | Rewired (< 0.3, low overlap) | Unexplained (< 0.3, high overlap) |
|---|---|---|---|
| Luminal A vs pooled BRCA | 653-654 of 768 | 3-6 | 0 |
| Luminal A vs basal | 392-411 of 782 | 63 | 0 |
| Luminal B vs basal | 378-395 of 806 | 90-93 | 0 |
| Cervical vs cervical squamous | 639-760 of 839 | 1 | 0 |
| Esophageal adeno vs squamous | 191 of 764 | 113 | 4 |

Almost no regulator disagrees for unexplained reasons: where two networks read a regulator differently, they also
give it different targets.

**4. ESR1 is the clearest rewired case** (`D_worked_ESR1.md`). Between the luminal A and basal breast networks ESR1
shares 1.4% of its targets (Jaccard 0.014; 207 and 145 targets), the two networks agree on its activity across drugs at
only 0.12 (early) and 0.17 (late), and 97% of the difference comes from unique targets. On the same anti-estrogen
conditions at 6 h, the luminal A network reads ESR1 as lowered by fulvestrant (-1.29) and clomifene (-2.35), the
expected direction, while the basal network reads both as raised (+1.52, +0.58); at 24 h both networks read ESR1 as
raised (inversion, as in Experiment B and D152). Early values rest on one condition per drug.

**5. Subtype networks differ from each other about as much as from another tissue.** The cross-tissue reference
(luminal A vs lung adenocarcinoma) agrees at 0.74-0.75 with overlap 0.11, the same range as luminal A vs basal
(0.71-0.73, overlap 0.10).

## Not resolved at current data depth

- Glioma (1 condition in its matched lines) and sarcoma (0) families: no drug conditions to read; not evaluable.
- Esophageal, endometrial (late) and PCPG (early): fewer than 20 conditions in that window.
- Directional perturbation controls (genetic knockdowns read through each network), which the spec asks for before a
  regulator's interpretation is called mechanistic, were not run here.

## What this supports

Cross-network differences in a regulator's reading are network rewiring, not noise: they track target overlap, sit on
unique targets, and reproduce on held-out drugs. For context-robust regulators an interpretation carries across
related networks; for rewired ones (ESR1 between breast subtypes) it is specific to the network it was read in, which
supports the per-cohort network design. No network was chosen or preferred on any label.

**Implementation note:** the first run stopped because the sarcoma family has no conditions (an empty array in the
exactness check); the rerun skips empty families and reports them in `run.json`. No setting changed.

**Files:** `D_pairs.tsv`, `D_flags.tsv`, `D_regulator_flags.tsv.gz`, `D_regulators.tsv.gz`, `D_worked_ESR1.md`,
`D_worked_ESR1.tsv`, `run.json`.
