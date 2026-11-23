# Experiment B — what changes between early and late, under an empirical noise model

Protocol `5cc2227` + B settings `2e55c43` (frozen before any B computation). Run on the 2026-11-23 inputs: 16 units
with drugs measured in both windows (14 cancers, 2 subgroups), 9,880 (cancer, drug, arm)
rows, 7,743 pseudo-replicate rows, 102 s. Observed windows only; nothing imputed. Score: gene-level reversal
`rev_gene` (production); `rev` beside it.

## What resolved

**1. The noise model is calibrated, so confident calls can be trusted.** Half-differences of the same drug within one
window fall inside 1.96 SE 95.3% of the time for `rev_gene` (95.8% in subgroups), 94.8% for `rev`, 97% for regulator
activity, against 95% expected. Per-condition noise (pooled within-drug SD): `rev_gene` 0.20 early, 0.14 late in
cancers (`B_noise.tsv`).

**2. Time changes the response beyond replicate noise.** Measured against the pseudo-replicate null (two halves of one
window, labelled early and late, classified the same way):

| Transition (cancers, `rev_gene`) | Matched-dose arm (primary; 3,712 rows, 779 drugs) | Production arm (4,128; 867) | High-confidence (1,224; 289) | Pseudo-replicate |
|---|---|---|---|---|
| Confidently active only at 24 h (late only) | **15.8%** (14.3-17.3) | **18.9%** (17.2-20.5) | **26.7%** (23.6-29.8) | 11.3% |
| Same confident sign at both times | **4.5%** (3.6-5.5) | **5.4%** (4.4-6.4) | 11.8% (9.3-14.6) | 11.5% |
| Confident opposite signs (sign reversal) | 0.30% (0.11-0.51) | **0.48%** (0.26-0.77) | **1.23%** (0.56-2.06) | 0.16% |
| Confidently active only at 6 h | 7.1% | 7.4% | 11.6% | 9.4% |

Delayed responses are the clearest temporal signal: a drug becoming confidently active only by 24 h is 1.4-2.4 times
the replicate rate in every arm. Persistence across 6 -> 24 h is about half of what two replicates of the same window
show, so the two windows genuinely measure different states (consistent with D183). Confident sign reversals are rare
but resolved above the false-flip rate in the production arm (3x) and the high-confidence arm (7.7x); in the
matched-dose arm the interval reaches down to the false rate (not resolved at that depth). Subgroups show the same
late-only excess (16-41% vs 11%).

The earlier 37% early/late sign discordance (D182) is therefore mostly a sign difference where one window is below
its own noise; confident reversals, when they happen, are real.

**3. Mechanism classes have distinct time profiles** (cancers, matched-dose arm; odds ratio against all other drugs;
classes fixed by the outcome-blind rule, 9 classes with >= 8 drugs):

| Class | Transition | Drugs (hits) | Odds ratio (95% CI) | Permutation p | FDR |
|---|---|---|---|---|---|
| PI3K inhibitors | late only | 10 (24) | **4.7** (1.9-11.6) | 0.001 | 0.023 |
| CDK inhibitors | persistent | 12 (15) | **5.0** (2.0-8.8) | 0.004 | 0.023 |
| HDAC inhibitors | early only | 18 (24) | **3.4** (1.7-5.6) | 0.003 | 0.023 |
| HSP inhibitors | early only | 8 (16) | **5.8** (2.1-11.9) | 0.004 | 0.023 |
| MEK inhibitors | late only | 11 (11) | 1.8 (1.05-3.0) | 0.10 | 0.40 |

All four clear both benchmarks (OR >= 1.5, FDR < 0.05). The production arm reproduces CDK persistent (OR 4.7,
2.4-7.2) and PI3K late only (OR 3.3, 1.4-8.5). The pattern is biologically coherent: chromatin and chaperone drugs act
on transcription fast, PI3K-pathway effects build by 24 h, and CDK inhibition holds.

**4. The change is spread across the network, not driven by a few regulators.** The effective number of regulators
carrying the early -> late change in `rev_gene` terms is about 357 (top five carry 4.8%); the ten largest
contributors at 6 h and at 24 h overlap at Jaccard 0.11, and shared top regulators never flip sign. Transitions are
program-wide replacement, not a handful of switches flipping.

**5. Worked cases** (`B_case_studies.tsv`, top 20 regulators by change). Clomifene in luminal A breast cancer is a
confident sign reversal (`rev_gene` +0.34 at 6 h, -0.40 at 24 h): the anti-estrogen reversal D152 saw in MCF7,
recovered here as a calibrated call. Tamoxifen and fulvestrant in basal breast become confidently active only at 24 h.
Erlotinib (lung adenocarcinoma) and trametinib (melanoma) are not measured in both windows in those cancers.

## Not resolved at current data depth

- Sign reversals in the matched-dose arm (0.30%; interval reaches the 0.16% false rate).
- High-confidence subgroup arm (64 rows, 32 drugs): intervals too wide.
- The original B4 confirmatory extras (class enrichment re-run on the pseudo-replicate split and on the high-confidence
  arm) were not run; the frozen config's rule (interval excludes 1, benchmarks reported) is what is applied here.

## Implementation note

After the smoke test (two cancers, before the full run) the odds ratio's +0.5 cell correction was removed: with rare
outcomes it produced odds ratios near 16 from zero hits. Odds ratios are now plain and marked not estimable when a
cell is zero. No threshold, class or rule changed.

**Files:** `B_summary.tsv`, `B_enrichment.tsv`, `B_noise.tsv`, `B_transitions.tsv.gz`, `B_pseudoreplicate.tsv.gz`,
`B_case_studies.tsv`, `run.json`.
