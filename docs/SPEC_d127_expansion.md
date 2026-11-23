# SPEC — Expanding the biomarker-expectation test (D207)

**Status: PRE-REGISTRATION. Every pair below is declared before any of them is scored.**
User, 2026-11-23: *"yes add these biomarkers."*

---

## 1. Why, and the one rule that makes it worth doing

D127 is the drug side's primary validation and rests on **49 instances across 7 primary pairs**
(38 correct, 78%). It is the only anchor in the project with a clear mechanistic signal, and it is
small. Everything else we have — approval recovery (D192: bimodal, 9 of 41 cohorts significantly
inverted), target role (D205: partly shares inputs with the lead score) — is weaker.

**The rule: the pairs are fixed here, before scoring.** D127 already separates `primary` (the
decision set) from `seen` (pairs whose values were known when prior decisions were taken) for
exactly this reason. If pairs were chosen, scored, and then kept on the basis of the result, the
expanded test would be worthless. So this file lists **all** candidates, including those expected
to fail, and the run happens once.

**The existing 7 primary pairs are not touched.** New pairs enter as the set `expansion`, scored
with the same pre-registered thresholds and reported separately; a pooled primary+expansion figure
is secondary.

## 2. Thresholds — unchanged from D127, not re-opened

```
SHARE_MIN = 0.70     share of instances with layered_pct(high) > layered_pct(low)
P_MAX     = 0.05     one-sided sign test
GAP_MIN   = 10.0     median gap, percentile points
SCORE     = layered_pct
```

## 3. The pairs, with the expectation declared

An instance is one compound scored in **both** cohorts of a pair.

### 3.1 Expansion-primary — clean single-target rationale

| class | target genes | high | low | expected | drugs |
|---|---|---|---|---|---|
| **AR** | `AR` | TCGA-PRAD | TCGA-LUSC | **pass** — androgen signalling drives prostate and is absent in squamous lung | ~14 |
| **PI3K** | `PIK3CA/CB/CD/CG` | BRCA-LumA | BRCA-Basal | **pass**, weakly — PIK3CA mutation is enriched in luminal breast | ~5 |
| **PI3K** | as above | BRCA-LumB | BRCA-Basal | **pass**, weakly | ~5 |

### 3.2 Expansion-secondary — multi-kinase, expected to be weaker

Declared **in advance** as weaker controls: `cabozantinib` is annotated against VEGFR, MET **and**
RET; `regorafenib` against VEGFR, KIT, RET and BCR-ABL. A pair built on promiscuous inhibitors
tests less than its drug count suggests, so these are reported and do **not** carry the result.

| class | target genes | high | low | expected | drugs |
|---|---|---|---|---|---|
| VEGFR | `KDR, FLT1, FLT4` | TCGA-KIRC | TCGA-LUSC | pass, weakly — VHL-loss renal is VEGF-driven | ~26 |
| VEGFR | as above | TCGA-KIRP | TCGA-LUSC | pass, weakly | ~26 |
| MET | `MET` | TCGA-KIRP | TCGA-LUSC | pass, weakly — papillary renal is MET-driven | ~16 |
| RET | `RET` | TCGA-THCA | TCGA-LUSC | pass, weakly — RET fusions in thyroid | ~9 |

### 3.3 Negative control — declared to FAIL

| class | target genes | high | low | expected |
|---|---|---|---|---|
| **HDAC** | `HDAC1/2/3` | TCGA-LUAD | TCGA-LUSC | **no difference (share ≈ 0.5)** — there is no HDAC biomarker distinguishing these two |

**This is the most informative pair in the document.** If the negative control "passes", the test
is detecting something other than biomarker biology — cohort-level score inflation, for instance —
and the positive pairs cannot be trusted either. It is a check on the instrument, not on a drug.

## 4. What is NOT added, and why

- **IDH.** The LGG-IDHmut-codel / IDHmut-noncodel / IDHwt subgroups exist and would be an ideal
  contrast, but **no IDH inhibitor is in the scored LINCS set** — ivosidenib and enasidenib
  postdate it. Recorded so it is not re-proposed.
- **KIT** (16 drugs): no GIST or mastocytosis cohort to serve as the high side.
- **SMO/Hedgehog** (4 drugs): no basal-cell or medulloblastoma cohort.
- **BTK** (2), **mTOR** (1): too few scored drugs to form an instance set.

## 5. Caveats registered in advance

- **TCGA-PRAD has only 7 in-key cell lines** (median across cohorts is 64). The AR pair is the
  strongest biology here and rests on the thinnest line matching; if it fails, that is a candidate
  explanation and must be reported as such rather than as a biology result.
- `TCGA-LUSC` is the low side of five pairs. A systematic artefact in LUSC's `layered_pct`
  distribution would move all of them together, so the per-pair shares are reported individually
  and a uniform result across LUSC pairs is treated as suspicious, not as replication.
- Instance counts above are from the scored-drug annotation and will shrink: an instance needs the
  compound scored, gate-passing and `layered_pct`-populated in **both** cohorts, and `layered_pct`
  is populated for only ~15% of rows.

## 6. Acceptance

Reported, per pair and pooled within set:

1. **expansion-primary** (AR + 2× PI3K) against the §2 thresholds.
2. **expansion-secondary**, same thresholds, reported — not decision-bearing.
3. **negative control**: share must be consistent with 0.5. A share ≥ 0.70 here **invalidates the
   expansion** and is written up as such.
4. **primary + expansion-primary pooled**, as a secondary figure, with the original 7 pairs' result
   restated unchanged alongside.

**Pre-committed:** no pair is removed after seeing its result; a pair with too few instances is
reported as underpowered, never dropped; the §2 thresholds are not re-opened.

## 7. Registration

The git commit adding this file is the registration. Later changes to §2–§3 go in a new dated
section, never an edit in place.
