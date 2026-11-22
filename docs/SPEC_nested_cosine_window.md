# SPEC — D150's nested cosine, re-run window-separated with the network mapping fixed (D191)

**Status: PRE-REGISTRATION. Written and committed before the re-run.**
User, 2026-11-22: *"let's re-run the D150 nested-cosine test window-separated with the
Her2-to-BRCA_BAS mapping bug fixed."*

---

## 1. What is being changed, and nothing else

Two changes to `scripts/systems_validation/integrate/nested_cosine.py`, re-implemented as
`scripts/79af_nested_cosine_window.py`. The cosine itself is **not** touched.

**(a) The network mapping bug.** `nested_cosine.py:37` hardcodes

```python
COHORTS = {"BRCA-LumA": "BRCA_LUMA", "TCGA-SKCM": "SKCM", "BRCA-Her2": "BRCA_BAS", ...}
```

There is **no Her2-specific ARACNe network**. `config/group_networks.tsv` carries overrides only
for BRCA-LumA / LumB / Basal, and production resolves an unlisted group through its TCGA cohort
(`TCGA-BRCA` → `BRCA`, the pooled breast network). So BRCA-Her2 was scored against **basal**
disease biology. This is the same defect D171 fixed in `window_reversal.py`, where group-ID string
parsing turned `BRCA-Her2` into `Her2` and silently dropped the cohort.

The fix is **not** a corrected literal: the new script imports `net_for()` from
`window_reversal.py`, the function that already encodes the production rule, so the mapping has one
source of truth and cannot drift again.

**(b) The window separation.** D150 used `drug_line_profiles_bing_specific.f32`, which `79c`
averages over doses **and times**. This uses `drug_condition_profiles_specific.f32` (same 10,112
gene space, verified) and splits by `pert_time`: **early ≤ 6 h, late ≥ 24 h**; `other` and the
unrecorded sentinel are dropped, not pooled (D180).

Aggregation **mirrors `window_reversal.py`, not `79e`** — conditions, in-key lines only,
line-weighted mean — because D182 established that following `79e` instead produces numbers that
silently disagree with the window table.

## 2. Why this is worth re-running at all

D150's controls failed: anti-estrogens scored 76th–96th opposing percentile in **basal** disease
against 0.9th–12.7th in LumA/LumB — the wrong subtype — and BRAF inhibitors sat at chance in
melanoma (vemurafenib 50.5th, dabrafenib 43.3rd).

Both failures have a candidate explanation that is **measured, not assumed**:

- the Her2 arm was scored against the wrong network (1a);
- D152 measured the ESR1 programme carrying **opposite signs** by exposure time in MCF7 — mean
  signed response **−0.056 at 6 h** (correct) and **+0.094 at 24 h** (inverted), and as cosine
  **−0.132 at 6 h vs +0.093 at 24 h**. Pooling averages a correct reading with an inverted one, and
  within MCF7 the anti-estrogens contribute **6 six-hour signatures against 59 at 24 h**, so the
  inverted window outvotes the correct one roughly 10:1.

**A caveat registered in advance.** 2026-11-22 measurement: the contribution cosine and the raw-z
cosine agree at Spearman 0.942 (0.968 against a w²-reweighted z cosine) with 98.5% sign agreement,
so reparameterising the vectors is not where any gain lives. **Only the time axis is a live
hypothesis here**, which is why it is the only thing being changed besides the bug.

## 3. Cohorts, drugs, controls

Same four cohorts as D150, so this is like-for-like: **BRCA-LumA, BRCA-Her2, TCGA-SKCM,
TCGA-LUAD**. `MIN_GENES = 15` per regulon, unchanged.

**Controls, named before the run:**

| control | cohort | D150 result (pooled) | prediction if timing is the explanation |
|---|---|---|---|
| anti-estrogens (fulvestrant, tamoxifen, raloxifene, clomifene, bazedoxifene, toremifene) | BRCA-LumA | 0.9th–12.7th opposing pct | **high early, low late** |
| BRAF inhibitors (vemurafenib, dabrafenib) | TCGA-SKCM | 50.5th / 43.3rd | no registered prediction — D150 gives no timing evidence for BRAF |
| anti-estrogens | BRCA-Her2 | 76th–96th (on the WRONG network) | no prediction; the point is that the old number was invalid |

`opposing percentile` = the within-(cohort, window) percentile of `−cos_weighted`, so higher means
the drug opposes the disease pattern more than other profiles scored in the same cell.

## 4. Acceptance — fixed here

This is a **development check on a known answer**, not a validation. It is registered so the result
cannot be reinterpreted afterwards.

1. **The stated success condition:** anti-estrogens in BRCA-LumA reach a **higher median opposing
   percentile in `early` than in `late`**, with the early median **above 50**.
2. Reported alongside, not gated: the same quantity for BRAF in SKCM, the corrected BRCA-Her2
   numbers, and the per-cohort shift from D150's pooled values.

**Pre-committed:** if (1) fails, the conclusion is that timing does **not** rescue D150's negative
result, and that is written up as prominently as a pass. Nothing about the drug set, the cohort
set, `MIN_GENES`, or the window boundaries is tuned afterwards. A cohort with no assayed matched
line is reported as **"insufficient matched evidence"**, which is not the same as a low score.

**What a pass would and would not mean.** It would mean the pooled D150 result was partly a timing
artefact. It would **not** establish that the nested cosine ranks drugs better than the production
reversal score — that is a separate, later test with its own registration.

## 5. Output

```
results/systems_validation/integrated/nested_cosine_window/
  per_drug.tsv.gz    cohort, network, window, pert_id, cmap_name, n_regulators,
                     cos_mean, cos_weighted, frac_opposing, opposing_pct
  controls.tsv       the section-3 controls, per cohort and window
  summary.tsv        per (cohort, window): medians, n drugs, n regulons usable
  run.json           networks resolved, seeds, input mtimes, acceptance verdict
```

Script: `scripts/79af_nested_cosine_window.py`.

## 6. Registration

The git commit adding this file is the registration. Later changes to sections 1–4 go in a new
dated section, never an edit in place.
