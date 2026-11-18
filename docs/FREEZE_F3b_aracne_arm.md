# F-3b PRE-REGISTRATION — full-ARACNe network as the alternative PKN

**Frozen:** 2026-11-18, before any CARNIVAL run on an ARACNe-derived PKN.
**Compares against:** F-3 on the curated PKN (D144: 0 of 4 genuine recoveries).
**Prompted by:** the user's observation that the *full* ARACNe regulator set contains kinases and
signalling proteins, not only TFs — which corrected my claim that ARACNe cannot supply upstream nodes.

## 1. What this arm tests

F-3 failed on the curated map. R-1 says that map is study-biased — only 33% of ARACNe regulators are
in it, and what separates the included from the excluded is curation effort, not activity. This arm
asks whether a **data-derived** map, with no curation bias, recovers the drivers that the curated
map missed.

## 2. The network

`data/processed/step5/regulons_aracne_project_v2.rds` — the **full** regulator set (TFs **and**
kinases and signalling proteins), 29 networks. Confirmed present as regulators with outgoing edges in
SKCM: BRAF (62 targets), ERBB2 (81), EGFR (133), AKT1, MAPK1, PIK3CA, MTOR, KRAS, NRAS, SRC, CDK4.

Edges: regulator → target, sign = `sign(tfmode)`. **`tfmode` is the sign of a co-expression
correlation across that cohort's tumours, not a measured direction of regulation**, and this arm is
labelled as such everywhere it is reported.

| cohort | ARACNe network | edges | drivers |
|---|---|---|---|
| TCGA-SKCM | SKCM | 291,570 | BRAF |
| BRCA-Her2 | BRCA (pooled; no Her2-specific full network) | 234,307 | ERBB2 |
| TCGA-LUAD | LUAD | 341,218 | EGFR |

BRCA-LumA/ESR1 is **excluded**: no `BRCA_LUMA` network exists in the full set, and ESR1's F-3 result
was circular anyway.

## 3. Measurements — unchanged, deliberately

The identical frozen F-1 TF activities used in F-3 (in curated PKN ∧ FDR < 0.05; 42 / 73 / 30 TFs).
**Only the map changes.** Keeping the measurements fixed is what makes this a controlled comparison;
all measured TFs are present in every ARACNe network tested (100% at every pruning level).

## 4. Pruning — none

Measured before freezing, per-regulator top-k by likelihood (k = 5, 10, 20, 50, all):

| k | SKCM edges | BRAF reaches | of 42 measured TFs |
|---|---|---|---|
| 5 | 22,790 | 6 nodes | **0** |
| 10 | 44,774 | 57 | **0** |
| 20 | 85,319 | 361 | **0** |
| 50 | 177,751 | 2,337 | 5 |
| **all** | **291,570** | **4,808** | **10** |

Global top-N-by-likelihood pruning was tested and **rejected**: matching the curated PKN's 15,646
edges leaves only 3 of 42 measured TFs in the network and strips EGFR of all outgoing edges — it
would rig the comparison.

**Any choice of k is a tuning choice**, so this arm uses the network as it is, unpruned. If the
solver does not return within the stated limit, the arm is reported **not evaluable**, not re-tuned.

## 5. The structural finding, recorded before running

In the ARACNe network the driver reaches the measured TFs **only when it also reaches nearly
everything**:

| map | driver reaches | measured TFs | selectivity |
|---|---|---|---|
| curated PKN | 561 nodes | 12 of 42 | **2.1%** |
| ARACNe, k=50 | 2,337 nodes | 5 of 42 | 0.21% |
| ARACNe, all | 4,808 nodes | 10 of 42 | 0.21% |

**The curated map is ~10× more selective.** There is no regime in which BRAF connects specifically to
the damaged TFs without also connecting to thousands of other nodes. This is predicted to show up as
**high ambiguity** — many equally good explanations — which the D141 enumerator will measure.

## 6. Settings, and what counts as recovery

Identical to F-2 except the network: CBC, `betaWeight` 0.2, `timelimit` 1800 s (raised from 600 for a
19× larger network; stated here, not after seeing a timeout), threads 1, enumeration capped at 25.
Recovery definitions (STRONG / WEAK / NOT RECOVERED / NOT EVALUABLE) and expected signs are
unchanged from `FREEZE_F2_driver_recovery.md` §7.

## 7. Pre-stated readings

| outcome | reading |
|---|---|
| ARACNe recovers drivers the curated map missed | the curated map's study bias (R-1) was the binding constraint |
| neither map recovers them | the failure is not about which map; it is about the approach |
| ARACNe "recovers" but with very high ambiguity | consistent with §5 — recovery by universal reachability, which is not evidence |
| solver does not return | **not evaluable**; do not prune to make it finish |

**Circularity note, carried into every result:** the TF activities fed in are computed by VIPER over
ARACNe regulons, so in this arm the map and the measurements share a data source in a way the curated
map did not. A recovery here is weaker evidence than the same recovery on the curated map.
