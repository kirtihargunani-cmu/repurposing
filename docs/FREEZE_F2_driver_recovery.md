# F-2 PRE-REGISTRATION — inverse-CARNIVAL driver recovery

**Frozen:** 2026-11-18, before any CARNIVAL run on real cohort data.
**Governs:** `docs/PLAN_systems_validation_implementation.md` Phase 1 step F-3 and F-4.
**Spec:** `SPEC_systems_theoretic_mechanism_validation.md` §4.5 (freeze forward predictions), §4.8
(F-MATH), §10.3 (no annotation leakage), §10.5 (freeze procedures).

**The commit containing this file is the freeze.** Any F-3 result produced from a run that predates
this commit, or that departs from the settings below, is void and must be re-run.

---

## 1. Pilot cohorts and the driver each test asks about

| cohort | driver | status | basis |
|---|---|---|---|
| TCGA-SKCM | **BRAF** | primary | D 99.99, A 100.0 — clean |
| BRCA-Her2 | **ERBB2** | primary | D 99.99, A 100.0 — clean |
| BRCA-LumA | **ESR1** | **seen development case** | spec §13 names it; its `o_source` is `NES`, the one place the regulon enters Layer 1, and that regulon reads anti-estrogens backwards in MCF7 (D83/D85) and is unstable at Jaccard 0.011 (D138.1). **Cannot count as independent evidence.** |
| TCGA-LUAD | **EGFR** | **deliberately weak case** | `dep_z = −0.025`, D percentile 47.7 — no dependency signal. Included to see what the method does when the driver evidence is thin. |

Structural feasibility, verified before freezing (so a failure cannot be blamed on unreachability):
all four drivers are PKN nodes with outgoing edges, and each reaches a share of its cohort's measured
TFs within 3 hops — BRAF 12/42, ERBB2 41/73, ESR1 23/74, EGFR 19/30.

## 2. Measurements (TF activities)

- Source: `results/systems_validation/fault_isolation/tf_activity_<cohort>.tsv` (F-1).
- **Selection rule, frozen:** regulator is in the PKN **and** `fdr < 0.05`. Counts: SKCM 42,
  Her2 73, LumA 74, LUAD 30.
- **Value passed:** the continuous VIPER `nes`, unaltered. CARNIVAL 2.20.0 discretises signs
  internally for its ILP constraints; we do not pre-binarise, and the installed version is recorded
  in every run manifest.
- Measurements are intersected with PKN nodes in F-1, never inside the run.

## 3. Candidate driver universe

`alt_share >= 0.05` **or** `dep_z > 2` in that cohort (from `gene_components.tsv.gz`), **and** present
in the PKN with at least one outgoing edge. Sizes: SKCM 74, Her2 548, LumA 358, LUAD 258.

## 4. Mode — and what may be claimed from each

- **PRIMARY: inverse mode** (`runInverseCarnival`). **No driver is supplied.** Only this arm may
  support a driver-recovery claim.
- **SECONDARY: vanilla mode** (`runVanillaCarnival`) with the driver supplied. Reported for
  downstream-prediction purposes only. **Supplying BRAF and recovering BRAF is not recovery** and
  will not be described as such.

## 5. Solver settings, frozen

| setting | value |
|---|---|
| solver | CBC 2.10.13, `solverPath` absolute (`Sys.which` fails — AUDIT C4.2) |
| `betaWeight` (complexity penalty) | **0.2** (CARNIVAL default) — primary |
| sensitivity arm | `betaWeight` ∈ {0.05, 0.5}, reported separately, never substituted for the primary |
| `timelimit` | 600 s per solve |
| seed / `threads` | threads 1 (determinism over speed) |
| CARNIVAL | 2.20.0, env `sysval`, `environment.sysval.lock.yml` |

The complexity penalty is frozen because the D140 smoke test returned **561 nodes and 2,292 edges
from 5 measurements**; tuning it after seeing results would let us choose how much mechanism to find.

## 6. Ambiguity

Every inverse run is followed by `enumerate_optima.py` (D141, validated 4/4) with
`--max-solutions 25`. Reported per cohort: the number of alternative optima, the **core nodes present
in every optimum**, the union, and the median pairwise Jaccard. **A run that hits the cap is reported
as `exhausted: false` and its count is a lower bound, never "the number of explanations".**

## 7. What counts as recovery — stated before any result is seen

| call | definition |
|---|---|
| **STRONG** | the driver is active, with the expected sign, in **every** alternative optimum |
| **WEAK** | the driver is active in **at least one but not all** optima |
| **NOT RECOVERED** | the driver is active in no optimum |
| **NOT EVALUABLE** | solver infeasible, timed out, or the cap was hit without exhaustion |

Expected sign: for an activating oncogenic driver (BRAF, ERBB2, EGFR) the inferred node activity is
**positive**; for ESR1 in a luminal tumour, **positive**. A recovery with the wrong sign is scored
**NOT RECOVERED** and reported as a sign failure.

## 8. Degree control — mandatory, because two pilots are hubs

ERBB2 and EGFR reach 1,908 and 2,382 nodes within 3 hops. For each cohort, draw **10 decoy
candidates** matched to the true driver on out-degree (±20%) and 3-hop reach (±20%), excluded from
the candidate universe for that draw, and run the identical inverse pipeline. **If decoys are
recovered at a rate comparable to true drivers, F reports that the method tracks degree, not
causality** — the §4.7 failure condition.

## 9. Pre-stated failure conditions for the whole F arm

F is reported as a **negative result**, not iterated, if any of these holds:

1. Degree-matched decoys are recovered as often as true drivers (§8).
2. The median number of alternative optima is so large that the core node set is empty or trivial —
   i.e. every explanation is different, so no explanation was found.
3. Recovery depends on `betaWeight`: the primary and sensitivity arms disagree on the call.
4. The F-4 negative controls close at a rate near D137's 3/5.

## 10. What F-3 may NOT read

The held-out perturbation outcomes used in F-4, any CTRP/PRISM/GDSC label, and
`results/lincs/d128/a4/*`. F-3 writes its predictions before F-4 reads any observation.

---

**Machine-readable copy:** `results/systems_validation/fault_isolation/f2_registration.json`
(same values; the JSON is what the run scripts read, so the two cannot drift).
