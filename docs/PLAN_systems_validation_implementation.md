# Implementation plan — SPEC_systems_theoretic_mechanism_validation

**Written:** 2026-11-18 · **Governs:** `docs/SPEC_systems_theoretic_mechanism_validation.md`
**Audit of record:** `docs/AUDIT_systems_theoretic_spec.md` · **Decisions:** D137–D142

This plan is ordered by **dependency and evaluability**, not by the spec's section numbers. Each step
states its inputs, its output artifact, and the **exit criterion that decides whether the next step
is allowed to start**. Steps marked ⛔ are not attempted, with the measurement that rules them out.

---

## Already done — do not re-plan

| | artifact | result |
|---|---|---|
| Dependencies | `environment.sysval.lock.yml` | CARNIVAL 2.20.0, CBC, lpSolve, sympy, scipy, networkx |
| Gate F0 | `gate_f0_carnival.json` | **4/4**, matches the package's own expected fixture at sign agreement 1.000 |
| Gate F1 | `gate_f1_adversarial.json` | 1 pass, 1 fail (ambiguity), 2 qualified — drove D141 |
| Gate F2 | `gate_f2_enum_*.json` | **4/4** — enumeration recovers known ambiguity counts (D141) |
| Gate O0 | `gate_o0_fixtures.json` | **10/10** |
| Gate S0 | `gate_s0_angulo.json` | **8/8** |
| Frozen PKN | `pkn/signed_pkn.tsv` | 15,646 edges, 5,264 nodes; SIGNOR separate (D140) |
| Entry set | `integrated/candidate_table.tsv` | 1,399 pairs, per-avenue coverage and limits (D142) |

**All mathematical replication gates pass.** The programme is cleared to attempt biology.

---

## ⛔ Not attempted, and why

| arm | measurement that rules it out |
|---|---|
| **S biological** (attractors, persistence) | 317 compound×line cells have ≥3 timepoints at production quality; `siginfo_beta` has **no washout/recovery field at all**. Persistence is not evaluable by construction (AUDIT §3.1). S-math is done and passes. |
| **O5 formal control** (Kalman on real data) | identifying `A`, `B`, `C` needs dense time-series interventions; same data limit. Spec's Gate O3 already forbids it. |
| **CPLEX solution pooling** | access refused (site number + business justification). Replaced by in-house enumeration, validated 4/4 (D141). |

---

## Phase 1 — Experiment F, biological pilot

The first step where the pipeline could learn something it does not already know.

### F-1. Build the frozen TF-activity inputs
**In:** `results/step6/aracne_project_v2_tf/group_mr.tsv` (disease TF activity per cohort)
**Out:** `results/systems_validation/fault_isolation/tf_activity_<cohort>.tsv` + manifest with sha256
**Do:** export per-cohort TF activities in CARNIVAL's `measurements` form; intersect with PKN nodes
**first** and record the dropped-TF count (AUDIT C3 — CARNIVAL aborts with a misleading "none of your
measurements are in PKN" when any node is missing).
**Exit:** every pilot cohort has ≥20 TFs present in the PKN, with the drop count recorded.
**Effort:** 0.5 day.

### F-2. Pre-register the inverse-mode driver-recovery test ← **FREEZE POINT**
**Do:** write and commit, *before any CARNIVAL biological run*: the candidate driver universe per
cohort, which drivers are supplied vs withheld, the ILP complexity weight (`betaWeight`), solver,
time limit, seed, and what counts as recovery. Commit with a checksum.
**Why it must be first:** the spec's §10.3 forbids a supplied driver counting as its own recovery,
and the smoke test (D140) returned **561 nodes and 2,292 edges from 5 measurements** — without a
frozen complexity weight the first result will be an impressively large network that explains
everything.
**Exit:** a committed prediction artifact. No CARNIVAL run on real data before this commit.
**Effort:** 1 day. **Needs a user decision:** the pilot cohort set (see Decisions below).

### F-3. Run inverse CARNIVAL + enumeration on the pilots
**In:** F-1 activities, frozen PKN, F-2 registration
**Out:** per cohort — the optimal subnetwork, **all alternative optima** (D141 enumerator), core
nodes shared by every explanation, and the median pairwise Jaccard.
**Exit:** for each pilot, a recovered-driver call *and* an ambiguity count. **A cohort whose
enumeration hits `--max-solutions` is reported as a lower bound, never as "the number".**
**Effort:** 1–2 days. Runtime grows with solutions found.

### F-4. Negative controls — **before any positive is reported**
**Do:** the spec's §4.6 set, with the D137 correction: randomise the **drug–target–line triple
together**, not just the cohort (my NC3 preserved a pharmacologically compatible line and closed).
Plus: genotype-incompatible context, degree-matched decoy targets, and the SIGNOR network as the
alternative-representation arm.
**Exit:** a false-positive rate on the controls. **If controls close at a rate near D137's 3/5, F
stops here and is reported as a negative.**
**Effort:** 1–2 days.

---

## Phase 2 — Experiment O, empirical arms

These have **written specs but were never run past S0 feasibility**. They are independent of F and
can proceed in parallel.

### O-1. D134 S1 — target recoverability
`docs/SPEC_target_recoverability.md` exists; `d134_s0_run.json` gives **71 testable units**.
**Open decisions carried from D134:** strength floor, chance pool, and the 16 guaranteed-retrieval
cells. **Exit:** chance-adjusted retrieval with the Morgan/Tanimoto stratification. **1–2 days.**

### O-2. D135 S1 — drug vs genetic perturbation
`docs/SPEC_drug_genetic_concordance.md` exists; **2,350 CRISPR-paired cells, 46–121 after stringent
exclusions**. **Open:** strength floor, positive-control essentiality, polypharmacology handling.
**Exit:** concordance with the asymmetry rule — a positive is strong evidence, a negative is weak.
**1–2 days.**

### O-3. D136 S0/S1 — dose response
`docs/SPEC_dose_response_bmd.md` exists, **never started**. 8,778 cells with ≥3 doses — by far the
best-powered evidence route in the project. **Exit:** BMD profiles across `k`, reported as profiles
not point values. **2–3 days.**

> **Sequencing note:** O-3 has ~10× the coverage of O-1 and O-2 and needs no new specification work.
> If only one O arm is run, run O-3.

---

## Phase 3 — Dossiers and integration

### P-1. Candidate dossiers
**In:** `integrated/candidate_table.tsv` (built), plus whatever F and O return
**Out:** `candidate_dossiers.jsonl` per SPEC §0.1
**Rule carried from D142:** the layered score and reversal are **uncorrelated** (ρ = −0.0002, top-50
overlap 0/50). The dossier reports both with their limits and does **not** combine them.
**Exit:** every dossier has a disposition — supported / preliminary / contradicted / not evaluable —
and no dossier claims a known use is novel (`labels.tsv` enforces this mechanically).

### P-2. Per-condition VIPER (only if O-3 runs)
`79c` collapses to one profile per (drug, line) at TAS ≥ 0.4, discarding dose and time. §0.1's atomic
record needs its own activity computation from Level 5. **Real cost, not in the spec's §9 tree.**

---

## Decisions needed from the user

1. **Pilot cohorts for F-2.** The spec proposes ESR1/BRCA-LumA, BRAF-V600/SKCM, EGFR/LUAD,
   ERBB2/BRCA-Her2. Audit findings: **EGFR/LUAD is weak** (`dep_z = −0.025`, D percentile 47.7 — no
   dependency signal), and **ESR1/BRCA-LumA inherits a known-broken input** (its `o_source` is `NES`,
   the one place the regulon enters Layer 1, and that regulon reads anti-estrogens backwards in MCF7
   per D83/D85 and is network-unstable at Jaccard 0.011 per D138.1). BRAF/SKCM and ERBB2/BRCA-Her2
   are clean (D ≥ 99.99, A = 100).
2. **Which O arm first** — recommendation above is O-3.
3. **What a dossier does to the ranking** — still open after D142; the measured answer is that they
   cannot be reconciled by weighting.

## Suggested order

**F-1 → F-2 (freeze) → F-3 → F-4**, with **O-3** in parallel from the start. Phase 3 only after at
least one of F-4 or an O arm returns. Nothing in Phase 1 or 2 depends on decision 3.

---

## Revisit after implementation — open caps and constraints

Recorded as we hit them, to be reopened once F, O and the dossiers are built. **None of these
blocks implementation; all of them bound what the finished thing can claim.**

| # | cap | measured | why revisit |
|---|---|---|---|
| **R-1** | **PKN coverage is study-biased** (D143, user flagged 2026-11-18) | only **33.0%** of 1,420 ARACNe regulators are PKN nodes; retained vs dropped are indistinguishable on activity (median \|NES\| 1.003 vs 0.970, significance 12.8% vs 11.2%) but differ entirely on curation (median effort 16 vs 0; **86% of dropped regulators have zero curation anywhere in OmniPath**). Only 35.7% of \|NES\| mass survives. | F can only reason about already-studied TFs, so a recovered driver is partly guaranteed by the PKN encoding curated biology. **Once F has run, revisit whether the result is explained by curation depth** — the fifth failure condition tests exactly this. Options to weigh then: restrict the candidate universe to well-curated genes so comparisons are like-for-like; add a non-curated network arm; or accept the cap and scope every F claim to the curated subnetwork. |
| **R-2** | reversal and the layered score are uncorrelated | ρ = −0.0002, top-50 overlap 0/50 (D142) | revisit once F and O give a third, independent read — that may indicate which of the two, if either, tracks mechanism |
| **R-3** | S biological is not evaluable | 317 cells with ≥3 timepoints; no washout field exists | revisit only if new longitudinal data arrives; nothing in LINCS can fix it |
| **R-4** | ambiguity enumeration is capped | `--max-solutions 25`; a cap hit is a lower bound | revisit the cap once typical optima counts on real cohorts are known |
