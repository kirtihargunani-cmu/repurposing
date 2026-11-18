# Audit: SPEC_systems_theoretic_mechanism_validation.md

**Auditor:** Claude (agent) · **Date:** 2026-11-17 · **Spec commits:** `d82661e`, `5f315db` (725 lines)
**Method:** every factual claim in §2 checked against the artifact it cites; every declared dependency
checked for installation; every data-feasibility premise recomputed from raw LINCS metadata.
**No experiment was run and no spec claim was adopted.**

---

## 1. Verdict in one paragraph

The spec is scientifically sound and its citations are **accurate** — all nine numeric claims in §2
reproduce exactly from the artifacts. Its statistical safeguards (§10) are the strongest in the
project and correctly bind the known failure modes, including one of my own (§4.6 on the NC3
confound). Three things block it as written: **Experiment S's biological arm is not merely gated but
structurally impossible on LINCS** (§3.1 below); **none of the five declared software dependencies is
installed** (§3.2); and **the spec never states that the production lead ranking is insensitive to
the regulon**, which undercuts its own framing (§3.3). Recommended sequencing in §6.

---

## 2. Claim-by-claim verification of §2

| § | Claim | Source checked | Result |
|---|---|---|---|
| 2.1 | 3 of 5 negative controls closed (60%) | `false_positive_rate.json` | **✓ exact** (3/5, 0.6) |
| 2.1 | Neither positive trace closed | `verdicts.json` | **✓** (both `chain_closes: false`) |
| 2.1 | ESR1 median pairwise Jaccard 0.011 over 41 networks | `network_stability.tsv` | **✓ exact** |
| 2.3 | 2,350 CRISPR-matched cells over 7 lines | `d135_s0_run.json` `bottom_line.trt_xpr` | **✓ exact** (`paired_cells` 2350, `lines` 7) |
| 2.3 | 46–121 cells survive stringent exclusions | same, `pert_nz2>=100` | **✓ exact** (46 single-target, 121 ≤2-target) |
| 2.3 | D134: 71 target units eligible | `d134_s0_run.json` | **✓ exact** (`units_testable: 71`) |
| 2.4 | 8,778 cells ≥3 doses; 2,736 ≥2 timepoints | `SPEC_dose_response_bmd.md` §7 | **✓ exact** |
| 2.5 | Production drug matrix = 73 cell lines | `d134`/`d135` `profiles.lines` | **✓ exact** |
| 2.5 | D129 inferred-gene reconstruction is strong | `settled_measurement_facts.md` | **✓** (median R² 0.875 vs 0.630 control) |

**No factual error found in §2.** This is unusual and worth stating plainly.

---

## 3. Blocking findings

### 3.1 Experiment S's biological arm is impossible, not merely gated — **the most important finding**

§6.6 and S-MATH-5 gate the biological S arm on longitudinal feasibility, and §6.8 makes persistence
conditional on washout data. Recomputed from `data/raw/lincs2020/siginfo_beta.txt` (720,216 `trt_cp`
signatures, 201,494 compound × line cells):

| timepoints per (compound, line) | all cells | at `is_hiq` & TAS ≥ 0.4 |
|---|---|---|
| ≥ 2 | 47,216 | 1,811 |
| **≥ 3** | **2,038 (1.0%)** | **317** |
| ≥ 4 | 226 | 28 |
| ≥ 5 | 24 | — |

97% of all signatures sit at 6 h or 24 h. **At production quality only 317 cells have three or more
timepoints**, and a continuous-time `A` cannot be identified from two. Separately, `siginfo_beta`
has **no washout, recovery or drug-removal field at all** — the columns are `pert_time`, `pert_dose`
and `pert_type` only. LINCS measures continuous exposure to a fixed endpoint.

**Consequence:** S-H3 (persistence) is **not evaluable by construction, not by sample size**, and no
amount of analysis changes that. The spec's own §6.8 rule then forces `persistence_status =
not_evaluable` for every compound. S-MATH-1/2 (the Angulo replication and Lie-bracket work) remain
fully doable — they are pure mathematics and need no data.

**Recommendation:** split S explicitly into **S-math (do it)** and **S-bio (declare not evaluable now,
with these numbers, rather than at Stage 4)**. Reaching this conclusion costs nothing today and saves
the entire Stage 4 audit. If temporal work is wanted, the realistic question is not attractors but
the 6 h vs 24 h divergence already flagged in `SPEC_dose_response_bmd.md` §8 — which has 1,811 usable
cells, not 317.

### 3.2 Every declared dependency is missing

| Dependency | Required by | `repurposing` env | `celligner` env |
|---|---|---|---|
| **CARNIVAL** | F-MATH, Gate F0 | not installed | not installed |
| **CBC solver** | F-MATH §2 (preferred) | not on PATH | not on PATH |
| **lpSolve** | F-MATH §2 (toys) | not installed | not installed |
| **SymPy** | S-MATH-2 symbolic brackets | not installed | not installed |
| **SciPy** | O-MATH, S-MATH | **not installed** | installed 1.9.1 |
| NetworkX | O-MATH-1 matching | installed | installed |

Note the env split this exposes: the drug-side Python is **`celligner`** (pandas 1.4.4) and R is
**`repurposing`** — `scipy` exists only in `celligner`, CARNIVAL would have to go in `repurposing`.
Any new script must declare which interpreter it runs under; the spec's §9 tree does not.

**Gate F0 cannot run until CARNIVAL and a solver are installed.** This is a half-day of environment
work, not research, and it should be Stage 0 item one.

### 3.3 The spec's framing omits today's most relevant result

The spec is built on the premise that regulatory-network mechanism is what the pipeline needs to
validate. Measured today (D139, ARACNe-only chain, merged 1:1 on `group_id` × `pert_id`):

| score | Spearman, ARACNe-only vs consensus |
|---|---|
| `rev` | 0.9905 |
| `spec_mech_pct` | 0.9866 |
| **`layered_pct` (the lead ranking since D122)** | **1.0000 — bit-identical, all 41 groups** |

`layered_pct` is built from DepMap dependency, TCGA alteration and ChEMBL binding; the regulon enters
only through the `O` component where `o_source == NES`, which is **4.9% of `gene_components` rows**.

**Consequence for the spec:** F, O and S all validate mechanism through regulatory networks, but the
production ranking barely consults them. Whatever F/O/S conclude, the current lead score will be
unchanged unless the spec also proposes how validated mechanism should *enter* the ranking. §7 and
§15 explicitly decline to produce a score — which is defensible — but then §8 Stage 6's "no
production ranking or tier changes without a separate validated proposal" leaves no path by which
this work ever reaches the output. **This should be stated as a known structural gap in §2.**

---

## 4. Non-blocking corrections and gaps

1. **§2.3 wording.** "D135's feasibility audit identified 2,350 CRISPR-matched cells" is right, but the
   artifact's own note says the clean-positive counts are an **upper bound** because limitation 2
   (non-catalytic scaffolding roles) is *not computable* — no curated list exists in the repo. The
   spec should carry that upper-bound caveat where it quotes 46–121.
2. **§5.3 / O2 sampling unit.** D135-S0 found the consensus-shRNA arm is the weakest-responding by ~6×
   (median `nz2` 14 vs 51 CRISPR, 88 raw shRNA); the spec correctly makes CRISPR primary but should
   note shRNA's floor must differ or that arm is empty.
3. **§4.2 F1 pilot feasibility.** Checked the four proposed pilots. Three are well-supported;
   **EGFR/LUAD is not**: `dep_z = -0.025`, D percentile 47.7 — no dependency signal, so an "upstream
   cause" there rests on alteration share (12.8%) alone. ESR1/BRCA-LumA D 99.98 / A 45.0 / O 94.3,
   BRAF/SKCM D 99.99 / A 100, ERBB2/BRCA-Her2 D 99.99 / A 100 are all strong. Consider replacing the
   EGFR pilot or declaring it the deliberately weak case.
4. **§4.2 ESR1 pilot caveat.** ESR1's `o_source` is `NES` — the one place the regulon enters Layer 1.
   So the ESR1 pilot is the single F case whose own disease-side evidence is regulon-dependent, and
   D138.1 showed that regulon is unstable (Jaccard 0.011) and reads anti-estrogens backwards in MCF7
   (D83/D85). Flag it, or the pilot inherits a known-broken input.
5. **§9 repository structure** omits the `run_carnival.R` adapter that F-MATH §1 requires, and gives no
   `tests/` entry for the math fixtures that Gates F1/O0/S0 depend on.
6. **§3.4 join grain.** "All joins must be audited at the correct grain" is right and today supplied a
   live example: the split tables are unique on `(group_id, pert_id)` but **not** on
   `(group_id, cmap_name)` — 3,762 duplicates, which silently inflated a 50,760-row comparison to
   60,120. Worth naming as the canonical key in the spec.
7. **Unsuffixed outputs.** `79t3` writes `layered_score_v2.tsv` with no arm tag, so any alternative arm
   silently overwrites production (fixed today under `SPLIT_SFX`). The spec's provenance requirements
   should extend to "every output file names its arm".

---

## 5. What is genuinely strong

- §10 safeguards 1–9 are correct and bind the project's actual historical failures: seen-vs-validation
  (D127/D128), dependence-preserving nulls (D130), annotation leakage, sampling unit, freeze
  discipline, label-blind calibration, negative-controls-first, no cross-platform arithmetic.
- §4.6 independently catches the NC3 confound I reported — the randomised cohort kept a
  pharmacologically compatible line — and requires randomising the **drug–target–line triple**.
- The source-method vs biological-adaptation split (§8 checkpoint table) is the right structure and
  directly prevents the "relabel D128 as CARNIVAL" failure.
- §11's edge-case register correctly encodes the inferred-gene rule: invalid for a gene's own
  knockout self-QC, not grounds to reject inferred genes generally.

---

## 6. Recommended sequencing and effort

| # | Work | Blocks | Effort |
|---|---|---|---|
| 1 | Install CARNIVAL + CBC (+ SymPy, SciPy in the right envs); record versions | Gate F0, all math gates | ~0.5 day |
| 2 | Declare **S-bio not evaluable** using §3.1's numbers; keep S-math | removes Stage 4 | ~0.5 day |
| 3 | Gate F0 — CARNIVAL's own standard + inverse toys | F entirely | 1–2 days |
| 4 | Gate S0 — Angulo invariant, symbolic + numeric | S-math | 1–2 days |
| 5 | Gate O0 — Kalman/matching/Gramian fixtures (table in O-MATH) | O-math | 1–2 days |
| 6 | Export frozen `signed_pkn.tsv` from `data/raw/omnipath/omnipath_signed_directed.tsv.gz` | F biological | 1 day |
| 7 | Stage 0 audit proper + resolve §3.3 (how mechanism reaches the ranking) | everything downstream | 2–3 days |

Items 3–5 are independent and parallelisable. **Nothing biological should start before item 7**, since
§3.3 determines whether the whole programme can change anything.

---

## 7. Single recommendation

Adopt the spec, with three amendments: **(a)** split S into S-math and S-bio and mark S-bio not
evaluable now, with the 317-cell and no-washout evidence; **(b)** add §3.3's layered-score finding to
§2 as a constraint, and state how validated mechanism is intended to reach the ranking; **(c)** make
Stage 0 item one the dependency installation, since no math gate can run today.

The spec's own standard — *"a method that fails its mathematical replication is not eligible for
biological inference"* — currently excludes all three experiments, because none of the three source
packages is installed. That is a fixable, half-day problem, and it is the only thing standing between
this document and its first real gate.

---

# Addendum — audit of the v2 amendment (`c440c52`, 2026-11-17 23:03)

**Change:** +94 / −15 lines. A new **§0** makes the primary deliverable *falsifiable new-indication
drug-repurposing proposals with explicit why/how*, adds a structured candidate dossier (§0.1), maps
F/O/S onto it (§0.2), adds a drug-side ARACNe arm (§0.3) and validation criteria (§0.4), plus
nomination clauses §4.7a / §5.7a / §6.9a and a rewritten Stage 6.

## A1. Effect on the three blocking findings

| Finding | Status after v2 |
|---|---|
| **§3.3 no path to output** | **Resolved.** The dossier *is* the output; §0.4 states "no production drug rankings are changed under this draft". My objection assumed the work had to reach `layered_pct`; v2 answers that it does not. The `layered_pct` insensitivity remains worth recording in §2 as context, but it is no longer a structural gap. **One question is now open instead:** the dossier and the production ranking coexist with no stated relationship — does a supported dossier supersede, annotate or ignore a `layered_pct` rank for the same pair? |
| **§3.1 S-bio impossible** | **De-blocked, not answered.** §0.2 makes S "conditional on longitudinal data and not a hard requirement for every candidate", and §6.9a marks S `not evaluable` "without automatically rejecting otherwise supported disease–drug mechanisms". That is the right structural fix. The measured numbers (317 quality cells with ≥3 timepoints; **no washout field exists in `siginfo_beta`**) should still be written into §6, so nobody re-audits this. Severity drops from blocking to documentation. |
| **§3.2 dependencies missing** | **Unchanged and still blocking.** CARNIVAL, CBC, lpSolve and SymPy are installed in neither env; SciPy only in `celligner`. No math gate can run. |

## A2. New finding — the new-indication requirement is feasible, and here is its denominator

§0 requires a drug with a *different* established indication, proposed for a *new* disease, with
verified source-line context. Computed from `results/lincs/known_use/labels.tsv` (7,570
cohort × FDA-ingredient rows, 41 cohorts, 391 ingredients), where `lab_hub` / `lab_chembl` mark
existing indications:

| filter | pairs remaining |
|---|---|
| all cohort × FDA-ingredient rows | 7,570 |
| minus known uses (Hub **or** ChEMBL) | 7,099 (−471) |
| and gate-eligible | 5,408 |
| **and source line matched, not borrowed** (§0.1) | **1,341** |
| and the drug is not already an anticancer | 618 |

**The 1,341 matched candidates span 266 drugs and 29 cohorts** — a workable denominator. The
618 non-anticancer subset (164 drugs, 20 cohorts) is the strongest repurposing case.

Two cautions the spec should carry:

- **Coverage is highly skewed:** median 5 candidates per cohort, range 1–132, and **12 of 41 cohorts
  have zero** matched non-known-use eligible pairs — SARC-DDLPS/-LMS/-UPS, TARGET-WT, TCGA-ACC,
  -HNSC, -KICH, -KIRC, -KIRP, -THCA, TGCT-NSGCT and -SEM. For those cohorts the honest §0 output is
  "no candidate meets the entry criteria", and that should be stated up front rather than discovered.
- **471 known-use pairs (97 drugs, 41 cohorts) must be excluded from novelty claims.** `labels.tsv`
  is the enforcement resource for §0's "avoid claiming a known use is novel"; name it in the spec so
  the check is mechanical rather than editorial.

## A3. New finding — the dossier's atomic grain is finer than the production pipeline's

§0.1 defines the atomic record as *new disease/subgroup × drug × target unit × observed source line,
**preserving relevant doses and times***. The production pipeline cannot supply that: `79c` keeps
**one consensus profile per (drug, line) at TAS ≥ 0.4**, discarding the dose and time structure
entirely (noted in `SPEC_dose_response_bmd.md` §7). §0.3 and §6.9a already say to compute VIPER per
condition from Level 5 rather than reuse aggregated `79d` output — correct, and it means the dossier
pipeline needs its **own** activity computation, not a join onto existing outputs. That is a real
implementation cost and is not reflected in §9's file tree.

## A4. Assessment of §0.3 (drug-side ARACNe arm)

Sound, and it states the right caveats. Two additions from D139, finished after the amendment:

- Switching the drug-side regulon from ARACNe+CollecTRI consensus to **ARACNe-only moves `rev` by
  Spearman 0.9905 and `spec_mech_pct` by 0.9866** — so the *set* choice barely matters. The
  instability that matters is across **cohort networks** (D138.1, Jaccard 0.011), which §0.3 already
  requires reporting. The sensitivity arm should therefore vary the *cohort network*, not the
  regulon set.
- The drug-side VIPER floor is now `MIN_TARGETS = 25` (D139, user decision), up from 10. Any dossier
  quoting regulator activity should record the floor, since `drug_viper_specific_v2_tf` (floor 10)
  and `..._ms25` (floor 25) both exist on disk.

## A5. Revised recommendation

The v2 amendment is an improvement: it gives the programme a deliverable that does not depend on
changing a ranking that is insensitive to the evidence being gathered, and it correctly demotes S
from blocker to conditional. **Two of my three blocking findings are now closed.** Remaining, in order:

1. **Install CARNIVAL + CBC, SymPy, and SciPy in the right envs** — still the only thing preventing
   any gate from running (~0.5 day).
2. **Write §A2's funnel into §0** as the entry criterion, naming `labels.tsv` as the known-use filter
   and listing the 12 zero-candidate cohorts.
3. **Record §3.1's temporal numbers in §6** so S-bio's status is evidence-backed rather than assumed.
4. **State the dossier↔production-ranking relationship** (§A1) — the one question v2 opens.
5. Add the per-condition VIPER computation (§A3) to §9's file tree and effort estimate.

---

# Addendum 2 — dependencies installed, Gate F0 **passed** (2026-11-18)

Blocking finding §3.2 is **closed**. All five missing dependencies are installed and the spec's first
real gate now runs.

## B1. Isolated environment, not the production ones

Installed into a **new `sysval` env**, not `repurposing` or `celligner`. Those two are pinned
(`environment.lock.yml`, `environment.celligner.lock.yml`) and production runs in them; CARNIVAL pulls
a large R dependency tree and CBC is a C++ solver, so adding them risked moving versions production
depends on. Lockfile committed as `environment.sysval.lock.yml` (247 packages).

| package | version | gate |
|---|---|---|
| `bioconductor-carnival` | **2.20.0** | F-MATH / Gate F0 |
| `coin-or-cbc` | 2.10.13 | F-MATH preferred solver |
| `r-lpsolve` | 5.6.23 | F-MATH toy solver |
| `sympy` | 1.14.0 | S-MATH-2 symbolic brackets |
| `scipy` | 1.17.1 | O-MATH, S-MATH-1 |
| `networkx` | 3.6.1 | O-MATH-1 matching |
| `r-base` / `python` | 4.5.3 / 3.11.18 | — |

Note the installed CARNIVAL is **2.20.0**, not the 2.8.0 suggested by the channel's version listing.
This is exactly why F-MATH §2 says to confirm the installed API — see B2.

## B2. Two findings about the package, for the spec

1. **The spec names only CARNIVAL's legacy API.** F-MATH §2 cites
   `runCARNIVAL(inputObj=, measObj=, netObj=)`. In 2.20.0 that function still exists, but the current
   entry points are **`runVanillaCarnival(perturbations, measurements, priorKnowledgeNetwork, weights,
   carnivalOptions)`** and **`runInverseCarnival(measurements, priorKnowledgeNetwork, weights,
   carnivalOptions)`**, with solver configuration through `defaultCbcSolveCarnivalOptions()` /
   `defaultLpSolveCarnivalOptions()`. The spec should name the modern pair; the gate script uses them.
2. **The package's own inverse fixture is stale relative to its own validator.** `toy_network_ex2.RData`
   ships with columns `Source/Effect/Target`, while 2.20.0's input check requires
   `source/interaction/target` — which `toy_network_ex1.RData` already uses. The inverse toy therefore
   **fails out of the box**. The gate renames the columns (format only; edges and signs untouched) and
   records the original names in the manifest. Anyone reproducing Gate F0 will hit this.

## B3. Gate F0 result — **PASS**

`scripts/systems_validation/fault_isolation/gate_f0_carnival_toys.R` →
`results/systems_validation/synthetic/gate_f0_carnival.json`

| run | result |
|---|---|
| vanilla (standard) / lpSolve | OK — 6 nodes, 6 edges, 0.2 s |
| vanilla (standard) / CBC | OK — 6 nodes, 6 edges, 0.1 s |
| inverse / lpSolve | OK — 9 nodes, 14 edges, 0.1 s |
| inverse / CBC | OK — 9 nodes, 14 edges, 0.1 s |

Crucially this is a **comparison, not just an execution check**: the package bundles
`toy_result_expected_ex1.RData`, and both vanilla runs match it on all 6 nodes with
**node-activity sign agreement 1.000**. Both solvers agree with each other and with the fixture.

**Gate F0 of the spec is satisfied.** F-MATH is now eligible to proceed to Gate F1 (synthetic
adversarial fixtures: the four-node signed chain, non-identifiable two-source case, contradictory
observations, disconnected TF, high-degree decoy, hidden cause).

## B4. Revised remaining blockers

| was | now |
|---|---|
| §3.2 dependencies missing | **CLOSED** — installed, isolated, locked, Gate F0 passed |
| §3.1 S-bio impossible | de-blocked by v2; record the numbers in §6 |
| §3.3 no path to output | closed by v2's §0; dossier↔ranking relationship still unstated |

Next, in order: **Gate F1** (F synthetic fixtures), **Gate O0** (Kalman/matching/Gramian table),
**Gate S0** (Angulo invariant). All three are now runnable in `sysval` and are independent of each
other, so they can proceed in parallel.

---

# Addendum 3 — implementation audit: tool calls, installed versions, edge cases (2026-11-18)

Scope: run **every hand-computable fixture the spec asserts**, against the actually installed
packages, and check each API call the spec specifies. A fixture whose stated answer the package does
not produce is a **spec defect**, not a code bug. Scripts under `scripts/systems_validation/`,
manifests under `results/systems_validation/synthetic/`.

## C1. Scoreboard

| gate | fixtures | result |
|---|---|---|
| **F0** CARNIVAL published toys | 4 runs + fixture comparison | **4/4 pass**, sign agreement 1.000 |
| **O0** O-MATH hand-computable table | 10 | **10/10 pass** |
| **S0** Angulo invariant + accessibility | 8 | **8/8 pass** |
| **F1** F-MATH adversarial table | 4 | **1 pass, 1 fail, 2 qualified** |

The mathematics in the spec is **correct wherever it can be checked**. The failures are all in
*behavioural claims about CARNIVAL*, which could only be found by running it.

## C2. The headline defect — ambiguity reporting requires CPLEX, which is not available

The spec repeatedly requires alternative-solution enumeration:

- F-MATH §1: "serialize **all competing feasible/optimal subnetwork solutions** permitted by the solver"
- F-MATH fixture: "Inverse mode reports **alternative explanations / non-identifiability**, not an invented unique driver"
- F endpoints (iv): "**causal ambiguity** (number and diversity of comparably feasible subnetworks)"
- §4.7 failure mode: "produces plausible explanations for negative controls"

Inspected the three option constructors in CARNIVAL 2.20.0:

| solver | pooling fields exposed |
|---|---|
| `defaultLpSolveCarnivalOptions()` | **none** — `lpFilename, outputFolder, workdir, betaWeight, cleanTmpFiles, keepLPFiles, solver` |
| `defaultCbcSolveCarnivalOptions()` | **`poolrelGap` only** |
| `defaultCplexCarnivalOptions()` | `poolrelGap, limitPop, poolCap, poolIntensity, poolReplace` — the full suite |

**Solution-pool enumeration is a CPLEX-only capability.** CPLEX is commercial and absent
(`sh: 1: NA/cplex: not found`). CBC exposes a pool gap but none of the population/capacity controls.

**Confirmed empirically (F1-b).** Fixture: two sources `S1 -> H`, `S2 -> H`, `H -> T1`, `H -> T2`,
with `T1 = T2 = +1` — a case where the driver is *provably* not identifiable. CARNIVAL inverse mode
returned **`S1 = 100`, `S2 = 0`**: it selected one of two indistinguishable explanations and reported
no ambiguity. The spec's required behaviour — "reports alternative explanations / non-identifiability,
not an invented unique driver" — **is not what the package does on the available solvers.**

**This matters because it is precisely the failure the spec was written to prevent.** D137's lesson was
that the procedure generated plausible stories for negative controls; F was designed to measure
causal ambiguity so that would be visible. On lpSolve/CBC that measurement is unavailable.

**Options:** (a) obtain a CPLEX academic licence (free for academics, and the only path to the spec as
written); (b) enumerate alternatives externally — re-solve with no-good cuts excluding each found
solution, which must be built and validated as new code; (c) amend the spec to state that ambiguity is
not measurable on open solvers, and treat every F result as "one of possibly many explanations".
**Not acceptable: silently reporting CARNIVAL's single returned solution as the explanation.**

## C3. The other three F1 fixtures

| fixture | spec's required behaviour | observed | verdict |
|---|---|---|---|
| **F1-a** signed chain `X->A(+), A-\|B, B->C(+)` | activating X gives `A+, B-, C-`; reversing X reverses | `X=+1 → A +100, B −100, C −100`; `X=−1 → A −100, B +100, C +100` | **exact pass** |
| **F1-c** contradictory observations | "report residual disagreement or infeasibility, **never invent a path**" | returned `A=+1, B=0` — no path invented, but **no disagreement reported**; B silently zeroed | **qualified** — the safe half holds, the reporting half does not. Any F pipeline must detect dropped/zeroed measurements itself |
| **F1-d** measured TF disconnected from PKN | "do not fabricate reachability" | **errors out**, refusing to run | **behaviour correct, diagnostic wrong** — see below |

**F1-d's error message is factually wrong and will waste someone's time.** Measuring `{A, Z}` where
`A` is in the PKN and `Z` is not produces: *"None of your measurements are in prior knowledge networks
(PKN)."* Verified that measuring `A` alone succeeds (4-node solution), so the claim "none" is false —
one matched, one did not. CARNIVAL refuses the whole call when **any** measurement is missing, which is
defensible, but it reports the cause incorrectly. **Implication for the dossier pipeline:** measurements
must be intersected with PKN nodes *before* calling CARNIVAL, and the count of dropped TFs recorded,
or runs will fail with a misleading message and the lost coverage will go unlogged.

## C4. API and version corrections for the spec

1. **Entry points** (already in Addendum 2): the spec names only the legacy
   `runCARNIVAL(inputObj=, measObj=, netObj=)`. Current: `runVanillaCarnival(perturbations,
   measurements, priorKnowledgeNetwork, weights, carnivalOptions)` and
   `runInverseCarnival(measurements, priorKnowledgeNetwork, weights, carnivalOptions)`.
2. **`solverPath` is mandatory for CBC and `Sys.which()` will not find it.** The env's `bin` is not on
   `Rscript`'s PATH when R is invoked by absolute path; `Sys.which("cbc")` returns `""` and CARNIVAL
   fails with *"Invalid path to solver provided"*. Pass the absolute path.
3. **The package's own inverse fixture is stale** (Addendum 2 B2): `toy_network_ex2.RData` ships
   `Source/Effect/Target`, the 2.20.0 validator requires `source/interaction/target`.
4. **Installed version is 2.20.0**, not the 2.8.0 the channel listing suggested. The spec's instruction
   to confirm the installed API is what caught items 1–3; keep it.
5. **Three-way environment split is now a real hazard.** `sysval` has **pandas 3.1.1 / numpy 2.4.6 /
   python 3.11**, `celligner` has **pandas 1.4.4 / python 3.10**, `repurposing` has R 4.x and no SciPy.
   Repo artifacts read fine under pandas 3 (spot-checked, no warnings), but any
   `scripts/systems_validation/` file must declare its interpreter in a shebang-style comment, as the
   spec's §9 tree does not. Today's `79i` failure (`ModuleNotFoundError: scipy` under the wrong Python)
   is the live example of this hazard.

## C5. What is confirmed sound

- **All 10 O-MATH fixtures** reproduce exactly: Kalman ranks (2/2, and 1 for the wrong actuator and
  wrong sensor), the structural-matching driver node for the directed path (`[1]`) and the perfect
  matching with no unmatched node for the 3-cycle, unreachability under a non-upstream actuator, the
  ill-conditioning case (full rank with a 1e-8 singular value), and a singular `Wc` that must not be
  inverted. The spec's bipartite construction and its stated answers are correct.
- **All 8 S-MATH fixtures** reproduce exactly. `d(x1·x2)/dt` simplifies **symbolically to 0** under
  x3-only control; with x1 actuation it is **exactly `x2·u1`**; numerically `|x1·x2 − 6| < 3.2e-10`
  across constant, pulse and sinusoidal controls; accessibility rank is **2** with one actuator and
  **3** with two at `x=(2,3,1)`; and `dξ` annihilates the one-actuator distribution. The spec's
  Lie-bracket convention and rank claims hold.

## C6. Revised recommendation

The spec's mathematics is sound and its gates F0, O0 and S0 now **pass**. One substantive amendment is
required before F proceeds to biology:

> **Amend F-MATH to state that causal-ambiguity enumeration requires CPLEX.** Either obtain an academic
> licence, or implement and validate external no-good-cut enumeration, or record that every F
> explanation is "one of possibly many" and drop endpoint (iv). The F1-b fixture as written cannot pass
> on lpSolve or CBC.

Secondary: add the measurement∩PKN pre-intersection (C3), pass `solverPath` explicitly (C4.2), name
the modern entry points (C4.1), and give each `systems_validation` script an interpreter declaration
(C4.5).
