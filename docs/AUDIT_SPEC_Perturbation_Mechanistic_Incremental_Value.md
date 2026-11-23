# Audit — `SPEC_Perturbation_Mechanistic_Incremental_Value.md`

**Audited:** 2026-11-23, against `main` at `9c22f17` and the outputs of the 2026-11-23 full run (13:45-15:17).
**Spec under audit:** commit `748038c` (2026-11-23 15:22), status "proposed pre-registration, not yet registered".
**Method:** every path, field, decision citation, formula and feasibility assumption in the spec was checked against
the repository, the decision log and the data on disk. Counts below are computed from current files, not copied
from the spec or from `DECISIONS.md`.

## Verdict

The protocol is careful and its guardrails are right: no imputation, observed windows only, target-unit
deduplication, quality-only and drug-class placebos, matched rows across models, an evidence-dependence graph, and
a stage that cannot silently disappear. **It is not ready to register as written.** It was drafted against the
repository before D210-D212 (committed four minutes after D211), so it misdescribes what production now is; the
primary outcome of Experiment A has already been examined against these scores; that outcome already has a
committed producer the spec does not know about; and Experiment B's sign rule needs an uncertainty estimate that
90% of conditions cannot supply. Five blocking amendments and five high-priority ones, listed below, should go into a
dated amendment section before the protocol is frozen. Stage 0, A, C, D and E are feasible on current data.

## Blocking — amend before registration

**1. The spec predates the production change (D210-D212).** It cites nothing after D209, describes `rev_gene` only
as a failed ranker (D203), and treats the layered lead score as production. Since D211, production *is* "filter on
target evidence (T >= 50), then rank by gene-level reversal `rev_gene` per timepoint". Consequences:
- Experiment A's question ("does perturbation add beyond target evidence?") is now the premise of the production
  ranking's second stage. The spec should say so, and its decision tree ("do not add reversal to the lead score")
  needs reframing: reversal is already in production by decision.
- D211 already measured one of A's endpoints. On the D127 biomarker pairs, `rev_gene` within T-eligible pairs orders
  2 of 8 target units correctly early and 1 of 8 late (`results/lincs/d127/d211_gene_summary.tsv`). The D127 pairs
  are therefore *seen* for exactly this feature-within-filter structure, and the D127-style endpoint in A can only be
  secondary or descriptive, never confirmatory.
- Section 2.3 says subgroups "have historically lacked window data". Stale: D211 built subgroup timepoint runs;
  8 of 14 subgroups have data (2,402 rows, `window_reversal_strata.tsv.gz`).
- D210 (the pipeline reproduces itself from `scripts/run_pipeline.sh`) is the provenance mechanism the spec's
  section 12 asks for; reference it rather than build a parallel one.

**2. Experiment A's outcome is not blind.** The spec requires committing the protocol "before looking at any new
outcome results". PRISM/GDSC selective sensitivity has already been tested against the P1 feature family: C16
(2026-10-12), C28 (2026-11-03) and `results/lincs/sensitivity_validation_v2/` (2026-11-08) report median
within-cohort Spearman of `rev`, `rev_mech` and related scores with selective killing of -0.06 to 0 (PRISM
secondary, 36 groups; `t1_summary.tsv`), and D206 looked again. The amendment must declare the outcome as
previously examined (as D127 did with its `seen` pairs) and confine confirmatory status to what is genuinely new:
the held-out *incremental* comparison B1 vs B1+P on matched rows. A marginal P1-outcome association is not new.

**3. The outcome producer already exists; do not rebuild it from scratch.** Section 2.3 and section 14 ask for the
D206 table to be rebuilt "from raw matrices with a committed producer". That producer is committed:
`scripts/79k_sensitivity_validation.py` reads PRISM 20Q2 (`data/raw/prism_20q2`) and GDSC 8.5 (`data/raw/gdsc_8.5`),
weights each unit's in-key lines by Celligner, and writes `obs_raw` (kill = 1 - AUC or -log fold change; sign
documented in the header) and `obs_sel` (the drug's kill z-scored across all lines, then weighted: selective
killing). This is the spec's `Y_cx` in A2. Its last build is 2026-11-08, before D208, so its *score* columns come from
the deleted consensus arm (the observed-sensitivity columns do not depend on the regulon arm). Required action:
rerun it on `_v2_tf` through `scripts/runlog.sh` (it was repointed in D208 but is not in the runner), and trace
whether D206's orphan table was derived from it.

**4. The primary outcome is not fixed.** Four datasets exist (PRISM primary, PRISM secondary, GDSC1, GDSC2) with
overlapping drugs and lines, and two measures (`obs_raw`, `obs_sel`). The spec names none as primary and does not
forbid pooling overlapping datasets. Fix one primary (dataset x measure), declare the rest secondary, never pool
them. `obs_sel` matches A2's within-drug standardization and removes general potency. Measured coverage of the real
A universe (target evidence present x observed window x outcome, counted in FDA ingredients):

| window | dataset | cohorts | cohorts with >= 20 ingredients | median ingredients |
|---|---|---|---|---|
| early | PRISM primary | 16 | 11 | 51.5 |
| early | PRISM secondary | 16 | 10 | 39.5 |
| early | GDSC1 / GDSC2 | 18 / 18 | 8 / 10 | 18 / 21 |
| late | PRISM primary | 24 | 12 | 20.5 |
| late | PRISM secondary | 24 | 11 | 11.5 |
| late | GDSC1 / GDSC2 | 21 / 26 | 10 / 12 | 16 / 17.5 |

The section A6 coverage gate (>= 5 cohorts with >= 20 compounds) passes for every dataset and window.

**5. Experiment B's sign rule cannot run as written.** It calls a response positive/negative only when "its two-sided
95% uncertainty interval excludes zero", from independent conditions. Measured: median `n_sig` per condition is 1
(10.1% of 26,130 conditions have >= 2 signatures), and only **400** drug x line pairs have >= 2 conditions in both
windows. The matched-dose arm is large (2,639 drug x line x dose combinations in both windows; 860 drugs, 27 lines,
in-key for 18 cohorts, 13 with >= 20 drugs) but has no per-condition uncertainty. The amendment must choose before
any outcome is seen: (a) confirmatory B restricted to the 400 pairs, stated as such; or (b) an empirical noise model
estimated from split halves or replicate signatures and applied to single-condition rows, with its false-flip rate
from the pseudo-replicate null fixed in advance.

## High priority

**6. Power for A's primary metric.** Within-cohort Spearman over 20-50 drugs has a standard error of roughly
0.15-0.23 per cohort, and the gate admits 10-12 cohorts. A median cohort delta of +0.05 is very likely below the
minimum detectable increment the spec itself requires (section 11.1), so the registered outcome would likely read
"underpowered". Decide before freezing: either make the pooled, cohort-stratified out-of-fold metric primary, or
accept "underpowered" as the expected result and state it.

**7. P feature families omit the production score.** P1 lists `rev`, `rev_mech` and magnitude but not `rev_gene`,
which is production since D211 and is available per regulator as `rev_gene_contrib` in
`regulator_contributions_by_window_all.tsv.gz`. Add it explicitly to P1 (or as P1b), early and late separately.

**8. The universe is FDA-only by construction.** The B1 baseline uses `layered_S` / `T`, which exist only for FDA,
gate-passing rows (14,300 rows over 55 units; T populated for 12,155; `layered_score_v2_v2_tf.tsv`). The schema of
`drug_disease.tsv` implies every LINCS compound. State the restriction; non-FDA compounds cannot enter A.

**9. Shared inputs missing from the independence section.** Section 1 and 11.4 warn about D and A. Two more shared
inputs belong in `evidence_dependence_graph.tsv`: (i) overactivity O is the disease VIPER NES of the target
regulator, and reversal multiplies the same disease NES vector, so B1's O and every P feature share disease input;
(ii) Celligner in-key weights aggregate both the perturbation features (window scorer) and the outcome Y (`79k`).
Neither invalidates A, but both are consistency links, not independent evidence.

**10. Process conflict to resolve.** Section 12.1 requires a new feature branch, a PR and no merge. This project's
working convention is commit and push to `main` after each step and merge into `celligner`. The owner should choose
one for this study before work starts.

## Medium

11. **Full edge decomposition is expensive.** `target_contributions_by_window.tsv.gz` holds the top 25 regulators
    per drug only (the spec notes this). Decision-bearing B and C need `79e3` with `TOP_REGS=0` or `QUERY`. Scope it to
    the U_B drugs and estimate size first; the all-regulator table alone is ~18.5M rows (D194).
12. **Parquet is not available.** `pyarrow` is not installed in the `celligner` environment. Use compressed TSV
    (the spec allows it) or add `pyarrow` to `environment.celligner.lock.yml` as a recorded environment change.
13. **Pin C's signaling network by file.** Name `data/raw/omnipath/omnipath_signed_directed.tsv.gz` (16,518 edges)
    and `signor_signed_directed.tsv.gz` (10,346), retrieved 2026-11-10 with checksums in `config/datasets.yaml`;
    forbid a re-pull (the registry calls a re-pull a new decision). These are `datasets=omnipath` signaling edges,
    not CollecTRI transcription regulons, so C3.2 stays within D208; say so explicitly.
14. **Thread budget.** Add `N_CORES` and BLAS thread limits (project rule: under half of the 64 threads in total).
15. **Citation slips.** D156 is the per-timepoint activity floor, not a proliferation-slope decision; D128 is the
    target-mechanism scoring proposal, not gene-dependence work; D158's primary-target enrichment is recorded as
    3.8x (the spec says 3.84x). `window_reversal.tsv.gz` has condition counts but no line counts (`n_lines` is in the
    regulator-contribution table).
16. **Naming collision.** Experiments A-F include "D" and "E", which read as decision numbers (`D_network_pairs.tsv`
    next to D127, D158). Prefix the experiments (for example `X-A` ... `X-F`) or name them.

## Verified correct

- All 24 paths in section 2.1 exist; field names of the window table, the regulator-contribution table and the
  target-contribution table match the spec.
- Section 3's reversal and per-regulator contribution match `docs/FORMULAS.md` and the window scorer.
- D137 (trace negative controls failed, 3 of 5 false positives), D154/D157 (per-window proliferation beta),
  D158 (primary targets enriched, secondary at chance), D183 (closed), D206 (orphan table), D208 and D209
  (16/20 units, p = 0.0059) are described accurately.
- The historical both-window set reproduces from the full run: 4,128 rows, 939 drugs, 14 cohorts.
- No unrecorded time sentinels exist in the condition index (0 negative times), so section 12.3's `-666` test passes
  trivially but should stay.
- Drug VIPER is computed per network over every condition (38 networks), so Experiment D's network swap
  (same conditions through BRCA LumA and Basal networks) is directly feasible.
- D134 outputs and `SPEC_target_recoverability.md` exist for C1.

## Feasibility by stage

| Stage | Feasible now? | Binding constraint |
|---|---|---|
| 0 Audit | Yes | Rerun `79k` on `_v2_tf` first (finding 3) |
| A Incremental value | Yes, with amendments 1-4, 6-9 | Outcome previously seen; power for +0.05 |
| B Temporal | Partly | Uncertainty for the sign rule (finding 5); full-edge cost (11) |
| C Target recovery | Yes | Pinned network files (13); secondary-target coverage threshold of 30 units |
| D Network sensitivity | Yes | None found |
| E Program access | Yes | Repeat conditions are sparse (same limit as B) |
| F Synthesis | After A-E | — |

## Recommended next step

Add a dated "Amendment 1 (2026-11-23)" section to the spec carrying findings 1-10 (and the medium items the owner
accepts), rerun `79k_sensitivity_validation.py` on `_v2_tf` under `runlog.sh`, then commit the amended spec plus
`config/mechanistic_incremental_v1.yaml` as the registration. Nothing in this audit computed any outcome of the
proposed experiments; the only quantities examined are coverage counts and previously published results.
