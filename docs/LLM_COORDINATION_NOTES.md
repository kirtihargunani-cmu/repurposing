# LLM coordination layer: requirements collected while building the pipeline

The project will add an LLM layer that coordinates the pipeline and reports its results (not yet designed). This file
collects requirements as they come up, so the design starts from them. Add an entry whenever a pipeline decision says
"the report must say X".

## Reporting requirements

### R1. Say what evidence each disease's drug list rests on (user 2026-10-28; D76-D78)
Every disease (and subgroup) result must state its evidence status, taken from `drug_reversal_summary*.tsv`
(`scripts/79e_drug_reversal.py`):

| Column | Meaning |
|---|---|
| `evidence_status` | `matched`: at least one tier 1-2 drug was measured in a cell line that resembles the disease; `borrowed only`: tier 1-2 drugs come only from other tissues' lines; `no adequate compound`: neither |
| `best_line_match`, `best_line` | the best Celligner match any LINCS-tested line in the disease's tissue reaches (ratio to the disease's top line; cutoff 0.1) |
| `n_t12_matched`, `n_t12_borrowed` | tier 1-2 drugs by evidence type |
| `evidence_note` | the sentence to use in reports |

- The layer must report the status and `evidence_note` with every drug list, never a ranked list alone.
- When the status is `no adequate compound`, it must say an adequate compound cannot be reported for that disease, not
  fall back to the top of the ranking.
- When it is `borrowed only`, every candidate must be labelled as borrowed: it rests on the drug acting consistently
  across cell lines (transfer confidence >= 0.5, D77), not on data from lines like the disease. Example: BRCA-Basal has
  no basal-like line with LINCS compound data (best BT-20, line match 0.018); anti-estrogens ranked 97-100th only
  because their data come from MCF7 (luminal).
- Per-drug flags behind this: `unmatched_lines` (direct drug, line match < 0.1; D76) and `low_transfer` (borrowed drug,
  transfer confidence < 0.5; D77). Flagged drugs stay in the ranking but never in tiers 1-2; if the layer shows them,
  it must show the flag.

### R2. Lead with disease specificity (user 2026-10-28; D80)
Drug lists are ordered and presented by disease specificity first: `spec_mech` (reversal of this disease relative to
other diseases, proliferation part removed; `79i` `drug_reversal_split*.tsv`), then `spec`. The overall reversal
score `rev` is shown as context, not as the ranking. Reason: by `rev` the same generic drugs top most diseases (32
drugs in the top 10% of >= 20 of 41 groups: statins, palbociclib, MEK inhibitors, anthracyclines); by `spec_mech` none.
A drug that ranks high on `rev` but not on specificity should be described as a broadly acting / growth-inhibiting
candidate, not as specific to the disease.
For subgroups, lead with `spec_sib_mech` (this subgroup vs its sibling, proliferation removed), then `spec_mech`
(`stratum_reversal_split_v2_consensus.tsv`). When siblings barely differ (LGG-noncodel, TGCT-SEM: `spec_sib` <= 0.6),
report one drug list for the parent and say the subgroups do not differ in predicted drug response.

### R3. Say where good-quality data run out (user 2026-11-01; C24)
If a step has no good-quality data for a disease, the report says so and names the step ("no good-quality data at
<step>: <reason>"), instead of showing results built on poor data. The rule applies at every step; where a step has a
data-quality check, its output gives the step and the reason.
- **Drug side (implemented):** `results/lincs/data_gaps_v2.tsv` and `data_gaps_strata_v2.tsv` (`79r_data_gaps.py`). When `good_quality_data` is False, the layer shows `data_gap_note` and no drug list. Steps:
  D1 no DepMap line of the disease's cell type; D2 no such line has good-quality LINCS compound data; D1b the lines with
  data do not resemble the disease (line match < 0.1); D2b the matched lines carry fewer than 10 compounds. Example:
  ESCA-SCC: "No good-quality data at drug-side Step 1, cell-line matching (Celligner resemblance): the squamous cell
  carcinoma lines with LINCS data do not resemble this disease (best HCC-15, line match 0.005 < 0.1)."
- Lineage-matched means the cell type matches (squamous, glandular, ...), whatever the organ
  (`config/group_lineages.tsv`); a line of the same organ but another cell type does not count.
- This sits above R1: with good-quality data, report the evidence status (R1) and the list (R2); without, only the gap.
- **Existing checks the layer should map to this rule** (step names to be given in each output as they are revisited):
  Step 4 signature `status` (not ok = no signature); Step 6/7 `low_mrs` (too few master regulators to rank drugs);
  Step 8 keep rule D71 (subgroup not kept); drug side `unmatched_lines` / `low_transfer` per drug.

### R4. Report performance against the noise ceiling (user 2026-11-09; D125 section A0)
Any reported accuracy of the drug-side scoring (known-use AUC, Spearman against response, recovery@k) must be shown
with the attainable ceiling beside it, not as a bare number: the split-half reliability `r` of the LINCS-derived
scores, the implied `AUC_max = (1 + sqrt(r)) / 2`, and attainment `D_observed / sqrt(r)` where `D = 2*AUC - 1`
(definitions and assumptions in [`SPEC_graph_metric_scoring.md`](SPEC_graph_metric_scoring.md) section A0).
- Raw and ceiling-relative numbers are both shown; attainment is never shown alone, because the ceiling treats the
  indication labels as noiseless and they are not.
- **Drug side (implemented 2026-11-09, A0a only so far):** `results/lincs/d125/a0/ceiling.tsv`, one row per reported
  known-use metric with auc, D_obs, r, r_estimator, r_optimistic, auc_max, attainment and a note; `ceiling_applies` is
  False for scores not built from LINCS (the lead layered_pct), which keep their raw AUC only. Both estimators ran (11-09):
  A0b covers 21 of 55 cohorts, the rest carry A0a's r and the flag `r_optimistic`; A0b exceeds A0a by 0.234, so the
  `r_differ_note` column is set on nearly every row and must be stated with the number.
- Where a comparison between pipeline arms is below the resolution the ceiling implies, the layer says the yardstick
  cannot resolve a change of that size (the existing finding in D90 / C28) instead of reporting the difference as a
  result.

## Related rules the layer should respect
- Thresholds and decision rules are fixed before results are looked at, and recorded in `docs/DECISIONS.md`; the layer
  reports them and does not tune them.
- Caveats per result live in `docs/DECISIONS.md` (disease side) and `docs/notes/drug_side_decisions_and_caveats.md`; a report about a disease should carry the caveats that apply to it.
