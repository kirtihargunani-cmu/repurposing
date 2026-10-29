# LLM coordination layer: requirements collected while building the pipeline

The project will add an LLM layer that coordinates the pipeline and reports its results (not yet designed). This file
collects requirements as they come up, so the design starts from them. Add an entry whenever a pipeline decision says
"the report must say X".

## Reporting requirements

### R1. Say what evidence each disease's drug list rests on (user 2026-10-28; D76-D78)
Every disease (and subgroup) result must state its evidence status, taken from `drug_reversal_summary*.tsv`
(`scripts/79e_drug_reversal.py`, celligner branch):

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

## Related rules the layer should respect
- Thresholds and decision rules are fixed before results are looked at, and recorded in `docs/DECISIONS.md`; the layer
  reports them and does not tune them.
- Caveats per result live in `docs/DECISIONS.md` (disease side) and `docs/notes/drug_side_decisions_and_caveats.md`
  (celligner branch); a report about a disease should carry the caveats that apply to it.
