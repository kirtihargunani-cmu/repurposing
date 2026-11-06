# D108 reading prompt (completed clinical trials)

You are extracting the outcome of completed clinical trials from the abstracts of their linked papers, for a
drug-repurposing pipeline. Work only from the text in the input file; do not use outside knowledge.

Input: one JSON object per line with `nct`, `title`, `drug` (the drug we care about), `condition`, `phase`,
`interventions`, and `papers` (each with `pmid`, `title`, `abstract`, `journal`, `year`).

For every input line write exactly one JSON line to the output file with these keys:

- `nct`: copied from the input.
- `outcome`, one of:
  - `positive`: the trial met its primary endpoint or the authors conclude the drug (or the regimen containing it)
    was effective / worth pursuing.
  - `negative`: it did not meet its primary endpoint, or the authors conclude no benefit, insufficient activity,
    unacceptable toxicity, or that it should not be pursued.
  - `mixed`: primary endpoint not met but a subgroup or secondary benefit is reported, or results are conflicting.
  - `not_reported`: the paper is about this trial but gives no efficacy conclusion (e.g. design paper, PK-only, safety-only with no verdict, biomarker substudy).
  - `not_this_trial`: none of the papers reports this trial's results (e.g. background review, a different study).
- `reason_class` (for `negative` and `mixed`; otherwise `not_applicable`): `lack_of_efficacy`, `toxicity`, `both`, or `other`.
- `reason`: one plain sentence saying why it failed or what it showed (empty for not_reported / not_this_trial).
- `quote`: a verbatim, contiguous excerpt (10-60 words) copied exactly from one abstract or paper title that supports
  the outcome. Copy characters exactly, no paraphrase, no ellipses. Empty for not_reported / not_this_trial.
- `pmid`: the PMID the quote comes from (empty if no quote).
- `drug_role`: `single_agent` if the drug is tested alone or vs placebo/control, `combination` if it is one part of a
  regimen and the result cannot be attributed to it alone, `background` if the drug is part of the standard backbone in
  all arms (so the trial does not test it).

Judge the outcome for the trial as a whole, not for the drug alone; `drug_role` records whether it is about the drug.
Phase 1 trials that only establish a dose with no efficacy verdict are `not_reported`.
Write valid JSON (one object per line, double quotes, escape internal quotes). Process every line; do not skip any.

## Amendment (2026-11-06, role pass): `drug_role` definitions, applied to every trial

- `single_agent`: the drug is the experimental intervention and the only difference between arms (alone, vs placebo /
  control, or added to a backbone and compared with the backbone alone); or a single-arm trial of the drug alone.
- `combination`: the drug is part of the experimental regimen together with other new drugs, so the result cannot be
  attributed to it alone (including single-arm trials of a multi-drug regimen).
- `comparator`: the drug is only in the control / standard arm against which another experimental treatment is tested
  (e.g. erlotinib as the control arm in FLAURA). The trial's outcome is about the other treatment.
- `background`: the drug is in the backbone of every arm, so the trial does not test it.
