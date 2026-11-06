# D109 Cmax extraction prompt (FDA labels)

You are extracting the peak plasma concentration (Cmax) of a drug from its FDA label text, to judge whether lab
concentrations are reachable in patients. Work only from the label text in the input file; do not use outside knowledge.

Input: one JSON object per line with `drug` and `label_text` (clinical pharmacology / pharmacokinetics / dosage text).

For every input line write exactly one JSON line to the output file with keys:

- `drug`: copied from the input.
- `cmax_value`: a number: the mean (or geometric mean / median if that is what is given) Cmax in plasma or serum of the
  parent drug, in adults, at the highest recommended dose for its approved systemic use, at steady state if reported,
  else after a single dose. If a range is given, use its upper end. Total (not unbound) concentration.
  If the label gives no usable Cmax (only AUC, only a metabolite, only animal data, only local / topical use), null.
- `cmax_unit`: the unit exactly as printed (e.g. `ng/mL`, `mcg/mL`, `µg/mL`, `mg/L`, `µM`, `nM`, `pg/mL`).
- `dose`: the dose and regimen the value refers to, as printed (e.g. `400 mg once daily, steady state`).
- `quote`: a verbatim, contiguous excerpt (10-60 words) copied exactly from `label_text` that contains the value.
  Copy characters exactly; no paraphrase, no ellipses. Empty if `cmax_value` is null.
- `note`: short note if anything was ambiguous (e.g. "pediatric only", "value for active metabolite", "prodrug: value for
  active moiety X"); empty otherwise. For prodrugs, report the active moiety's Cmax if given, and say so here.

Write valid JSON, one object per line. Process every line.
