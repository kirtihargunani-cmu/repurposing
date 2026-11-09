# D123 reading prompt: free fraction in plasma (and Cmax where missing)

You are reading FDA label / review passages to record (1) how much of a drug is free (not bound to plasma proteins)
in human plasma and, when asked, (2) its peak plasma concentration (Cmax). Work only from the passages in the input
line; do not use outside knowledge, even if you know the value.

Input: one JSON object per line with `drug`, `binding_passages` (list of text), `need_cmax` (true / false) and
`cmax_passages` (list of text).

For every input line write exactly one JSON line to the output file with keys:

- `drug`: copied from the input.
- `fu`: fraction unbound in human plasma (0-1) of the parent drug, a number, or null.
  - "X% bound" -> fu = 1 - X/100 (e.g. "95% bound" -> 0.05).
  - A range of % bound ("95 to 99%") -> use the midpoint of the bound range (97% -> fu 0.03).
  - "> X%" or "more than X%" bound -> fu = 1 - X/100 and put `bound_limit` in `fu_flags`.
  - Binding that depends on concentration -> the value at therapeutic concentrations if the text says which;
    otherwise the midpoint, and put `concentration_dependent` in `fu_flags`.
  - Prodrugs: use the active moiety if the text gives it, and say so in `note`.
  - Only animal data, only a metabolite, or no statement about plasma / serum protein binding -> null.
  - "low" / "negligible" / "not significantly bound" without a number -> null, and say so in `note`.
- `fu_text`: the binding statement as printed (e.g. "approximately 95% bound to plasma proteins").
- `fu_quote`: a verbatim, contiguous excerpt (10-60 words) copied exactly from one of `binding_passages` that contains
  the value. Copy characters exactly; no paraphrase, no ellipses. Empty if `fu` is null.
- `fu_flags`: `bound_limit`, `concentration_dependent`, both separated by `|`, or empty.
- If `need_cmax` is true, also apply the D109 Cmax rules below to `cmax_passages`; otherwise leave these empty / null:
  - `cmax_value`: mean (or geometric mean / median if that is what is given) Cmax in plasma or serum of the parent
    drug, in adults, at the highest recommended dose for its approved systemic use, at steady state if reported, else
    after a single dose; a range -> its upper end; total (not unbound) concentration; null if no usable Cmax (only
    AUC, only a metabolite, only animal data, only local / topical use).
  - `cmax_unit`: the unit exactly as printed (e.g. `ng/mL`, `mcg/mL`, `µg/mL`, `mg/L`, `µM`, `nM`).
  - `cmax_dose`: the dose and regimen the value refers to, as printed.
  - `cmax_quote`: verbatim excerpt (10-60 words) from one of `cmax_passages` containing the value; empty if null.
- `note`: short note if anything was ambiguous; empty otherwise.

Write valid JSON, one object per line. Process every line of the input.
