# ambiguous action

**What this case is.** ChEMBL gives no unambiguous direction for the drug on the target (MODULATOR, mixed), and P11 forbids assuming one for M.

**Why it matters.** a real annotation gap, not a model failure. R is still computable for these rows: it carries no action sign.

**Rows:** 207 of 10131 drug-level rows (2.0%)

**Largest cases** (by M where defined):

| cohort | pert_id | cmap_name | M_abstain_reason | R | layered_pct |
|---|---|---|---|---|---|
| UCEC-SER | BRD-A70858459 | estramustine | ambiguous_action:direction_assumed;no_footprint | -0.0381 | 28.18 |
| UCEC-SER | BRD-K04218075 | clomifene | ambiguous_action:direction_assumed | -0.008118 | 79.56 |
| UCEC-SER | BRD-K29950728 | clomifene | ambiguous_action:direction_assumed | -0.02951 | 79.56 |
| UCEC-SER | BRD-K63828191 | raloxifene | ambiguous_action:direction_assumed | -0.004145 | 56.35 |
| UCEC-SER | BRD-K90195324 | bazedoxifene | ambiguous_action:direction_assumed | 0.01794 | 31.49 |
| UCEC-EEC | BRD-A70858459 | estramustine | ambiguous_action:direction_assumed;no_footprint | -0.04026 | 89.39 |
| UCEC-EEC | BRD-K04218075 | clomifene | ambiguous_action:direction_assumed | -0.01576 | 100 |
| UCEC-EEC | BRD-K29950728 | clomifene | ambiguous_action:direction_assumed | -0.01994 | 100 |
| UCEC-EEC | BRD-K63828191 | raloxifene | ambiguous_action:direction_assumed | 0.003193 | 92.18 |
| UCEC-EEC | BRD-K90195324 | bazedoxifene | ambiguous_action:direction_assumed | 0.01608 | 90.5 |
| THCA-RAS-like-follicular | BRD-A70858459 | estramustine | ambiguous_action:direction_assumed;no_footprint | -0.02303 | 39.89 |
| THCA-RAS-like-follicular | BRD-K04218075 | clomifene | ambiguous_action:direction_assumed | 0.009183 | 84.15 |
