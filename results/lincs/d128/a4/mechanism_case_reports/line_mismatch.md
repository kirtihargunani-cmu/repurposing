# line mismatch

**What this case is.** the drug has no profiled line inside the disease's tissue key, so its response is borrowed from other tissues (or the best matched line is a poor match).

**Why it matters.** the mechanism may hold in the profiled line and not in the disease's tissue; `source` and `line_match` travel with every row so this is never invisible.

**Rows:** 4278 of 10131 drug-level rows (42.2%)

**Largest cases** (by M where defined):

| cohort | pert_id | cmap_name | source | line_match | transfer_conf | M | M_pct | layered_pct |
|---|---|---|---|---|---|---|---|---|
| TCGA-LUAD | BRD-K57886322 | fluocinonide | direct | 0.0503 | 1 | 0.1581 | 100 | 98.34 |
| TCGA-LUSC | BRD-K57886322 | fluocinonide | direct | 0.0087 | 1 | 0.1581 | 100 | 100 |
| TCGA-LUAD | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0503 | 1 | 0.1516 | 99.09 | 98.34 |
| TCGA-LUSC | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0087 | 1 | 0.1516 | 99.09 | 100 |
| TCGA-LUSC | BRD-A27887842 | prednisolone | direct | 0.0157 | 1 | 0.148 | 98.18 | 100 |
| CESC-SCC | BRD-K62310379 | fluticasone-propionate | direct | 0.0494 | 1 | 0.1461 | 100 | 64.93 |
| CESC-SCC | BRD-A88138582 | halcinonide | direct | 0.0494 | 1 | 0.1448 | 98.98 | 64.93 |
| BRCA-Basal | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0053 | 1 | 0.144 | 100 | 75.33 |
| BRCA-LumB | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0151 | 1 | 0.144 | 100 | 100 |
| BRCA-LumA | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0169 | 1 | 0.144 | 100 | 100 |
| BRCA-Her2 | BRD-A15297126 | fluocinolone-acetonide | direct | 0.0233 | 1 | 0.144 | 100 | 100 |
| TCGA-LUAD | BRD-K56515112 | medrysone | direct | 0.0503 | 1 | 0.1391 | 97.27 | 98.34 |
