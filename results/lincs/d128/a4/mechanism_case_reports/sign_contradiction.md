# sign contradiction

**What this case is.** at least half of the genes the target reaches are excluded from the signed numerator because the routes to them disagree.

**Why it matters.** M is computed on the minority of genes whose sign survives; the coverage denominator records it, and a high value here means the prior is ambiguous about this target, not that the drug failed.

**Rows:** 1723

**Largest cases** (by M where defined):

| cohort | pert_id | cmap_name | unit_genes | M | M_n_genes | gene_ambiguity_fraction | n_conflicted |
|---|---|---|---|---|---|---|---|
| CESC-SCC | BRD-K23984367 | sorafenib | RAF1 | 0.1392 | 156 | 0.832 | 10 |
| CESC-AD | BRD-K23984367 | sorafenib | RAF1 | 0.1392 | 156 | 0.832 | 10 |
| BRCA-LumB | BRD-K16730910 | regorafenib | RAF1 | 0.07433 | 156 | 0.832 | 10 |
| BRCA-LumA | BRD-K16730910 | regorafenib | RAF1 | 0.07433 | 156 | 0.832 | 10 |
| BRCA-Basal | BRD-K16730910 | regorafenib | RAF1 | 0.07433 | 156 | 0.832 | 10 |
| BRCA-Her2 | BRD-K16730910 | regorafenib | RAF1 | 0.07433 | 156 | 0.832 | 10 |
| TCGA-LUSC | BRD-K77625799 | vandetanib | EGFR|ERBB2|ERBB3|ERBB4 | 0.06337 | 1660 | 0.549 | 775 |
| TCGA-LUAD | BRD-K77625799 | vandetanib | EGFR|ERBB2|ERBB3|ERBB4 | 0.06337 | 1660 | 0.549 | 775 |
| TCGA-LUSC | BRD-K19796430 | sonidegib | SMO | 0.04617 | 180 | 0.909 | 3 |
| TCGA-LUAD | BRD-K19796430 | sonidegib | SMO | 0.04617 | 180 | 0.909 | 3 |
| BRCA-LumB | BRD-K77625799 | vandetanib | EGFR|ERBB2|ERBB3|ERBB4 | 0.04561 | 1660 | 0.549 | 775 |
| BRCA-LumA | BRD-K77625799 | vandetanib | EGFR|ERBB2|ERBB3|ERBB4 | 0.0451 | 1660 | 0.549 | 775 |
