# polypharmacology

**What this case is.** the drug's own annotated target units disagree: at least one gives a positive M and another a negative one, and P13 averages them.

**Why it matters.** the averaged M hides a real disagreement about which target the response belongs to; the per-unit values are kept in `units` so the average is never the only number.

**Rows:** 1232 of 10131 drug-level rows (12.2%)

**Largest cases** (by M where defined):

| cohort | pert_id | cmap_name | n_units_M | M | M_unit_min | M_unit_max | units |
|---|---|---|---|---|---|---|---|
| TCGA-PRAD | BRD-K49328571 | dasatinib | 6 | 0.06654 | -0.01234 | 0.1852 | ABL1=0.1852|PDGFRB=0.0115|KIT=-0.0123|EPHA2=0.0340|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=-0.0045|ABL1=0.1852 |
| BRCA-Her2 | BRD-K49328571 | dasatinib | 6 | 0.06581 | -0.01545 | 0.2115 | ABL1=0.2115|PDGFRB=-0.0105|KIT=-0.0071|EPHA2=0.0050|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=-0.0155|ABL1=0.2115 |
| BRCA-Basal | BRD-K49328571 | dasatinib | 6 | 0.06013 | -0.01239 | 0.188 | ABL1=0.1880|PDGFRB=-0.0060|KIT=-0.0068|EPHA2=0.0100|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=-0.0124|ABL1=0.1880 |
| BRCA-Her2 | BRD-K92723993 | imatinib | 4 | 0.04838 | -0.001002 | 0.0934 | ABL1=0.0934|PDGFRB=-0.0010|KIT=0.0077|ABL1=0.0934 |
| BRCA-Basal | BRD-K92723993 | imatinib | 4 | 0.04838 | -0.001002 | 0.0934 | ABL1=0.0934|PDGFRB=-0.0010|KIT=0.0077|ABL1=0.0934 |
| BRCA-LumB | BRD-K92723993 | imatinib | 4 | 0.04838 | -0.001002 | 0.0934 | ABL1=0.0934|PDGFRB=-0.0010|KIT=0.0077|ABL1=0.0934 |
| BRCA-LumA | BRD-K92723993 | imatinib | 4 | 0.04838 | -0.001002 | 0.0934 | ABL1=0.0934|PDGFRB=-0.0010|KIT=0.0077|ABL1=0.0934 |
| TCGA-LUAD | BRD-K49328571 | dasatinib | 6 | 0.04632 | -0.01695 | 0.1305 | ABL1=0.1305|PDGFRB=0.0092|KIT=-0.0170|EPHA2=0.0250|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=-0.0004|ABL1=0.1305 |
| TCGA-LUSC | BRD-K49328571 | dasatinib | 6 | 0.0431 | -0.01678 | 0.12 | ABL1=0.1200|PDGFRB=0.0102|KIT=-0.0168|EPHA2=0.0233|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=0.0019|ABL1=0.1200 |
| BRCA-Her2 | BRD-K81418486 | vorinostat | 4 | 0.04119 | -0.01167 | 0.1527 | HDAC1=0.0032|HDAC3=0.1527|HDAC6=-0.0117|HDAC2=0.0206 |
| TCGA-PAAD | BRD-K49328571 | dasatinib | 6 | 0.03989 | -0.02739 | 0.08625 | ABL1=0.0702|PDGFRB=0.0315|KIT=-0.0274|EPHA2=0.0863|BLK|FGR|FRK|FYN|HCK|LCK|LYN|SRC|SRMS|YES1=0.0085|ABL1=0.0702 |
| BRCA-LumB | BRD-K66175015 | afatinib | 3 | 0.03607 | -0.001927 | 0.06713 | EGFR=0.0430|ERBB2=-0.0019|ERBB4=0.0671 |
