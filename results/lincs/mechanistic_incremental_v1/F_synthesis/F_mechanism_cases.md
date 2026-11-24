# Worked mechanism cases (chosen by the frozen category rules, not by outcome)
One case per category and window: the highest gene-level percentile in the category, ties by ingredient name. Statuses come from the stage outputs; nothing here re-ranks drugs.

## i positive control (D127 primary pair, high cohort)
- **early: LAPATINIB DITOSYLATE** in BRCA-Her2 (target unit ERBB2; T 99.8; gene-level percentile 98; 1 conditions). Evidence: disease supported · perturbation supported · temporal supported · target unknown · network not_tested · program unknown.
- **late: OSIMERTINIB MESYLATE** in TCGA-LUAD (target unit EGFR; T 65.8; gene-level percentile 100; 7 conditions). Evidence: disease supported · perturbation supported · temporal not_tested · target unknown · network not_tested · program unknown.

## ii high target evidence, low gene-level reversal
- **early: NERATINIB MALEATE** in ESCA-SCC (target unit EGFR; T 96.8; gene-level percentile 25; 1 conditions). Evidence: disease supported · perturbation conflicting · temporal not_tested · target unknown · network not_tested · program unknown.
- **late: BOSUTINIB MONOHYDRATE** in CESC-AD (target unit HCK; T 79.8; gene-level percentile 25; 2 conditions). Evidence: disease supported · perturbation conflicting · temporal not_tested · target unknown · network not_tested · program supported.

## iii lower target evidence, high gene-level reversal
- **early: AZACITIDINE** in TCGA-SKCM (target unit DNMT3A; T 50.1; gene-level percentile 100; 1 conditions). Evidence: disease supported · perturbation supported · temporal not_tested · target unknown · network not_tested · program unknown.
- **late: AFATINIB** in TCGA-OV (target unit ERBB4; T 58.5; gene-level percentile 100; 2 conditions). Evidence: disease supported · perturbation supported · temporal not_tested · target unknown · network not_tested · program unknown.

## iv confident time transition
- **early: ALPELISIB** in BRCA-Her2 (target unit PIK3CA; T 75.9; gene-level percentile 100; 1 conditions). Evidence: disease supported · perturbation supported · temporal supported · target unknown · network not_tested · program not_tested.
- **late: LAPATINIB DITOSYLATE** in BRCA-Basal (target unit EGFR; T 55.3; gene-level percentile 100; 5 conditions). Evidence: disease supported · perturbation supported · temporal supported · target unknown · network not_tested · program unknown.

## v secondary-target drug (measured binding)
- **early: CERITINIB** in TCGA-LUAD (target unit ALK; T 60.8; gene-level percentile 100; 2 conditions). Evidence: disease supported · perturbation supported · temporal supported · target unknown · network not_tested · program not_tested.
- **late: AFATINIB** in TCGA-OV (target unit ERBB4; T 58.5; gene-level percentile 100; 2 conditions). Evidence: disease supported · perturbation supported · temporal not_tested · target unknown · network not_tested · program unknown.

## vi filtered out (T < 50)
- **early: ALPELISIB** in BRCA-Basal (target unit PIK3CA; T 45.0; gene-level percentile 100; 1 conditions).
- **late: BELINOSTAT** in TARGET-NBL (target unit HDAC1|HDAC10|HDAC11|HDAC2|HDAC3|HDAC4|HDAC5|HDAC6|HDAC7|HDAC8|HDAC9; T 47.7; gene-level percentile 100; 3 conditions).

## vii chemotherapy table top
- **early: AURANOFIN** in BRCA-Basal (target unit nan; T nan; gene-level percentile 100; 2 conditions).
- **late: ATORVASTATIN CALCIUM** in BRCA-Basal (target unit HMGCR; T 61.6; gene-level percentile 100; 5 conditions).

## viii top novel, >= 2 conditions
- **early: CERITINIB** in TCGA-LUSC (target unit ALK; T 66.7; gene-level percentile 100; 2 conditions). Evidence: disease supported · perturbation supported · temporal supported · target unknown · network not_tested · program not_tested.
- **late: AFATINIB** in TCGA-OV (target unit ERBB4; T 58.5; gene-level percentile 100; 2 conditions). Evidence: disease supported · perturbation supported · temporal not_tested · target unknown · network not_tested · program unknown.
