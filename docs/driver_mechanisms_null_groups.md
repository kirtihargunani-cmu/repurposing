# Driver mechanisms of groups with no CollecTRI regulator at FDR < 0.10

> **Status 2026-10-31: v1 context** (written for the v1 CollecTRI run, 2026-10-03). The literature summary still holds; the 'no regulator at FDR < 0.10' framing refers to v1.

Compiled 2026-10-03 from the literature (background research agent; PMIDs checked against PubMed titles, journals and years). Frequencies are approximate (about ±5 points) and quoted from the TCGA marker papers and Sanchez-Vega et al., Cell 2018 (PMID 29625050). Lines marked **[Inference]** are interpretation, not cited fact. Our own patients' alteration frequencies are in `results/step7/driver_classes.tsv` (`scripts/60_driver_classes.py`).

Driver classes: (a) kinase / RTK / RAS-MAPK / PI3K; (b) transcription factor; (c) tumor suppressor / chromatin; (d) viral. "mRNA-visible" asks whether the driving event shows up in target-gene mRNA, which is what VIPER reads.

| Cohort | Primary driver class | Key events (approx. freq.) | Known subtypes | mRNA-visible driver? |
|---|---|---|---|---|
| LUAD | (a) + (c) | KRAS 33%, EGFR 14%, BRAF 10%, STK11 17%, KEAP1 17%, NKX2-1 amp ~12% (TCGA Nature 2014, PMID 25079552) | TRU / PP / PI | Mostly no; NRF2 targets in ~20% |
| LUSC | (b) + (c) | TP53 90%, CDKN2A 72%, SOX2/TP63 3q amp ~40%, NRF2 pathway 34% (PMID 22960745) | classical / basal / secretory / primitive (Wilkerson 2010, PMID 20643781) | Yes (SOX2, TP63, NRF2), subtype-restricted |
| KIRC | (c) | VHL ~85%, 3p loss > 90%, PBRM1 41%, SETD2 12%, BAP1 10% (PMID 23792563) | ccA / ccB (Brannon 2010, PMID 20871783); TCGA m1-m4 | Yes (HIF targets); HIF mRNA itself unchanged |
| HNSC | (c) + (d) + (a) | HPV ~15%, TP53 84% (HPV-), PIK3CA 34-56%, 3q amp, EGFR amp ~15% (PMID 25631445) | atypical / basal / classical / mesenchymal (Walter 2013, PMID 23451093) | Partial (E2F, TP63/SOX2) |
| COAD / READ | (c) + (a) | APC 80% (WNT > 90%), KRAS 43%, TP53 60%, BRAF 10% (PMID 22810696) | CMS1-4 (Guinney 2015, PMID 26457759) | WNT targets yes, but shared with normal crypts |
| PRAD | (b) | ERG fusion 46%, other ETS 7%, SPOP 11%, PTEN 17% (PMID 26544944) | 7 classes: ERG, ETV1, ETV4, FLI1, SPOP, FOXA1, IDH1 | Yes (ERG), in about half of tumors |
| THCA | (a) | BRAF V600E 60%, RAS 13%, fusions 15%; MAPK driver in ~97% (PMID 25417114) | BRAF-like / RAS-like; thyroid differentiation score | No (protein-level); only dedifferentiation is visible |
| BRCA LumA | (a) + (b) lineage | PIK3CA 45%, GATA3 14%, CDH1, MAP3K1 (PMID 23000897) | single PAM50 class (PMID 19204204) | ER/GATA3 yes, but shared with normal breast |
| SKCM | (a) + (c) | BRAF 52%, NRAS 28%, NF1 14%, CDKN2A 40% (PMID 26091043) | BRAF / RAS / NF1 / triple-WT; immune / keratin / MITF-low | No (MAPK); MITF variable |
| STAD | mixed | CIN 50%, MSI 22%, GS 20%, EBV 9%; RTK amplifications in CIN (PMID 25079317) | CIN / GS / MSI / EBV; ACRG (PMID 25894828) | Partial, subtype-specific |
| LIHC | (c) + (d) | TERT 44%, TP53 31%, CTNNB1 27%, HBV 30% (PMID 28622513) | iCluster 1-3; Hoshida S1-S3 (PMID 19723656) | Partial (WNT targets in ~27%) |
| BLCA | (c) + (b) | TP53 48%, KMT2D 28%, KDM6A 26%, FGFR3 14%, PPARG amp (PMID 28988769) | LumP / LumU / LumNS / stroma-rich / Ba-Sq / NE (Kamoun 2020, PMID 31563503) | Yes, but opposite directions by subtype |
| PAAD | (a) + (c) | KRAS 93%, TP53 72%, CDKN2A ~90% inactivated, SMAD4 32%; low purity (PMID 28810144) | classical / basal-like (Moffitt 2015, PMID 26343385); Bailey 4 classes (PMID 26909576) | No (KRAS); low purity |

## Per cohort notes
- **LUAD**: RTK/RAS/RAF altered in ~75%. TRU (NKX2-1-high), PP (KRAS/STK11) and PI (NF1/TP53) each make up about a third and differ in TF programs. [Inference] Kinase drivers plus three subtypes with opposing TF programs dilute any consensus signal.
- **LUSC**: SOX2/TP63 amplification and NRF2 activation are TF-level and mRNA-visible, concentrated in the classical subtype (~35-40%). [Inference] The most surprising null; normal lung basal cells also express TP63.
- **KIRC**: VHL loss stabilizes HIF protein; HIF targets (VEGFA, CA9, EGLN3) are strongly mRNA-visible. [Inference] A miss here points to the method (small or mixed-sign HIF regulons in CollecTRI, conservative null, GTEx vs TCGA differences). Use it as a sensitivity positive control.
- **HNSC**: HPV E6/E7 inactivate p53 and RB at the protein level; E2F targets are mRNA-visible. [Inference] Viral vs non-viral tumors and different anatomic sites dilute the signal.
- **COAD / READ**: APC loss stabilizes beta-catenin; TCF7L2 targets are visible but normal crypts also have high WNT activity. CMS2 is WNT/MYC-driven; CMS4 is mesenchymal and stromal.
- **PRAD**: TMPRSS2-ERG strongly raises ERG mRNA, but only in about half of tumors. AR, FOXA1 and NKX3-1 lineage activity is shared with normal prostate.
- **THCA**: low mutation burden and kinase-driven; much of normal thyroid identity is retained. BRAF-like vs RAS-like tumors may move TFs in opposite directions.
- **BRCA LumA**: the least transformed, lowest-proliferation breast subtype; ER/FOXA1/GATA3 are also active in normal luminal cells; PI3K is not TF-encoded. LumB and Basal pass on proliferation (E2F/MYB/FOXM1).
- **SKCM**: GTEx skin is keratinocyte-dominated; most TCGA samples are metastases; MITF phenotype switching may cancel.
- **STAD**: four etiologically different subtypes; GATA4/GATA6 amplification in CIN is TF-level.
- **LIHC**: the hepatocyte-like subtype keeps HNF4A/HNF1A/CEBP identity while others lose it, so the directions cancel against normal liver.
- **BLCA**: near-symmetric luminal (PPARG/GATA3/FOXA1) vs basal-squamous (TP63/KRT5/STAT3) split.
- **PAAD**: low purity, desmoplastic stroma, and an acinar-dominated GTEx reference make the contrast mostly tissue composition.

## Pattern [Inference]
1. Kinase-driven tumors whose drivers act on protein, not mRNA: THCA, PAAD, SKCM, LUAD.
2. Lineage TFs shared with the matched normal tissue: LumA, PRAD, COAD, THCA.
3. Subtypes with opposing TF programs that cancel when pooled: BLCA, PAAD, LIHC, HNSC, LUSC, LUAD, colorectal.
Groups that passed are TF-defined or strongly proliferative and differ sharply from their GTEx comparator (seminoma, Wilms, gliomas, LumB, Basal, sarcomas). TF-activity benchmarking across TCGA: Garcia-Alonso, Genome Res 2019 (PMID 31340985).
