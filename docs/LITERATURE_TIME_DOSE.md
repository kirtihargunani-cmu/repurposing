# Time and dose in drug signatures: literature check

Written 2026-11-03 to check the drug-side time / dose recommendations (D88-D90) against published work. Abstracts and
publisher pages were read (PMC full texts blocked automated access); claims marked "abstract" rest on the abstract only.

## What the literature supports

| Claim used in our scoring | Source | What it says | Checked |
|---|---|---|---|
| Early profiles capture mechanism | Lamb et al. 2006, *Science* 313:1929 (Connectivity Map; doi:10.1126/science.1132939) | The original CMap profiled MCF7 and other lines after about 6 h of treatment against 6 h vehicle controls | records of the paper and of its GEO data (GSE5258) |
| 24 h profiles record phenotype more than mechanism | De Wolf et al. 2016, *Assay Drug Dev Technol* 14:252 (PMID 27187605) | Public LINCS compound data at 6 h and 24 h; "the longer time point is likely to result in phenotypic rather than mechanistic effects being recorded"; two time points add information | abstract |
| Direct and secondary estrogen targets differ; secondary ones are proliferation genes | Bourdeau et al. 2008, *Nucleic Acids Res* 36:76 (PMID 17986456) | MCF7, estradiol with / without cycloheximide: >700 primary and ~1,300 secondary targets; EREs enriched in up-regulated primary targets; secondary targets enriched for E2F sites | abstract |
| "Early = direct" is only an approximation | Tullai et al. 2007, *J Biol Chem* 282:23981 (PMID 17575275) | Of 133 genes induced within 4 h, 37% immediate-early, 44% delayed primary (no protein synthesis needed), 19% secondary | abstract |
| Most estradiol-regulated genes in MCF7 go down | Frasor et al. 2003, *Endocrinology* 144:4562 (PMID 12959972) | >400 robustly regulated genes, ~70% down-regulated. The abstract does **not** say whether the down-regulated genes respond later | abstract |
| Low doses avoid cytotoxic secondary effects when reading master regulators | Mundi et al. 2023, *Cancer Discovery* 13:1386 (OncoTarget / OncoTreat; doi:10.1158/2159-8290.CD-22-1020) | VIPER with context-matched networks; drugs profiled at **24 h** at sublethal doses (EC20 and EC20 / 10); 15 of 18 predicted drugs inverted the MR module in vivo | abstract, figure 2 |
| Time and dose can be modelled jointly per gene | Serra et al. 2020, *GigaScience* 9:giaa055 (TinderMIX; PMID 32449777) | Integrated time x dose models per gene, effect maps and an integrated point of departure (Open TG-GATEs) | abstract |

## What it does not support (corrections to our earlier statements)

- **"The late luminal rise under anti-estrogens is a compensation linked to resistance."** Not supported. In
  endocrine-resistant MCF7 derivatives luminal genes are *lost*, GATA3 is methylated, and restoring GATA3 resensitizes
  the cells (Hu et al., GATA3 / AP1 enhancer reprogramming, PMC7737911). MSigDB's set of genes raised in MCF7 by
  4-hydroxytamoxifen, raloxifene or fulvestrant (Frasor et al. 2004, PMID 14973112) contains none of GATA3, FOXA1, ERBB4,
  ANXA9, HEBP1, TBX3, SCCPDH. The late luminal rise in LINCS is our observation, not an established mechanism; flags
  describe it without interpreting it.
- **"Frasor 2003 shows the late changes are secondary."** The abstract gives no early / late comparison.
- **"Our approach should use early profiles because the field does."** The closest pipeline (OncoTreat) reads drugs
  at 24 h and controls secondary effects with low doses. No published head-to-head of 6 h vs 24 h for mechanism
  calls on L1000 was found.

## Methods for later per-gene / per-regulator modelling

| Method | Use | Fit for LINCS |
|---|---|---|
| BMDExpress 2 (Phillips et al. 2019, *Bioinformatics* 35:1780; doi:10.1093/bioinformatics/bty878) | benchmark dose per gene / pathway | yes, where >= 3-4 doses |
| DRomics (Larras et al. 2018, *Environ Sci Technol*; CRAN) | dose-response model choice, response shape, benchmark dose | yes, where >= 3-4 doses |
| ImpulseDE2 (Fischer, Theis, Yosef 2018, *Nucleic Acids Res* 46:e119) | impulse model of a time course | no: needs several time points; LINCS has 2 |
| DREM (Ernst et al. 2007, *Mol Syst Biol* 3:74) | when TF-driven gene groups split over time | no: needs several time points |
| TimeDelay-ARACNE (Zoppoli et al. 2010, *BMC Bioinformatics* 11:154) | ARACNe with time lags, directed links | not as published (one long time series, small networks); a lagged variant across drugs is explored in `scripts/explore/` |

## Links

- Lamb 2006: https://golublab.broadinstitute.org/publications/connectivity-map-using-gene-expression-signatures-connect-small-molecules ; GSE5258 re-analysis https://link.springer.com/content/pdf/10.1186/1471-2105-9-258.pdf
- Subramanian 2017 (L1000): https://pubmed.ncbi.nlm.nih.gov/29195078/
- De Wolf 2016: https://www.ncbi.nlm.nih.gov/pubmed/27187605
- Bourdeau 2008: https://pmc.ncbi.nlm.nih.gov/articles/PMC2248750
- Tullai 2007: doi:10.1074/jbc.m702044200 (PMID 17575275)
- Frasor 2003: https://pubmed.ncbi.nlm.nih.gov/12959972/
- Frasor 2004 gene set: https://www.gsea-msigdb.org/gsea/msigdb/cards/FRASOR_RESPONSE_TO_SERM_OR_FULVESTRANT_UP
- Hu et al. (GATA3 / AP1, endocrine resistance): https://pmc.ncbi.nlm.nih.gov/articles/PMC7737911
- Mundi 2023 (OncoTreat): https://pmc.ncbi.nlm.nih.gov/articles/PMC10239356 ; figure 2 https://pmc.ncbi.nlm.nih.gov/articles/PMC10239356/figure/F2
- Alvarez 2018 (OncoTreat, *Nature Genetics*): https://www.altmetric.com/details/43842407
- TinderMIX: https://pubmed.ncbi.nlm.nih.gov/32449777/
- BMDExpress: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2198920/ ; v2 https://agris.fao.org/search/en/records/65df51ec0f3e94b9e5d758ac
- DRomics: https://cran.r-project.org/web/packages/DRomics/index.html
- ImpulseDE2: https://pmc.ncbi.nlm.nih.gov/articles/PMC6237758
- TimeDelay-ARACNE: https://pubmed.ncbi.nlm.nih.gov/20338053/
- DREM: https://ernstlab.github.io/DREM
- Sirota 2011 (signature reversal): https://pmc.ncbi.nlm.nih.gov/articles/PMC3502016
