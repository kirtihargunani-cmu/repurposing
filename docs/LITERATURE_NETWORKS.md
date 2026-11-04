# Network building and drug-signature processing: literature check

Written 2026-11-03 to check the network configuration question (D92-D94: pooled cohort networks vs cell-type /
subtype networks) and related processing against the most-cited work. Methods text was read for Paull 2021 and
Alvarez 2016 (PDF text); other entries rest on abstracts and publisher / Bioconductor pages.

| Question | Source | What it says | Fit with our pipeline |
|---|---|---|---|
| One network per what? | Alvarez et al. 2016, *Nat Genet* 48:838 (VIPER; doi:10.1038/ng.3593) | Networks reverse engineered by ARACNe (100 bootstraps) from 17 TCGA tumor data sets plus B-cell and glioma sets: one per tissue / tumor type | Matches our pooled per-cohort networks |
| | Paull et al. 2021, *Cell* 184:334 (MOMA; PMID 33434495) | 9,738 tumors from 20 TCGA cohorts "with n >= 100 samples"; the minimum "reflected the need to generate accurate regulatory network models" with ARACNe; VIPER with tissue-matched regulons; the 112 subtypes are found afterwards by clustering MR activity, not by building per-subtype networks | Matches CUR; subtype networks below 100 tumors depart from it |
| | `aracne.networks` (Bioconductor; Giorgi / Lachmann 2016) | 25 ARACNe-AP networks, one per TCGA tumor type | Same |
| Does tissue matching matter? | Alvarez 2016, Fig. 2a | Accuracy drops with a non-tissue-matched interactome (B-cell or glioma network on other tissues) | Supports context-specific networks in principle |
| How much network noise does VIPER tolerate? | Alvarez 2016, Fig. 2b | Performance degrades noticeably only when > 60% of regulon interactions are randomized; accurate as long as > ~28% of a regulon's targets are real | Smaller, noisier subtype networks may still give usable activity; matching the cell context matters more than edge precision |
| Sample size for ARACNe | Margolin et al. 2006 (ARACNE; arXiv q-bio/0410037); Lachmann 2016 (ARACNe-AP; PMC4937200) | > 100 samples "generally sufficient" for pairwise MI; with fewer samples true edges decay gracefully while false positives stay low; ARACNe-AP's MI threshold depends on sample number | Matches what D94 measured (smaller networks less reproducible) |
| Networks per cell type when cells differ fundamentally | PISCES (Obradovic, Vlahos et al.; Califano lab, bioRxiv 2021; github califano-lab/PISCES) | Single-cell pipeline: cluster first, then build lineage / cluster-specific ARACNe networks | Supports the user's "split fundamentally different cell types" idea, for compartments that are different cell types |
| Mixed or uncertain context | Ding et al. 2018, *Nat Commun* 9:1471 (metaVIPER; PMC5902599) | Score with several interactomes and keep, per regulator, the regulon that best recapitulates its targets in the signature | A middle path: pooled + cell-type networks together instead of choosing one |
| Drug side: viability / proliferation confound | Szalai et al. 2019, *Nucleic Acids Res* (bioRxiv 454348) | Cell viability is a major factor in L1000 perturbation signatures and makes toxic compounds look alike; removing viability-correlated genes reduces the confound | Matches our proliferation split (C3 / C4, 79i) |
| Reversal scoring and line choice | Chen et al. 2017, *Nat Commun* 8:16022 (RGES) | Reversal score correlates with IC50 on LINCS; lineage-matched cell lines with the most compounds chosen per cancer | Matches lineage matching (Celligner weights, D76) |
| Time / dose | Mundi et al. 2023, *Cancer Discov* (OncoTreat); see LITERATURE_TIME_DOSE.md | 24 h at sublethal doses in context-matched lines | D88-D90 |

## Bottom line for D94 / D95

- The best-cited practice is our current production: one tissue-matched network per cohort, >= 100 tumors, subtypes
  defined downstream from master-regulator activity.
- Cell-type networks have precedent (PISCES) when compartments are different cell types; the published minimum for
  bulk TCGA networks is 100 samples (Paull 2021). Breast LumA 560, LumB 209, Basal 190 meet it; Her2 (81), and most
  of the D93 subtypes (ESCA 89 / 92, TGCT 62 / 70, SARC-DDLPS 57 / UPS 76, UCEC-SER 60, LGG-IDHwt 94) do not.
- metaVIPER is the published answer to "which network fits this sample": use both and let each regulator pick.

## Links

- Alvarez 2016: https://cris.unibo.it/bitstream/11585/657580.5/2/ng.3593.pdf ; https://scholars.mssm.edu/en/publications/functional-characterization-of-somatic-mutations-in-cancer-using--2/
- Paull 2021: https://pubmed.ncbi.nlm.nih.gov/33434495/ ; PDF https://cris.unibo.it/retrieve/1e7ce56b-a234-4d1e-b84b-5c3ff1dcb0d1/A%20modular%20master%20regulator%20landscape%20controls%20cancer%20transcriptional%20identity_compressed.pdf
- aracne.networks: https://bioconductor.org/packages/release/data/experiment/html/aracne.networks.html
- ARACNe-AP: https://pmc.ncbi.nlm.nih.gov/articles/PMC4937200
- ARACNE (original): https://arxiv.org/pdf/q-bio/0410037
- PISCES: https://github.com/Aleksobrad/PISCES ; https://rdrr.io/github/califano-lab/PISCES/f/README.md
- metaVIPER: https://pmc.ncbi.nlm.nih.gov/articles/PMC5902599
- Szalai 2019: https://www.biorxiv.org/content/10.1101/454348v2.full-text
- Chen 2017 (RGES): https://pmc.ncbi.nlm.nih.gov/articles/PMC5510182
- Mundi 2023: https://pmc.ncbi.nlm.nih.gov/articles/PMC10239356
