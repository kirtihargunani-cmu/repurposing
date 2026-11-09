# Drug-disease scoring methods: literature check (gene-expression and activity space)

Written 2026-11-03. Abstracts, PMC / publisher pages and reviews were read; no full-text benchmark tables were
re-derived. Our scores for comparison: `rev` (weighted reversal of tumor-vs-normal master regulators by drug-induced
regulator activity), `spec` / `spec_mech` (robust z of rev minus the mean over reference diseases; mech = without
proliferation-linked regulators, activity adjusted for proliferation), tiers, and the D98 selectivity flag.

## Methods in use

| Space | Method | Score | Source |
|---|---|---|---|
| Genes | CMap 1.0 connectivity score (CS) | Kolmogorov-Smirnov enrichment of the disease up / down genes in the drug's ranked profile (GSEA-like), permutation-normalized | Lamb 2006; definitions in Samart 2021 |
| Genes | CMap 2.0 WTCS / NCS / tau | WTCS = (ES_up - ES_down) / 2 when signs differ (weighted KS); NCS = WTCS normalized within cell line and perturbagen type; **tau** = signed percentile of a connection vs the query's connections to all reference perturbagens (touchstone), -100..100 | Subramanian 2017; Samart 2021 |
| Genes | XSum, XCos, Zhang (sscMap) | sum / cosine over the most changed genes; Zhang's rank-weighted score | reviewed in Samart 2021; Frontiers Syst Biol 2022 benchmark favours Zhang for drug-disease |
| Genes | RGES / sRGES | reversal score from the KS enrichments, summarized over cell lines and doses; correlates with IC50 | Chen 2017 |
| Activity | OncoTreat / OncoTarget | VIPER (aREA) NES of the patient's top master regulators (25 up / 25 down) in each drug's protein-activity signature, from drug screens in context-matched lines at sublethal doses, 24 h; unmatched models combined by Fisher's method | Alvarez 2018; Mundi 2023 |
| Activity | DeMAND | mechanism of action: regulon edges whose expression dependence a compound disrupts (KL divergence), cell-context network, >= 6 profiles | Woo 2015 |
| Activity | Footprints (DoRothEA / PROGENy) | TF / pathway activity from target genes, used to compare perturbations and diseases | Holland / Schubert et al.; cross-species study 2019 |

## Large issues and gaps

1. **Reversal mostly reflects generic anti-proliferation.** Koudijs 2023 validated transcriptome signature reversion
   against PRISM: predictive value was much lower than reported and disappeared once drug signatures were adjusted
   for general anti-proliferative (viability) effects. Szalai 2019: viability is a major factor in L1000 signatures.
   Our D85 / D94 dominance results and the C28 / C29 validations (known-use AUC ~0.6, specificity tests null, PRISM
   correlations ~0) are consistent with this. The proliferation split and `spec_mech` address it; whether any
   disease-specific signal remains beyond that is not yet shown by our validations either.
2. **Reproducibility of the drug signatures is limited.** Lim & Pavlidis 2021: CMap 1 signatures recovered the same
   compound in CMap 2 only 17% of the time; reproducibility tracks differential-expression strength, dose and
   cell-line responsiveness. Our TAS >= 0.4 activity filter follows this; weak signatures should not be scored.
3. **The scores disagree and none is established as best.** Samart 2021 unifies CS, WTCS, NCS, tau, XSum, XCos, RGES
   and finds they are built differently enough to disagree; head-to-head benchmarks are few and use different ground
   truths (IC50, approved indications, PRISM).
4. **Activity-space reversal has little public benchmarking.** OncoTreat is the main validated activity-space method
   (prospective disease control, in vivo MR inversion for 15 of 18 drugs), but with purpose-built drug screens in
   matched lines and small cohorts. No large benchmark of VIPER-based reversal on L1000 was found; the closest public
   ground truth is PANACEA (below), which is purpose-built for MR-based inference but is a few dozen drugs, not
   L1000-scale. The identity
   confound we traced (tumor-built regulons reading late identity shifts as regulator activity) is not discussed
   in the papers read.
5. **Normalization against the drug's behaviour elsewhere is standard in CMap.** Tau ranks a connection against the
   same query's connections to every reference perturbagen. Our `spec` normalizes the other way (a drug's reversal
   of this disease vs reference diseases), and the D98 margin (a drug's percentile here vs its own median) is
   closer to tau's idea for the drug axis.
6. **Ground truth is weak everywhere.** Validations use IC50, PRISM viability, or indication lists; viability-based
   truths reward exactly the generic anti-proliferation that issue 1 warns about.

## Checked, not in use (2026-11-09)

Two options noted after a repo review; neither is implemented and neither has a pre-registered rule yet. Citations
below were not fetched in-session and need checking against the papers and the current releases.

**PANACEA** — perturbational profiles of clinical oncology drugs across a panel of cancer cell lines, assembled for
the DREAM drug mechanism-of-action challenge (Douglass et al.). Relevance: it is ground truth for the
activity-space step specifically (does VIPER on a drug profile recover the drug's known mechanism), and it is
independent of GDSC, so it does not inherit the circularity that rules GDSC out as a validation set for us. Scale is
tens of drugs in ~10 lines, so it can check the regulator reading, not the ranking. Would address issue 4 above and
sits alongside the D84 / D87 knockdown and overexpression benchmarks rather than replacing them.

**CARNIVAL** (and PHONEMeS, CausalPath) — takes TF activities plus a signed prior signaling network and infers which
upstream nodes could produce the observed TF state. Relevance: it consumes the disease side's output natively
(regulator NES) and returns candidate upstream druggable proteins, i.e. a route from the regulon arm into target
selection rather than a second score competing with it. Since D122 moved the lead ranking to target engagement
(`layered_pct`) and D121 fitted the reversal term to zero weight (v = 0), this is the way the regulon arm could still
inform the lead score instead of being bypassed by it. Open question before any work: whether the signaling prior is
any less literature-biased than CollecTRI, which is the limit that already bites us.

Also noted on the scoring function itself: our `rev` (normalized signed dot product over every tested regulator, no
MR cutoff, D73) is not the statistic OncoTreat validated (two-tail aREA of a fixed top 25 up / 25 down MR set in the
drug's rank-sorted activity profile, models combined by Fisher's). D73 ruled out an FDR-cut MR list, which is a
different object from a fixed top-k. A fixed-k aREA arm is a small code change and would give a direct comparison
against the one validated precedent; it needs its own rule if it is to be more than reported.

## Links

- Samart et al. 2021, Brief Bioinform 22:bbab161, Reconciling multiple connectivity scores: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8597919/ ; https://arxiv.org/abs/2009.09317
- Lim & Pavlidis 2021, Sci Rep 11:17624: https://pmc.ncbi.nlm.nih.gov/articles/PMC8413422
- Koudijs, Böhringer, Guchelaar 2023, Brief Bioinform 24:bbac490 (PMID 36445193): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9851289/
- Fisher et al. 2024, FEBS Open Bio 14:803 (signature reversion, three signature methods): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11073506/
- Frontiers in Systems Biology 2022 (KS vs Zhang vs XSum): https://www.frontiersin.org/journals/systems-biology/articles/10.3389/fsysb.2022.1050730/full
- Chen et al. 2017, Nat Commun 8:16022 (RGES): https://pmc.ncbi.nlm.nih.gov/articles/PMC5510182
- Subramanian et al. 2017, Cell (L1000, tau): https://pubmed.ncbi.nlm.nih.gov/29195078/
- Alvarez et al. 2018, Nat Genet (OncoTreat, GEP-NET): https://pmc.ncbi.nlm.nih.gov/articles/PMC6421579
- Mundi et al. 2023, Cancer Discov (OncoTarget / OncoTreat): https://pmc.ncbi.nlm.nih.gov/articles/PMC10239356 ; figure 2 https://pmc.ncbi.nlm.nih.gov/articles/PMC10239356/figure/F2
- Woo et al. 2015, Cell 162:441 (DeMAND): https://pmc.ncbi.nlm.nih.gov/articles/PMC4506491
- Szalai et al. 2019, Nucleic Acids Res (viability signature): https://www.biorxiv.org/content/10.1101/454348v2.full-text
- DoRothEA / PROGENy cross-species footprints: https://www.biorxiv.org/content/10.1101/532739.full.pdf
- Douglass et al. (PANACEA / DREAM MoA challenge), drug mechanism inference from perturbational profiles - citation
  and dataset release to verify: https://www.synapse.org/ (challenge portal); search "PANACEA perturbational
  profiles oncology"
- Liu et al. 2019, npj Syst Biol Appl (CARNIVAL, causal network contextualization from TF activity):
  https://saezlab.github.io/CARNIVAL/ ; PHONEMeS: https://saezlab.github.io/PHONEMeS/
