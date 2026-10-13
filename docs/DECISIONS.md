# Disease-side decisions and assumptions (running list)

What we decided, why, and what we are assuming. One line per item, in plain language. Add new items at the end of each
section; never delete. If a decision is replaced, strike it through (`~~...~~`) and add the new one with a pointer.
Details and numbers: [`PROGRESS.md`](PROGRESS.md). Drug-side decisions are kept in the `celligner` worktree.

"PI" = decided by the PI / user. "Claude" = chosen by Claude after the user delegated the choice.

## 1. Decisions

### Scope and disease groups (Step 1)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D1 | 10-01 | Use the Toil TCGA + TARGET + GTEx data (all on GENCODE v23) as the main dataset | Tumors and healthy tissue processed by one pipeline | PI |
| D2 | 10-01 | Split mixed cohorts into subtypes (BRCA by PAM50, LGG by IDH/codel, ESCA, CESC, UCEC, PCPG, TGCT, SARC; ALL by lineage) | Different subtypes are different diseases with different drugs | PI |
| D3 | 10-01 | A group needs >= 30 tumor patients; 49 groups included | Too few patients give unstable signatures | PI |
| D4 | 10-01 | GTEx tissue choices reviewed per cancer (e.g. HNSC -> esophagus mucosa; SARC -> muscle + fat) | Best available healthy match for each cancer | PI |
| D5 | 10-01 | Steps 5-7 for every group with a signature; Steps 8-10 only for breast and colorectal | Focus deep work where data are richest | PI |
| D6 | 10-01 | Blood cancers (AML, ALL, DLBC) de-prioritized; their normal-reference search parked | Data more limited and variable; no good healthy reference | PI |

### Covariates and signatures (Steps 3-4)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D7 | 10-01 | One sample per patient (diagnosis sample first) and one per GTEx donor | Repeat samples would over-count patients | PI |
| D8 | 10-01 | Tumor vs healthy signature with limma-voom: ~ tumor/normal + age + sex + ESTIMATE purity + sample-type indicators; no plate term | Standard method; plate cannot be separated from tumor vs GTEx | PI |
| D9 | 10-01 | Protein-coding genes for drug matching; extended gene set only for mechanism work | Drug signatures measure protein-coding genes | PI |
| D10 | 10-01 | Keep every gene's score (signed z); no significance cut-off | Drug matching uses the whole ranked list | PI |
| D11 | 10-01 | CESC-SCC normals = cervix + vagina; CESC-AD = cervix only (minimum lowered to 5 normals) | Same tissue type; too few cervix donors otherwise | PI |
| D12 | 10-01 | Breast is female-only throughout | Male breast cancer is biologically distinct | PI |
| D13 | 10-01 | BLCA compared with adjacent normal instead of GTEx | GTEx has only 9 bladder samples | PI |

### Regulator networks and activity (Steps 5-6)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D14 | 10-02 | Build our own ARACNe networks per cancer, with the full regulator set (TFs, kinases, signaling; ~6,200), no target cap | Drivers may be kinases or signaling proteins; published networks failed our checks | PI |
| D15 | 10-02 | Significance from a "normal vs normal" null (1,000 random splits of the healthy samples); the faster analytic null is never reported | The analytic null assumes genes are independent and overstates significance | PI |
| D16 | 10-02 | Regulator kept only if expressed in the tumors (median logCPM >= 1) and >= 10 targets / >= 50% of targets measured | Avoid scoring regulators that are absent or barely measured | PI |
| D17 | 10-08 | ARACNe is the primary network; CollecTRI and DoRothEA are supporting evidence | ARACNe passes known-biology checks; curated sets do not fit tissue context | PI |
| D18 | 10-08 | Confirmed regulator = ARACNe FDR < 0.05, or FDR < 0.10 with a same-direction FDR < 0.10 call in CollecTRI or DoRothEA; FDR < 0.10 for exploratory selection | Balance discovery and confidence | PI |
| D19 | 10-08 | Failing known-biology checks are reported as limitations; no tuning to make them pass | Avoid fitting the method to expectations | PI |
| D20 | 10-08 | Small/weak groups (CESC-AD, ESCA-SCC, UCS) kept for drug matching | See what they match; flagged | PI |
| D21 | 10-08 | ARACNe consolidation threshold kept as built | Networks pass confounding checks | PI |

### Patient-level corrections (Step 6, input to clustering)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D22 | 10-07 | Remove tumor purity from patient scores (ESTIMATE; ABSOLUTE for LGG from 10-08) | Purity was the main pattern between patients, not mechanism | PI |
| D23 | 10-08 | Gene-level correction for purity + RNA quality (long-minus-short gene index) before patient VIPER | Plate pattern traced to RNA quality differences between batches | PI |
| D24 | 10-09 | CQN normalization tested and not adopted | Did not fix the recurring module; mixed side effects | PI |
| D25 | 10-09 | Mappability filter adopted as standard (Saha & Battle 2018): genes with mappability < 0.8 and cross-mappable edges removed | Removes multi-mapping artifacts at no cost to biology | PI |
| D26 | 10-09 | Recurring PILRB/CCNL2 module kept in clustering | PI reads it as growth / mitochondrial / immune-rescue biology (our data show no proliferation or immune link and a negative OXPHOS link) | PI |
| D27 | 10-08 | TGCT-SEM immune question tabled | Seminomas are lymphocyte-rich by definition | PI |

### Clustering and replication (Steps 8-9)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D28 | 10-09 | Cluster each of BRCA-Basal/Her2/LumA/LumB, COAD, READ separately | Subtypes and colon/rectum differ too much to pool | PI |
| D29 | 10-09 | Features = one score per regulator module (>= 10 co-varying regulators); all modules as a sensitivity check | Thousands of regulators collapse into a few hundred groups | PI |
| D30 | 10-09 | M3C consensus clustering; rule fixed before running (best K significant on both p-values); "no clusters" is an accepted outcome | Avoid finding clusters in continuous variation | PI |
| D31 | 10-08 | Two replication cohorts per cancer (breast: METABRIC + CPTAC-BRCA; colorectal: GSE39582 + CPTAC-COAD) | Independent platforms and populations | PI |
| D32 | 10-09 | In replication cohorts each patient is compared with the cohort's average tumor | No healthy samples in these cohorts; clustering compares patients with each other | PI |
| D33 | 10-12 | Colorectal replication by three methods (transfer TCGA clusters, cluster independently, compare module differences) + CMS comparison | Each method answers a different question | PI |

### Strata and Step 10 signatures
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D34 | 10-12 | Colon split named COAD-CMS1/3-like (MSI-like, 92) and COAD-CMS2/4-like (195); used as the colon strata | Replicated in both cohorts; matches CMS 1+3 vs 2+4 | PI |
| D35 | 10-12 | The three READ clusters kept as separate subgroups, each matched to drugs; flagged tentative | No rectal cohort to replicate; plate-associated | PI |
| D36 | 10-12 | Breast: no clusters, so continuous axes = principal components of the module scores within each subtype | Few independent axes, no hand-picking, transferable to other cohorts | Claude (delegated) |
| D37 | 10-12 | Axis kept if above a shuffled-data null, stable under resampling, and not technical (plate above chance < 0.2; RNA quality and purity rho < 0.5); 13 axes kept | Rule fixed before use; plate measured above chance because small groups have high chance overlap | Claude (delegated) |
| D38 | 10-12 | Every stratum signature = stratum tumors vs healthy GTEx tissue, same model as Step 4 | All signatures mean the same thing ("disease state to reverse") | PI |
| D39 | 10-12 | Breast axis signatures: top third ("high") and bottom third ("low") each vs GTEx (26 signatures), plus a companion list of genes that follow each axis (explanation only) | Consistent with D38; companion shows what drives the axis | PI |
| D40 | 10-12 | Breast axes checked in METABRIC and CPTAC-BRCA by projecting TCGA axis weights (rule fixed before running) | Axes should be used as fixed TCGA directions | PI / Claude |
| D41 | 10-12 | Her2 PC4 kept, flagged weak | Did not replicate in METABRIC | PI |
| D42 | 10-12 | Stratum master regulators built exactly as Step 6, with the regulator-expressed filter judged on the stratum's own tumors | Drug side scores reversal of master regulators | PI |
| D43 | 10-12 | Drugs for strata scored on reversal + specificity (specificity vs sibling strata proposed) | Stratum signatures are 0.97-0.99 similar to the parent; specificity is what separates strata | PI |
| D44 | 10-12 | For strata, cell-line matching uses the average of the stratum's own patients; target gate judged on the stratum's own patients | Strata may match different lines and express different targets | PI (drug side implements) |

## 2. Assumptions

What the results rely on. "Impact" = how much results could change if the assumption is wrong. "Check" = what we
did or could do about it.

| # | Assumption | Impact | Check |
|---|---|---|---|
| A1 | Purity (stroma/immune content) is a nuisance to remove, but it may carry real mechanism (immune evasion, fibrosis) | High | Sensitivity run without ESTIMATE adjustment: proposed, not run |
| A2 | limma treats tumors and normals as equally variable; tumors are more variable | High | voomWithQualityWeights / group variances: proposed, not run |
| A3 | Gene-gene correlation in healthy tissue represents "chance" for tumors (normal-split null); the null is unadjusted for covariates while the signature is adjusted | High | Conservative in practice (GTEx donor structure widens it) |
| A4 | GTEx is a fair healthy reference despite post-mortem collection and source being collinear with tumor vs normal; adjacent normals carry field effects | High | GTEx vs adjacent agreement median ~0.6 (18 groups) |
| A5 | Partial GTEx tissue matches add cell-type differences to signatures (pancreas acinar vs ductal, uterus myometrium, adult adrenal for NBL) | High | Flagged per group |
| A6 | Regulator activity can be read from its targets' mRNA; misses non-transcriptional activation (kinases) | Medium | PROGENy / CPTAC phospho proposed |
| A7 | The expressed-regulator filter drops active but low-mRNA regulators | Medium | None yet |
| A8 | CollecTRI is context-free; unsigned edges (~1.5k of 64.5k) are treated as activation | Medium | CollecTRI is supporting only |
| A9 | ARACNe: co-expression = regulation; one network per cancer, shared by its subtypes and strata | Medium | Edge confounding check: < 3% of edges lose > 50% strength after purity/quality/plate |
| A10 | Subtype labels (PanCanAtlas PAM50, IDH/codel, marker-based ALL lineage; METABRIC CLAUDIN_SUBTYPE; CPTAC PAM50) are correct | Medium | Taken as given |
| A11 | RSEM expected counts behave like counts for voom | Low | Standard practice |
| A12 | Patients with missing covariates (< 2%) are missing at random | Low | Dropped |
| A13 | GTEx age = middle of the 10-year bin | Low | |
| A14 | The diagnosis sample represents the patient; metastasis/recurrence handled by an indicator | Low | |
| A15 | BH FDR is valid under positive dependence between regulators | Low | |
| A16 | Filter thresholds (>= 10 targets, >= 50% present, logCPM >= 1, >= 30 patients, >= 10 normals) are judgment calls | Low | Sensitivity analyses advisable |
| A17 | The long-minus-short gene index measures RNA quality | Medium | Agrees only partly with GTEx RIN / ischemic time; treated as a plate-linked technical axis |
| A18 | The recurring module is biology (PI reading), not a technical axis | Medium | Plate-linked above chance in 39 groups; RNA quality, GC bias and multi-mapping ruled out as causes; cause unknown |
| A19 | Regulator modules (correlation >= 0.7) are meaningful units, and module scores capture patient differences | Medium | Mean-member and PC1 scores agree (r 0.998) |
| A20 | Comparing replication patients with their cohort's average tumor is equivalent to the TCGA comparison with healthy tissue once scores are standardized | Medium | Colon: TCGA modules replicate (rho 0.89-0.94); breast axes replicate in gene content |
| A21 | Array (GSE39582, METABRIC) and RNA-seq scores are comparable after standardizing within each cohort | Medium | Same patterns found on both platforms |
| A22 | CMScaller calls on TCGA/replication data are correct, and the CMS1/3 vs 2/4 naming describes the colon split | Low | Consistent across three cohorts |
| A23 | READ clusters reflect biology despite their plate association | High (for READ) | Not testable without a rectal cohort; flagged tentative |
| A24 | Breast axes are continuous biology, not technical; principal-component signs are arbitrary (high/low are labels only) | Medium | Technical screen (D37); replication (D40); low ends of Basal PC3 and Her2 PC4 look like housekeeping genes |
| A25 | Cutting axes into thirds captures the axis ends; the cut point is arbitrary and patients appear in several axis strata | Low | Ends agree with the companion gene list (rho 0.82-0.93) |
| A26 | One healthy reference (the parent's GTEx tissue) suits every stratum of a cancer | Low | Same as A4 |
| A27 | Stratum signatures are dominated by the shared cancer-vs-normal signal (0.97-0.99 similar to the parent); stratum differences are a thin layer | High (for drug ranking) | Addressed by specificity scoring (D43) |

## 3. Open / parked
| Item | Status |
|---|---|
| Sensitivity checks: no ESTIMATE adjustment (A1), voomWithQualityWeights (A2), PROGENy (A6) | Proposed |
| Residual plate in BRCA-Basal and STAD | Noted |
| Blood cancers' healthy references; recount3 non-cancer diseases | Parked (PI) |
| Exact formula for sibling specificity | Proposed in `docs/notes/subgrouping_for_drug_side.md`; to confirm |
