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

### Step 4 rebuild (corrections carried to signatures and networks)
| # | Date | Decision | Why | By |
|---|---|---|---|---|
| D45 | 10-18 | ESTIMATE purity wherever possible, including LGG patient-level correction (replaces the LGG ABSOLUTE part of D22) | ESTIMATE exists for tumors and normals alike, so it is the only purity measure usable in tumor-vs-normal models; one measure everywhere. Note: ESTIMATE agrees poorly with ABSOLUTE in LGG (rho -0.28 to -0.33) | PI |
| D46 | 10-18 | Genes with mappability < 0.8 (Saha & Battle 2018) and readthrough loci removed before Step 4 and before rebuilding the ARACNe networks (1,824 of 19,168 Step 4 genes); cross-mappable edges removed | Until now these filters reached only patient-level scores, not the signatures or networks the drug side uses | PI |
| D47 | 10-18 | Before the rebuild, measure the TCGA-vs-GTEx lab difference (adjacent vs GTEx normals), GTEx post-mortem effects (ischemic time, ventilator death) and the RNA-quality index across tumors and normals together (`scripts/31_reference_diagnostics.R`) | Decide the corrections from measurements, not by default | PI |
| D48 | 10-18 | ~~Remove the shared TCGA-vs-GTEx lab difference from every GTEx-referenced signature: per gene, median adjacent-vs-GTEx logFC over 10 organ units, leave-one-unit-out, added to the GTEx normals' log-expression (gene corrected only if estimated in >= 6 units)~~ Replaced by D57 (user 10-18: the gap is not lab alone). | Healthy-vs-healthy difference is as large as a tumor signature (median 3,844 genes), shared across organs (rho 0.61) and overlaps the tumor signatures (median rho 0.62); GC/length explain only ~10% so CQN alone is not enough | PI |
| D49 | 10-18 | ~~GTEx-only covariates in Step 4 (and the lab-unit models): ischemic time (standardized within GTEx normals) and ventilator death (Hardy 0, centered within GTEx normals; only if >= 5 donors in each class); 0 for tumors~~ Replaced by D57 (ischemic time + death circumstance, corrected toward best-preserved donors). | Post-mortem effects are large in lung, esophagus/HNSC, OV and align with the tumor signature there | PI |
| D50 | 10-18 | RNA-quality index NOT added to Step 4; kept only for within-tumor / within-normal patient-level correction | Computed jointly, TCGA adjacent normals score higher than GTEx (d +0.22) while tumors score lower (d -0.70), and the index barely tracks GTEx RIN (rho 0.13): the tumor-vs-GTEx gap is not RNA integrity, so adjusting would remove tumor biology | PI |
| D51 | 10-18 | Rebuild all 28 ARACNe networks on the v2 gene set (input otherwise as v1, no covariate regression), remove cross-mappable edges after the build | Networks and signatures must use the same genes; v1 edge checks showed covariates barely affect edges (D21) | PI |
| D52 | 10-18 | ~~Patient-level scores (Step 6b) use the same lab shift on the normals, ESTIMATE purity in every group, and ischemic time + ventilator in the normals' correction~~ Lab-shift part replaced by D57; ESTIMATE everywhere stays (D45). | Patient scores are also tumor vs GTEx normals | PI |
| D53 | 10-18 | All v2 outputs in `_v2` paths; v1 kept for comparison; acceptance checks C1-C10 fixed before running and reported whatever they show | Comparable before/after; no tuning (D19) | PI |
| D54 | 10-18 | Genes without a mappability value (71) are kept and flagged | No evidence of a mapping problem | Claude (delegated) |
| D55 | 10-18 | VIPER normal-split null kept exactly as v1 (unadjusted Welch t on normals' logCPM, v2 genes) | Change one thing at a time; the lab shift cancels in the null | Claude (delegated) |
| D56 | 10-18 | Steps 8-10 and the drug side are rerun on v2 only after the checks are reviewed; the recurring-module decision (D26) is revisited first | The recurring module is now linked to the lab axis | PI |
| D57 | 10-18 | GTEx stays the healthy reference. GTEx normals are corrected from within GTEx only: per tissue, each gene modelled on log ischemic time + death circumstance (ventilator / fast / slow; Hardy 0 / 1-2 / 3-4) with age, sex, ESTIMATE; all donors kept and shifted so their average equals the best-preserved donors (ventilator, ischemic time in the tissue's shortest third, autolysis <= 1): the observed best-preserved average where the tissue has >= 45 such donors, the model's prediction where fewer; tissues with < 10 such donors (brain cortex 0, kidney 7) corrected to the average donor only, flagged. No subtraction of the TCGA-adjacent vs GTEx gap | Ischemic time + death circumstance lower held-out donor variance in about half the tissues (to ~6%); autolysis / fixative / RIN / batch add nothing. The shifted reference is closer to hidden best-preserved donors than the average-donor reference in 20 of 20 tissues (gap closed median 0.13); the plain best-preserved average is closer still where enough donors exist (median 0.18). The adjacent-vs-GTEx gap mixes lab, surgery vs death, donor population and tumor-adjacent biology, and GTEx-internal collection variation explains only 0-10% of it (scripts 33, 34) | PI |
| D58 | 10-18 | TCGA adjacent normals: check only for the disease signature (never subtracted); candidate secondary use in drug ranking under discussion | GTEx gives equally or more reproducible signatures across disjoint patient halves in 8 of 11 organs; adjacent tissue is not healthy (field effect, tumor cells; Aran 2017) | PI |
| D59 | 10-18 | GTEx colon (Colon - Transverse) and kidney (Kidney - Cortex): drop samples from a different part of the organ before any use. Colon: keep samples with mean log-CPM of 8 epithelial markers (CDX1, CDX2, VIL1, KRT20, CDH17, LGALS4, CEACAM5, EPCAM) >= 7, the valley between a mucosa mode and a muscle-wall mode and the lowest value seen in TCGA adjacent colon (165 -> 109 samples). Kidney: keep samples with mean log-CPM of 8 proximal-tubule markers (SLC34A1, LRP2, CUBN, SLC22A6, SLC22A8, SLC5A12, ALDOB, GATM) >= 5 (28 -> 21; removed: 2 near-empty samples, 5 medulla-dominant) | Colon GTEx contains muscle-wall-only samples; kidney GTEx contains medulla-dominant and failed samples. After removal GTEx variability matches TCGA adjacent: colon BCV 0.67 -> 0.41 (adjacent 0.44), PC spread 92 -> 54 (56); kidney BCV 0.78 -> 0.47 (0.50), spread 87 -> 71 (70) | PI (rule: Claude) |
| D60 | 10-18 | Kidney flagged for the PI: after D59 only 21 GTEx donors (7 best-preserved, so no best-preserved shift); KICH arises from collecting-duct cells, which the proximal-tubule rule de-emphasizes | TCGA adjacent kidney gave more reproducible signatures than GTEx (0.92 vs 0.87) before D59 | Open (PI) |
| D61 | 10-18 | Measure tissue composition per cohort (tumors vs v2 GTEx reference vs TCGA adjacent): MCP-counter immune / stromal scores + Human Protein Atlas single-cell marker panels; share of each signature explained by cell-type membership. Spec: [`SPEC_tissue_composition.md`](SPEC_tissue_composition.md); not executed | Signature genes that only reflect a different cell mix in normal tissue (fat in breast, muscle in uterus) cannot be reversed by a drug acting on tumor cells | PI |
| D62 | 10-18 | Project all v2 GTEx reference donors (and TCGA adjacent normals) through Celligner per cohort: line weights, lineage composition, shared weight with the tumors, organ-driven vs cancer-specific lines. Spec on the `celligner` branch: `docs/notes/SPEC_celligner_normals.md`; not executed | Shows how much of each disease's cell-line match reflects the organ rather than the cancer | PI |
| D63 | 10-18 | Network build questions settled before the v2 build, prompted by drug-side findings (proliferation carries most known-drug signal, C15; ESR1 regulon = luminal identity, C3/C14) | A regulator list or sample pooling cannot be changed after the build without rebuilding | PI |
| D64 | 10-18 | Build two v2 network sets: the full regulator set as D14 (5,789 candidates) and a TF-only set (Lambert TFs, ~1,150-1,360 per cohort), same inputs and settings | Restricting a full build to TFs afterwards is not the same as a TF-only build (DPI makes regulators compete for targets); TF-only is comparable with CollecTRI and avoids non-TF proliferation 'regulators' (MKI67, TOP2A, CKS1B); costs ~0.7 extra days | PI |
| D65 | 10-18 | Add a luminal-only breast network (BRCA_LUM: LumA + LumB, 769 tumors) next to the pooled BRCA network | In the pooled network the ESR1 regulon reflects luminal-vs-basal identity, so anti-estrogens read wrongly on the drug side; within luminal tumors ESR1 targets should reflect ER signaling. Her2 (81) and Basal (190) stay pooled (too small alone) | PI |
| D66 | 10-18 | No proliferation regression of the network input; keep the v2 gene filter (D46) although it removes 28 of 978 L1000 landmark genes (e.g. FOXO3, TFDP1, CHEK2, HSPA1A) and 412 candidate regulators (100 TFs; none of the known-biology controls) | E2F / FOXM1 regulons are real biology and pass controls; proliferation is handled at drug scoring (79i). The lost regulators are mostly artefacts (AGAP4/5/6, readthrough) or proliferation genes (CKS1B, CKS2, UHRF1) | PI |
| D67 | 10-18 | Composition marker panels from cell-type families: 15 fixed families of closely related HPA types (e.g. enterocytes + Paneth, keratinocytes + squamous epithelium, alveolar type 1 + 2, blood + lymphatic endothelium, fibroblasts + endometrial / ovarian stroma) scored as one unit, markers >= 4x every type outside the family; positive controls renamed to the families | Under the single-type rule 9 of 18 controls had no panel; with families 51 of 56 units get a panel. Smooth muscle still has none (its genes are shared with HPA Sertoli, skeletal muscle and myoepithelial cells) | User |
| D68 | 10-20 | Subtype networks: only the breast luminal test (BRCA_LUM vs pooled BRCA master regulators, cutoff for a material change fixed before the comparison); other subtypes stay on the shared cohort network. If luminal regulators change materially, build networks for the other subtypes with >= 100 tumors (BRCA LumA / LumB / Basal, CESC-SCC, LGG IDHmut-codel / -noncodel, PCPG-PHEO, UCEC-EEC, SARC-LMS). Strata found with the networks (COAD / READ) tested only after v2 Steps 8-10, against size- and split-matched control networks. No synthetic patients to reach the ARACNe sample size | Subtypes are assumed to share wiring and differ in regulator use, which VIPER measures on a shared network (the colon strata were found from activity differences on the shared network). A network built on part of a cohort differs even with identical wiring (less variation, fewer patients), so a difference alone is not rewiring. Resampled or generated patients repeat the real ones and inflate mutual information | User |
| D69 | 10-20 | Step 4 v2 uses voomByGroup (separate mean-variance trends for tumors and normals) instead of voom; the voom version is kept as check C12 and the 13-group test results are kept (`results/step4/voom_group_variance/`) | Tumors are 1.2-1.8x more variable than normals (A2 wrong in 12 of 13 groups); split-half reproducibility of z and NES slightly higher (median +0.001), controls unchanged. Whole rankings barely move (rho >= 0.996, below the pre-set threshold for a practical effect), but the top-500 genes and the FDR-cut regulator lists that the drug side uses shift 10-30% (READ 157 -> 230, THCA 1 -> 55), which the user judged significant | User |
| D70 | 10-20 | Step 4 v2 signatures accepted as built. The larger GC tilt (C1 median -0.24 -> -0.30) needs no correction; C2 still reports whether it reaches the regulators. Pediatric groups (TARGET-WT, NBL) are recorded as more variable (A36) and left as they are | A GC difference between healthy donors and tumors may be partly real and is not by itself a problem. Corrections are made only when they apply to all cohorts: the goal is an accurate, reproducible workflow, not fixes for single edge cases (WT genes changed 2,801 -> 6,916 when kidney normals went 28 -> 21) | User |

Full rebuild spec: [`SPEC_rebuild_v2.md`](SPEC_rebuild_v2.md).

## 2. Assumptions

What the results rely on. "Impact" = how much results could change if the assumption is wrong. "Check" = what we
did or could do about it.

| # | Assumption | Impact | Check |
|---|---|---|---|
| A1 | Purity (stroma/immune content) is a nuisance to remove, but it may carry real mechanism (immune evasion, fibrosis) | High | Sensitivity run without ESTIMATE adjustment: proposed, not run |
| ~~A2~~ | ~~limma treats tumors and normals as equally variable; tumors are more variable~~ | Resolved (D69) | Tested 2026-10-20: tumors 1.2-1.8x more variable; v2 models the arms separately (voomByGroup) |
| A3 | Gene-gene correlation in healthy tissue represents "chance" for tumors (normal-split null); the null is unadjusted for covariates while the signature is adjusted | High | Conservative in practice (GTEx donor structure widens it) |
| A4 | GTEx is a fair healthy reference despite post-mortem collection and source being collinear with tumor vs normal; adjacent normals carry field effects | High | GTEx vs adjacent agreement median ~0.6 (18 groups) |
| A5 | Partial GTEx tissue matches add cell-type differences to signatures (pancreas acinar vs ductal, uterus myometrium, adult adrenal for NBL) | High | Flagged per group |
| A6 | Regulator activity can be read from its targets' mRNA; misses non-transcriptional activation (kinases) | Medium | PROGENy / CPTAC phospho proposed |
| A7 | The expressed-regulator filter drops active but low-mRNA regulators | Medium | None yet |
| A8 | CollecTRI is context-free; unsigned edges (~1.5k of 64.5k) are treated as activation | Medium | CollecTRI is supporting only |
| A9 | ARACNe: co-expression = regulation; one network per cancer, shared by its subtypes and strata | Medium | Edge confounding check: < 3% of edges lose > 50% strength after purity/quality/plate; BRCA_LUM vs BRCA tests the subtype part (D68) |
| A10 | Subtype labels (PanCanAtlas PAM50, IDH/codel, marker-based ALL lineage; METABRIC CLAUDIN_SUBTYPE; CPTAC PAM50) are correct | Medium | Taken as given |
| A11 | RSEM expected counts behave like counts for voom | Low | Standard practice |
| A12 | Patients with missing covariates (< 2%) are missing at random | Low | Dropped |
| A13 | GTEx age = middle of the 10-year bin | Low | |
| A14 | The diagnosis sample represents the patient; metastasis/recurrence handled by an indicator | Low | |
| A15 | BH FDR is valid under positive dependence between regulators | Low | |
| A16 | Filter thresholds (>= 10 targets, >= 50% present, logCPM >= 1, >= 30 patients, >= 10 normals) are judgment calls | Low | Sensitivity analyses advisable |
| A17 | The long-minus-short gene index measures RNA quality | Medium | Agrees only partly with GTEx RIN / ischemic time; treated as a plate-linked technical axis. 10-18: computed jointly, TCGA normals score higher than GTEx and tumors lower; not used in Step 4 (D50) |
| A18 | The recurring module is biology (PI reading), not a technical axis | Medium | Plate-linked above chance in 39 groups; RNA quality, GC bias and multi-mapping ruled out as causes; cause unknown. 10-18: 89% of its 555 members (all 17 core) are higher in GTEx than in TCGA normals, i.e. it matches the TCGA-vs-GTEx lab axis; revisit before Step 8 v2 (D56) |
| A19 | Regulator modules (correlation >= 0.7) are meaningful units, and module scores capture patient differences | Medium | Mean-member and PC1 scores agree (r 0.998) |
| A20 | Comparing replication patients with their cohort's average tumor is equivalent to the TCGA comparison with healthy tissue once scores are standardized | Medium | Colon: TCGA modules replicate (rho 0.89-0.94); breast axes replicate in gene content |
| A21 | Array (GSE39582, METABRIC) and RNA-seq scores are comparable after standardizing within each cohort | Medium | Same patterns found on both platforms |
| A22 | CMScaller calls on TCGA/replication data are correct, and the CMS1/3 vs 2/4 naming describes the colon split | Low | Consistent across three cohorts |
| A23 | READ clusters reflect biology despite their plate association | High (for READ) | Not testable without a rectal cohort; flagged tentative |
| A24 | Breast axes are continuous biology, not technical; principal-component signs are arbitrary (high/low are labels only) | Medium | Technical screen (D37); replication (D40); low ends of Basal PC3 and Her2 PC4 look like housekeeping genes |
| A25 | Cutting axes into thirds captures the axis ends; the cut point is arbitrary and patients appear in several axis strata | Low | Ends agree with the companion gene list (rho 0.82-0.93) |
| A26 | One healthy reference (the parent's GTEx tissue) suits every stratum of a cancer | Low | Same as A4 |
| A27 | Stratum signatures are dominated by the shared cancer-vs-normal signal (0.97-0.99 similar to the parent); stratum differences are a thin layer | High (for drug ranking) | Addressed by specificity scoring (D43) |
| A28 | ~~(v2) The part of the adjacent-vs-GTEx difference shared across organs is lab/protocol + average post-mortem effect; some may be real biology of tumor-adjacent tissue (field effect, tumor cells, wound response) and is removed with it~~ (lab subtraction dropped, D57) | High | Acceptance checks C5 (known biology), C6 (proliferation retained) |
| A29 | ~~(v2) One lab correction per gene (median of 10 organ units) fits organs without adjacent normals; organ-specific lab effects remain~~ (lab subtraction dropped, D57) | Medium | C4 (held-out unit) |
| A30 | ~~(v2) Uncertainty of the lab correction is not carried into the moderated t, so z is slightly over-confident~~ (lab subtraction dropped, D57) | Low | z used as ranks |
| A31 | ~~(v2) Ischemic time and ventilator death capture GTEx post-mortem variation linearly; other agonal factors not modelled~~ (superseded by A33) | Medium | Ischemic effect agrees poorly across tissues (rho 0.13) |
| A32 | (v2) The 71 genes without a mappability value behave like normally mappable genes | Low | Flagged |
| A33 | (v2) GTEx collection effects are additive per gene in log ischemic time and death circumstance; other agonal / handling factors are unmodelled | Medium | Held-out variance (script 34): autolysis, fixative time, RIN and batch add nothing beyond M1 |
| A34 | (v2) Best-preserved GTEx donors (ventilator, short ischemia) are the closest available stand-in for living tissue; they are still post-mortem, and surgical tissue has its own warm-ischemia changes that TCGA did not record | High | Moving GTEx toward them closes only 0-10% of the TCGA-adjacent vs GTEx gap; remaining gap = lab + surgery vs death + donor population + tumor-adjacent biology, not separable with available data |
| A35 | (v2) Marker-score thresholds pick the right part of the organ in GTEx colon (epithelium >= 7) and kidney (proximal tubule >= 5) | Low | Thresholds sit at the valley of a bimodal distribution; variability after removal matches TCGA adjacent |
| A36 | (v2) Pediatric tumors (TARGET-WT, NBL) can be compared with adult GTEx donors; age is nearly collinear with tumor vs normal, so the age adjustment and the condition estimate trade off | Medium (for WT / NBL) | Recorded 10-20 (D70): WT genes changed 2,801 -> 6,916 between v1 and v2; no cohort-wide correction needed, so no action |

## 3. Open / parked
| Item | Status |
|---|---|
| Sensitivity checks: no ESTIMATE adjustment (A1), PROGENy (A6) | Proposed (A2 tested 2026-10-20, D69) |
| Residual plate in BRCA-Basal and STAD | Noted |
| Blood cancers' healthy references; recount3 non-cancer diseases | Parked (PI) |
| Exact formula for sibling specificity | Proposed in `docs/notes/subgrouping_for_drug_side.md`; to confirm |
| ~~Step 4 rebuild: how to remove the TCGA-vs-GTEx lab component, whether to adjust GTEx normals for ischemic time / ventilator death, RNA-quality index in Step 4~~ | Decided 10-18: D48-D55; spec in SPEC_rebuild_v2.md; not yet executed |
| Recurring module in clustering (D26) given its link to the lab axis | Open (PI), before Step 8 v2 |
| ~~Composition marker rule: 9 of 18 positive controls not evaluable; proposed: families~~ | Decided 10-18 (D67). Remaining: smooth muscle has no panel (5 controls); curated panel proposed | Open (user) |
| Using TCGA adjacent normals in drug ranking | User 10-18: develop (1) a safety filter (target higher in tumor than in the patients' own adjacent normal tissue; normal-tissue impact) and (2) the inflammation around tumors as a secondary target (with surgery / immediate-early genes separated out). Sketch in the drug-side notes (`celligner` branch); not designed yet |
| ~~D48 premise under review (user 10-18): the adjacent-vs-GTEx difference is not only lab. It mixes lab/protocol, collection (surgery vs post-mortem), donor population (cancer patients, ~12 years older) and tissue state next to a tumor (field effect, inflammation, tumor cells). Removing it removes all four. Reference choice pending the variance / split-half test (`scripts/33_normal_variance.R`) ~~ | Resolved 10-18: D57-D60 |
