# Glossary: terms added or changed in v2

Status 2026-11-02. The full glossary (every term and acronym of Steps 1-7, grouped by topic) is the Glossary tab of the
Steps 1-7 specification doc. This file lists what v2 changed: entries that are now out of date, and terms v2 introduced.
How each step runs: [`PIPELINE_V2.md`](PIPELINE_V2.md); decisions: [`DECISIONS.md`](DECISIONS.md).

## Entries that changed

| Term | Old entry said | Now (v2) |
|---|---|---|
| voomWithQualityWeights | Proposed fix for unequal tumor vs normal variability (A2) | Not used; A2 was resolved with **voomByGroup** (D69) |
| Master regulator (MR) | Median ~600 per group | v1 figure. v2 TF-only arm, corrected null: median 151 at FDR < 0.10, 116 at FDR < 0.05 |
| Exploratory / confirmed (D18) | Confirmed: FDR < 0.05, or < 0.10 plus a same-direction CollecTRI / DoRothEA call | Rule unchanged, but CollecTRI and DoRothEA rarely reach significance (median 0-1 per group), so confirmed is in practice ARACNe FDR < 0.05 (median 117 per group). Whether to keep the rule is open (Step 7 question) |
| Normal-split null (D15) | 1,000 random halvings of the healthy samples | Same idea; since D73 the halvings use exactly the Step 4a GTEx samples the signature used. A covariate-adjusted version was tested and not adopted |
| Group level | Group scores feed the drug side | Every tested regulator feeds the drug side, weighted by its group NES; no MR cutoff (D74) |
| Module | 37-112 modules of >= 10 per group | v1 (full regulator set). v2 TF-only arm: median 18 per group, 15 regulators each |
| Recurring module (D26) | In 39 of 41 groups; open | v2: in 36 of 41 groups; subgroups are defined from the clustering without it (D79); PI to confirm D26 |
| Known-biology check | 6 of 16 tissue-specific checks pass | v1. v2 TF arm (all controls incl. proliferation): 77 pass, 55 right direction below cutoff, 4 wrong direction (all near zero) |
| Network | 28 in v1, 29 per set in v2 | v2 adds BRCA_LUM; the active build uses the TF-only set |
| z-score / NES | NES is on the same scale as z | Roughly: NES is normalized against the null, but its spread differs between regulon sets (CollecTRI's is narrower) |

## Terms introduced in v2

### Data preparation and signatures
| Term | Meaning |
|---|---|
| v2 gene set | 17,784 protein-coding genes: genes with mappability < 0.8 and readthrough loci removed (D46) |
| Readthrough locus | A transcript that runs across two neighboring genes; removed because its reads mix both genes |
| Step 4a reference | The prepared GTEx healthy samples: wrong-part-of-organ samples removed (D59), collection covariates added (D57) |
| Anchor | Per GTEx tissue, the reference point the signature compares tumors with: the best-preserved donors, modelled or observed, or the average donor where too few exist |
| Collection covariates (`coll_isch`, `coll_fast`, `coll_slow`, `coll_nonbp`) | Ischemic time, death circumstance and non-best-preserved status, centered at the anchor, 0 for tumors; in the Step 4 model |
| voomByGroup | voom with a separate mean-variance trend for tumors and normals (D69) |
| Build | A named set of paths for one pipeline version (`config/builds.yaml`): `v1`, `v2_full_xmap`, `v2_tf_xmap` (active) |

### Networks and regulator activity
| Term | Meaning |
|---|---|
| TF-only set | ARACNe networks built with Lambert TFs as the only candidate regulators (D64); the active arm |
| Identity / signaling edge (D82) | Identity: an edge that exists only because both genes differ between subtypes (keeps < 50% of its correlation once subtype is removed); signaling: it remains within subtypes (`48`) |
| Subtype-specific edge | Strong in one subtype (absolute rho >= 0.3) and absent or reversed in another: evidence that wiring differs (A9) |
| Signaling-only lists: base / sig1 / sig2 (D83) | TF-only target lists for a comparison: base = current; sig1 drops identity-only edges under histology labels; sig2 also under molecular subtypes. Rule met, effect small; not adopted (open) |
| Lineage-matched line (C24) | A cell line of the disease's cell type (squamous, glandular, ...), any organ (`config/group_lineages.tsv`, celligner) |
| Data gap (R3) | The first step where good-quality lineage-matched data run out; reported instead of a drug list (`79r`) |
| Cross-mappable edge filter | Edges between genes that share reads (cross-mappability > 0) removed before building regulons, 11-15% of edges |
| Fix B / adjusted null (C) | The two null variants tested in `64`: B uses the Step 4a normals (adopted), C also removes covariates (not adopted) (D73) |
| Prediction strength | How well clusters found in one half of the patients predict those found in the other half (Tibshirani & Walther 2005); keep threshold 0.75 (D71) |
| Keep rule (D71) | M3C significant, prediction strength >= 0.75, smallest cluster >= 30, feature sets agree (ARI >= 0.3), plate excess < 0.1, purity / quality R^2 < 0.1 |
| `modules10` / `all` | Step 8 feature sets: one score per module of >= 10 regulators / every module and single regulator |
| Subgroup (stratum) | A patient cluster kept by D71 / D79 that gets its own signature, VIPER and drug list; 14 in v2 |
| Coherence | How closely a patient's regulator activity matches their subgroup's average (`86`) |
| Small program | A module of 2-9 regulators or a single regulator, present only in the `all` feature set (`87`) |

### Drug side
| Term | Meaning |
|---|---|
| Direct / borrowed drug | Tested in the disease's tissue / tested only in >= 3 other cancer lines |
| Consensus regulons (D75) | Per regulator, the mean of the ARACNe TF and CollecTRI activity, over the sets scored on both the disease and drug side |
| `rev` (reversal) | How strongly a drug pushes the disease's regulators back toward healthy: minus the NES-weighted sum of drug activity (D74) |
| `spec` (specificity) | `rev` relative to the same drug scored against diseases from other tissues |
| `rev_mech`, `spec_mech` | The same with proliferation-linked regulators removed and activity adjusted for the drug's proliferation effect; reports lead with `spec_mech` (D80) |
| `spec_sib`, `spec_sib_mech` | Subgroups: specificity against the sibling subgroup |
| Line match / `unmatched_lines` (D76) | Celligner weight of the drug's best tested line relative to the disease's best line; < 0.1 = flagged |
| Transfer confidence / `low_transfer` (D77) | Predicted accuracy of a borrowed drug's effect relative to a replicate; < 0.5 = flagged |
| Evidence status (D78) | Per disease: matched / borrowed only / no adequate compound, with a sentence for reports |
| Knockdown benchmark | LINCS shRNA / CRISPR TF knockdowns used to test which regulon set reads inhibition correctly; reported, not used (C19) |
| Drugs@FDA | FDA's database of approved drugs; source of FDA status and first approval date (`79p`) |
