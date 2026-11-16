"""D128 (target-mechanism scoring): the one config block for every D128 script.

Values follow docs/IMPLEMENTATION_target_mechanism_scoring.md Section 0 (P1-P14). **FROZEN 2026-11-15** by the
`D128-A0 freeze` entry in docs/DECISIONS.md (tag `d128-a0-freeze`), after the label-blind A0 audit, the A0.1 rule
correction and the user's D0.1 / D0.2 decisions. From that tag onward a change to any value here is a dated amendment
row in DECISIONS.md, never an edit of the frozen entry. Each script imports the names it uses, so the values live
here only.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
FROZEN = True                       # frozen by the D128-A0 freeze entry / tag d128-a0-freeze (2026-11-15)

# --- inputs (read-only) -------------------------------------------------------------------------------------------
LINCS = ROOT / "data/processed/lincs"
PROFILES_SPECIFIC = LINCS / "drug_line_profiles_bing_specific.f32"   # P5 primary
PROFILES_RAW = LINCS / "drug_line_profiles_bing.f32"                 # P5 sensitivity
PROFILE_INDEX = LINCS / "drug_line_profiles_index.tsv"               # rows: pert_id, ModelID, n_sig
PROFILE_GENES = LINCS / "drug_line_profiles_genes.tsv"               # cols: gene_id, symbol, feature_space
DRUG_ANNOT = ROOT / "results/lincs/drug_annotations.tsv"
LAYERED = ROOT / "results/lincs/layered_score_v2.tsv"                # production lead rows (population)
SPLIT_GROUPS = ROOT / "results/lincs/drug_reversal_split_v2_consensus.tsv"
SPLIT_STRATA = ROOT / "results/lincs/stratum_reversal_split_v2_consensus.tsv"
UNITS = ROOT / "results/lincs/target_role/d121/units.tsv"           # D121 target units (79t2 chembl_units)
LINE_W_GROUPS = ROOT / "results/celligner/disease_line_weights_long.tsv"
LINE_W_STRATA = ROOT / "results/celligner/stratum_line_weights_long_v2.tsv"
COLLECTRI = LINCS / "precursor_regulons/collectri.tsv.gz"           # production regulons_core.rds collectri export
OMNIPATH = ROOT / "data/raw/omnipath/omnipath_signed_directed.tsv.gz"  # D126 cache, primary signaling prior (P6)
SIGNOR = ROOT / "data/raw/omnipath/signor_signed_directed.tsv.gz"      # P6 sensitivity
DE_GROUPS = (ROOT / "data/processed/step4_v2", ROOT / "results/step4_v2/signature_summary.tsv")
DE_STRATA = (ROOT / "data/processed/step10_v2", ROOT / "results/step10_v2/signature_summary.tsv")
DEPMAP = ROOT / "data/raw/depmap_24q4"
TCGA = ROOT / "data/raw/tcga_pancan"

# --- Section 0: ALL ROWS FROZEN (D128-A0 freeze, 2026-11-15) ---------------------------------------------------
DE_SIGN_COL, DE_WEIGHT_COL, DE_P_COL = "logFC", "z", "adj.P.Val"   # P1 settled (limma + voomByGroup)
CONTRAST_ORDER = ("gtex", "adjacent")                               # P2 settled: 79t rule
NETWORK = "collectri"                                               # P3
FEATURE_SPACE = "all"                                               # P4: all BING genes; landmark-only arm
PROFILE_ARM = "specific"                                            # P5
MAX_HOPS = 2                                                        # P6 (sensitivity 0, 1, 3)
HOP_ARMS = (0, 1, 2, 3)
WM = "unweighted"                                                   # P7
R_P_MAX = 0.05                                                      # P8: R gene set adj.P.Val < 0.05 (median 7,930 genes/cohort)
F_DIVISOR = 2.0                                                     # P9: f(z) = tanh(z/2); A0 kept it (saturation 1.3%)
FLOOR_GENES = 10                                                    # P10: resolved genes, per subset (D128-A0.1)
FLOOR_TFS = 2                                                       # P10: resolved TFs, signaling subset only
FLOOR_MULTIPLIERS = (0.5, 1.0, 2.0)                                 # P10 reported arms
M_AMBIGUOUS_FLAGS = ("direction_assumed", "mixed")                  # P11: no default sign for M

# --- D0.1 / D0.2 (user, 2026-11-15) --------------------------------------------------------------------------------
R_BREAST_REFERENCE = "production_P2"        # D0.1(a): no special breast reference; (b) adjacent and (c) LumA/LumB-vs-
                                            # Basal are reported-only arms. Breast R predicted uninformative on
                                            # CollecTRI, inverted on ARACNe; the control is gated on M alone.
ESR1_PRIMARY = ("fulvestrant",)             # D0.2(a): the primary positive set
ESR1_SERM_SET = ("tamoxifen", "raloxifene", "bazedoxifene", "toremifene", "clomifene")   # reported separately
SERM_LINEAGE = "Breast"                     # SERMs count as ESR1 antagonists in breast-lineage lines only
ESR1_EXCLUDED_UNITS = ("ESRRG",)            # tamoxifen's ESRRG unit must not stand in for ESR1

# --- A0 audit -------------------------------------------------------------------------------------------------------
PANEL = ("BRAF", "EGFR", "KRAS", "NRAS", "ERBB2", "ESR1", "PIK3CA", "RB1", "BRCA1", "BRCA2")
Z_QUANTILES = (0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
Z_THRESHOLDS = (1.0, 2.0, 4.0)
CONTROLS = {   # named in the plan (Sections 4-5); availability only, no outcome
    "ESR1": ["fulvestrant", "tamoxifen", "raloxifene", "bazedoxifene", "toremifene", "clomifene", "estradiol"],
    "ERBB2": ["lapatinib", "neratinib", "afatinib", "canertinib", "tucatinib"],
    "BRAF": ["vemurafenib", "dabrafenib", "encorafenib", "PLX-4720", "SB-590885"],
    "MEK": ["trametinib", "selumetinib", "cobimetinib", "binimetinib", "PD-0325901"],
    "EGFR": ["gefitinib", "erlotinib", "afatinib", "osimertinib", "dacomitinib"],
    "CDK4/6": ["palbociclib", "ribociclib", "abemaciclib"],
    "PARP": ["olaparib", "niraparib", "rucaparib", "talazoparib", "veliparib"],
}
OUT = ROOT / "results/lincs/d128"
SEED = 128
THREADS = 8                                                         # < half of the 64 threads
