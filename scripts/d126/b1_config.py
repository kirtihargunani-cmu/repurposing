#!/usr/bin/env python
"""D126 run-order B1/B3 config: every pre-registered constant in one block.

From the spec (docs/SPEC_graph_metric_scoring.md B3, B5) and from D126-B1. Fixed before running; never adjusted
because a result looks wrong. Label note: the spec's run order is B0-B3 while its sections are B1-B6; run-order B1 =
section B3's primary method, run-order B2 = the breast CARNIVAL arm, run-order B3 = section B5's criteria.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# --- spec B3: the quantity ----------------------------------------------------------------------------------------
TOP_K = 25                  # top regulators by |activity| in the drug's scoring-line consensus profile
# no hop cap: the spec sets none and inventing one would add a threshold (D126-B1 (4))

# --- D126-B1 (6): the flag. Neither number is tuned -----------------------------------------------------------------
FLAG_COVERAGE_MIN = 0.4     # = the B5.3 floor, reused so the flag cannot fire on near-empty evidence
FLAG_CONSISTENT_MAX = 0.5   # the no-information point of a sign comparison

# --- spec B5: the criteria ------------------------------------------------------------------------------------------
B51_CONTROL_MIN = 0.5       # D126-B1 (7): the flag must fire on > 50% of scorable anti-estrogen rows in MCF7 breast groups
B52_RHO_MAX = 0.3           # |Spearman(flag rate, main-target out-degree)| at or above this triggers B4
B53_COVERAGE_MIN = 0.4      # mechanism_coverage floor ...
B53_SHARE_MIN = 0.6         # ... met by at least this share of gate-passing FDA drugs

ANTI_ESTROGEN_TARGETS = {"ESR1", "ESR2", "CYP19A1"}   # the class as D127 defines it
MCF7_LINE = "MCF7"

# --- paths ------------------------------------------------------------------------------------------------------------
PRIOR = ROOT / "data/raw/omnipath/omnipath_signed_directed.tsv.gz"
PRIOR_SIGNOR = ROOT / "data/raw/omnipath/signor_signed_directed.tsv.gz"
ACT_DIR = ROOT / "data/processed/lincs/d126_act"      # 79e D126_ACT_DIR dump
OUT = ROOT / "results/lincs/d126"
