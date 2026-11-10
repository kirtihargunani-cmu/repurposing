#!/usr/bin/env python
"""D125 A1 config: every pre-registered constant of the lambda sweep, in one block.

Thresholds come from the spec (docs/SPEC_graph_metric_scoring.md A2, A3) and from D125-A1. They are fixed before
running and are never adjusted because a result looks wrong. If one looks unworkable, stop and report.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# --- spec A3: the grid and the early stop -----------------------------------------------------------------------
LAMBDAS = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
MIN_GAIN_PTS = 1.0          # a lambda must beat the lambda = inf arm by >= this many percentage points, else D125 is
                            # recorded as not met and A2 does not run

# --- spec A2: the graph ------------------------------------------------------------------------------------------
EDGE_MIN = 0.1              # Szymkiewicz-Simpson overlap below this is dropped
# no self-loops; sign ignored; L = I - D^-1/2 A D^-1/2; isolated nodes keep D_ii = 0 (D125-A1 (5))

# --- D125-A1: what is read, and on which slice --------------------------------------------------------------------
SETS = ["aracne_v2_tf", "collectri"]    # the D75 consensus; dorothea_ABC is scored by 79n but is not in the score
PRIMARY_ANALYSIS = "main"               # D84 main analysis (own_ok): the selection slice
PRIMARY_READ = "consensus"              # mean of the two sets where both score the TF, else the one (D85 criterion 1)

# --- identity check (D125-A1.1): the lambda = inf arm must reproduce 79n's own output on the CURRENT inputs ---------
IDENTITY_TOL = 1e-4         # max |diff| in NES and in percentile vs the d84_signatures table of the same 79n run
# Reported only, not a gate: the published D84 main-analysis numbers (results/lincs/d84_summary.tsv). They predate the
# 2026-11-03 rebuild of regulons_aracne_project_v2_tf.rds, which added the 12 subtype networks, so ARACNe no longer
# matches them (kd 0.597 -> 0.594, oe 0.690 -> 0.694); CollecTRI, unchanged since 10-01, still does.
D84_PUBLISHED = {("aracne_v2_tf", "kd"): 0.597, ("aracne_v2_tf", "oe"): 0.690,
                 ("collectri", "kd"): 0.604, ("collectri", "oe"): 0.701}

# --- paths and compute --------------------------------------------------------------------------------------------
NES_DIR = ROOT / "data/processed/lincs/d125_a1"     # 79n D125_NES_DIR dump; nothing at the top of results/lincs
OUT = ROOT / "results/lincs/d125/a1"
THREADS = 8
