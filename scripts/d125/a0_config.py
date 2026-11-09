"""D125 section A0 (noise ceiling): the one config block for every A0 script (docs/SPEC_graph_metric_scoring.md, D125-A0).

Every threshold below is pre-registered in the spec or in D125-A0 (docs/DECISIONS.md). Do not tune them after results
are visible; a threshold that looks unworkable is reported, not changed. Each A0 script imports the names it uses at its
top, so the values live here only.
"""

# --- signature quality (D125-A0 (1); the signatures 79c averages) -------------------------------------------------
TAS_MIN = 0.4                 # active = is_hiq == 1 and TAS >= 0.4 (79c, user 2026-10-09)

# --- A0.0 eligible-pair gate (spec A0.0) ----------------------------------------------------------------------------
A00_MIN_REPLICATES = 6        # active level-5 signatures per drug x line pair (D125-A0 (1))
A00_MIN_PAIRS = 500           # fewer qualifying pairs -> A0b is a line-level bound
A00_MIN_LINES = 3             # qualifying pairs in fewer lines -> A0b is a line-level bound

# --- A0a paired-signature reliability (spec A0a, D125-A0 (2)-(3)) ---------------------------------------------------
A0A_STRATUM = ("pert_dose", "pert_dose_unit", "pert_time")   # independent stratum = (dose, time) cell
A0A_MIN_STRATA = 2            # pairs with >= 2 cells are split into pass A / pass B
A0A_SEED = 125                # draws the two cells per pair

# --- A0b split-half reliability (spec A0b) --------------------------------------------------------------------------
A0B_BATCH = 10                # splits per batch
A0B_FLOOR = 20                # minimum splits
A0B_CAP = 100                 # maximum splits; reached without stabilising -> ceiling indicative only
A0B_SE_STABLE = None          # NOT pre-registered: the spec says "until the SE of the median r_c stabilises" /
                              # "SE small relative to 0.05" without a number. Execution refuses while None; set it only
                              # through a recorded decision (D125-A0 amendment), never after seeing results.
A0B_SEED = 1250               # split s uses seed A0B_SEED + s

# --- A0c / A0d reporting --------------------------------------------------------------------------------------------
R_POOLED_MIN_DRUGS = 50       # fewer drugs in a cohort's r_c -> corrected with the median, flagged r_pooled
R_DIFF_NOTE = 0.1             # |r_A0b - r_A0a| >= this -> stated wherever a ceiling-relative number is reported
CI_LEVEL = 0.95               # r_c CI: Fisher z, Bonett-Wright variance (A0a); split percentiles (A0b)

# --- harness checks (D125-A0 (4); implementation, not analysis thresholds) -------------------------------------------
PROFILE_IDENTITY_TOL = 1e-4   # rebuilt production profiles vs drug_line_profiles_bing[_specific].f32 (z units)
VIPER_IDENTITY_TOL = 1e-3     # 79d on rebuilt columns vs production VIPER columns (NES units)
SCORE_IDENTITY_TOL = 0.1      # sandbox scoring of production inputs vs production spec_mech_pct_pregate (pct points)

# --- compute (user: < half of the 64 threads) -----------------------------------------------------------------------
THREADS = 8                   # total budget
PASS_PARALLEL = 2             # A0b passes run concurrently; each gets THREADS // PASS_PARALLEL threads
