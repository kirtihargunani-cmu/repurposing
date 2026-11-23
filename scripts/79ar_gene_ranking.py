#!/usr/bin/env python
"""Production ranking (D211): filter on target evidence, rank by gene-level reversal, per timepoint.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-23: "use the dependency, alteration, and overactivity as a filter and then calculate the
gene-level reversal score on the eligible (drug, disease) pairs. make sure that everything is
time-resolved." Choices made by the user the same day: filter = combined T >= 50; chemotherapy /
essential-gene / DNA-acting rows reported in a separate table; adopted as production directly (D127 is
reported on it as a description, not a gate); subgroups get their own window runs.

STAGE 1 - FILTER (tumour-side, so it has no timepoint: dependency, alteration and overactivity describe
the disease, not the drug's time course). A (unit, compound) row is ELIGIBLE when
    layered_T >= T_MIN (50)  and  layered_n_axes >= 1
layered_T (79t3) is the D/A/O mean, orientation-aware (D177: dependency and alteration are dropped for
`raise` units, overactivity only counts for them when it is a measured NES) and shrunk toward 50 by the
number of axes present (D177c). Rows scored through the essential-unit or DNA-acting path have no D/A/O
(layered_n_axes is empty) and cannot pass a D/A/O filter by construction: they go to the separate table.

STAGE 2 - RANK, per (unit, window), window in {early (<= 6 h), late (>= 24 h)} (D180). The score is
    rev_gene = - sum_i sum_j d_ij * C_ij / ||s||          (docs/FORMULAS.md s3, D203)
summed from 79e3's per-regulator `rev_gene_contrib` (groups: ALL=1 table; subgroups: STRATA=1 ALL=1).
gene_pct = percentile of rev_gene among the ELIGIBLE compounds of that (unit, window); higher = the drug's
gene-level response opposes the disease pattern more. Early and late are ranked separately, never pooled.

ONE COMPOUND PER INGREDIENT, per (unit, window): an FDA ingredient often has several LINCS compounds. The
representative is the compound with the most conditions in that window (ties: pert_id order) - NOT the
best-scoring one, which would select on the outcome.

NO TIMEPOINT -> `unknown` with the D180 reason carried from 79j (no_window_run / no_condition_profile /
no_condition_in_lines). Those rows stay in the table, unranked; nothing is pooled across windows.

Outputs (results/lincs/):
  gene_ranking<SFX>.tsv        eligible + ineligible D/A/O rows, one row per (unit, ingredient, window)
  gene_ranking_chemo<SFX>.tsv  essential-unit / DNA-acting rows, same windowed columns, ranked among themselves
  gene_ranking_run<SFX>.json
Env: SFX (_v2_tf) | T_MIN (50)
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402

RES = ROOT / "results/lincs"
SFX = os.environ.get("SFX", "_v2_tf")
if "consensus" in SFX:
    raise SystemExit("D208: the consensus arm has been removed")
T_MIN = float(os.environ.get("T_MIN", "50"))
CONTRIB = RES / f"regulon_contributions{SFX}"
SOURCES = {"group_id": CONTRIB / "regulator_contributions_by_window_all.tsv.gz",
           "stratum": CONTRIB / "regulator_contributions_by_window_all_strata.tsv.gz"}
SPLITS = {"group_id": RES / f"drug_reversal_split{SFX}.tsv", "stratum": RES / f"stratum_reversal_split{SFX}.tsv"}
WINDOWS = ("early", "late")


def window_scores(path):
    """(unit, pert_id, window) -> rev, rev_gene, n_conditions, n_lines, summed from the regulator table."""
    if not path.exists():
        return pd.DataFrame(columns=["unit", "pert_id", "window", "rev", "rev_gene", "n_conditions", "n_lines"])
    acc = []
    cols = ["cohort", "pert_id", "window", "n_conditions", "n_lines", "rev_contrib", "rev_gene_contrib"]
    for ch in pd.read_csv(path, sep="\t", usecols=cols, chunksize=3_000_000, low_memory=False):
        acc.append(ch.groupby(["cohort", "pert_id", "window"], sort=False)
                   .agg(rev=("rev_contrib", "sum"), rev_gene=("rev_gene_contrib", "sum"),
                        n_conditions=("n_conditions", "first"), n_lines=("n_lines", "first")).reset_index())
    S = pd.concat(acc, ignore_index=True)
    # a (cohort, drug, window) can straddle two chunks: re-sum the additive terms, keep the counts
    S = S.groupby(["cohort", "pert_id", "window"], as_index=False).agg(
        rev=("rev", "sum"), rev_gene=("rev_gene", "sum"), n_conditions=("n_conditions", "first"),
        n_lines=("n_lines", "first"))
    return S.rename(columns={"cohort": "unit"})


def pct(x):
    return 100.0 * x.rank(method="average") / x.notna().sum() if x.notna().any() else x


def main():
    L = pd.read_csv(RES / f"layered_score_v2{SFX}.tsv", sep="\t", low_memory=False)
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    L["fda_ingredient"] = L.pert_id.map(fm)
    L["dao_path"] = L.layered_n_axes.notna() & (L.layered_n_axes >= 1)
    L["eligible"] = L.dao_path & (L.layered_T >= T_MIN)
    L["path"] = np.where(L.dao_path, "D/A/O",
                np.where(L.get("dna_acting", False).map(truthy) if "dna_acting" in L else False, "DNA-acting",
                np.where(L.get("essential_unit", False).map(truthy) if "essential_unit" in L else False,
                         "essential-unit", "other")))

    W, reason = [], []
    for unit, path in SOURCES.items():
        w = window_scores(path)
        w["level"] = unit
        W.append(w)
        sp = pd.read_csv(SPLITS[unit], sep="\t", low_memory=False,
                         usecols=lambda c: c in (unit, "pert_id", "window_unknown_reason"))
        reason.append(sp.rename(columns={unit: "unit"}).assign(level=unit))
    W = pd.concat(W, ignore_index=True)
    R = pd.concat(reason, ignore_index=True).drop_duplicates(["unit", "pert_id", "level"])

    base = L[["unit", "level", "pert_id", "cmap_name", "fda_ingredient", "path", "eligible", "layered_T",
              "layered_D", "layered_A", "layered_O", "layered_n_axes", "layered_orient", "layered_genes",
              "layered_pct"]]
    X = base.merge(W, on=["unit", "level", "pert_id"], how="left")
    has = X.window.isin(WINDOWS)
    # rows with no timepoint at all: one `unknown` row per (unit, compound), with the D180 reason
    U = base.merge(X[has][["unit", "level", "pert_id"]].drop_duplicates().assign(_w=1),
                   on=["unit", "level", "pert_id"], how="left")
    U = U[U._w.isna()].drop(columns="_w").assign(window="unknown")
    U = U.merge(R, on=["unit", "level", "pert_id"], how="left")
    U["window_unknown_reason"] = U.window_unknown_reason.fillna("no_window_run")
    X = pd.concat([X[has], U], ignore_index=True)

    # one compound per FDA ingredient per (unit, window): most conditions, then pert_id (not best score)
    X["_n"] = X.n_conditions.fillna(-1)
    X["_k"] = X.fda_ingredient.fillna(X.cmap_name).fillna(X.pert_id)
    X = (X.sort_values(["unit", "window", "_k", "_n", "pert_id"], ascending=[True, True, True, False, True])
          .drop_duplicates(["unit", "level", "window", "_k"]).drop(columns=["_n", "_k"]))

    main_t = X[X.path == "D/A/O"].copy()
    chemo = X[X.path != "D/A/O"].copy()
    for df, elig_col in ((main_t, "eligible"), (chemo, None)):
        df["gene_pct"] = np.nan
        df["gene_rank"] = np.nan
        for (u, w), g in df[df.window.isin(WINDOWS)].groupby(["unit", "window"]):
            gg = g[g[elig_col]] if elig_col else g
            gg = gg[gg.rev_gene.notna()]
            if not len(gg):
                continue
            df.loc[gg.index, "gene_pct"] = pct(gg.rev_gene).round(2)
            df.loc[gg.index, "gene_rank"] = gg.rev_gene.rank(ascending=False, method="min")
        df.sort_values(["level", "unit", "window", "gene_rank"], inplace=True)

    main_t.to_csv(RES / f"gene_ranking{SFX}.tsv", sep="\t", index=False, float_format="%.6g")
    chemo.to_csv(RES / f"gene_ranking_chemo{SFX}.tsv", sep="\t", index=False, float_format="%.6g")

    e = main_t[main_t.eligible]
    run = dict(script="79ar_gene_ranking.py", decision="D211", t_min=T_MIN, score="rev_gene (FORMULAS s3)",
               filter="layered_T >= t_min and layered_n_axes >= 1 (D/A/O path)",
               representative="most conditions in the window, then pert_id",
               rows_main=int(len(main_t)), rows_chemo=int(len(chemo)),
               eligible_unit_ingredients=int(e[["unit", "level", "fda_ingredient"]].drop_duplicates().shape[0]),
               units_with_ranked_rows={w: int(e[(e.window == w) & e.gene_pct.notna()].unit.nunique()) for w in WINDOWS},
               ranked_rows={w: int(((e.window == w) & e.gene_pct.notna()).sum()) for w in WINDOWS},
               unknown_by_reason={k: int(v) for k, v in
                                  e[e.window == "unknown"].window_unknown_reason.value_counts().items()},
               by_level={lv: {w: int(((e.level == lv) & (e.window == w) & e.gene_pct.notna()).sum()) for w in WINDOWS}
                         for lv in ("group_id", "stratum")})
    (RES / f"gene_ranking_run{SFX}.json").write_text(json.dumps(run, indent=1))
    print(json.dumps(run, indent=1))


if __name__ == "__main__":
    main()
