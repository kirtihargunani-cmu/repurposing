"""D128 A3: aggregation from per-line, per-unit scores to one row per (cohort, drug) (plan 3.4; D128-A3).

Two steps, in the order the decision row fixes:

  1. source lines -> target unit   weighted mean over the drug's source lines, with the 79e / P14 direct-or-borrowed
                                   population and Celligner weights carried unchanged from A0. Weights are
                                   renormalised over the lines that actually scored; abstained lines are dropped,
                                   never counted as 0. If every scored line carries weight 0 the weights fall back to
                                   equal, exactly as `79e_drug_reversal.py` does (`ww if ww.sum() > 0 else ones`).
  2. target units -> drug          unweighted mean over scoreable units (P13). The fraction of annotated units that
                                   stayed unscorable travels with the row instead of being hidden.

Nothing here reads a label, a baseline column or a production score.
"""
import numpy as np
import pandas as pd


def _w(weights):
    """79e's weighting rule: the given weights, or equal weights if they sum to zero."""
    w = np.asarray(weights, dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    return w if w.sum() > 0 else np.ones(len(w))


def weighted_stats(values, weights):
    """Weighted mean and the spread around it. n = 1 gives sd / range NA, not 0."""
    v = np.asarray(values, dtype=float)
    ok = np.isfinite(v)
    v, w = v[ok], _w(weights)[ok]
    if not len(v):
        return dict(n=0, mean=np.nan, sd=np.nan, min=np.nan, max=np.nan, range=np.nan, unweighted_mean=np.nan)
    w = w if w.sum() > 0 else np.ones(len(v))
    mu = float(np.average(v, weights=w))
    sd = float(np.sqrt(np.average((v - mu) ** 2, weights=w))) if len(v) > 1 else np.nan
    return dict(n=len(v), mean=mu, sd=sd, min=float(v.min()), max=float(v.max()),
                range=float(v.max() - v.min()) if len(v) > 1 else np.nan, unweighted_mean=float(v.mean()))


def lines_to_unit(S):
    """Step 1: one row per (cohort, level, pert_id, unit_id), aggregating over source lines."""
    out = []
    keys = ["cohort", "level", "pert_id", "cmap_name", "unit_id"]
    for k, d in S.groupby(keys, sort=False):
        m = weighted_stats(np.where(d.M_status == "scored", d.M, np.nan), d.line_weight)
        r = weighted_stats(np.where(d.R_status == "scored", d.R, np.nan), d.line_weight)
        sc_m = d[d.M_status == "scored"]
        row = dict(zip(keys, k))
        row.update(
            unit_genes=d.unit_genes.iloc[0], unit_kind=d.unit_kind.iloc[0], action=d.action.iloc[0],
            action_flag=d.action_flag.iloc[0], source=d.source.iloc[0], contrast=d.contrast.iloc[0],
            n_lines_total=len(d), n_lines_M=m["n"], n_lines_R=r["n"],
            M=m["mean"], M_sd=m["sd"], M_min=m["min"], M_max=m["max"], M_range=m["range"],
            M_unweighted=m["unweighted_mean"],
            R=r["mean"], R_sd=r["sd"], R_min=r["min"], R_max=r["max"], R_range=r["range"],
            R_unweighted=r["unweighted_mean"],
            M_n_genes=float(sc_m.M_n_genes.mean()) if len(sc_m) else np.nan,
            M_n_landmark=float(sc_m.M_n_landmark.mean()) if len(sc_m) else np.nan,
            M_n_inferred=float(sc_m.M_n_inferred.mean()) if len(sc_m) else np.nan,
            M_n_excluded_unresolved=float(sc_m.M_n_excluded_unresolved.mean()) if len(sc_m) else np.nan,
            R_n_genes=float(d.loc[d.R_status == "scored", "R_n_genes"].mean()) if (d.R_status == "scored").any() else np.nan,
            hub_fraction_structural=d.hub_fraction_structural.iloc[0], hub_tf=d.hub_tf.iloc[0],
            M_abstain_reason="" if m["n"] else d.M_abstain_reason.mode().iloc[0] if len(d) else "",
            R_abstain_reason="" if r["n"] else d.R_abstain_reason.mode().iloc[0] if len(d) else "",
            lines="|".join(f"{a}:{b:.4f}" for a, b in zip(sc_m.ModelID, sc_m.M)) if len(sc_m) else "")
        out.append(row)
    return pd.DataFrame(out)


def units_to_drug(U, n_units_annotated):
    """Step 2: one row per (cohort, level, pert_id), unweighted mean over scoreable target units (P13)."""
    out = []
    keys = ["cohort", "level", "pert_id", "cmap_name"]
    for k, d in U.groupby(keys, sort=False):
        okM, okR = d[d.M.notna()], d[d.R.notna()]
        row = dict(zip(keys, k))
        n_ann = n_units_annotated.get(d.cmap_name.iloc[0], len(d))
        row.update(
            source=d.source.iloc[0], contrast=d.contrast.iloc[0],
            n_units_annotated=n_ann, n_units_in_run=len(d), n_units_M=len(okM), n_units_R=len(okR),
            unit_scored_fraction=round(len(okM) / n_ann, 4) if n_ann else np.nan,
            M=float(okM.M.mean()) if len(okM) else np.nan,
            M_unit_sd=float(okM.M.std(ddof=0)) if len(okM) > 1 else np.nan,
            M_unit_min=float(okM.M.min()) if len(okM) else np.nan,
            M_unit_max=float(okM.M.max()) if len(okM) else np.nan,
            R=float(okR.R.mean()) if len(okR) else np.nan,
            R_unit_sd=float(okR.R.std(ddof=0)) if len(okR) > 1 else np.nan,
            M_line_sd_max=float(okM.M_sd.max()) if len(okM) else np.nan,
            n_lines_M=int(okM.n_lines_M.max()) if len(okM) else 0,
            n_lines_R=int(okR.n_lines_R.max()) if len(okR) else 0,
            M_n_genes=float(okM.M_n_genes.mean()) if len(okM) else np.nan,
            M_n_landmark=float(okM.M_n_landmark.mean()) if len(okM) else np.nan,
            R_n_genes=float(okR.R_n_genes.mean()) if len(okR) else np.nan,
            hub_fraction_structural=float(okM.hub_fraction_structural.max()) if len(okM) else np.nan,
            hub_tf=okM.hub_tf.iloc[0] if len(okM) else "",
            units="|".join(f"{a}={b:.4f}" for a, b in zip(okM.unit_genes, okM.M)) if len(okM) else "",
            unit_kinds="|".join(sorted(set(d.unit_kind.dropna()))),
            M_abstain_reason="" if len(okM) else ";".join(sorted(set(d.M_abstain_reason) - {""}))[:120],
            R_abstain_reason="" if len(okR) else ";".join(sorted(set(d.R_abstain_reason) - {""}))[:120])
        out.append(row)
    return pd.DataFrame(out)


def add_percentiles(D, cols=("M", "R"), within=("cohort", "level"), mask=None):
    """Percentile WITHIN the declared comparison population, never across cohorts (plan A3 action 5)."""
    D = D.copy()
    keep = pd.Series(True, index=D.index) if mask is None else mask
    for c in cols:
        D[f"{c}_pct"] = np.nan
        sub = D[keep & D[c].notna()]
        if len(sub):
            D.loc[sub.index, f"{c}_pct"] = (sub.groupby(list(within))[c].rank(pct=True) * 100).round(3)
        D[f"n_{c}_ranked"] = D.groupby(list(within))[f"{c}_pct"].transform("count")
    return D
