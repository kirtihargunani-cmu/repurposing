"""D128-A3 landmark sensitivity (amendment 2026-11-16, user; post-freeze, reported only).

P4 freezes the primary gene space as all 10,112 BING genes. Most of those are not measured: L1000 measures 978
landmarks and the rest are inferred from them by a fixed model. The question this module answers is not whether the
inferred genes improve coverage (A0 already measured that) but whether they change the **scores**. Every unit is
scored twice through the identical functions - `all_bing` and `landmark_only` - and the three criteria the decision
row fixed **before** the result was seen are evaluated here:

    (i)  Spearman rho(M_all, M_landmark) >= 0.90
    (ii) sign-flip rate <= 5% of units scored in both arms
    (iii) median |delta percentile| <= 5 points

All three are reported whichever way they fall. `all_bing` remains the primary arm unless the user decides
otherwise; this module never changes a score.
"""
import numpy as np
import pandas as pd

def spearman(x, y):
    """Spearman rho with average ranks for ties. scipy is not in this environment; pandas ranks are enough."""
    x, y = pd.Series(np.asarray(x, float)), pd.Series(np.asarray(y, float))
    ok = x.notna() & y.notna()
    if int(ok.sum()) < 3:
        return np.nan
    rx, ry = x[ok].rank(), y[ok].rank()
    if rx.std(ddof=0) == 0 or ry.std(ddof=0) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


KEY_UNIT = ["cohort", "level", "pert_id", "cmap_name", "unit_id"]
KEY_DRUG = ["cohort", "level", "pert_id", "cmap_name"]


def compare(un_all, un_lm, dr_all, dr_lm):
    """Unit-level table with both arms side by side, plus the drug-level percentile of each arm."""
    a = un_all[KEY_UNIT + ["unit_genes", "unit_kind", "M", "R", "M_n_genes", "M_n_landmark", "R_n_genes",
                           "M_abstain_reason", "R_abstain_reason", "n_lines_M"]]
    b = un_lm[KEY_UNIT + ["M", "R", "M_n_genes", "R_n_genes", "M_abstain_reason", "R_abstain_reason"]]
    T = a.merge(b, on=KEY_UNIT, how="outer", suffixes=("_all", "_landmark"), validate="one_to_one")
    T["landmark_fraction"] = (T.M_n_landmark / T.M_n_genes_all).round(4)
    T["M_delta"] = T.M_landmark - T.M_all
    T["R_delta"] = T.R_landmark - T.R_all
    T["M_sign_flip"] = (np.sign(T.M_all) * np.sign(T.M_landmark) < 0)
    T["R_sign_flip"] = (np.sign(T.R_all) * np.sign(T.R_landmark) < 0)
    T["floor_lost_in_landmark"] = T.M_all.notna() & T.M_landmark.isna()
    T["floor_gained_in_landmark"] = T.M_all.isna() & T.M_landmark.notna()

    pa = dr_all[KEY_DRUG + ["M_pct", "R_pct"]].rename(columns={"M_pct": "M_pct_all", "R_pct": "R_pct_all"})
    pb = dr_lm[KEY_DRUG + ["M_pct", "R_pct"]].rename(columns={"M_pct": "M_pct_landmark", "R_pct": "R_pct_landmark"})
    P = pa.merge(pb, on=KEY_DRUG, how="outer", validate="one_to_one")
    P["M_pct_delta"] = P.M_pct_landmark - P.M_pct_all
    P["R_pct_delta"] = P.R_pct_landmark - P.R_pct_all
    T = T.merge(P, on=KEY_DRUG, how="left")
    return T, P


def _criteria(T, P, score, subset_name, mask=None):
    t = T if mask is None else T[mask]
    p = P if mask is None else P[P.set_index(KEY_DRUG).index.isin(t.set_index(KEY_DRUG).index)]
    both = t[t[f"{score}_all"].notna() & t[f"{score}_landmark"].notna()]
    rho = spearman(both[f"{score}_all"], both[f"{score}_landmark"]) if len(both) > 2 else np.nan
    flips = float(both[f"{score}_sign_flip"].mean() * 100) if len(both) else np.nan
    dp = p[f"{score}_pct_delta"].abs().dropna()
    med = float(dp.median()) if len(dp) else np.nan
    met = bool(np.isfinite(rho) and rho >= 0.90 and np.isfinite(flips) and flips <= 5.0
               and np.isfinite(med) and med <= 5.0)
    return dict(
        subset=subset_name, score=score, n_units_both_arms=len(both),
        spearman_rho=round(float(rho), 4) if np.isfinite(rho) else np.nan,
        crit_i_rho_ge_090=bool(np.isfinite(rho) and rho >= 0.90),
        sign_flip_pct=round(flips, 3) if np.isfinite(flips) else np.nan,
        crit_ii_flip_le_5pct=bool(np.isfinite(flips) and flips <= 5.0),
        n_drug_rows_pct=int(len(dp)), median_abs_delta_pct=round(med, 3) if np.isfinite(med) else np.nan,
        crit_iii_median_delta_le_5=bool(np.isfinite(med) and med <= 5.0),
        median_delta_raw=round(float(both[f"{score}_delta"].median()), 5) if len(both) else np.nan,
        interchangeable=met,
        n_floor_lost_in_landmark=int(t.floor_lost_in_landmark.sum()),
        n_floor_gained_in_landmark=int(t.floor_gained_in_landmark.sum()))


def summary(T, P, low_landmark_quantile=0.25):
    """The three criteria for M and R, over the whole population and over the low-landmark-fraction subset."""
    rows = [_criteria(T, P, s, "all units") for s in ("M", "R")]
    lf = T.landmark_fraction
    if lf.notna().any():
        cut = float(lf.quantile(low_landmark_quantile))
        m = lf <= cut
        for s in ("M", "R"):
            r = _criteria(T, P, s, f"landmark fraction <= {cut:.3f} (lowest {int(low_landmark_quantile*100)}%)", m)
            rows.append(r)
        m2 = lf < 0.15                     # the A0 median landmark share
        for s in ("M", "R"):
            rows.append(_criteria(T, P, s, "landmark fraction < 0.150 (below the A0 median)", m2))
    return pd.DataFrame(rows)
