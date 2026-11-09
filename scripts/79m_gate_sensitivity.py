#!/usr/bin/env python
"""Measured drug sensitivity vs the gated v2 rankings: receptor gate, evidence, stand-ins (protocol C28, `celligner`).

C16 tested the v1 scores against PRISM / GDSC before the receptor gate (79w), the FDA filter (79p), the matched /
borrowed split and the stand-ins (79y) existed. This asks whether those pieces carry signal.
Observed selective killing as in 79k: per dataset, the drug's kill is z-scored across all screened cancer lines
(removes general potency), then averaged over the cohort's own in-key lines with the Celligner weights, renormalized
over screened lines; a cohort x drug needs >= MIN_LINES lines and >= MIN_WCOV of its in-key weight screened.
Tests (all on the published FDA rows, diseases; pass rule C28: same direction, p < 0.05 or CI clear of 0, in >= 2 of
the 4 datasets, PRISM secondary primary):
  B1 gate       per cohort: median obs_sel of drugs with >= 1 passing receptor vs drugs whose receptors all fail;
                difference, median over cohorts, permutation of the pass label within cohort
  B2 evidence   per cohort: Spearman(lead score, obs_sel) for matched rows and for borrowed rows; medians and a
                paired bootstrap over cohorts (only cohorts with both)
  B3 stand-ins  rows with a stand-in (79y2): Spearman(obs_sel, stand-in average) vs Spearman(obs_sel, own score),
                pooled over cohorts, paired bootstrap over rows
Outputs (results/lincs/gate_sensitivity/): obs.tsv.gz, b1.tsv, b2.tsv, b3.tsv, summary.tsv, run.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fda_rows import fda_rows, truthy

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SENS = ROOT / "data/processed/sensitivity"
OUT = RES / "gate_sensitivity"
DATASETS = ["prism_sec", "prism_pri", "gdsc1", "gdsc2"]
LEAD = "layered_pct"   # D122 lead (was spec_mech_pct)
MIN_LINES, MIN_WCOV, MIN_DRUGS, MIN_SIDE, N_PERM, N_BOOT, SEED = 2, 0.25, 20, 5, 1000, 1000, 0


def spear(a, b):
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 5:
        return np.nan
    return np.corrcoef(rankdata(a[ok]), rankdata(b[ok]))[0, 1]


def boot_ci(x, rng, n=N_BOOT):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan, np.nan
    b = [np.median(rng.choice(x, len(x))) for _ in range(n)]
    return round(np.percentile(b, 2.5), 4), round(np.percentile(b, 97.5), 4)


def unit_obs(X, Z, w):
    w = w.reindex(X.columns).fillna(0).values
    m = X.notna().values
    den = m @ w
    nl = m.astype(int) @ (w > 0).astype(int)
    sel = np.nan_to_num(Z.values) @ w / np.where(den > 0, den, np.nan)
    ok = (nl >= MIN_LINES) & (den >= MIN_WCOV * w.sum())
    return np.where(ok, sel, np.nan), nl


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    D, _ = fda_rows()
    D["evidence"] = np.where(D.source == "direct", "matched", "borrowed")

    # receptor gate -> per cohort x drug: any receptor passing / all failing (drugs with no receptor row are excluded)
    g = pd.read_csv(RES / "receptor_gate_v2.tsv", sep="\t", low_memory=False)
    g = g[g.level == "disease"]
    gp = g.assign(p=g.receptor_gate.eq("pass")).groupby(["group_id", "pert_id"]).p.any().rename("gate_pass").reset_index()
    D = D.merge(gp, on=["group_id", "pert_id"], how="left")

    # stand-in rows (79y2 table already pairs the drug's own borrowed percentile with the stand-in average)
    sb = pd.read_csv(RES / "standin_vs_borrowed_v2.tsv", sep="\t")
    sb = sb[sb.level == "disease"] if "level" in sb else sb

    W = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    W = {u: d.set_index("ModelID").weight for u, d in W[W.in_key].groupby("group_id")}

    obs_rows, b1, b2, b3, summ = [], [], [], [], []
    for ds in DATASETS:
        X = pd.read_csv(SENS / f"{ds}.tsv.gz", sep="\t", index_col=0)
        Z = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
        per = {}
        for u, d in D.groupby("group_id"):
            if u not in W:
                continue
            sel, nl = unit_obs(X, Z, W[u])
            o = pd.DataFrame({"pert_id": X.index, "obs_sel": sel, "n_lines": nl}).dropna(subset=["obs_sel"])
            d = d.merge(o, on="pert_id")
            if len(d) >= MIN_DRUGS:
                per[u] = d
                obs_rows.append(d.assign(dataset=ds)[["dataset", "group_id", "pert_id", "drug", "evidence",
                                                      "gate_pass", "obs_sel", "n_lines", LEAD, "rev_pct"]])
        # B1 receptor gate
        diffs, nulls = [], []
        for u, d in per.items():
            d = d.dropna(subset=["gate_pass"])
            p, f = d.obs_sel[d.gate_pass.astype(bool)], d.obs_sel[~d.gate_pass.astype(bool)]
            if len(p) < MIN_SIDE or len(f) < MIN_SIDE:
                continue
            diffs.append(p.median() - f.median())
            lab = d.gate_pass.astype(bool).values
            v = d.obs_sel.values
            nulls.append([np.median(v[s]) - np.median(v[~s]) for s in (rng.permutation(lab) for _ in range(N_PERM))])
            b1.append(dict(dataset=ds, cohort=u, n_pass=len(p), n_fail=len(f), median_pass=round(p.median(), 4),
                           median_fail=round(f.median(), 4), diff=round(p.median() - f.median(), 4)))
        if diffs:
            obs, null = np.median(diffs), np.median(np.array(nulls), axis=0)
            lo, hi = boot_ci(diffs, rng)
            summ.append(dict(test="B1_gate_pass_vs_fail", dataset=ds, n_cohorts=len(diffs), estimate=round(obs, 4),
                             ci_lo=lo, ci_hi=hi, p=(1 + (null >= obs).sum()) / (1 + len(null)),
                             detail="median obs_sel difference, drugs with a passing receptor minus the rest"))
        # B2 matched vs borrowed
        pair = []
        for u, d in per.items():
            r = {}
            for e in ("matched", "borrowed"):
                x = d[d.evidence == e]
                r[e] = spear(x[LEAD].values.astype(float), x.obs_sel.values) if len(x) >= MIN_DRUGS else np.nan
            b2.append(dict(dataset=ds, cohort=u, n_matched=int((d.evidence == "matched").sum()),
                           n_borrowed=int((d.evidence == "borrowed").sum()),
                           rho_matched=round(r["matched"], 4) if np.isfinite(r["matched"]) else np.nan,
                           rho_borrowed=round(r["borrowed"], 4) if np.isfinite(r["borrowed"]) else np.nan))
            if np.isfinite(r["matched"]) and np.isfinite(r["borrowed"]):
                pair.append((r["matched"], r["borrowed"]))
        for name, vals in (("B2_rho_matched", [a for a, b in pair]), ("B2_rho_borrowed", [b for a, b in pair]),
                           ("B2_matched_minus_borrowed", [a - b for a, b in pair])):
            if len(vals) >= 3:
                lo, hi = boot_ci(vals, rng)
                summ.append(dict(test=name, dataset=ds, n_cohorts=len(vals), estimate=round(np.median(vals), 4),
                                 ci_lo=lo, ci_hi=hi, p=np.nan,
                                 detail="within-cohort Spearman of the lead score with obs_sel (paired cohorts)"))
        # B3 stand-ins: does the stand-in average track obs_sel better than the drug's own score?
        O = pd.concat(per.values()) if per else pd.DataFrame()
        if len(O):
            m = sb.merge(O[["group_id", "pert_id", "obs_sel"]].rename(columns={"pert_id": "fda_pert"}),
                         on=["group_id", "fda_pert"])
            if len(m) >= MIN_DRUGS:
                r_own, r_avg = spear(m.borrowed_pct.values.astype(float), m.obs_sel.values), \
                               spear(m.average.values.astype(float), m.obs_sel.values)
                r_si = spear(m.standin.values.astype(float), m.obs_sel.values)
                d_b = []
                for _ in range(N_BOOT):
                    i = rng.integers(0, len(m), len(m))
                    d_b.append(spear(m.average.values[i].astype(float), m.obs_sel.values[i]) -
                               spear(m.borrowed_pct.values[i].astype(float), m.obs_sel.values[i]))
                b3.append(dict(dataset=ds, n_rows=len(m), rho_own=round(r_own, 4), rho_standin=round(r_si, 4),
                               rho_average=round(r_avg, 4)))
                summ.append(dict(test="B3_average_minus_own", dataset=ds, n_cohorts=int(m.group_id.nunique()),
                                 estimate=round(r_avg - r_own, 4), ci_lo=round(np.nanpercentile(d_b, 2.5), 4),
                                 ci_hi=round(np.nanpercentile(d_b, 97.5), 4), p=np.nan,
                                 detail=f"Spearman with obs_sel: average {r_avg:.3f} vs own borrowed {r_own:.3f}, {len(m)} rows"))
    pd.concat(obs_rows).to_csv(OUT / "obs.tsv.gz", sep="\t", index=False)
    for n, t in (("b1", b1), ("b2", b2), ("b3", b3), ("summary", summ)):
        pd.DataFrame(t).to_csv(OUT / f"{n}.tsv", sep="\t", index=False)
    json.dump(dict(datasets=DATASETS, lead=LEAD, min_lines=MIN_LINES, min_wcov=MIN_WCOV, min_drugs=MIN_DRUGS,
                   min_side=MIN_SIDE, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED,
                   fda_rows=len(D), gate_known=int(D.gate_pass.notna().sum())), open(OUT / "run.json", "w"), indent=1)
    print(pd.DataFrame(summ).to_string(index=False))
    print(pd.DataFrame(b3).to_string(index=False))


if __name__ == "__main__":
    main()
