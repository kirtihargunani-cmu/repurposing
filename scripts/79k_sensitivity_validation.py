#!/usr/bin/env python
"""Validate the drug scores against measured drug sensitivity (PRISM, GDSC) (runs in `celligner`).

Question: do drugs that score high for a cancer kill that cancer's cell lines more, beyond being generally potent?
Independent of LINCS (viability, not expression). Protocol fixed before looking at results (user 2026-10-12):
Observed sensitivity per unit (41 groups and 31 Step 10 strata) x drug, per dataset (scripts/79j):
  lines    the unit's own in-key lines (strict tissue key), weighted by its Celligner weights (groups: scripts/74,
           strata: 74b), renormalized over the lines screened for that drug; >= MIN_LINES lines and >= MIN_WCOV of the
           unit's in-key weight on screened lines, else missing
  obs_raw  weighted mean kill (1 - AUC, or - log fold change)
  obs_sel  weighted mean of the drug's z-scored kill across all cancer lines screened (selective sensitivity: does it
           hit these lines more than lines in general; removes general toxicity)
  potency  the drug's mean kill over all lines (baseline; uses the sensitivity data itself, so an oracle-like bar)
Scores tested (scripts/79i): rev, rev_prolif, rev_nonprolif, rev_mech, prolif_effect, spec, spec_mech
(strata also spec_sib, spec_sib_mech); plus target_frac (fraction of the drug's targets expressed in the tumors,
scripts/79b) as a non-reversal baseline.
Tests:
  T1 within unit   Spearman(score, obs_sel) over the unit's drugs (primary endpoint); Spearman(score, obs_raw) and
                   the partial Spearman given potency (rank residuals); enrichment = share of the top TOP_PCT% drugs
                   by score that are in the top TOP_PCT% by obs_sel, / TOP_PCT%. Per dataset: median over units with
                   >= MIN_DRUGS drugs; p = permutation of obs within each unit (N_PERM), one-sided for median > 0.
                   Units split into two fixed halves (seeded, strata follow their parent) to check consistency.
                   Low-MR groups (and strata of them) excluded; direct-only rerun (drugs measured in the tissue).
  T2 across units  per drug scored in >= MIN_UNITS of the 41 groups: Spearman across groups of the score percentile vs
                   obs_sel (does the drug score highest where it kills selectively); median over drugs, permutation of
                   groups within drug.
  T3 siblings      strata: Spearman over drugs of spec_sib (and spec_sib_mech, rev) vs obs_sel minus the mean of the
                   siblings' obs_sel; with the overlap of sibling weights on screened lines (1 = same lines, no
                   contrast possible)
  Controls         obs_sel percentile of the known-drug controls in their unit
Outputs (results/lincs/sensitivity_validation/): unit_drug_obs.tsv.gz, t1_per_unit.tsv, t1_summary.tsv,
         t2_across_units.tsv, t3_siblings.tsv, controls.tsv, summary.tsv (headline numbers)
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SENS = ROOT / "data/processed/sensitivity"
VAL = RES / "sensitivity_validation"
DATASETS = ["prism_sec", "prism_pri", "gdsc1", "gdsc2"]
SCORES = ["rev", "rev_prolif", "rev_nonprolif", "rev_mech", "prolif_effect", "spec", "spec_mech", "target_frac"]
SSCORES = SCORES + ["spec_sib", "spec_sib_mech"]
MIN_LINES, MIN_WCOV, MIN_DRUGS, MIN_UNITS, TOP_PCT, N_PERM, SEED = 2, 0.25, 30, 10, 10, 1000, 0
CONTROLS = [("BRCA-LumA", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("BRCA-LumB", ["fulvestrant", "tamoxifen", "raloxifene", "palbociclib"]),
            ("BRCA-Her2", ["lapatinib", "neratinib", "afatinib"]),
            ("TCGA-PRAD", ["bicalutamide", "enzalutamide", "flutamide", "nilutamide"]),
            ("TCGA-COAD", ["fluorouracil", "irinotecan", "oxaliplatin"]),
            ("TCGA-LUAD", ["erlotinib", "gefitinib", "afatinib", "osimertinib"]),
            ("TCGA-SKCM", ["vemurafenib", "dabrafenib", "trametinib", "cobimetinib"]),
            ("BRCA-Basal", ["fulvestrant", "tamoxifen", "raloxifene"])]


def srank(x):
    return rankdata(x) / len(x)


def spear(a, b):
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 5:
        return np.nan
    return np.corrcoef(rankdata(a[ok]), rankdata(b[ok]))[0, 1]


def partial_spear(a, b, c):
    """Spearman of a and b given c: Pearson of rank residuals."""
    ok = ~(np.isnan(a) | np.isnan(b) | np.isnan(c))
    if ok.sum() < 5:
        return np.nan
    ra, rb, rc = (rankdata(v[ok]) for v in (a, b, c))
    res = lambda y: y - np.polyval(np.polyfit(rc, y, 1), rc)
    return np.corrcoef(res(ra), res(rb))[0, 1]


def unit_obs(X, Z, w):
    """Weighted mean of X and Z over screened lines per drug; coverage checks."""
    w = w.reindex(X.columns).fillna(0).values
    m = X.notna().values
    den = m @ w
    nl = m.astype(int) @ (w > 0).astype(int)
    raw = np.nan_to_num(X.values) @ w / np.where(den > 0, den, np.nan)
    sel = np.nan_to_num(Z.values) @ w / np.where(den > 0, den, np.nan)
    ok = (nl >= MIN_LINES) & (den >= MIN_WCOV * w.sum())
    return np.where(ok, raw, np.nan), np.where(ok, sel, np.nan), nl, den / max(w.sum(), 1e-12)


def main():
    rng = np.random.default_rng(SEED)
    VAL.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(RES / "drug_reversal_split.tsv", sep="\t", low_memory=False)
    s = pd.read_csv(RES / "stratum_reversal_split.tsv", sep="\t", low_memory=False)
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_targets", "targets_expressed"])
    sd = pd.read_csv(RES / "stratum_drugs.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_targets", "targets_expressed"])
    nexp = lambda d: d.targets_expressed.fillna("").map(lambda x: len([t for t in str(x).split("|") if t]))
    for d in (dd, sd):
        d["target_frac"] = np.where(d.n_targets > 0, nexp(d) / d.n_targets.where(d.n_targets > 0), np.nan)
    g = g.merge(dd.drop_duplicates(["group_id", "pert_id"])[["group_id", "pert_id", "target_frac"]], how="left")
    s = s.merge(sd.drop_duplicates(["group_id", "pert_id"]).rename(columns={"group_id": "stratum"})
                [["stratum", "pert_id", "target_frac"]], how="left")
    g["unit"], g["parent"], g["kind"] = g.group_id, g.group_id, "group"
    s["unit"], s["parent"], s["kind"] = s.stratum, s.parent_group, "stratum"
    low = set(g.group_id[g.low_mrs])
    sc = pd.concat([g, s], ignore_index=True)
    sc = sc[~sc.parent.isin(low)]

    gw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    sw = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long.tsv", sep="\t")
    W = {u: d.set_index("ModelID").weight for u, d in gw[gw.in_key].groupby("group_id")}
    W.update({u: d.set_index("ModelID").weight for u, d in sw[sw.in_key].groupby("stratum")})
    parents = sorted(sc.parent.unique())
    half = dict(zip(parents, rng.permutation(np.arange(len(parents)) % 2)))

    obs_rows, t1, t1sum, t2, t3, ctl = [], [], [], [], [], []
    for ds in DATASETS:
        X = pd.read_csv(SENS / f"{ds}.tsv.gz", sep="\t", index_col=0)
        Z = X.sub(X.mean(1), axis=0).div(X.std(1), axis=0)
        pot = X.mean(1)
        per = []
        for u, d in sc.groupby("unit"):
            if u not in W:
                continue
            raw, sel, nl, cov = unit_obs(X, Z, W[u])
            o = pd.DataFrame({"pert_id": X.index, "obs_raw": raw, "obs_sel": sel, "n_lines": nl, "wcov": cov,
                              "potency": pot.values}).dropna(subset=["obs_sel"])
            d = d.merge(o, on="pert_id")
            if len(d) < MIN_DRUGS:
                continue
            d["dataset"] = ds
            obs_rows.append(d[["dataset", "unit", "kind", "parent", "pert_id", "cmap_name", "source", "obs_raw",
                               "obs_sel", "potency", "n_lines", "wcov"] + [c for c in SSCORES if c in d]])
            row = dict(dataset=ds, unit=u, kind=d.kind.iloc[0], parent=d.parent.iloc[0], half=half[d.parent.iloc[0]],
                       n_drugs=len(d), n_direct=int((d.source == "direct").sum()),
                       n_lines_median=float(d.n_lines.median()), rho_potency_raw=spear(d.potency.values,
                                                                                        d.obs_raw.values))
            hit = d.obs_sel.values >= np.nanpercentile(d.obs_sel.values, 100 - TOP_PCT)
            for c in SSCORES:
                if c not in d or d[c].isna().all():
                    continue
                v = d[c].values.astype(float)
                row[f"rho_sel_{c}"] = spear(v, d.obs_sel.values)
                row[f"rho_raw_{c}"] = spear(v, d.obs_raw.values)
                row[f"prho_raw_{c}"] = partial_spear(v, d.obs_raw.values, d.potency.values)
                top = v >= np.nanpercentile(v, 100 - TOP_PCT)
                row[f"enr_{c}"] = (hit & top).sum() / max(top.sum(), 1) / (TOP_PCT / 100)
                dr = d.source.values == "direct"
                row[f"rho_sel_direct_{c}"] = spear(v[dr], d.obs_sel.values[dr]) if dr.sum() >= MIN_DRUGS else np.nan
            per.append((row, d))
            t1.append(row)
        # permutation null for the median within-unit Spearman with obs_sel (groups and strata separately)
        for kind in ["group", "stratum"]:
            units = [(r, d) for r, d in per if r["kind"] == kind]
            if not units:
                continue
            for c in SSCORES:
                data = [(srank(x[c].values.astype(float)), srank(x.obs_sel.values)) for x in
                        (d.dropna(subset=[c]) for r, d in units if c in d) if len(x) >= MIN_DRUGS]
                if not data:
                    continue
                obs = np.median([np.corrcoef(a, b)[0, 1] for a, b in data])
                null = np.median([[np.corrcoef(a, rng.permutation(b))[0, 1] for a, b in data] for _ in range(N_PERM)],
                                 axis=1)
                rr = pd.DataFrame([r for r, d in units])
                rec = dict(dataset=ds, kind=kind, score=c, n_units=len(data), median_rho_sel=round(obs, 4),
                           p_perm=(1 + (null >= obs).sum()) / (1 + N_PERM),
                           frac_units_rho_pos=round((rr[f"rho_sel_{c}"] > 0).mean(), 3),
                           median_rho_sel_half0=round(rr.loc[rr.half == 0, f"rho_sel_{c}"].median(), 4),
                           median_rho_sel_half1=round(rr.loc[rr.half == 1, f"rho_sel_{c}"].median(), 4),
                           median_rho_sel_direct=round(rr[f"rho_sel_direct_{c}"].median(), 4),
                           median_rho_raw=round(rr[f"rho_raw_{c}"].median(), 4),
                           median_prho_raw_given_potency=round(rr[f"prho_raw_{c}"].median(), 4),
                           median_rho_potency_raw=round(rr.rho_potency_raw.median(), 4),
                           median_enrichment_top10=round(rr[f"enr_{c}"].median(), 3))
                t1sum.append(rec)

        # T2 across groups per drug
        gd = pd.concat([d for r, d in per if r["kind"] == "group"]) if per else pd.DataFrame()
        for c in SCORES:
            if gd.empty or c not in gd:
                continue
            gd[c + "_u"] = gd.groupby("unit")[c].rank(pct=True)
            rhos, nulls = [], []
            for pid, x in gd.groupby("pert_id"):
                x = x.dropna(subset=[c + "_u", "obs_sel"])
                if len(x) < MIN_UNITS:
                    continue
                a, b = srank(x[c + "_u"].values), srank(x.obs_sel.values)
                rhos.append(np.corrcoef(a, b)[0, 1])
                nulls.append([np.corrcoef(a, rng.permutation(b))[0, 1] for _ in range(200)])
            if rhos:
                obs, null = np.nanmedian(rhos), np.nanmedian(np.array(nulls), axis=0)
                t2.append(dict(dataset=ds, score=c, n_drugs=len(rhos), median_rho=round(obs, 4),
                               frac_rho_pos=round(np.mean(np.array(rhos) > 0), 3),
                               p_perm=(1 + (null >= obs).sum()) / (1 + len(null))))

        # T3 sibling contrasts
        sd_ = {r["unit"]: d.set_index("pert_id") for r, d in per if r["kind"] == "stratum"}
        for u, d in sd_.items():
            sibs = [t for t in str(d.siblings.iloc[0]).split(";") if t in sd_]
            if not sibs:
                continue
            ob = pd.concat([sd_[t].obs_sel.groupby(level=0).mean() for t in sibs], axis=1).mean(1)
            x = d[~d.index.duplicated()].join(ob.rename("sib_obs"), how="inner")
            x["d_obs"] = x.obs_sel - x.sib_obs
            wa = W[u]
            ov = np.mean([1 - 0.5 * (wa / wa.sum()).sub(W[t] / W[t].sum(), fill_value=0).abs().sum() for t in sibs])
            rec = dict(dataset=ds, stratum=u, parent=d.parent.iloc[0], siblings=";".join(sibs), n_drugs=len(x),
                       weight_overlap_siblings=round(ov, 3), sd_d_obs=round(x.d_obs.std(), 3))
            for c in ["spec_sib", "spec_sib_mech", "rev", "rev_mech", "rev_prolif"]:
                rec[f"rho_{c}"] = round(spear(x[c].values.astype(float), x.d_obs.values), 4)
            t3.append(rec)

        # controls
        for u, names in CONTROLS:
            for r, d in per:
                if r["unit"] != u:
                    continue
                pc = d.obs_sel.rank(pct=True) * 100
                for nm in names:
                    m = d.cmap_name == nm
                    ctl.append(dict(dataset=ds, unit=u, drug=nm, in_data=bool(m.any()),
                                    obs_sel_pct=round(pc[m].max(), 1) if m.any() else np.nan,
                                    obs_sel=round(d.obs_sel[m].max(), 3) if m.any() else np.nan,
                                    rev_pct=d.rev_pct[m].max() if m.any() else np.nan,
                                    rev_mech_pct=d.rev_mech_pct[m].max() if m.any() else np.nan))
        print(f"{ds}: units {len(per)}", flush=True)

    pd.concat(obs_rows).round(4).to_csv(VAL / "unit_drug_obs.tsv.gz", sep="\t", index=False)
    pd.DataFrame(t1).round(4).to_csv(VAL / "t1_per_unit.tsv", sep="\t", index=False)
    t1s = pd.DataFrame(t1sum)
    t1s.to_csv(VAL / "t1_summary.tsv", sep="\t", index=False)
    pd.DataFrame(t2).to_csv(VAL / "t2_across_units.tsv", sep="\t", index=False)
    pd.DataFrame(t3).to_csv(VAL / "t3_siblings.tsv", sep="\t", index=False)
    pd.DataFrame(ctl).to_csv(VAL / "controls.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    print(t1s.to_string(index=False))
    print(pd.DataFrame(t2).to_string(index=False))


if __name__ == "__main__":
    main()
