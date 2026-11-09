#!/usr/bin/env python
"""Do the cohort's cell lines depend on the genes our drugs target? (protocol C29, runs in `celligner`).

DepMap 24Q4 CRISPR Chronos gene effect; each gene z-scored across all screened cancer lines (common essentials carry
no signal), dep = -z so higher = this line needs the gene more than the typical cancer line. Per cohort x gene: the
weighted mean of dep over the cohort's in-key lines (Celligner weights, as 79e / 79k), needing >= MIN_LINES screened
lines and >= MIN_WCOV of the in-key weight. Per cohort x drug: dep_max over the drug's tumor-expressed targets
(dep_mean reported alongside).
Tests (C29): yardstick of known driver-target pairs (runs first); D1 Spearman(lead score, dep_max) per cohort;
D2 receptor-gate pass vs fail on the receptor's own dep; D3 known-use labelled drugs (C28) vs unlabelled.
Outputs (results/lincs/dependency/): gene_dep.tsv.gz, drug_dep.tsv.gz, yardstick.tsv, d1.tsv, d2.tsv, d3.tsv,
         summary.tsv, run.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fda_rows import fda_rows

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "dependency"
LEAD = "layered_pct"   # D122 lead (was spec_mech_pct)
MIN_LINES, MIN_WCOV, MIN_DRUGS, MIN_SIDE, N_PERM, N_BOOT, SEED = 2, 0.25, 20, 5, 1000, 1000, 0
YARDSTICK = [("TCGA-LUAD", "EGFR"), ("BRCA-Her2", "ERBB2"), ("BRCA-LumA", "ESR1"), ("BRCA-LumB", "ESR1"),
             ("TCGA-PRAD", "AR"), ("TCGA-SKCM", "BRAF"), ("TCGA-STAD", "KIT"), ("TCGA-LIHC", "CTNNB1")]


def spear(a, b):
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() < 5:
        return np.nan
    return np.corrcoef(rankdata(a[ok]), rankdata(b[ok]))[0, 1]


def boot_ci(x, rng):
    x = np.asarray(x, float)[np.isfinite(np.asarray(x, float))]
    if len(x) < 3:
        return np.nan, np.nan
    b = [np.median(rng.choice(x, len(x))) for _ in range(N_BOOT)]
    return round(np.percentile(b, 2.5), 4), round(np.percentile(b, 97.5), 4)


def cohort_gene_dep(Dz, W):
    """{cohort: Series(gene -> weighted mean dep)} over in-key screened lines."""
    out = {}
    for u, w in W.items():
        w = w.reindex(Dz.index).fillna(0).values
        if w.sum() <= 0:
            continue
        m = Dz.notna().values
        den = m.T @ w
        nl = m.astype(int).T @ (w > 0).astype(int)
        v = np.nan_to_num(Dz.values).T @ w / np.where(den > 0, den, np.nan)
        ok = (nl >= MIN_LINES) & (den >= MIN_WCOV * w.sum())
        out[u] = pd.Series(np.where(ok, v, np.nan), index=Dz.columns)
    return out


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    G = pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRGeneEffect.csv", index_col=0)
    G.columns = [c.split(" (")[0] for c in G.columns]
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv")
    cancer = set(models.ModelID[models.OncotreePrimaryDisease != "Non-Cancerous"])
    G = G.loc[G.index.isin(cancer)]
    Dz = -G.sub(G.mean(0), axis=1).div(G.std(0).replace(0, np.nan), axis=1)  # dep: higher = more dependent than typical

    W = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    W = {u: d.set_index("ModelID").weight for u, d in W[W.in_key].groupby("group_id")}
    gd = cohort_gene_dep(Dz, W)
    GD = pd.DataFrame(gd)
    GD.to_csv(OUT / "gene_dep.tsv.gz", sep="\t")

    yard = []
    for coh, gene in YARDSTICK:
        s = gd.get(coh)
        v = s.get(gene, np.nan) if s is not None else np.nan
        pct = float((s.dropna() < v).mean() * 100) if s is not None and np.isfinite(v) else np.nan
        yard.append(dict(cohort=coh, gene=gene, dep=round(v, 3) if np.isfinite(v) else np.nan,
                         pct_in_cohort=round(pct, 1) if np.isfinite(pct) else np.nan,
                         n_genes=int(s.notna().sum()) if s is not None else 0))
    Y = pd.DataFrame(yard)
    Y.to_csv(OUT / "yardstick.tsv", sep="\t", index=False)
    n_ok = int((Y.pct_in_cohort >= 80).sum())
    print(Y.to_string(index=False))
    print(f"yardstick: {n_ok} of {len(Y)} known pairs above the 80th percentile (C29 needs >= 4)")

    D, _ = fda_rows()
    D["evidence"] = np.where(D.source == "direct", "matched", "borrowed")
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "targets_expressed"]).drop_duplicates(["group_id", "pert_id"])
    D = D.merge(dd, on=["group_id", "pert_id"], how="left")
    dmax, dmean, ntg = [], [], []
    for coh, tx in zip(D.group_id, D.targets_expressed):
        s = gd.get(coh)
        vs = [s[t] for t in str(tx).split("|") if t and s is not None and t in s.index and np.isfinite(s[t])] if s is not None else []
        dmax.append(max(vs) if vs else np.nan)
        dmean.append(float(np.mean(vs)) if vs else np.nan)
        ntg.append(len(vs))
    D["dep_max"], D["dep_mean"], D["n_targets_scored"] = dmax, dmean, ntg
    D[["group_id", "drug", "pert_id", "evidence", LEAD, "rev_pct", "dep_max", "dep_mean",
       "n_targets_scored"]].to_csv(OUT / "drug_dep.tsv.gz", sep="\t", index=False)

    summ, d1, d2, d3 = [], [], [], []
    # D1 lead score vs dep_max
    for name, sub in (("all", D), ("matched", D[D.evidence == "matched"]), ("borrowed", D[D.evidence == "borrowed"])):
        rhos, nulls = [], []
        for coh, d in sub.groupby("group_id"):
            d = d.dropna(subset=["dep_max", LEAD])
            if len(d) < MIN_DRUGS:
                continue
            r = spear(d[LEAD].values.astype(float), d.dep_max.values)
            rhos.append(r)
            a, b = d[LEAD].values.astype(float), d.dep_max.values
            nulls.append([spear(rng.permutation(a), b) for _ in range(N_PERM)])
            if name == "all":
                d1.append(dict(cohort=coh, n_drugs=len(d), rho=round(r, 4),
                               median_dep_top10=round(d.dep_max[d[LEAD] >= np.nanpercentile(d[LEAD], 90)].median(), 4),
                               median_dep_rest=round(d.dep_max[d[LEAD] < np.nanpercentile(d[LEAD], 90)].median(), 4)))
        if len(rhos) >= 3:
            obs, null = np.median(rhos), np.median(np.array(nulls), axis=0)
            lo, hi = boot_ci(rhos, rng)
            summ.append(dict(test="D1_lead_vs_dependency", subset=name, n_cohorts=len(rhos), estimate=round(obs, 4),
                             ci_lo=lo, ci_hi=hi, p=(1 + (null >= obs).sum()) / (1 + len(null)),
                             detail="within-cohort Spearman of the lead score with the drug's most depended-on target"))
    # D2 receptor gate pass vs fail, on the receptor gene's own dependency
    g = pd.read_csv(RES / "receptor_gate_v2.tsv", sep="\t", low_memory=False)
    g = g[g.level == "disease"].copy()
    g["dep"] = [gd[c][t] if c in gd and t in gd[c].index else np.nan for c, t in zip(g.group_id, g.target)]
    g["pass"] = g.receptor_gate.eq("pass")
    diffs, nulls = [], []
    for coh, d in g.dropna(subset=["dep"]).groupby("group_id"):
        p, f = d.dep[d["pass"]], d.dep[~d["pass"]]
        if len(p) < MIN_SIDE or len(f) < MIN_SIDE:
            continue
        diffs.append(p.median() - f.median())
        v, lab = d.dep.values, d["pass"].values
        nulls.append([np.median(v[s]) - np.median(v[~s]) for s in (rng.permutation(lab) for _ in range(N_PERM))])
        d2.append(dict(cohort=coh, n_pass=len(p), n_fail=len(f), median_dep_pass=round(p.median(), 4),
                       median_dep_fail=round(f.median(), 4), diff=round(p.median() - f.median(), 4)))
    if len(diffs) >= 3:
        obs, null = np.median(diffs), np.median(np.array(nulls), axis=0)
        lo, hi = boot_ci(diffs, rng)
        summ.append(dict(test="D2_gate_pass_vs_fail", subset="all", n_cohorts=len(diffs), estimate=round(obs, 4),
                         ci_lo=lo, ci_hi=hi, p=(1 + (null >= obs).sum()) / (1 + len(null)),
                         detail="median dependency of the receptor gene, gate-passing pairs minus failing pairs"))
    # D3 known-use labelled drugs vs unlabelled (C28 labels)
    lab = RES / "known_use/labels.tsv"
    if lab.exists():
        L = pd.read_csv(lab, sep="\t")[["cohort", "fda_ingredient", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl"]]   # D118: drug level
        M = D.merge(L, left_on=["group_id", "fda_ingredient"], right_on=["cohort", "fda_ingredient"], how="left")
        for src in ("hub", "chembl"):
            diffs, nulls = [], []
            for coh, d in M[~M[f"excl_{src}"].fillna(False).astype(bool)].dropna(subset=["dep_max"]).groupby("group_id"):
                k = d[f"lab_{src}"].fillna(False).astype(bool)
                if k.sum() < 3 or (~k).sum() < MIN_DRUGS:
                    continue
                diffs.append(d.dep_max[k].median() - d.dep_max[~k].median())
                v = d.dep_max.values
                nulls.append([np.median(v[s]) - np.median(v[~s]) for s in (rng.permutation(k.values) for _ in range(N_PERM))])
                d3.append(dict(source=src, cohort=coh, n_labelled=int(k.sum()), diff=round(diffs[-1], 4)))
            if len(diffs) >= 3:
                obs, null = np.median(diffs), np.median(np.array(nulls), axis=0)
                lo, hi = boot_ci(diffs, rng)
                summ.append(dict(test="D3_known_use_dependency", subset=src, n_cohorts=len(diffs),
                                 estimate=round(obs, 4), ci_lo=lo, ci_hi=hi,
                                 p=(1 + (null >= obs).sum()) / (1 + len(null)),
                                 detail="dep_max of drugs approved for this cancer minus the rest"))
    for n, t in (("d1", d1), ("d2", d2), ("d3", d3), ("summary", summ)):
        pd.DataFrame(t).to_csv(OUT / f"{n}.tsv", sep="\t", index=False)
    json.dump(dict(lines_screened=int(len(Dz)), genes=int(Dz.shape[1]), cohorts_scored=len(gd),
                   yardstick_above_80pct=n_ok, yardstick_n=len(Y), fda_rows=len(D),
                   rows_with_dep=int(D.dep_max.notna().sum()), min_lines=MIN_LINES, min_wcov=MIN_WCOV,
                   min_drugs=MIN_DRUGS, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED), open(OUT / "run.json", "w"), indent=1)
    print(pd.DataFrame(summ).to_string(index=False))


if __name__ == "__main__":
    main()
