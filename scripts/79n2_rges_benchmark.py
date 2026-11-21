#!/usr/bin/env python
"""RGES benchmark: gene-level reversal (Chen et al.) against our regulon-level reversal (D174).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
SPEC: docs/SPEC_rges_benchmark.md — written and committed BEFORE this was run.

Everything but the level of reversal is held identical to 79e: same disease contrast, same
`drug_line_profiles_bing_specific.f32` matrix, same Celligner in-key lines and weights, same
weighted-mean aggregation. The only difference is genes vs master regulators.

RGES per the spec: top 150 up / bottom 150 down by `z` among adj.P < 0.05 genes present in the
LINCS universe; Lamb/CMap KS enrichment of each set against the drug's descending gene ranking;
RGES = ES_up - ES_down, forced to 0 when the two share a sign. Reported as `rges_rev = -RGES`
so that, like our `rev`, higher means stronger reversal.

Output: results/lincs/rges_benchmark.tsv (+ _run.json)
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
S4 = ROOT / "data/processed/step4_v2"
RES = ROOT / "results/lincs"
OUT = RES / "rges_benchmark.tsv"
N_SET, MIN_SET, FDR = 150, 50, 0.05
PROFILE = os.environ.get("PROFILE", "specific")
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "6")


def ks_es(ranks, n_genes):
    """Lamb/CMap KS enrichment of a gene set, given its 1-based ranks in the drug profile.

    ranks: (n_profiles, n_set) int array, already sorted ascending along axis 1.
    Returns a (n_profiles,) array: the signed maximum deviation, positive if the set is enriched
    at the TOP of the ranking (genes the drug raises).
    """
    n_prof, t = ranks.shape
    j = np.arange(1, t + 1, dtype=np.float64)[None, :]
    a = (j / t - ranks / n_genes).max(axis=1)          # set runs ahead of uniform
    b = (ranks / n_genes - (j - 1) / t).max(axis=1)    # set lags behind uniform
    return np.where(a > b, a, -b)


def main():
    t0 = time.time()
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    f = PROC / (f"drug_line_profiles_bing{'_specific' if PROFILE == 'specific' else ''}.f32")
    P = np.fromfile(f, dtype=np.float32).reshape(len(idx), len(genes))
    print(f"profiles ({PROFILE}): {P.shape[0]:,} drug x line x {P.shape[1]:,} genes", flush=True)

    # rank every gene within every profile, once. rank 1 = most raised by the drug.
    order = np.argsort(-P, axis=1, kind="stable")
    rank = np.empty_like(order)
    np.put_along_axis(rank, order, np.arange(1, P.shape[1] + 1, dtype=order.dtype)[None, :], axis=1)
    del order, P
    gpos = {g: i for i, g in enumerate(genes.gene_id)}
    n_genes = rank.shape[1]

    summ = pd.read_csv(ROOT / "results/step4_v2/signature_summary.tsv", sep="\t")
    summ = summ[summ.status == "ok"]
    lw = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    lw = lw[lw.in_key]
    rows, skipped = [], []

    for _, s in summ.iterrows():
        g, contrast = s.group_id, s.contrast
        sf = S4 / f"{g}__{contrast}.tsv.gz"
        if not sf.exists():
            skipped.append((g, "no signature file"))
            continue
        sig = pd.read_csv(sf, sep="\t", usecols=["gene_id", "z", "adj.P.Val"])
        sig = sig[(sig["adj.P.Val"] < FDR) & sig.gene_id.isin(gpos)].dropna(subset=["z"])
        up = sig.nlargest(N_SET, "z").gene_id.map(gpos).values
        dn = sig.nsmallest(N_SET, "z").gene_id.map(gpos).values
        if len(up) < MIN_SET or len(dn) < MIN_SET:
            skipped.append((g, f"sets too small ({len(up)} up, {len(dn)} down)"))
            continue
        wg = lw[lw.group_id == g].set_index("ModelID").weight
        sel = idx.index[idx.ModelID.isin(wg.index)].values
        if not len(sel):
            skipped.append((g, "no in-key lines with profiles"))
            continue
        ru = np.sort(rank[np.ix_(sel, up)].astype(np.float64), axis=1)
        rd = np.sort(rank[np.ix_(sel, dn)].astype(np.float64), axis=1)
        es_u, es_d = ks_es(ru, n_genes), ks_es(rd, n_genes)
        rges = np.where(np.sign(es_u) == np.sign(es_d), 0.0, es_u - es_d)
        sub = idx.loc[sel, ["pert_id", "ModelID"]].copy()
        sub["rges_rev"] = -rges
        sub["w"] = sub.ModelID.map(wg).fillna(0).values
        agg = (sub.assign(wx=sub.rges_rev * sub.w).groupby("pert_id")[["wx", "w"]].sum())
        agg = agg[agg.w > 0]
        rows.append(pd.DataFrame({"group_id": g, "pert_id": agg.index,
                                  "rges_rev": (agg.wx / agg.w).values,
                                  "n_lines": sub.groupby("pert_id").size().reindex(agg.index).values}))
        print(f"  {g:<22} {len(agg):>5,} drugs   {len(sel):>4} line-profiles", flush=True)

    d = pd.concat(rows, ignore_index=True)
    d["rges_pct"] = (d.groupby("group_id").rges_rev
                     .transform(lambda v: v.rank(method="max") / len(v) * 100).round(2))
    d.to_csv(OUT, sep="\t", index=False)
    run = dict(spec="docs/SPEC_rges_benchmark.md", profile=PROFILE, n_set=N_SET, fdr=FDR,
               n_rows=int(len(d)), n_cohorts=int(d.group_id.nunique()),
               skipped=[f"{a}: {b}" for a, b in skipped], runtime_s=round(time.time() - t0, 1))
    (RES / "rges_benchmark_run.json").write_text(json.dumps(run, indent=1))
    print(f"\n{len(d):,} (cohort, drug) RGES scores over {d.group_id.nunique()} cohorts "
          f"in {run['runtime_s']}s; skipped {len(skipped)}")
    for a, b in skipped:
        print(f"  skipped {a}: {b}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
