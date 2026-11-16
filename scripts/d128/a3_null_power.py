#!/usr/bin/env python
"""D128-A3: how much power the degree-matched-target null actually has, across the whole population.

The A4 ESR1 control returned a blank z for this null, which means only one candidate survived. That is a property
of the candidate pool, not of the drug, so it can be measured for every unit WITHOUT scoring anything: the pool is
the 497 D121 target units, the matching rule is degree (footprint size and TF count within +-25%), the shared-target
exclusion removes candidates that share a target gene with the observation, and the D128-A3 de-duplication removes
candidates whose resolved footprint is identical to one already drawn.

This runs no null and reads no drug response. It reports, per unit and weighted by how many population rows use
that unit, how many candidates survive each step.

Output: results/lincs/d128/a3/null_power.tsv (per unit) and null_power_summary.tsv.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import nulls as NU  # noqa: E402

A3 = C.OUT / "a3"


def main():
    cov = pd.read_csv(C.OUT / "a0/coverage.tsv", sep="\t")
    F = pd.read_parquet(C.OUT / "a1/target_gene_footprints.parquet")
    S = pd.read_parquet(C.OUT / "a2/scores_by_line.parquet")
    S = S[(S.gene_arm == "all_bing") & (S.M_status == "scored")] if "gene_arm" in S.columns else S[S.M_status == "scored"]
    use = S.groupby("unit_id").size().rename("n_population_rows")

    sig = {}
    for uid, d in F[F.status == "resolved"].groupby("unit_id"):
        sig[uid] = hash(tuple(sorted(zip(d.gene_id.tolist(), d.expected_sign_pre_action.tolist()))))

    rows = []
    for r in cov.itertuples(index=False):
        if r.M_unit_status != "eligible":
            continue
        n_g, n_t = int(r.n_sign_resolved), int(r.n_resolved_tfs)
        lo, hi = n_g * 0.75, n_g * 1.25
        tlo, thi = max(0, n_t * 0.75), n_t * 1.25 + 1
        pool = cov[(cov.unit_id != r.unit_id) & (cov.n_sign_resolved >= lo) & (cov.n_sign_resolved <= hi)
                   & (cov.n_resolved_tfs >= tlo) & (cov.n_resolved_tfs <= thi)]
        after = NU.degree_matched_units(r.unit_id, cov, n_g, n_t, k=10 ** 6,
                                        exclude_genes=str(r.genes).split("|"))
        distinct = {sig[u] for u in after.unit_id if u in sig}
        rows.append(dict(unit_id=r.unit_id, cmap_name=r.cmap_name, genes=r.genes, unit_kind=r.unit_kind,
                         n_sign_resolved=n_g, n_resolved_tfs=n_t,
                         n_degree_matched=len(pool), n_after_shared_target_exclusion=len(after),
                         n_distinct_footprints=len(distinct),
                         n_drawn_at_k20=min(len(distinct), 20),
                         n_population_rows=int(use.get(r.unit_id, 0))))
    D = pd.DataFrame(rows)
    D.to_csv(A3 / "null_power.tsv", sep="\t", index=False)

    w = D.n_population_rows
    def share(mask, weights=None):
        ww = (np.ones(len(D)) if weights is None else weights.to_numpy(float))
        return round(float((ww * mask.to_numpy(float)).sum() / ww.sum()), 4) if ww.sum() else np.nan

    SU = pd.DataFrame([
        dict(level="per unit (M-eligible units)", n=len(D),
             median_distinct=float(D.n_distinct_footprints.median()),
             share_lt_10=share(D.n_distinct_footprints < 10), share_lt_3=share(D.n_distinct_footprints < 3),
             share_eq_0=share(D.n_distinct_footprints == 0)),
        dict(level="weighted by population rows that use the unit", n=int(w.sum()),
             median_distinct=float(np.median(np.repeat(D.n_distinct_footprints.values,
                                                       D.n_population_rows.values)) if w.sum() else np.nan),
             share_lt_10=share(D.n_distinct_footprints < 10, w), share_lt_3=share(D.n_distinct_footprints < 3, w),
             share_eq_0=share(D.n_distinct_footprints == 0, w)),
    ])
    SU.to_csv(A3 / "null_power_summary.tsv", sep="\t", index=False)
    print(SU.to_string())
    print("\nunits with the thinnest degree-matched null (most used first):")
    print(D[D.n_distinct_footprints < 10].sort_values("n_population_rows", ascending=False)
          .head(20)[["cmap_name", "genes", "n_sign_resolved", "n_resolved_tfs", "n_degree_matched",
                     "n_after_shared_target_exclusion", "n_distinct_footprints", "n_population_rows"]].to_string())


if __name__ == "__main__":
    main()
