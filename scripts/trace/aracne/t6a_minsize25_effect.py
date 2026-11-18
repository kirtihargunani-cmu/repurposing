#!/usr/bin/env python
"""D139-S0: do small regulons actually affect the drug score, and what would minsize=25 do?

Exact, not a re-run. Two facts make that possible, both verified first:
  * VIPER's aREA NES for a regulator is computed independently of the other regulators, EXCEPT
    through `eset.filter=TRUE`, which restricts the gene universe to genes appearing in the
    regulon set. Measured: raising the floor 10 -> 25 in BRCA_LUMA shrinks that universe only
    9,146 -> 9,133 genes (99.9% retained); median |dNES| 2.0e-3 and only 0.02% of
    regulator x profile cells move by more than 0.1. With eset.filter=FALSE the NES is exactly
    invariant (max difference 0.000).
  * `rev = -(D . s)/||s||` is LINEAR in the activity D, and a drug's D is a weighted mean of its
    profiles' activities, so per-profile rev averages to drug-level rev. Correlating per profile
    is therefore faithful to the drug score.

So raising minsize is, to within 0.02% of cells, exactly dropping the small regulons from the sum.
This script computes that drop for every group under SCORE=weighted (the production setting:
every tested regulator, weight = disease-side group NES).

Output: results/lincs/trace_aracne/minsize25_effect.tsv
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd


ROOT = Path("/home/kirti/projects/repurposing")
PROC = ROOT / "data/processed/lincs"
OUT = ROOT / "results/lincs/trace_aracne"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
FLOOR = 25


def net_for(cohort, override):
    if cohort in override:
        return override[cohort]
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def main():
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    mr = pd.read_csv(ROOT / "results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t")
    dreg = pd.read_csv(ROOT / "results/lincs/drug_viper_regulators_v2_tf.tsv", sep="\t")
    nlincs = {(r.net, r.regulator): r.n_lincs for r in dreg.itertuples()}
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    grp = pd.read_csv(ROOT / "results/lincs/fda_ranking_v2_consensus.tsv", sep="\t")[
        ["group_id"]].drop_duplicates()

    rows = []
    for g in sorted(mr.group_id.unique()):
        net = net_for(g, ovr)
        f = PROC / "drug_viper_specific_v2_tf" / f"{net}.f32"
        if not f.exists():
            continue
        regs = pd.Index((PROC / "drug_viper_specific_v2_tf" / f"{net}.regulators.txt").read_text().split())
        A = np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
        m = mr[mr.group_id == g].copy()
        m["ci"] = regs.get_indexer(m.regulator)
        m = m[m.ci >= 0]
        m["n_lincs"] = [nlincs.get((net, r), 0) for r in m.regulator]
        if len(m) < 2:
            continue
        s = m.nes.values
        rev_all = -(A[:, m.ci.values] @ s) / np.sqrt((s ** 2).sum())
        # the per-regulator contribution mass, averaged over profiles
        contrib = np.abs(A[:, m.ci.values] * s).mean(0)
        for side, col in (("disease", "size"), ("drug", "n_lincs"), ("either", None)):
            keep = ((m["size"] >= FLOOR) & (m.n_lincs >= FLOOR)) if col is None else (m[col] >= FLOOR)
            if keep.sum() < 2:
                continue
            s2 = s[keep.values]
            rev25 = -(A[:, m.ci.values[keep.values]] @ s2) / np.sqrt((s2 ** 2).sum())
            rows.append(dict(
                group_id=g, network=net, side=side, n_tested=len(m), n_kept=int(keep.sum()),
                pct_regulators_dropped=round(100 * float((~keep).mean()), 2),
                pct_weight_mass_dropped=round(100 * float((s[~keep.values] ** 2).sum() / (s ** 2).sum()), 2),
                pct_contribution_dropped=round(100 * float(contrib[~keep.values].sum() / contrib.sum()), 2),
                spearman_rev=round(float(np.corrcoef(pd.Series(rev_all).rank(), pd.Series(rev25).rank())[0, 1]), 5),
                pearson_rev=round(float(np.corrcoef(rev_all, rev25)[0, 1]), 5)))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "minsize25_effect.tsv", sep="\t", index=False)
    for side in ("disease", "drug", "either"):
        d = df[df.side == side]
        print(f"\n=== floor {FLOOR} applied on the {side} side ({len(d)} groups) ===")
        print(f"  regulators dropped   median {d.pct_regulators_dropped.median():5.1f}%  "
              f"range {d.pct_regulators_dropped.min():.1f}-{d.pct_regulators_dropped.max():.1f}%")
        print(f"  weight mass dropped  median {d.pct_weight_mass_dropped.median():5.1f}%  "
              f"range {d.pct_weight_mass_dropped.min():.1f}-{d.pct_weight_mass_dropped.max():.1f}%")
        print(f"  contribution dropped median {d.pct_contribution_dropped.median():5.1f}%  "
              f"range {d.pct_contribution_dropped.min():.1f}-{d.pct_contribution_dropped.max():.1f}%")
        print(f"  Spearman rev(10) vs rev(25) across profiles: median {d.spearman_rev.median():.4f}  "
              f"min {d.spearman_rev.min():.4f}  "
              f"groups below 0.99: {int((d.spearman_rev < 0.99).sum())} of {len(d)}")
    print("\nwrote", OUT / "minsize25_effect.tsv")


if __name__ == "__main__":
    sys.exit(main())
