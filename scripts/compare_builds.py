#!/usr/bin/env python
"""Decompose what changed between two versions of the same result table.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

    python scripts/compare_builds.py --old OLD.tsv --new NEW.tsv \
        --key group_id,pert_id --value rev [--by group_id] [--top 10]

Written for the D156 rebuild (D164). The question a rebuild has to answer is never
"did the number move" but "*what moved it*". Three mechanisms produce a different
score and they need different responses, so the report separates them:

  COMPOSITION  rows the new build has that the old did not (and vice versa). If the
               new rows score differently from the shared ones, the headline number
               moved because the POPULATION changed, not because any drug-disease
               pair was re-scored. This is the D157 NEW-stratum dilution: a worse
               aggregate that is not a worse method.
  RE-SCORING   on rows present in BOTH, how far and which way the value moved.
               Correlation, signed median shift, and how many cross zero.
  MOVERS       the individual pairs that moved most, so a shift can be traced to
               specific drugs rather than described only as a distribution.

Reports only. It makes no accept/reject call -- that is a judgement for the decision
log, not for a script.
"""
import argparse, sys
import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr, pearsonr
except ImportError:
    spearmanr = pearsonr = None


def pct(n, d):
    return f"{n:,} ({100 * n / d:.1f}%)" if d else f"{n:,} (-)"


def dist(v, label):
    v = pd.Series(v).dropna()
    if not len(v):
        return f"  {label:22} (empty)"
    return (f"  {label:22} n={len(v):>7,}  mean {v.mean():+.3f}  median {v.median():+.3f}  "
            f"sd {v.std():.3f}  p10 {v.quantile(.1):+.3f}  p90 {v.quantile(.9):+.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--key", required=True, help="comma-separated key columns")
    ap.add_argument("--value", required=True, help="comma-separated value columns")
    ap.add_argument("--by", default="", help="optional column for a per-stratum breakdown")
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--label-old", default="OLD")
    ap.add_argument("--label-new", default="NEW")
    a = ap.parse_args()

    keys = a.key.split(",")
    vals = a.value.split(",")
    o = pd.read_csv(a.old, sep="\t", low_memory=False)
    n = pd.read_csv(a.new, sep="\t", low_memory=False)

    print(f"\n{'=' * 78}\n{a.label_old}: {a.old}\n{a.label_new}: {a.new}\n{'=' * 78}")
    missing = [c for c in keys + vals if c not in o.columns or c not in n.columns]
    if missing:
        print(f"!! columns absent from one side, cannot compare: {missing}")
        print(f"   old has: {list(o.columns)[:40]}")
        return 2
    print(f"rows      {len(o):,} -> {len(n):,}   ({len(n) - len(o):+,})")
    print(f"columns   {len(o.columns)} -> {len(n.columns)}"
          + (f"   new cols: {sorted(set(n.columns) - set(o.columns))}" if set(n.columns) - set(o.columns) else "")
          + (f"   dropped: {sorted(set(o.columns) - set(n.columns))}" if set(o.columns) - set(n.columns) else ""))

    ok = o.set_index(keys)
    nk = n.set_index(keys)
    if ok.index.has_duplicates or nk.index.has_duplicates:
        print(f"!! duplicate keys (old {ok.index.duplicated().sum():,}, new {nk.index.duplicated().sum():,}) "
              f"-- keeping first of each; the key may be incomplete")
        ok, nk = ok[~ok.index.duplicated()], nk[~nk.index.duplicated()]
    shared = ok.index.intersection(nk.index)
    gained = nk.index.difference(ok.index)
    lost = ok.index.difference(nk.index)

    print(f"\n-- COMPOSITION  (did the population change?)")
    print(f"  shared  {pct(len(shared), len(nk))} of new")
    print(f"  gained  {pct(len(gained), len(nk))} of new")
    print(f"  lost    {pct(len(lost), len(ok))} of old")

    for v in vals:
        print(f"\n-- {v}")
        so, sn = pd.to_numeric(ok.loc[shared, v], errors="coerce"), pd.to_numeric(nk.loc[shared, v], errors="coerce")
        gv = pd.to_numeric(nk.loc[gained, v], errors="coerce") if len(gained) else pd.Series(dtype=float)
        lv = pd.to_numeric(ok.loc[lost, v], errors="coerce") if len(lost) else pd.Series(dtype=float)
        print(dist(pd.to_numeric(ok[v], errors="coerce"), f"{a.label_old} all"))
        print(dist(pd.to_numeric(nk[v], errors="coerce"), f"{a.label_new} all"))
        print(dist(sn, f"{a.label_new} shared"))
        print(dist(gv, f"{a.label_new} gained"))
        if len(lv):
            print(dist(lv, f"{a.label_old} lost"))

        # is the aggregate shift composition or re-scoring?
        if len(gv.dropna()) and len(sn.dropna()):
            d_comp = gv.mean() - sn.mean()
            print(f"  >> gained rows score {d_comp:+.3f} vs shared rows in the same build.")
            print(f"     {'DILUTION: the aggregate moves because the population changed.' if abs(d_comp) > 0.05 else 'gained rows look like the shared ones; composition is not driving the aggregate.'}")

        m = so.notna() & sn.notna()
        print(f"\n  RE-SCORING on {m.sum():,} shared rows with both values:")
        if m.sum() >= 3:
            d = (sn - so)[m]
            if spearmanr is not None:
                rho = spearmanr(so[m], sn[m]).correlation
                r = pearsonr(so[m], sn[m])[0]
                print(f"    Spearman {rho:+.3f}   Pearson {r:+.3f}")
            print(f"    delta: median {d.median():+.3f}  mean {d.mean():+.3f}  "
                  f"p10 {d.quantile(.1):+.3f}  p90 {d.quantile(.9):+.3f}")
            print(f"    |delta| > 0.5: {pct((d.abs() > .5).sum(), len(d))}    "
                  f"sign flips: {pct(((so[m] > 0) != (sn[m] > 0)).sum(), int(m.sum()))}")
            print(f"    unchanged (|d|<1e-9): {pct((d.abs() < 1e-9).sum(), len(d))}")
            big = d.reindex(d.abs().sort_values(ascending=False).index)[:a.top]
            print(f"    biggest movers:")
            for k, dv in big.items():
                nm = nk.loc[k].get("cmap_name", "") if "cmap_name" in nk.columns else ""
                print(f"      {str(k):<46} {so.loc[k]:+7.3f} -> {sn.loc[k]:+7.3f}  ({dv:+.3f})  {nm}")
        else:
            print("    too few shared rows to compare")

        if a.by and a.by in nk.reset_index().columns and m.sum() >= 3:
            print(f"\n  by {a.by} (shared rows, median delta, top/bottom 5):")
            t = pd.DataFrame({"by": nk.loc[shared].reset_index()[a.by].values if a.by not in keys
                                    else shared.get_level_values(a.by),
                              "d": (sn - so).values})
            g = t.dropna().groupby("by").d.agg(["median", "size"]).sort_values("median")
            for side in (g.head(5), g.tail(5)):
                for k, row in side.iterrows():
                    print(f"      {str(k):<32} {row['median']:+7.3f}  (n={int(row['size']):,})")
                print("      ...")
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
