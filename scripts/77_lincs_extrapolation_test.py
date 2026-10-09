#!/usr/bin/env python
"""Can a drug's LINCS response in one cell line be predicted from its response in other (neighboring) lines?

Runs in the `celligner` env. Inputs: data/processed/lincs/level5_landmark.{npy,_sigs.tsv.gz} (scripts/76),
Celligner line x line distances (data/processed/celligner/dist_line_line.tsv.gz, scripts/74).
1. Active signatures only: is_hiq == 1 and TAS >= 0.4 (user, 2026-10-09). One profile per drug x line = mean of its
   active signatures over doses / times (978 landmark genes).
2. Drug-specific effect: strong signatures share a generic response (random strong pairs r ~0.2). Per line, the
   generic response g_l = mean active profile of that line over all its drugs; specific = profile minus its
   projection on g_l. Checked by correlation of same drug / different line vs random pairs, before and after.
3. Leave-one-line-out: for each deep line L (>= MIN_DRUGS active drugs, cancer, in the Celligner distances) and each
   drug d active in L and in >= 1 other line, predict d's specific profile in L from the other lines without L:
     nearest   the closest other line (Celligner distance) with d
     weighted  softmax(-distance / tau) average over the other lines with d, tau in TAUS
     mean      unweighted mean over the other lines with d (no distance information)
     null      a random other drug's profile in the nearest line (chance level for this line pair)
   Score: Pearson r of predicted vs observed over the 978 genes; also by distance of the nearest line with d.
   Ceiling: r between two single active signatures of the same drug in the same line (other dose / time), specific.
Outputs: results/lincs/extrapolation_{generic_check, lolo_pairs, lolo_summary, lolo_by_distance, lolo_by_line}.tsv
         data/processed/lincs/active_profiles.npy (+ _index.tsv.gz: pert_id, ModelID), specific_profiles.npy
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
TAS_MIN, MIN_DRUGS, SEED = 0.4, 50, 0
TAUS = [10, 25, 50]
DIST_BINS = [0, 75, 100, 125, 150, 200, 1000]


def rowcorr(a, b):
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))


def pair_check(P, idx, rng, n=5000):
    """median r: same drug / different line vs random pairs"""
    by = idx.groupby("pert_id").indices
    multi = [v for v in by.values() if len(v) >= 2]
    sd = np.array([rng.choice(v, 2, replace=False) for v in (multi[i] for i in rng.integers(0, len(multi), n))])
    rd = rng.integers(0, len(idx), (n, 2))
    rd = rd[idx.pert_id.values[rd[:, 0]] != idx.pert_id.values[rd[:, 1]]]
    return float(np.median(rowcorr(P[sd[:, 0]], P[sd[:, 1]]))), float(np.median(rowcorr(P[rd[:, 0]], P[rd[:, 1]])))


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    X = np.load(PROC / "level5_landmark.npy", mmap_mode="r")
    sig = pd.read_csv(PROC / "level5_landmark_sigs.tsv.gz", sep="\t")
    dist = pd.read_csv(ROOT / "data/processed/celligner/dist_line_line.tsv.gz", sep="\t", index_col=0)

    act = sig[(sig.is_hiq == 1) & (sig.tas >= TAS_MIN) & sig.ModelID.isin(dist.index)]
    print(f"active signatures: {len(act)} ({act.pert_id.nunique()} drugs, {act.ModelID.nunique()} lines)")
    grp = act.groupby(["pert_id", "ModelID"]).indices
    idx = pd.DataFrame(list(grp.keys()), columns=["pert_id", "ModelID"])
    rows = act.index.values
    P = np.stack([np.asarray(X[np.sort(rows[v])]).mean(axis=0) for v in grp.values()]).astype(np.float64)
    idx["n_sig"] = [len(v) for v in grp.values()]
    np.save(PROC / "active_profiles.npy", P.astype(np.float32))
    idx.to_csv(PROC / "active_profiles_index.tsv.gz", sep="\t", index=False)

    # generic response per line; specific = residual after projecting it out
    S = P.copy()
    G = {}
    for m, ix in idx.groupby("ModelID").indices.items():
        g = P[ix].mean(axis=0)
        G[m] = g = g / np.linalg.norm(g)
        S[ix] = P[ix] - np.outer(P[ix] @ g, g)
    np.save(PROC / "specific_profiles.npy", S.astype(np.float32))
    chk = []
    for name, M in [("raw", P), ("specific", S)]:
        same, rand = pair_check(M, idx, rng)
        chk.append(dict(profiles=name, same_drug_diff_line_median_r=round(same, 3), random_pair_median_r=round(rand, 3),
                        gap=round(same - rand, 3)))
    # ceiling: two single active signatures, same drug and line
    ceil = []
    for (d, m), ix in act.groupby(["pert_id", "ModelID"]).indices.items():
        if len(ix) >= 2:
            a, b = rng.choice(act.index.values[ix], 2, replace=False)
            va, vb = np.asarray(X[a], float), np.asarray(X[b], float)
            va, vb = va - (va @ G[m]) * G[m], vb - (vb @ G[m]) * G[m]
            ceil.append(rowcorr(va[None], vb[None])[0])
    chk.append(dict(profiles="ceiling: same drug, same line, 2 single signatures (specific)",
                    same_drug_diff_line_median_r=round(float(np.median(ceil)), 3), random_pair_median_r=np.nan,
                    gap=np.nan))
    chk = pd.DataFrame(chk)
    chk.to_csv(RES / "extrapolation_generic_check.tsv", sep="\t", index=False)
    print(chk.to_string(index=False))

    # leave-one-line-out on the specific profiles
    per_line = idx.groupby("ModelID").size()
    deep = per_line[per_line >= MIN_DRUGS].index
    drug_rows = idx.groupby("pert_id").indices
    out = []
    for L in deep:
        for i in np.where(idx.ModelID.values == L)[0]:
            d = idx.pert_id.values[i]
            others = [j for j in drug_rows[d] if j != i]
            if not others:
                continue
            om = idx.ModelID.values[others]
            dd = dist.loc[L, om].values
            o = np.argsort(dd)
            j_near = others[o[0]]
            rec = dict(ModelID=L, pert_id=d, n_other_lines=len(others), nearest_line=om[o[0]],
                       dist_nearest=float(dd[o[0]]), dist_mean=float(dd.mean()))
            rec["r_nearest"] = rowcorr(S[[i]], S[[j_near]])[0]
            rec["r_mean"] = rowcorr(S[[i]], S[others].mean(axis=0, keepdims=True))[0]
            for tau in TAUS:
                w = np.exp(-(dd - dd.min()) / tau)
                rec[f"r_weighted_tau{tau}"] = rowcorr(S[[i]], (w[:, None] * S[others]).sum(0, keepdims=True) / w.sum())[0]
            pool = np.where(idx.ModelID.values == om[o[0]])[0]
            pool = pool[idx.pert_id.values[pool] != d]
            rec["r_null"] = rowcorr(S[[i]], S[[rng.choice(pool)]])[0] if len(pool) else np.nan
            out.append(rec)
    res = pd.DataFrame(out)
    res.to_csv(RES / "extrapolation_lolo_pairs.tsv", sep="\t", index=False, float_format="%.4f")

    meth = ["r_nearest", "r_mean"] + [f"r_weighted_tau{t}" for t in TAUS] + ["r_null"]
    summ = pd.DataFrame({m: [res[m].median(), (res[m] > 0.3).mean()] for m in meth},
                        index=["median_r", "frac_r_gt_0.3"]).T.round(3)
    summ["n_pairs"] = len(res)
    summ.to_csv(RES / "extrapolation_lolo_summary.tsv", sep="\t", index_label="method")
    res["dist_bin"] = pd.cut(res.dist_nearest, DIST_BINS)
    bd = res.groupby("dist_bin").agg(n=("pert_id", "size"), r_nearest=("r_nearest", "median"),
                                     r_null=("r_null", "median"), frac_gt_0_3=("r_nearest", lambda x: (x > 0.3).mean()))
    bd.round(3).to_csv(RES / "extrapolation_lolo_by_distance.tsv", sep="\t")
    bl = res.groupby("ModelID").agg(n_drugs=("pert_id", "size"), median_dist_nearest=("dist_nearest", "median"),
                                    r_nearest=("r_nearest", "median"), r_mean=("r_mean", "median"),
                                    r_null=("r_null", "median")).round(3).sort_values("n_drugs", ascending=False)
    bl.to_csv(RES / "extrapolation_lolo_by_line.tsv", sep="\t")

    pd.set_option("display.width", 200)
    print(f"deep lines: {len(deep)}; drug x held-out-line pairs: {len(res)}")
    print(summ.to_string())
    print(bd.round(3).to_string())
    print(bl.to_string())
    print(f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
