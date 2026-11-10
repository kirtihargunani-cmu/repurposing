#!/usr/bin/env python
"""D125 A1: the lambda sweep on the D84 perturbation benchmark (docs/SPEC_graph_metric_scoring.md A2/A3; D125-A1).

What it asks: does reading a regulator through the graph metric K = (L + lambda I)^-1 read TF perturbations more
often in the expected direction than reading it directly? `rev_L` needs a disease vector; the benchmark has none, so
the read is `rev_L` at d = e_i, the unit vector on regulator i (D125-A1 (1)):

    s_i = (e_i' K g) / sqrt(e_i' K e_i) = (K g)_i / sqrt(K_ii),   and at K = I/lambda, s_i = g_i / sqrt(lambda)

so the lambda -> inf arm is production exactly. That arm is computed at K = I and must reproduce the published D84
numbers before any swept arm is read.

Graph (A2): nodes = the regulators scored in that network, edge = Szymkiewicz-Simpson overlap of the pruned target
sets, edges < EDGE_MIN dropped, no self-loops, sign ignored; L = I - D^-1/2 A D^-1/2, isolated nodes keep D_ii = 0.
Built per network and per regulon set; K is applied inside a network, then the read is averaged over the signature's
lineage-matched networks exactly as 79n averages NES.

Selection quantity (D125-A1 (3)): the consensus read -- mean of ARACNe v2 TF and CollecTRI where both score the TF,
else the one that does -- share of D84 main-analysis signatures with the expected sign. Per-set rates, the
knockdown-only subset and the no-own-check sensitivity are reported beside it.

Early stop (A3): if no lambda beats the lambda = inf arm by >= MIN_GAIN_PTS percentage points, D125 is recorded as
not met and A2 does not run. D84 is the selection set and is never an adoption criterion.

Outputs: results/lincs/d125/a1/{reads,rates,per_tf,graph}.tsv, run.json.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parent))
from a1_config import (D84_PUBLISHED, EDGE_MIN, IDENTITY_TOL, LAMBDAS, MIN_GAIN_PTS, NES_DIR, OUT, PRIMARY_ANALYSIS,
                       PRIMARY_READ, SETS, THREADS)

INF = "inf"   # the lambda = inf arm, computed at K = I


def load_tag(tag):
    regs = (NES_DIR / f"{tag}.regulators.txt").read_text().split()
    sigs = (NES_DIR / f"{tag}.sigs.txt").read_text().split()
    g = np.fromfile(NES_DIR / f"{tag}.nes.f32", dtype=np.float32).reshape((len(regs), len(sigs)), order="F")
    tg = pd.read_csv(NES_DIR / f"{tag}.targets.tsv.gz", sep="\t")
    return regs, sigs, g.astype(np.float64), tg


def laplacian(regs, tg):
    """A2: Szymkiewicz-Simpson overlap of the pruned target sets, thresholded; symmetric normalized Laplacian."""
    ri = {r: i for i, r in enumerate(regs)}
    tg = tg[tg.regulator.isin(ri)]
    genes = {g: j for j, g in enumerate(pd.unique(tg.target))}
    M = csr_matrix((np.ones(len(tg)), (tg.regulator.map(ri).to_numpy(), tg.target.map(genes).to_numpy())),
                   shape=(len(regs), len(genes)))
    inter = (M @ M.T).toarray()
    size = np.asarray(M.sum(1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        A = np.where(np.minimum.outer(size, size) > 0, inter / np.minimum.outer(size, size), 0.0)
    np.fill_diagonal(A, 0.0)
    A[A < EDGE_MIN] = 0.0
    deg = A.sum(1)
    dinv = np.where(deg > 0, 1.0 / np.sqrt(np.where(deg > 0, deg, 1.0)), 0.0)   # isolated: D_ii = 0 (D125-A1 (5))
    L = np.eye(len(regs)) - dinv[:, None] * A * dinv[None, :]
    stats = dict(regulators=len(regs), isolated=int((deg == 0).sum()),
                 edges=int((A > 0).sum() // 2), mean_degree=float((A > 0).sum(1).mean()),
                 mean_weight=float(A[A > 0].mean()) if (A > 0).any() else np.nan)
    return L, stats


def read_one(args):
    """One (set, network): the perturbed TF's transformed read and its percentile, for every lambda."""
    tag, tf_of_sig = args
    regs, sigs, G, tg = load_tag(tag)
    L, stats = laplacian(regs, tg)
    ri = {r: i for i, r in enumerate(regs)}
    rows = [(k, ri[tf], j) for j, k in enumerate(sigs) if (tf := tf_of_sig.get(k)) in ri]
    out = []
    for lam in [INF] + LAMBDAS:
        if lam == INF:
            S = G                                   # K = I: production, exactly
        else:
            K = cho_solve(cho_factor(L + lam * np.eye(len(regs)), lower=True), np.eye(len(regs)))
            S = (K @ G) / np.sqrt(np.diag(K))[:, None]
        pct = (S.argsort(0).argsort(0)) / (len(regs) - 1)   # rank percentile among all regulators, as 79n
        for sid, i, j in rows:
            out.append((tag, lam, sid, float(S[i, j]), float(pct[i, j])))
    return pd.DataFrame(out, columns=["tag", "lam", "sig_id", "nes", "pct"]), dict(tag=tag, **stats)


def combine(r, sg):
    """79n's averaging: K inside a network, then mean over the signature's matched networks; then the D75 consensus."""
    net_of = r.tag.str.extract(r"^aracne_v2_tf__(.+)$", expand=False)
    r = r.assign(set=np.where(net_of.notna(), "aracne_v2_tf", r.tag), net=net_of)
    use = (sg[sg.set == "aracne_v2_tf"][["sig_id", "nets"]].dropna()
           .assign(net=lambda d: d.nets.str.split(";")).explode("net"))
    a = (r[r.set == "aracne_v2_tf"].merge(use[["sig_id", "net"]], on=["sig_id", "net"])
         .groupby(["lam", "sig_id"], as_index=False)[["nes", "pct"]].mean().assign(set="aracne_v2_tf"))
    c = r[r.set == "collectri"][["lam", "sig_id", "nes", "pct"]].assign(set="collectri")
    per_set = pd.concat([a, c], ignore_index=True)
    cons = per_set.groupby(["lam", "sig_id"], as_index=False)[["nes", "pct"]].mean().assign(set="consensus")
    return pd.concat([per_set, cons], ignore_index=True)


def rates(reads):
    """D84 / 79v reading: correct = NES sign as expected; fpct = pct (knockdown) or 1 - pct (overexpression)."""
    out = []
    for analysis, d in [("main", reads[reads.own_ok]), ("no_own_check", reads)]:
        for kind, x in [("all", d), ("kd", d[d.direction < 0]), ("oe", d[d.direction > 0])]:
            for (st, lam), y in x.groupby(["set", "lam"]):
                out.append(dict(analysis=analysis, kind=kind, read=st, lam=lam, n=len(y),
                                n_tfs=int(y.cmap_name.nunique()), correct=round(float(y.correct.mean()), 4),
                                auc=round(1 - float(y.fpct.mean()), 4), median_fpct=round(float(y.fpct.median()), 4)))
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sg = pd.read_csv(NES_DIR / "d84_signatures.tsv.gz", sep="\t")
    sg["own_ok"] = sg.own_ok.astype(bool)
    meta = sg.drop_duplicates("sig_id")[["sig_id", "cmap_name", "direction", "own_ok", "technology", "cell_iname"]]
    tf_of_sig = dict(zip(meta.sig_id, meta.cmap_name))
    tags = ["collectri"] + sorted(p.name[: -len(".regulators.txt")] for p in NES_DIR.glob("aracne_v2_tf__*.regulators.txt"))
    print(f"{len(tags)} tags, {len(meta)} signatures, {len(LAMBDAS) + 1} arms")

    with Pool(THREADS) as pool:
        res = pool.map(read_one, [(t, tf_of_sig) for t in tags])
    r = pd.concat([x for x, _ in res], ignore_index=True)
    graph = pd.DataFrame([s for _, s in res])
    graph.to_csv(OUT / "graph.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"graph: regulators {graph.regulators.min()}-{graph.regulators.max()}, isolated "
          f"{graph.isolated.min()}-{graph.isolated.max()}, mean degree {graph.mean_degree.mean():.1f}")

    reads = combine(r, sg).merge(meta, on="sig_id")
    reads["correct"] = np.sign(reads.nes) == reads.direction
    reads["fpct"] = np.where(reads.direction < 0, reads.pct, 1 - reads.pct)
    reads.to_csv(OUT / "reads.tsv", sep="\t", index=False, float_format="%.5g")

    # identity: the lambda = inf arm must reproduce 79n's own per-set numbers and the published D84 summary
    ident = {}
    base = reads[reads.lam == INF]
    chk = sg.merge(base[["sig_id", "set", "nes", "pct"]], on=["sig_id", "set"], suffixes=("_79n", "_a1"))
    ident["n_rows_checked"] = int(len(chk))
    ident["max_abs_diff_nes"] = float((chk.nes_79n - chk.nes_a1).abs().max())
    ident["max_abs_diff_pct"] = float((chk.pct_79n - chk.pct_a1).abs().max())
    rt = rates(reads)
    for (st, kind), want in D84_PUBLISHED.items():
        got = float(rt[(rt.analysis == "main") & (rt.kind == kind) & (rt.read == st) & (rt.lam == INF)].correct.iloc[0])
        ident[f"published_{st}_{kind}"] = dict(published=want, recomputed=round(got, 4), diff=round(got - want, 4))
    ident["ok"] = (ident["max_abs_diff_nes"] < 1e-4 and ident["max_abs_diff_pct"] < 1e-4
                   and all(abs(v["diff"]) <= IDENTITY_TOL for k, v in ident.items() if k.startswith("published_")))
    print("identity:", json.dumps(ident, indent=1))
    if not ident["ok"]:
        (OUT / "run.json").write_text(json.dumps(dict(identity=ident, verdict="IDENTITY FAILED"), indent=1))
        sys.exit("lambda = inf arm does not reproduce production; stopping before any swept arm is read")

    rt.to_csv(OUT / "rates.tsv", sep="\t", index=False)
    per_tf = (reads[reads.own_ok & (reads.set == PRIMARY_READ)]
              .groupby(["lam", "cmap_name"], as_index=False)
              .agg(n=("sig_id", "size"), correct=("correct", "mean"), median_fpct=("fpct", "median")))
    per_tf.to_csv(OUT / "per_tf.tsv", sep="\t", index=False, float_format="%.4g")

    # the pre-registered decision
    p = rt[(rt.analysis == PRIMARY_ANALYSIS) & (rt.kind == "all") & (rt.read == PRIMARY_READ)].set_index("lam")
    base_rate = float(p.loc[INF, "correct"])
    gains = {str(l): round(100 * (float(p.loc[l, "correct"]) - base_rate), 2) for l in LAMBDAS}
    best = max(gains, key=gains.get)
    met = gains[best] >= MIN_GAIN_PTS
    run = dict(n_signatures=int(len(meta)), n_tags=len(tags), identity=ident,
               primary=dict(analysis=PRIMARY_ANALYSIS, read=PRIMARY_READ, kind="all",
                            n=int(p.loc[INF, "n"]), n_tfs=int(p.loc[INF, "n_tfs"])),
               base_rate_inf=round(base_rate, 4), rate_by_lambda={str(l): round(float(p.loc[l, "correct"]), 4) for l in LAMBDAS},
               gain_points=gains, best_lambda=float(best), best_gain_points=gains[best],
               min_gain_points=MIN_GAIN_PTS, criterion_met=bool(met),
               verdict=("lambda = %s adopted for A2" % best) if met else
                       "EARLY STOP: D125 recorded as not met; A2 does not run")
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 220)
    print("\nPrimary (%s, %s, all signatures):" % (PRIMARY_ANALYSIS, PRIMARY_READ))
    print(p[["n", "n_tfs", "correct", "auc", "median_fpct"]].to_string())
    print("\nGain over lambda = inf, percentage points:", gains)
    print("\n" + run["verdict"])
    print("\nAll arms:")
    print(rt[rt.analysis == "main"].pivot_table(index=["read", "kind"], columns="lam", values="correct").to_string())


if __name__ == "__main__":
    main()
