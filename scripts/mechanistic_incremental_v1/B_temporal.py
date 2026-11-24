#!/usr/bin/env python
"""Experiment B (A1.6 + config block B): early vs late transitions under an empirical noise model.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227 + amendment 1b (2e55c43). Observed windows only; nothing is imputed.

Per (unit): condition-level regulator activity (N.A), rev_c = -(nes . a_c)/||s||, and the per-regulator gene-level
terms of rev_gene, q_c = -sum_{j in r} d_rj C_cj / ||s|| (rev_gene_c = sum_r q_c), for every condition of a drug that
has both windows or >= 2 conditions in a window.
Noise: pooled within-drug variance across single conditions, per window (scalar for rev / rev_gene, per regulator for
activity); SE of a weighted estimate = sigma * sqrt(sum w^2) / sum w.
Arms: matched_dose (conditions whose (line, dose) the drug has in both windows), production (all in-key conditions),
high_confidence (>= 2 conditions in both windows).
Calls: |estimate / SE| > 1.96 -> positive / negative, else uncertain; transition classes as in the config.
Pseudo-replicate null: one window's conditions split in halves, labelled early / late, same calls.
Outputs: results/lincs/mechanistic_incremental_v1/B_temporal/{B_transitions.tsv.gz, B_regulator_summary.tsv.gz,
  B_noise.tsv, B_pseudoreplicate.tsv.gz, B_summary.tsv, B_enrichment.tsv, B_case_studies.tsv, run.json}
"""
import importlib.util
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("A", HERE / "A_resolvable.py")
A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(A)
CFG, IN, RES, PROC = A.CFG, A.IN, A.RES, A.PROC
BC = CFG["B"]
OUT = A.WT / CFG["output_root"] / "B_temporal"
OUT.mkdir(parents=True, exist_ok=True)
Z = 1.96
CLASSES = BC["transition_classes"]


def call(est, se):
    z = np.divide(est, se, out=np.zeros_like(np.asarray(est, float)), where=np.asarray(se) > 0)
    return np.where(z > Z, 1, np.where(z < -Z, -1, 0))


def transition(ce, cl):
    ce, cl = np.asarray(ce), np.asarray(cl)
    return np.select([(ce * cl) == -1, (ce * cl) == 1, (ce != 0) & (cl == 0), (ce == 0) & (cl != 0)],
                     ["sign_reversal", "persistent_same_sign", "early_only", "late_only"], "both_uncertain")


def setup(level, unit, net):
    P = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype=np.float32, mode="r", shape=(len(A.idx), A.NG))
    edges = A.cached("edges", lambda: pd.read_csv(IN / "data/processed/step5/regulons_aracne_project_v2_tf_edges.tsv.gz", sep="\t"))
    N = A.Net(net, edges)
    mrd = "aracne_project_v2_tf" + ("_strata" if level == "stratum" else "")
    mr = A.cached(("mr", mrd), lambda: pd.read_csv(RES / f"step6/{mrd}/group_mr.tsv", sep="\t"))
    m = mr[(mr.group_id == unit) & mr.regulator.isin(N.regs)]
    ci = N.regs.get_indexer(m.regulator)
    nes = np.zeros(len(N.regs)); nes[ci] = m.nes.values
    tested = np.zeros(len(N.regs), bool); tested[ci] = True
    norm = np.linalg.norm(m.nes.values)
    dzall = A.cached(("dz", mrd), lambda: pd.read_csv(RES / f"step6/{mrd}/regulon_contributions.tsv.gz", sep="\t",
                                                       usecols=["group_id", "regulator", "target", "contrib_nes"]))
    dz = dzall[dzall.group_id == unit].drop_duplicates(["regulator", "target"]).set_index(["regulator", "target"]).contrib_nes
    dvec = np.nan_to_num(dz.reindex(pd.MultiIndex.from_arrays([N.regs[N.reg], N.target])).to_numpy(dtype=np.float64))
    if level == "stratum":
        lw = pd.read_csv(IN / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
        lw["in_key"] = lw.in_key.astype(str).str.lower().eq("true")
    else:
        lw = pd.read_csv(IN / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    wg = lw[(lw.group_id == unit) & lw.in_key].set_index("ModelID").weight
    sub = A.idx[A.idx.ModelID.isin(wg.index) & A.idx.window.isin(("early", "late"))].copy()
    sub["w"] = wg.reindex(sub.ModelID).fillna(0).values
    return P, N, nes, tested, norm, dvec, sub


def est(vals, w):
    """weighted mean and the factor sqrt(sum w^2)/sum w (multiplies sigma to give the SE)."""
    w = np.asarray(w, float); w = w if w.sum() > 0 else np.ones(len(w))
    return (w / w.sum()) @ vals, float(np.sqrt((w ** 2).sum()) / w.sum())


def work(task):
    level, unit, net, cases = task
    P, N, nes, tested, norm, dvec, sub = setup(level, unit, net)
    nt = nes[tested]
    cnt = sub.groupby(["pert_id", "window"]).size().unstack(fill_value=0)
    both = set(cnt.index[(cnt.get("early", 0) > 0) & (cnt.get("late", 0) > 0)])
    multi = {w: set(cnt.index[cnt.get(w, 0) >= 2]) for w in ("early", "late")}
    need = sub[sub.pert_id.isin(both | multi["early"] | multi["late"])]
    cols = need.col.values
    act = np.asarray(N.A[cols], dtype=np.float64)[:, tested]
    Q = np.zeros((len(cols), int(tested.sum())))
    for s in range(0, len(cols), 100):                     # per-condition gene-level regulator terms
        C = N.contrib(P, cols[s:s + 100])
        Q[s:s + 100] = (-np.add.reduceat(dvec * C, N.starts, axis=1) / norm)[:, tested]
    rev = -(act @ nt) / norm
    revg = Q.sum(1)
    pos = pd.Series(np.arange(len(cols)), index=need.index)
    # noise: pooled within-drug variance across single conditions, per window
    noise = {}
    for w in ("early", "late"):
        ss_s, ss_g, ss_a, dof = 0.0, 0.0, np.zeros(act.shape[1]), 0
        for p in multi[w]:
            r = pos[need[(need.pert_id == p) & (need.window == w)].index].values
            ss_s += ((rev[r] - rev[r].mean()) ** 2).sum(); ss_g += ((revg[r] - revg[r].mean()) ** 2).sum()
            ss_a += ((act[r] - act[r].mean(0)) ** 2).sum(0); dof += len(r) - 1
        noise[w] = dict(sd_rev=np.sqrt(ss_s / dof) if dof else np.nan, sd_revg=np.sqrt(ss_g / dof) if dof else np.nan,
                        sd_act=np.sqrt(ss_a / dof) if dof else np.full(act.shape[1], np.nan), dof=dof)
    NO = [dict(level=level, unit=unit, window=w, dof=v["dof"], sd_rev=v["sd_rev"], sd_rev_gene=v["sd_revg"],
               sd_act_median=float(np.nanmedian(v["sd_act"]))) for w, v in noise.items()]

    def block(p, w, mask=None):
        g = need[(need.pert_id == p) & (need.window == w)]
        if mask is not None:
            g = g[mask(g)]
        if not len(g):
            return None
        r = pos[g.index].values
        f = est(rev[r], g.w.values)[1]
        return dict(rev=est(rev[r], g.w.values)[0], revg=est(revg[r], g.w.values)[0], act=est(act[r], g.w.values)[0],
                    q=est(Q[r], g.w.values)[0], f=f, n=len(g))

    T, RS, CS = [], [], []
    regs = N.regs[tested]
    for p in sorted(both):
        de = need[(need.pert_id == p) & (need.window == "early")]; dl = need[(need.pert_id == p) & (need.window == "late")]
        shared = set(zip(de.ModelID, de.pert_dose)) & set(zip(dl.ModelID, dl.pert_dose))
        arms = {"production": None,
                "matched_dose": (lambda g, s=shared: pd.Series([(a, b) in s for a, b in zip(g.ModelID, g.pert_dose)], index=g.index)) if shared else False,
                "high_confidence": None if (cnt.loc[p, "early"] >= 2 and cnt.loc[p, "late"] >= 2) else False}
        for arm, mask in arms.items():
            if mask is False:
                continue
            e, l = block(p, "early", mask), block(p, "late", mask)
            if e is None or l is None:
                continue
            se = lambda b, w, k: noise[w][k] * b["f"]
            row = dict(level=level, unit=unit, pert_id=p, arm=arm, n_early=e["n"], n_late=l["n"])
            for sc, k in (("rev_gene", "sd_revg"), ("rev", "sd_rev")):
                ve, vl = e["revg" if sc == "rev_gene" else "rev"], l["revg" if sc == "rev_gene" else "rev"]
                ce, cl = call(np.array([ve]), np.array([se(e, "early", k)]))[0], call(np.array([vl]), np.array([se(l, "late", k)]))[0]
                row.update({f"{sc}_early": ve, f"{sc}_late": vl, f"{sc}_se_early": se(e, "early", k),
                            f"{sc}_se_late": se(l, "late", k), f"{sc}_class": str(transition(ce, cl))})
            # regulator level (activity)
            ca = call(e["act"], noise["early"]["sd_act"] * e["f"]); cb = call(l["act"], noise["late"]["sd_act"] * l["f"])
            tr = transition(ca, cb)
            for c in CLASSES:
                row[f"reg_share_{c}"] = float(np.mean(tr == c))
            d = l["q"] - e["q"]; ad = np.abs(d)
            row["dq_top5_share"] = float(np.sort(ad)[::-1][:5].sum() / ad.sum()) if ad.sum() > 0 else np.nan
            row["dq_neff"] = float(ad.sum() ** 2 / (d ** 2).sum()) if (d ** 2).sum() > 0 else np.nan
            te, tl = set(np.argsort(-np.abs(e["q"]))[:10]), set(np.argsort(-np.abs(l["q"]))[:10])
            row["top10_jaccard"] = len(te & tl) / len(te | tl)
            sh = list(te & tl)
            row["top10_shared_flip"] = float(np.mean(e["q"][sh] * l["q"][sh] < 0)) if sh else np.nan
            lm = np.abs(l["q"]); row["late_mass_outside_early_top10"] = float(lm[[i for i in range(len(lm)) if i not in te]].sum() / lm.sum()) if lm.sum() > 0 else np.nan
            T.append(row)
            if arm == "production" and p in cases:
                top = np.argsort(-ad)[:20]
                for i in top:
                    CS.append(dict(level=level, unit=unit, pert_id=p, regulator=regs[i], disease_nes=nt[i],
                                   act_early=e["act"][i], act_late=l["act"][i], q_early=e["q"][i], q_late=l["q"][i], dq=d[i]))
    # pseudo-replicate null: split one window into halves labelled early / late
    PS = []
    for w in ("early", "late"):
        for p in sorted(multi[w]):
            g = need[(need.pert_id == p) & (need.window == w)]
            rng = np.random.default_rng(A.rseed("pseudo", level, unit, p, w))
            perm = rng.permutation(len(g)); ha, hb = g.iloc[perm[: len(g) // 2]], g.iloc[perm[len(g) // 2:]]
            ra, rb = pos[ha.index].values, pos[hb.index].values
            row = dict(level=level, unit=unit, pert_id=p, window=w, n=len(g))
            for sc, vals, k in (("rev_gene", revg, "sd_revg"), ("rev", rev, "sd_rev")):
                (va, fa), (vb, fb) = est(vals[ra], ha.w.values), est(vals[rb], hb.w.values)
                sa, sb = noise[w][k] * fa, noise[w][k] * fb
                row[f"{sc}_class"] = str(transition(call(np.array([va]), np.array([sa]))[0], call(np.array([vb]), np.array([sb]))[0]))
                row[f"{sc}_z_halfdiff"] = (va - vb) / np.sqrt(sa ** 2 + sb ** 2) if sa > 0 else np.nan
            (aa, fa), (ab, fb) = est(act[ra], ha.w.values), est(act[rb], hb.w.values)
            sd = noise[w]["sd_act"]
            tr = transition(call(aa, sd * fa), call(ab, sd * fb))
            row["reg_share_sign_reversal"] = float(np.mean(tr == "sign_reversal"))
            zz = (aa - ab) / (sd * np.sqrt(fa ** 2 + fb ** 2))
            row["act_share_within_1.96"] = float(np.mean(np.abs(zz[np.isfinite(zz)]) < Z))
            PS.append(row)
    return pd.DataFrame(T), pd.DataFrame(PS), pd.DataFrame(NO), pd.DataFrame(CS)


def main():
    t0 = time.time()
    ovr = pd.read_csv(IN / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = pd.read_csv(RES / "step1/disease_groups_toil.tsv", sep="\t").drop_duplicates("group_id").set_index("group_id").cohort.to_dict()
    parent_of = pd.read_csv(RES / "celligner/stratum_line_weights_long_v2.tsv", sep="\t").drop_duplicates("stratum").set_index("stratum").parent.to_dict()
    avail = {p.stem for p in A.VDIR.glob("*.f32")}
    ann = pd.read_csv(RES / "lincs/drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    name2pert = ann.reset_index().groupby("cmap_name").pert_id.apply(set).to_dict()
    tasks = []
    for level, f in (("group_id", "window_reversal.tsv.gz"), ("stratum", "window_reversal_strata.tsv.gz")):
        w = pd.read_csv(RES / "systems_validation/integrated" / f, sep="\t")
        for unit in sorted(w[w.window_coverage == "both"].cohort.unique()):
            cs = set().union(*[name2pert.get(n, set()) for n in BC["cases"].get(unit, [])])
            tasks.append((level, unit, A.net_for(parent_of.get(unit, unit), ovr, avail, cohort_of), cs))
    only = [x for x in os.environ.get("ONLY_UNITS", "").split(",") if x]   # smoke tests only; never in a reported run
    if only:
        tasks = [t for t in tasks if t[1] in only]
        global OUT
        OUT = OUT.parent / "B_smoke"; OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(tasks)} units with both-window drugs", flush=True)
    Ts, Ps, Ns, Cs = [], [], [], []
    with ProcessPoolExecutor(A.NC) as ex:
        for i, (t, p, n, c) in enumerate(ex.map(work, tasks), 1):
            Ts.append(t); Ps.append(p); Ns.append(n); Cs.append(c)
            print(f"  [{i}/{len(tasks)}] {time.time() - t0:.0f}s", flush=True)
    T, PS, NO, CS = (pd.concat(x, ignore_index=True) for x in (Ts, Ps, Ns, Cs))
    fm = pd.read_csv(RES / "lincs/fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    T["cluster"] = T.pert_id.map(lambda p: fm.get(p, ann.cmap_name.get(p, p)))
    T["moa"] = T.pert_id.map(ann.hub_moa)
    reg_cols = [c for c in T.columns if c.startswith("reg_share_")]
    T.to_csv(OUT / "B_transitions.tsv.gz", sep="\t", index=False, float_format="%.6g")
    PS.to_csv(OUT / "B_pseudoreplicate.tsv.gz", sep="\t", index=False, float_format="%.6g")
    NO.to_csv(OUT / "B_noise.tsv", sep="\t", index=False, float_format="%.6g")
    CS.to_csv(OUT / "B_case_studies.tsv", sep="\t", index=False, float_format="%.5g")

    # summary: transition-class shares, observed vs pseudo-replicate, with cluster bootstrap intervals
    S = []
    for (lv, arm), d in T.groupby(["level", "arm"]):
        ps = PS[PS.level == lv]
        for sc in ("rev_gene", "rev"):
            for c in CLASSES:
                v = (d[f"{sc}_class"] == c).astype(float).values
                lo, hi = A.boot_ci(v, d.cluster.values, stat=np.mean)
                nullv = float((ps[f"{sc}_class"] == c).mean()) if len(ps) else np.nan
                S.append(dict(level=lv, arm=arm, score=sc, cls=c, n_rows=len(d), n_clusters=d.cluster.nunique(),
                              share=float(v.mean()), ci_lo=lo, ci_hi=hi, pseudo_replicate_share=nullv,
                              resolved_vs_pseudo=bool(lo > nullv) if c == "sign_reversal" and np.isfinite(nullv) else None))
        for c in reg_cols:
            lo, hi = A.boot_ci(d[c].values, d.cluster.values)
            S.append(dict(level=lv, arm=arm, score="regulator activity", cls=c.replace("reg_share_", ""), n_rows=len(d),
                          n_clusters=d.cluster.nunique(), share=float(d[c].median()), ci_lo=lo, ci_hi=hi,
                          pseudo_replicate_share=float(ps.reg_share_sign_reversal.median()) if c.endswith("sign_reversal") and len(ps) else np.nan))
    SUM = pd.DataFrame(S); SUM.to_csv(OUT / "B_summary.tsv", sep="\t", index=False, float_format="%.4f")
    CAL = PS.groupby("level").agg(rev_gene_within=("rev_gene_z_halfdiff", lambda z: float(np.mean(np.abs(z.dropna()) < Z))),
                                  rev_within=("rev_z_halfdiff", lambda z: float(np.mean(np.abs(z.dropna()) < Z))),
                                  act_within=("act_share_within_1.96", "median"), n=("pert_id", "size")).reset_index()

    # mechanism-class enrichment (cancers; matched_dose primary, production secondary)
    E = []
    base = T[(T.level == "group_id") & (T.arm == "production")]
    classes = [m for m, k in base.groupby("moa").cluster.nunique().items() if k >= 8]
    for arm in ("matched_dose", "production"):
        d = T[(T.level == "group_id") & (T.arm == arm)].copy()
        cl_moa = d.drop_duplicates("cluster").set_index("cluster").moa
        for m in classes:
            for c in ("sign_reversal", "persistent_same_sign", "early_only", "late_only"):
                x = (d.moa == m).values; y = (d.rev_gene_class == c).values
                def orr(x, y):
                    # plain odds ratio; not estimable (nan) when any cell is zero. (A +0.5 correction was
                    # removed after the smoke test: with rare outcomes it produced OR ~16 from zero hits.)
                    a, b = (x & y).sum(), (x & ~y).sum(); cc, dd = (~x & y).sum(), (~x & ~y).sum()
                    return (a * dd) / (b * cc) if min(a, b, cc, dd) > 0 else np.nan
                o = orr(x, y)
                rng = np.random.default_rng(A.rseed("enrich", arm, m, c))
                null = []
                for _ in range(A.NP):
                    sh = dict(zip(cl_moa.index, rng.permutation(cl_moa.values)))
                    null.append(orr((d.cluster.map(sh) == m).values, y))
                nl = np.array(null, float); nl = nl[np.isfinite(nl)]
                p = float((np.sum(nl >= o) + 1) / (len(nl) + 1)) if np.isfinite(o) else np.nan
                # cluster bootstrap of the OR
                groups = [np.flatnonzero(d.cluster.values == k) for k in d.cluster.unique()]
                rb = np.random.default_rng(A.rseed("enrichboot", arm, m, c)); bs = []
                for _ in range(A.NB):
                    ii = np.concatenate([groups[k] for k in rb.integers(len(groups), size=len(groups))])
                    bs.append(orr(x[ii], y[ii]))
                bs = np.array(bs, float); bs = bs[np.isfinite(bs)]
                ok = np.isfinite(o) and len(bs) >= 0.9 * A.NB
                E.append(dict(arm=arm, moa=m, transition=c, n_class_rows=int(x.sum()), n_class_drugs=int(d[x].cluster.nunique()),
                              n_hits=int((x & y).sum()), odds_ratio=o, ci_lo=float(np.percentile(bs, 2.5)) if ok else np.nan,
                              ci_hi=float(np.percentile(bs, 97.5)) if ok else np.nan, p_perm=p, estimable=bool(ok)))
    EN = pd.DataFrame(E)
    if len(EN):
        for arm, g in EN.groupby("arm"):
            g = g[g.p_perm.notna()]
            if not len(g):
                continue
            pv = g.p_perm.values; o = np.argsort(pv); q = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
            q = np.minimum.accumulate(q[::-1])[::-1]; qq = np.empty_like(q); qq[o] = np.minimum(q, 1)
            EN.loc[g.index, "fdr"] = qq
        EN["resolved"] = EN.estimable & ((EN.ci_lo > 1) | (EN.ci_hi < 1))
        EN["benchmark_or15_fdr05"] = (EN.odds_ratio >= 1.5) & (EN.fdr < 0.05)
    EN.to_csv(OUT / "B_enrichment.tsv", sep="\t", index=False, float_format="%.4f")
    run = dict(stage="B_temporal", protocol="5cc2227+2e55c43", units=len(tasks), rows=len(T), pseudo_rows=len(PS),
               classes=classes, calibration=CAL.to_dict("records"), runtime_s=round(time.time() - t0, 1),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1, default=float))
    pd.set_option("display.width", 250)
    print(CAL.to_string(index=False)); print(SUM[SUM.score != "regulator activity"].to_string(index=False))
    print(EN.sort_values("p_perm").head(25).to_string(index=False) if len(EN) else "no enrichment classes")


if __name__ == "__main__":
    sys.exit(main())
