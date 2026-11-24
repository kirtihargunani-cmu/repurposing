#!/usr/bin/env python
"""Experiment E (config block E, amendment 1e): which disease programs can drugs reproducibly perturb?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227 + 6583b74. Programs = D198 modules rebuilt in the study namespace (E0_programs.py), D199 flags kept.
Per (cancer, window): window activity per drug (in-key weighted, as production), program opposition
O = -sum_{r in P n S} d_r a_r / sqrt(sum d_r^2) (n_used >= 3), split halves seeded as in Experiment A, Spearman-Brown
reliability per program with a size- and |NES|-tertile-matched pseudo-program null (200 sets), target-link test
(primary targets and their 1-hop OmniPath targets among members), amplitude control, flag comparison, line
selectivity. Explanatory only (D201): no ranking.
Outputs: results/lincs/mechanistic_incremental_v1/E_program_access/{E_program_response.tsv.gz, E_repeatability.tsv,
  E_target_link.tsv, E_summary.tsv, E_line_selectivity.tsv, run.json}
"""
import importlib.util
import json
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("A", HERE / "A_resolvable.py")
A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(A)
CFG, IN, RES = A.CFG, A.IN, A.RES
OUT = A.WT / CFG["output_root"] / "E_program_access"
PROG = OUT / "programs"
N_PSEUDO, MIN_DRUGS = 200, 10


def sb(r):
    return 2 * r / (1 + r) if np.isfinite(r) and r > -1 else np.nan


def opposition(Amat, d, idx_sets):
    """drugs x programs opposition; NaN where fewer than 3 members are scored."""
    out = np.full((Amat.shape[0], len(idx_sets)), np.nan)
    for k, ii in enumerate(idx_sets):
        if len(ii) >= 3:
            dn = np.sqrt((d[ii] ** 2).sum())
            if dn > 0:
                out[:, k] = -(Amat[:, ii] @ d[ii]) / dn
    return out


def work(task):
    unit, net, progs, links = task
    regs = pd.Index((A.VDIR / f"{net}.regulators.txt").read_text().split())
    Av = np.memmap(A.VDIR / f"{net}.f32", dtype=np.float32, mode="r", shape=(len(A.idx), len(regs)))
    mr = A.cached("mrE", lambda: pd.read_csv(RES / "step6/aracne_project_v2_tf/group_mr.tsv", sep="\t"))
    m = mr[(mr.group_id == unit) & mr.regulator.isin(regs)]
    S = regs.get_indexer(m.regulator); d = np.zeros(len(regs)); d[S] = m.nes.values
    lw = pd.read_csv(IN / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    wg = lw[(lw.group_id == unit) & lw.in_key].set_index("ModelID").weight
    sub = A.idx[A.idx.ModelID.isin(wg.index) & A.idx.window.isin(("early", "late"))].copy()
    sub["w"] = wg.reindex(sub.ModelID).fillna(0).values
    scored = set(m.regulator)
    psets = [np.array([regs.get_loc(r) for r in mem if r in scored], dtype=int) for mem in progs.members_list]
    tert = pd.qcut(np.abs(m.nes.values), 3, labels=False, duplicates="drop")
    tert_of = dict(zip(S, tert))
    pool = {t: np.array([s for s in S if tert_of[s] == t]) for t in set(tert)}
    rng0 = np.random.default_rng(A.rseed("Epseudo", unit))
    pseudo = []                                   # per program: N_PSEUDO matched index sets
    for ii in psets:
        comp = pd.Series([tert_of[i] for i in ii]).value_counts().to_dict() if len(ii) else {}
        pseudo.append([np.concatenate([rng0.choice(pool[t], size=min(c, len(pool[t])), replace=False) for t, c in comp.items()])
                       if comp else np.array([], int) for _ in range(N_PSEUDO)])
    resp, rel, lines_rows, link_rows = [], [], [], []
    for w in ("early", "late"):
        g = sub[sub.window == w]
        if not len(g):
            continue
        drugs = sorted(g.pert_id.unique())
        full, ha, hb, nmag = [], [], [], []
        split_ok = []
        for p in drugs:
            gp = g[g.pert_id == p]
            X = np.asarray(Av[gp.col.values], dtype=np.float64)
            full.append(A.wmean(X, gp.w.values)); nmag.append(np.linalg.norm(full[-1][S]))
            if len(gp) >= 2:
                perm = np.random.default_rng(A.rseed("group_id", unit, p, w)).permutation(len(gp))
                a, b = perm[: len(gp) // 2], perm[len(gp) // 2:]
                ha.append(A.wmean(X[a], gp.w.values[a])); hb.append(A.wmean(X[b], gp.w.values[b])); split_ok.append(True)
            else:
                ha.append(np.full(len(regs), np.nan)); hb.append(np.full(len(regs), np.nan)); split_ok.append(False)
        F, HA, HB = np.array(full), np.array(ha), np.array(hb)
        ok = np.array(split_ok)
        O = opposition(F, d, psets)
        Oa, Ob = opposition(HA[ok], d, psets), opposition(HB[ok], d, psets)
        for k, prow in enumerate(progs.itertuples()):
            cov = len(psets[k]) / max(1, len(prow.members_list))
            for i, p in enumerate(drugs):
                if np.isfinite(O[i, k]):
                    resp.append((unit, w, p, prow.program, cov, len(psets[k]), O[i, k], nmag[i]))
            r_obs = spearmanr(Oa[:, k], Ob[:, k]).correlation if ok.sum() >= MIN_DRUGS and np.isfinite(Oa[:, k]).all() else np.nan
            if not np.isfinite(r_obs):
                rel.append(dict(unit=unit, window=w, program=prow.program, n_drugs=int(ok.sum()), reliability=np.nan))
                continue
            rb = np.random.default_rng(A.rseed("Eboot", unit, w, prow.program))
            n = ok.sum(); bs = []
            for _ in range(A.NB):
                ii = rb.integers(n, size=n); bs.append(sb(spearmanr(Oa[ii, k], Ob[ii, k]).correlation))
            nulls = []
            for ps in pseudo[k]:
                if len(ps) >= 3:
                    oa = opposition(HA[ok], d, [ps])[:, 0]; ob = opposition(HB[ok], d, [ps])[:, 0]
                    nulls.append(sb(spearmanr(oa, ob).correlation))
            q95 = float(np.nanpercentile(nulls, 95)) if nulls else np.nan
            lo, hi = np.nanpercentile(bs, [2.5, 97.5])
            rel.append(dict(unit=unit, window=w, program=prow.program, n_drugs=int(n), n_members_scored=len(psets[k]),
                            reliability=sb(r_obs), ci_lo=float(lo), ci_hi=float(hi), pseudo_median=float(np.nanmedian(nulls)) if nulls else np.nan,
                            pseudo_q95=q95, perturbable=bool(lo > q95) if np.isfinite(q95) else False,
                            flag=prow.flag))
            # target link: drugs whose primary target (or its 1-hop targets) is a member
            mem = set(prow.members_list)
            lk = np.array([bool(links.get(p, set()) & mem) for p in drugs])
            col = O[:, k]; good = np.isfinite(col)
            if lk[good].sum() >= 3 and (~lk[good]).sum() >= 3:
                pct = pd.Series(np.abs(col[good])).rank(pct=True).to_numpy()
                mag = np.array(nmag)[good]
                res = np.abs(col[good]) - np.polyval(np.polyfit(mag, np.abs(col[good]), 1), mag)
                pr = pd.Series(res).rank(pct=True).to_numpy()
                link_rows.append(dict(unit=unit, window=w, program=prow.program, n_linked=int(lk[good].sum()),
                                      n_drugs=int(good.sum()), linked=lk[good], pct=pct, pct_resid=pr, flag=prow.flag))
        # line selectivity: share of opposition variance (after drug means) explained by line, programs pooled
        multi = g.groupby("pert_id").ModelID.nunique()
        gm = g[g.pert_id.isin(multi[multi >= 2].index)]
        if len(gm) >= 20:
            Xc = np.asarray(Av[gm.col.values], dtype=np.float64)
            Oc = opposition(Xc, d, psets)
            for k, prow in enumerate(progs.itertuples()):
                y = Oc[:, k]
                if not np.isfinite(y).all():
                    continue
                df = pd.DataFrame({"y": y, "drug": gm.pert_id.values, "line": gm.ModelID.values})
                df["r"] = df.y - df.groupby("drug").y.transform("mean")
                ss = (df.r ** 2).sum()
                ssl = (df.groupby("line").r.transform("mean") ** 2).sum()
                lines_rows.append(dict(unit=unit, window=w, program=prow.program, n_conditions=len(df),
                                       n_lines=df.line.nunique(), line_share=float(ssl / ss) if ss > 0 else np.nan, flag=prow.flag))
    R = pd.DataFrame(resp, columns=["unit", "window", "pert_id", "program", "coverage", "n_used", "opposition", "drug_magnitude"])
    return R, pd.DataFrame(rel), link_rows, pd.DataFrame(lines_rows)


def main():
    t0 = time.time()
    P = pd.read_csv(PROG / "programs.tsv", sep="\t")
    Fl = pd.read_csv(PROG / "program_filter.tsv", sep="\t")
    P = P.merge(Fl[["cohort", "program", "fail_plate", "fail_quality", "fail_proliferation", "fail_generic"]], on=["cohort", "program"], how="left")
    P["flag"] = np.select([P.fail_plate.fillna(False), P.fail_quality.fillna(False), P.fail_proliferation.fillna(False), P.fail_generic.fillna(False)],
                          ["plate", "quality", "proliferation", "generic"], "clean")
    P["members_list"] = P.members.str.split(",")
    # primary targets + 1-hop OmniPath, per drug
    sp = pd.read_csv(RES / "lincs/drug_reversal_split_v2_tf.tsv", sep="\t", low_memory=False,
                     usecols=["pert_id", "evidence_targets", "target_source"]).drop_duplicates("pert_id")
    om = pd.read_csv(IN / "data/raw/omnipath/omnipath_signed_directed.tsv.gz", sep="\t", usecols=["source_genesymbol", "target_genesymbol"])
    hop = om.groupby("source_genesymbol").target_genesymbol.apply(set).to_dict()
    links = {}
    for p, e, s in zip(sp.pert_id, sp.evidence_targets, sp.target_source):
        if s in ("chembl_mechanism", "confirmed_annotation") and isinstance(e, str):
            t = {x.split(":")[0].strip() for x in e.replace(";", "|").split("|") if x.strip()}
            links[p] = t | set().union(*[hop.get(g, set()) for g in t])
    ovr = pd.read_csv(IN / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = pd.read_csv(RES / "step1/disease_groups_toil.tsv", sep="\t").drop_duplicates("group_id").set_index("group_id").cohort.to_dict()
    avail = {p.stem for p in A.VDIR.glob("*.f32")}
    win_units = set(pd.read_csv(RES / "systems_validation/integrated/window_reversal.tsv.gz", sep="\t").cohort)
    tasks = [(u, A.net_for(u, ovr, avail, cohort_of), g.reset_index(drop=True), links)
             for u, g in P.groupby("cohort") if u in win_units]
    print(f"{len(tasks)} cancers with programs and window data", flush=True)
    Rs, Ls, Ks, Ss = [], [], [], []
    with ProcessPoolExecutor(A.NC) as ex:
        for i, (r, l, k, s) in enumerate(ex.map(work, tasks), 1):
            Rs.append(r); Ls.append(l); Ks.extend(k); Ss.append(s)
            print(f"  [{i}/{len(tasks)}] {time.time() - t0:.0f}s", flush=True)
    R = pd.concat(Rs, ignore_index=True); REL = pd.concat(Ls, ignore_index=True); LS = pd.concat(Ss, ignore_index=True)
    R.to_csv(OUT / "E_program_response.tsv.gz", sep="\t", index=False, float_format="%.5g")
    REL.to_csv(OUT / "E_repeatability.tsv", sep="\t", index=False, float_format="%.4f")
    LS.to_csv(OUT / "E_line_selectivity.tsv", sep="\t", index=False, float_format="%.4f")
    # pooled target-link test per window (and per flag), permutation of links within program
    T = []
    for w in ("early", "late"):
        for flag in ("all", "clean"):
            ks = [k for k in Ks if k["window"] == w and (flag == "all" or k["flag"] == "clean")]
            if len(ks) < 5:
                continue
            for col in ("pct", "pct_resid"):
                diffs = np.array([k[col][k["linked"]].mean() - k[col][~k["linked"]].mean() for k in ks])
                obs = float(diffs.mean())
                rng = np.random.default_rng(A.rseed("Elink", w, flag, col)); null = []
                for _ in range(A.NP):
                    s = 0.0
                    for k in ks:
                        lk = rng.permutation(k["linked"]); s += k[col][lk].mean() - k[col][~lk].mean()
                    null.append(s / len(ks))
                lo, hi = A.boot_ci(diffs, np.arange(len(diffs)), stat=np.mean)
                nm = float(np.median(null))
                T.append(dict(window=w, programs=flag, measure="|O| percentile" if col == "pct" else "|O| residual on drug magnitude",
                              n_programs=len(ks), n_cancers=len({k["unit"] for k in ks}), linked_minus_unlinked=obs, null=nm,
                              ci_lo=lo, ci_hi=hi, p_perm=float((np.sum(np.array(null) >= obs) + 1) / (A.NP + 1)),
                              resolved=bool(lo > nm or hi < nm)))
    TL = pd.DataFrame(T); TL.to_csv(OUT / "E_target_link.tsv", sep="\t", index=False, float_format="%.4f")
    rv = REL.dropna(subset=["reliability"])
    SUM = rv.groupby(["window", "flag"]).agg(programs=("program", "size"), cancers=("unit", "nunique"),
                                             median_reliability=("reliability", "median"),
                                             median_pseudo_q95=("pseudo_q95", "median"),
                                             perturbable=("perturbable", "sum")).reset_index()
    SUM["perturbable_share"] = SUM.perturbable / SUM.programs
    SUM.to_csv(OUT / "E_summary.tsv", sep="\t", index=False, float_format="%.4f")
    run = dict(stage="E_program_access", protocol="5cc2227+6583b74", cancers=len(tasks), programs=int(P.shape[0]),
               response_rows=len(R), reliability_rows=len(rv), link_programs=len(Ks), runtime_s=round(time.time() - t0, 1),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 250)
    print(SUM.to_string(index=False)); print(TL.to_string(index=False))
    print(LS.groupby(["window", "flag"]).line_share.median().round(3).to_string())


if __name__ == "__main__":
    sys.exit(main())
