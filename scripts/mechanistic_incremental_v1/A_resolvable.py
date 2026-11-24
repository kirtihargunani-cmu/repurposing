#!/usr/bin/env python
"""Experiment A (Amendment 1, A1.5): how much drug-mechanism signal each processing step resolves, per window.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227; config/mechanistic_incremental_v1.yaml. Reads the 2026-11-23 run read-only from `input_root`.

PART 1 - reproducibility (rows with >= 2 in-key conditions; conditions permuted with crc32(key) ^ seed; half A =
first floor(n/2)). Per row, on each half: S1 gene profile (10,112 genes), S3 window regulator activity, S4 rev,
S5 rev_gene; S2 = two single conditions' activity. Vector steps: same-drug cosine (half A vs half B) against a
matched different-drug cosine (same cell, same condition-count stratum, seeded choice). Scalar steps: per cell with
>= 10 drugs, Spearman(half A, half B) across drugs, Spearman-Brown reliability; S6 = S5 within the eligible set.

PART 2 - mechanism labels on full-window values (the production numbers):
 (a) primary target percentile in the drug's own regulator response at S3 (signed activity), S4 (rev_contrib),
     S5 (rev_gene_contrib), and S5 within the eligible set (S6); null = target labels permuted among drugs of the cell
 (b) PROGENy on-target pathway (D99 labels) on the window gene profile (S1); null = labels permuted within cell
 (c) drugs sharing a target unit vs matched non-sharing pairs: cosine at S1, S3 and S5 (rev_gene_contrib vector);
     null = target-unit labels permuted within cell
Uncertainty: cluster bootstrap over FDA ingredient (else pert_id) or target unit, n_boot from config.
Resolved = the 95% interval excludes the matched null (A1.7).
Outputs (results/lincs/mechanistic_incremental_v1/A_resolvable/): A_halves.tsv.gz, A_vector_pairs.tsv.gz,
  A_reliability_cells.tsv, A_progeny.tsv.gz, A_target_rank.tsv.gz, A_shared_target_pairs.tsv.gz, A_curve.tsv,
  A_labels.tsv, run.json
"""
import json
import os
import re
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
import yaml
from scipy.special import ndtri
from scipy.stats import rankdata, spearmanr

WT = Path(__file__).resolve().parents[2]
CFG = yaml.safe_load((WT / "config/mechanistic_incremental_v1.yaml").read_text())
IN = Path(CFG["input_root"])
OUT = WT / CFG["output_root"] / "A_resolvable"
OUT.mkdir(parents=True, exist_ok=True)
SEED, NB, NP, NC = int(CFG["seed"]), int(CFG["n_boot"]), int(CFG["n_perm"]), int(CFG["n_cores"])
MIN_DRUGS = int(CFG["reproducibility"]["min_drugs_per_cell"])
PROC, RES = IN / "data/processed/lincs", IN / "results"
VDIR = PROC / "drug_viper_specific_v2_tf_cond"


def rseed(*key):
    return (zlib.crc32("|".join(map(str, key)).encode()) ^ SEED) & 0xFFFFFFFF


def stratum(n):
    return "2" if n == 2 else ("3-4" if n <= 4 else "5+")


genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
NG = len(genes)
idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
idx["col"] = np.arange(len(idx))
_t = idx.pert_time.where(idx.pert_time >= 0)
idx["window"] = np.where(_t.isna(), "unknown", np.where(_t <= 6, "early", np.where(_t >= 24, "late", "other")))
gpos = pd.Series(np.arange(NG), index=genes.gene_id)
_CACHE = {}


def cached(key, loader):
    """Per-process cache: each worker reads the large shared inputs once."""
    if key not in _CACHE:
        _CACHE[key] = loader()
    return _CACHE[key]


class Net:
    """Exact copy of 79e3's Net (aREA per-edge contributions; sum over a regulator's edges = its activity)."""

    def __init__(self, net, edges):
        self.regs = pd.Index((VDIR / f"{net}.regulators.txt").read_text().split())
        e = edges[edges.net == net]
        e = e[e.regulator.isin(set(self.regs)) & e.target.isin(gpos.index)]
        e = e.assign(_k=self.regs.get_indexer(e.regulator)).sort_values("_k", kind="stable")
        self.U = np.sort(gpos[e.target.unique()].values)
        upos = pd.Series(np.arange(len(self.U)), index=gpos.index[self.U])
        self.reg = e._k.values
        self.gene = upos[e.target].values
        self.target = e.target.values
        self.m, self.l = e.tfmode.values, e.likelihood.values
        self.starts = np.r_[0, np.flatnonzero(np.diff(self.reg)) + 1]
        L = np.add.reduceat(self.l, self.starts)
        lmax = np.maximum.reduceat(self.l, self.starts)
        self.w = np.sqrt(np.add.reduceat((self.l / lmax[self.reg]) ** 2, self.starts))
        self.k = (self.w / L)[self.reg]
        self.A = np.memmap(VDIR / f"{net}.f32", dtype=np.float32, mode="r", shape=(len(idx), len(self.regs)))

    def contrib(self, P, cols):
        X = np.asarray(P[cols][:, self.U], dtype=np.float64)
        r = rankdata(X, axis=1) / (X.shape[1] + 1)
        t1 = np.abs(r - 0.5) * 2
        t1 = ndtri(t1 + (1 - t1.max(1, keepdims=True)) / 2)
        t2 = ndtri(r)
        a = self.m * self.l * t2[:, self.gene]
        b = (1 - np.abs(self.m)) * self.l * t1[:, self.gene]
        S1 = np.add.reduceat(a, self.starts, axis=1)
        S2 = np.add.reduceat(b, self.starts, axis=1)
        s = np.where(S1 == 0, 1.0, np.sign(S1))
        return (a + (s * (S2 > 0))[:, self.reg] * b) * self.k


def cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else np.nan


def wmean(M, w):
    w = w if w.sum() > 0 else np.ones(len(w))
    return (w / w.sum()) @ M


def work(task):
    """One (level, unit): halves, vector pairs, full-window gene profiles -> PROGENy and shared-target pairs (S1, S3)."""
    level, unit, net, rows, prog_w, prog_genes, share = task
    P = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype=np.float32, mode="r", shape=(len(idx), NG))
    edges = cached("edges", lambda: pd.read_csv(IN / "data/processed/step5/regulons_aracne_project_v2_tf_edges.tsv.gz", sep="\t"))
    N = Net(net, edges)
    mrd = "aracne_project_v2_tf" + ("_strata" if level == "stratum" else "")
    mr = cached(("mr", mrd), lambda: pd.read_csv(RES / f"step6/{mrd}/group_mr.tsv", sep="\t"))
    m = mr[(mr.group_id == unit) & mr.regulator.isin(N.regs)]
    ci = N.regs.get_indexer(m.regulator)
    nes = np.zeros(len(N.regs)); nes[ci] = m.nes.values
    tested = np.zeros(len(N.regs), bool); tested[ci] = True
    norm = np.linalg.norm(m.nes.values)
    dzall = cached(("dz", mrd), lambda: pd.read_csv(RES / f"step6/{mrd}/regulon_contributions.tsv.gz", sep="\t",
                                                     usecols=["group_id", "regulator", "target", "contrib_nes"]))
    dz = dzall[dzall.group_id == unit].drop_duplicates(["regulator", "target"]).set_index(["regulator", "target"]).contrib_nes
    dvec = np.nan_to_num(dz.reindex(pd.MultiIndex.from_arrays([N.regs[N.reg], N.target])).to_numpy(dtype=np.float64))
    if level == "stratum":
        lw = pd.read_csv(IN / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
        lw["in_key"] = lw.in_key.astype(str).str.lower().eq("true")
    else:
        lw = pd.read_csv(IN / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    wg = lw[(lw.group_id == unit) & lw.in_key].set_index("ModelID").weight
    sub = idx[idx.ModelID.isin(wg.index) & idx.window.isin(("early", "late"))]
    pg = gpos.reindex(prog_genes).values
    H, FULL = [], {}
    for pert, win in rows:
        g = sub[(sub.pert_id == pert) & (sub.window == win)]
        if not len(g):
            continue
        ww = wg.reindex(g.ModelID).fillna(0).values
        prof = wmean(np.asarray(P[g.col.values], dtype=np.float64), ww)
        act = wmean(np.asarray(N.A[g.col.values], dtype=np.float64), ww)[tested]
        FULL[(pert, win)] = (prof, act, len(g))
        if len(g) < 2:
            continue
        rng = np.random.default_rng(rseed(level, unit, pert, win))
        perm = rng.permutation(len(g)); ha, hb = perm[: len(g) // 2], perm[len(g) // 2:]
        cols = g.col.values
        out = dict(level=level, unit=unit, pert_id=pert, window=win, n_conditions=len(g), stratum=stratum(len(g)))
        half = {}
        for tag, h in (("a", ha), ("b", hb)):
            pa = wmean(np.asarray(P[cols[h]], dtype=np.float64), ww[h])
            aa = wmean(np.asarray(N.A[cols[h]], dtype=np.float64), ww[h])
            C = wmean(N.contrib(P, cols[h]), ww[h])
            half[tag] = (pa, aa[tested])
            out[f"rev_{tag}"] = float(-(nes * aa)[tested].sum() / norm)
            out[f"rev_gene_{tag}"] = float(-(dvec * C).sum() / norm)
        out["S1_cos_same"] = cos(half["a"][0], half["b"][0])
        out["S3_cos_same"] = cos(half["a"][1], half["b"][1])
        c0, c1 = cols[perm[0]], cols[perm[1]]
        out["S2_cos_same"] = cos(np.asarray(N.A[c0], dtype=np.float64)[tested], np.asarray(N.A[c1], dtype=np.float64)[tested])
        H.append((out, half, (c0, c1)))
    # matched different-drug null: same cell (window) and condition-count stratum, seeded choice of partner
    rowsH, pairs = [], []
    by = {}
    for o, half, cc in H:
        by.setdefault((o["window"], o["stratum"]), []).append((o, half, cc))
    for key, lst in by.items():
        for i, (o, half, cc) in enumerate(lst):
            others = [j for j in range(len(lst)) if lst[j][0]["pert_id"] != o["pert_id"]]
            if others:
                j = others[np.random.default_rng(rseed("null", level, unit, o["pert_id"], o["window"])).integers(len(others))]
                o2, half2, cc2 = lst[j]
                o["null_partner"] = o2["pert_id"]
                o["S1_cos_null"] = cos(half["a"][0], half2["b"][0])
                o["S3_cos_null"] = cos(half["a"][1], half2["b"][1])
                o["S2_cos_null"] = cos(np.asarray(N.A[cc[0]], dtype=np.float64)[tested], np.asarray(N.A[cc2[1]], dtype=np.float64)[tested])
            rowsH.append(o)
    # PROGENy on full-window profiles (S1)
    PR = []
    for (pert, win), (prof, act, n) in FULL.items():
        z = prof[pg[~np.isnan(pg)].astype(int)]
        sc = {pw: float(z @ wv) for pw, wv in prog_w.items()}
        PR.append(dict(level=level, unit=unit, pert_id=pert, window=win, n_conditions=n, **sc))
    # shared-target pairs (S1, S3) among drugs with a target unit in this cell
    SP = []
    for win in ("early", "late"):
        keys = [(p, share[p]) for (p, w) in FULL if w == win and p in share]
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                (p1, u1), (p2, u2) = keys[i], keys[j]
                f1, f2 = FULL[(p1, win)], FULL[(p2, win)]
                SP.append(dict(level=level, unit=unit, window=win, p1=p1, p2=p2, unit1=u1, unit2=u2, same=u1 == u2,
                               n_cond_min=int(min(f1[2], f2[2])),
                               S1_cos=cos(f1[0], f2[0]), S3_cos=cos(f1[1], f2[1])))
    return pd.DataFrame(rowsH), pd.DataFrame(PR), pd.DataFrame(SP)


def net_for(cohort, ovr, available, cohort_of):
    cand = [ovr[cohort]] if cohort in ovr else []
    if cohort in ("TCGA-COAD", "TCGA-READ"):
        cand.append("COADREAD")
    if cohort in cohort_of:
        tc = cohort_of[cohort]
        cand.append("COADREAD" if tc in ("TCGA-COAD", "TCGA-READ") else tc.split("-", 1)[-1])
    cand += [cohort.split("-", 1)[1], cohort.split("-", 1)[1].split("_")[0]]
    for c in cand:
        if c and c in available:
            return c
    return None


def boot_ci(values, clusters, stat=np.nanmedian, nb=NB, seed=SEED):
    """Cluster bootstrap 95% interval of stat(values)."""
    v = np.asarray(values, float); c = np.asarray(clusters)
    ok = ~np.isnan(v); v, c = v[ok], c[ok]
    if len(v) < 3:
        return np.nan, np.nan
    uc, inv = np.unique(c, return_inverse=True)
    groups = [np.flatnonzero(inv == k) for k in range(len(uc))]
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(nb):
        pick = rng.integers(len(groups), size=len(groups))
        bs.append(stat(v[np.concatenate([groups[k] for k in pick])]))
    return float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))


def main():
    t0 = time.time()
    fm = pd.read_csv(RES / "lincs/fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    clus = lambda p: fm.get(p, p)
    L = pd.read_csv(RES / "lincs/layered_score_v2_v2_tf.tsv", sep="\t", low_memory=False)
    L["eligible"] = (L.layered_n_axes >= 1) & (L.layered_T >= 50)
    ovr = pd.read_csv(IN / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = pd.read_csv(RES / "step1/disease_groups_toil.tsv", sep="\t").drop_duplicates("group_id").set_index("group_id").cohort.to_dict()
    parent_of = pd.read_csv(RES / "celligner/stratum_line_weights_long_v2.tsv", sep="\t").drop_duplicates("stratum").set_index("stratum").parent.to_dict()
    avail = {p.stem for p in VDIR.glob("*.f32")}
    # PROGENy weights and D99 labels
    pw = pd.read_csv(IN / "data/raw/progeny/progeny_human_top100.tsv", sep="\t")
    sym2id = genes.drop_duplicates("symbol").set_index("symbol").gene_id
    pw = pw[pw.target.isin(sym2id.index)]
    prog_genes = list(dict.fromkeys(sym2id[pw.target].values))
    gi = {g: i for i, g in enumerate(prog_genes)}
    prog_w = {}
    for p, d in pw.groupby("source"):
        v = np.zeros(len(prog_genes)); v[[gi[sym2id[t]] for t in d.target]] = d.weight.values; prog_w[p] = v
    tasks = []
    for level, f in (("group_id", "window_reversal.tsv.gz"), ("stratum", "window_reversal_strata.tsv.gz")):
        w = pd.read_csv(RES / "systems_validation/integrated" / f, sep="\t")
        for unit, d in w.groupby("cohort"):
            net = net_for(parent_of.get(unit, unit), ovr, avail, cohort_of)
            rows = [(r.pert_id, win) for r in d.itertuples() for win in ("early", "late") if pd.notna(getattr(r, f"rev_{win}"))]
            lu = L[(L.unit == unit) & L.layered_genes.notna()]
            share = dict(zip(lu.pert_id, lu.layered_genes))
            tasks.append((level, unit, net, rows, prog_w, prog_genes, share))
    only = [x for x in os.environ.get("ONLY_UNITS", "").split(",") if x]   # smoke tests only; never in a reported run
    if only:
        tasks = [t for t in tasks if t[1] in only]
        global OUT
        OUT = OUT.parent / "A_smoke"; OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(tasks)} units, {sum(len(t[3]) for t in tasks):,} rows; {NC} workers", flush=True)
    Hs, Ps, Ss = [], [], []
    with ProcessPoolExecutor(NC) as ex:
        for i, (h, p, s) in enumerate(ex.map(work, tasks), 1):
            Hs.append(h); Ps.append(p); Ss.append(s)
            print(f"  [{i}/{len(tasks)}] {time.time() - t0:.0f}s", flush=True)
    HA = pd.concat(Hs, ignore_index=True); PR = pd.concat(Ps, ignore_index=True); SP = pd.concat(Ss, ignore_index=True)
    HA["cluster"] = HA.pert_id.map(clus)
    elig = set(zip(L[L.eligible].unit, L[L.eligible].pert_id))
    HA["eligible"] = [(u, p) in elig for u, p in zip(HA.unit, HA.pert_id)]
    HA.to_csv(OUT / "A_halves.tsv.gz", sep="\t", index=False, float_format="%.6g")
    PR.to_csv(OUT / "A_progeny.tsv.gz", sep="\t", index=False, float_format="%.6g")

    curve, cells = [], []
    # vector steps: same-drug minus matched different-drug cosine, cluster bootstrap over ingredients
    for (lv, win), d in HA.groupby(["level", "window"]):
        for st in ("S1", "S2", "S3"):
            x = d.dropna(subset=[f"{st}_cos_same", f"{st}_cos_null"])
            diff = (x[f"{st}_cos_same"] - x[f"{st}_cos_null"]).values
            lo, hi = boot_ci(diff, x.cluster.values)
            curve.append(dict(level=lv, window=win, step=st, measure="same-drug minus different-drug cosine", n_rows=len(x),
                              n_clusters=x.cluster.nunique(), same=float(x[f"{st}_cos_same"].median()),
                              null=float(x[f"{st}_cos_null"].median()), signal=float(np.median(diff)), ci_lo=lo, ci_hi=hi,
                              resolved=bool(lo > 0)))
        # scalar steps: Spearman-Brown reliability across drugs within each cell
        for st, a, b, sel in (("S4", "rev_a", "rev_b", None), ("S5", "rev_gene_a", "rev_gene_b", None),
                              ("S6", "rev_gene_a", "rev_gene_b", "eligible")):
            rel = []
            for unit, c in d.groupby("unit"):
                c = c[c[sel]] if sel else c
                c = c.dropna(subset=[a, b])
                if len(c) < MIN_DRUGS:
                    continue
                r = spearmanr(c[a], c[b]).correlation
                rng = np.random.default_rng(rseed("perm", lv, win, unit, st))
                nullr = [spearmanr(c[a], rng.permutation(c[b].values)).correlation for _ in range(200)]
                sb = 2 * r / (1 + r) if r > -1 else np.nan
                cells.append(dict(level=lv, window=win, step=st, unit=unit, n_drugs=len(c), spearman_half=r,
                                  reliability=sb, ceiling=np.sqrt(sb) if sb > 0 else 0.0,
                                  p_perm=float((np.sum(np.array(nullr) >= r) + 1) / 201)))
                rel.append((unit, sb))
            if rel:
                u, v = zip(*rel)
                lo, hi = boot_ci(np.array(v), np.array(u))
                curve.append(dict(level=lv, window=win, step=st, measure="split-half reliability across drugs (Spearman-Brown)",
                                  n_rows=len(v), n_clusters=len(v), same=float(np.median(v)), null=0.0,
                                  signal=float(np.median(v)), ci_lo=lo, ci_hi=hi, resolved=bool(lo > 0)))
            else:
                curve.append(dict(level=lv, window=win, step=st, measure="split-half reliability across drugs (Spearman-Brown)",
                                  n_rows=0, n_clusters=0, signal=np.nan, resolved=False, note=f"no cell with >= {MIN_DRUGS} drugs"))
    CELLS = pd.DataFrame(cells); CELLS.to_csv(OUT / "A_reliability_cells.tsv", sep="\t", index=False, float_format="%.4f")

    # ---------------- PART 2: mechanism labels ----------------
    labels = []
    # (b) PROGENy on-target direction, S1
    ot = pd.read_csv(IN / "config/progeny_on_target.tsv", sep="\t", comment="#")
    ann = pd.read_csv(RES / "lincs/drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    moa = ann.hub_moa.fillna(ann.get("lincs_moa", pd.Series(dtype=str))).fillna("")
    lab = []
    for p in PR.pert_id.unique():
        mm = moa.get(p, "")
        for r in ot.itertuples():
            if mm and re.search(r.moa_regex, mm, flags=re.I):
                lab.append((p, r.pathway, int(r.direction), r.receptor)); break
    LB = pd.DataFrame(lab, columns=["pert_id", "pathway", "direction", "receptor"])
    X = PR.merge(LB, on="pert_id")
    paths = list(prog_w)
    for (lv, win), d in X.groupby(["level", "window"]):
        allp = PR[(PR.level == lv) & (PR.window == win)]
        # standardise each pathway within cell so cells and pathways are comparable
        z = allp.groupby("unit")[paths].transform(lambda s: (s - s.mean()) / (s.std() or 1))
        zz = pd.concat([allp[["unit", "pert_id"]], z], axis=1)
        dd = d[["unit", "pert_id", "pathway", "direction"]].merge(zz, on=["unit", "pert_id"])
        if not len(dd):
            continue
        Z = dd[paths].to_numpy(dtype=float)
        pidx = np.array([paths.index(x) for x in dd.pathway]); dirs = dd.direction.to_numpy(dtype=float)
        rows_i = np.arange(len(dd))
        signed = dirs * Z[rows_i, pidx]
        share = float(np.mean(signed > 0))
        rng = np.random.default_rng(rseed("progeny", lv, win))
        null = []
        for _ in range(NP):
            k = rng.permutation(len(dd))          # labels (pathway, direction) move together
            null.append(float(np.mean(dirs[k] * Z[rows_i, pidx[k]] > 0)))
        lo, hi = boot_ci((signed > 0).astype(float), dd.pert_id.map(clus).values, stat=np.mean)
        labels.append(dict(label="b PROGENy on-target direction", step="S1", level=lv, window=win, n=len(dd),
                           n_clusters=dd.pert_id.map(clus).nunique(), value=share, null=float(np.median(null)),
                           ci_lo=lo, ci_hi=hi, resolved=bool(lo > float(np.median(null))),   # frozen rule: interval excludes the matched null
                           p_perm=float((np.sum(np.array(null) >= share) + 1) / (NP + 1))))
    # (c) shared target at S1 and S3
    SP.to_csv(OUT / "A_shared_target_pairs.tsv.gz", sep="\t", index=False, float_format="%.5g")
    for (lv, win), d in SP.groupby(["level", "window"]):
        s, n = d[d.same], d[~d.same]
        if len(s) < 5:
            continue
        for st in ("S1", "S3"):
            # matched null: each sharing pair against non-sharing pairs of the same cell
            diffs, cl = [], []
            for (unit, u1), g in s.groupby(["unit", "unit1"]):
                base = n[n.unit == unit][f"{st}_cos"].median()
                diffs.append(g[f"{st}_cos"].median() - base); cl.append(f"{unit}|{u1}")
            lo, hi = boot_ci(np.array(diffs), np.array([c.split("|")[1] for c in cl]))
            labels.append(dict(label="c shared target unit vs other pairs", step=st, level=lv, window=win,
                               n=len(diffs), n_clusters=len(set(c.split("|")[1] for c in cl)),
                               value=float(np.nanmedian(diffs)), null=0.0, ci_lo=lo, ci_hi=hi, resolved=bool(lo > 0)))
    # (a) primary target rank and (c) at S5, from the production regulator table (FDA drugs only; chunked)
    sp = pd.concat([pd.read_csv(RES / f"lincs/{f}", sep="\t", low_memory=False,
                                usecols=lambda c: c in (k, "pert_id", "evidence_targets", "target_source")).rename(columns={k: "unit"})
                    for f, k in (("drug_reversal_split_v2_tf.tsv", "group_id"), ("stratum_reversal_split_v2_tf.tsv", "stratum"))])
    sp = sp[sp.target_source.isin(["chembl_mechanism", "confirmed_annotation"]) & sp.evidence_targets.notna()]
    tgt = {(u, p): {x.split(":")[0].strip() for x in re.split(r"[|;]", e) if x.strip()} for u, p, e in zip(sp.unit, sp.pert_id, sp.evidence_targets)}
    orient = dict(zip(zip(L.unit, L.pert_id), L.layered_orient))
    shareu = dict(zip(zip(L.unit, L.pert_id), L.layered_genes))
    fda = set(L.pert_id)
    TR, S5V, PAIRS5 = [], {}, []
    for lv, f in (("group_id", "regulator_contributions_by_window_all.tsv.gz"), ("stratum", "regulator_contributions_by_window_all_strata.tsv.gz")):
        parts = [ch[ch.pert_id.isin(fda)] for ch in pd.read_csv(
            RES / "lincs/regulon_contributions_v2_tf" / f, sep="\t", chunksize=3_000_000,
            usecols=["cohort", "pert_id", "window", "regulator", "drug_activity", "rev_contrib", "rev_gene_contrib"])]
        allf = pd.concat(parts, ignore_index=True)       # whole groups: no group is split across chunks
        if True:
            for (u, p, w), g in allf.groupby(["cohort", "pert_id", "window"], sort=False):
                if (u, p) in shareu and isinstance(shareu[(u, p)], str):
                    S5V.setdefault((lv, u, w), {}).setdefault(p, g.set_index("regulator").rev_gene_contrib)
                T = tgt.get((u, p))
                if not T:
                    continue
                hit = g[g.regulator.isin(T)]
                if not len(hit):
                    continue
                o = orient.get((u, p), "")
                sa = -g.drug_activity if o == "lower" else (g.drug_activity if o == "raise" else g.drug_activity.abs())
                pct = lambda s: s.rank(pct=True)
                ps3, ps4, ps5 = pct(sa), pct(g.rev_contrib), pct(g.rev_gene_contrib)
                for i in hit.index:
                    TR.append(dict(level=lv, unit=u, pert_id=p, window=w, target=g.regulator[i], n_regs=len(g),
                                   eligible=(u, p) in elig, S3=ps3[i], S4=ps4[i], S5=ps5[i]))
    TRD = pd.DataFrame(TR); TRD.to_csv(OUT / "A_target_rank.tsv.gz", sep="\t", index=False, float_format="%.5g")
    for (lv, win), d in TRD.groupby(["level", "window"]):
        for st, dd in (("S3", d), ("S4", d), ("S5", d), ("S6", d[d.eligible])):
            col = "S5" if st == "S6" else st
            if len(dd) < 5:
                continue
            v = dd[col].values
            lo, hi = boot_ci(v, dd.target.values, stat=np.mean)
            labels.append(dict(label="a primary target percentile in own response", step=st, level=lv, window=win,
                               n=len(dd), n_clusters=dd.target.nunique(), value=float(np.mean(v)), null=0.5,
                               ci_lo=lo, ci_hi=hi, resolved=bool(lo > 0.5),
                               top5pct=float(np.mean(v >= 0.95)), chance_top5pct=0.05))
    # (c) at S5: cosine of rev_gene_contrib vectors, sharing vs non-sharing pairs
    for (lv, u, w), vecs in S5V.items():
        ps = list(vecs)
        if len(ps) < 3:
            continue
        M = pd.concat([vecs[p] for p in ps], axis=1).fillna(0).values.T
        nrm = np.linalg.norm(M, axis=1); nrm[nrm == 0] = np.nan
        Cm = (M @ M.T) / np.outer(nrm, nrm)
        units = [shareu[(u, p)] for p in ps]
        for i in range(len(ps)):
            for j in range(i + 1, len(ps)):
                PAIRS5.append(dict(level=lv, unit=u, window=w, unit1=units[i], same=units[i] == units[j], S5_cos=Cm[i, j]))
    P5 = pd.DataFrame(PAIRS5)
    for (lv, win), d in P5.groupby(["level", "window"]) if len(P5) else []:
        s, n = d[d.same], d[~d.same]
        diffs, cl = [], []
        for (unit, u1), g in s.groupby(["unit", "unit1"]):
            diffs.append(g.S5_cos.median() - n[n.unit == unit].S5_cos.median()); cl.append(u1)
        if len(diffs) < 5:
            continue
        lo, hi = boot_ci(np.array(diffs), np.array(cl))
        labels.append(dict(label="c shared target unit vs other pairs", step="S5", level=lv, window=win, n=len(diffs),
                           n_clusters=len(set(cl)), value=float(np.nanmedian(diffs)), null=0.0, ci_lo=lo, ci_hi=hi,
                           resolved=bool(lo > 0)))
    LAB = pd.DataFrame(labels); LAB.to_csv(OUT / "A_labels.tsv", sep="\t", index=False, float_format="%.4f")
    CUR = pd.DataFrame(curve); CUR.to_csv(OUT / "A_curve.tsv", sep="\t", index=False, float_format="%.4f")
    run = dict(stage="A_resolvable", protocol="5cc2227", rows_split=len(HA), runtime_s=round(time.time() - t0, 1),
               seed=SEED, n_boot=NB, n_perm=NP, finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 250)
    print(CUR.to_string(index=False)); print(LAB.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
