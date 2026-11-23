#!/usr/bin/env python
"""Per-window breakdown of the drug reversal score: window -> regulator -> regulon target (D182, 2026-11-22).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

The timepoint-resolved twin of `79e2_regulon_contributions.py`. 79e2 splits the POOLED score (79c averages
every dose and time per drug x line), so every one of its rows carries `window = "unknown"` (D180) and it
cannot say WHEN a drug acts. This one splits `rev_early` and `rev_late` - the D155 per-window scores built on
the condition-level profiles - so the drivers come out separately for each window.

It re-scores nothing. It rebuilds the exact input of the window scorer and splits it additively.

WHAT IT MIRRORS. `systems_validation/integrate/window_reversal.py`, not 79e: that script is what produced
`rev_early` / `rev_late`, and it differs from 79e in ways that change the answer -
  * conditions, not drug x line profiles (`drug_condition_profiles_specific.f32`, 26,130 x 10,112);
  * in-key lines only (`disease_line_weights_long.in_key`), so there is no direct / borrowed split;
  * one score per (cohort, drug, window), window = early (<= 6 h) / late (>= 24 h), `other` and
    unrecorded times (D180) dropped;
  * `rev = -(D . s) / ||s||` over the cohort's regulators, as 79e:285.
Following 79e instead would silently produce numbers that do not match the window table.

THE SPLIT, as in 79e2. Condition-level VIPER (`drug_viper_specific_v2_tf_cond`) has no null model, so
NES = ES x w and the per-target split is EXACT. Universe U = condition-profile genes that are a target of
any scored regulator of the network (viper eset.filter); per condition rank over U,
t2 = qnorm(rank/(|U|+1)), t1 = qnorm(|2 rank/(|U|+1) - 1| + (1 - max)/2); for regulator i, target j
(mode m_j, likelihood l_j): a_j = m_j l_j t2_j, b_j = (1-|m_j|) l_j t1_j, S1 = sum a, S2 = sum b,
s = sign(S1) (0 -> +1), and A_i = sum_j c_ij with c_ij = w (a_j + s b_j [S2>0]) / sum l. The window mean over
conditions is linear, so the terms carry through to `rev_<window>` exactly.

CHECKS (the run stops if either fails):
  1. activity recomputed for CHECK_N random conditions per network vs the production
     `drug_viper_specific_v2_tf_cond/<net>.f32`;
  2. the regulator terms summed vs `rev_early` / `rev_late` in `window_reversal.tsv.gz`, for every row written
     (that table stores rev to 5 dp).

SCOPE. Only the 29 cohorts with condition-level data have windows at all; no subgroup does (D180). Written for
every (cohort, drug, window) of the FDA-ranking drugs by default:
  regulator_contributions_by_window.tsv.gz  every tested regulator:
      cohort pert_id cmap_name window n_conditions n_lines regulator disease_nes drug_activity rev_contrib
      rev rank_in_drug
  target_contributions_by_window.tsv.gz     the TOP_REGS (default 25) regulators per row, every target:
      ... regulator target target_symbol tfmode likelihood drug_z activity_contrib rev_contrib
      disease_z disease_contrib_nes (the same edge on the disease side, 50b)
  window_shift.tsv  per (cohort, drug) with BOTH windows: the regulators whose contribution moves most
      between early and late (`delta` = early - late) - the point of the whole exercise.
ALL=1 scores every drug with window data; QUERY=<cohort>:<drug>[,...] writes the full array for those rows
only (query_window_<cohort>_<drug>.tsv) and skips the bulk tables.

Env: SFX (_v2_tf) | TOP_REGS | CHECK_N (200) | N_CORES (6) | ONLY_COHORTS | ALL | QUERY | SHIFT_TOP (15)
Needs: 50b (disease side, for the joined columns) and the D155 window run.
Output dir: results/lincs/regulon_contributions<SFX>/
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # N_CORES workers x 1 thread (BLAS would otherwise spawn its own)

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
PROC, RES = ROOT / "data/processed/lincs", ROOT / "results/lincs"
SFX = os.environ.get("SFX", "_v2_tf")
VDIR = PROC / f"drug_viper_specific{SFX}_cond"
WIN = ROOT / "results/systems_validation/integrated/window_reversal.tsv.gz"
MR_DIR = ROOT / f"results/step6/aracne_project{SFX}"
OUT = RES / f"regulon_contributions{SFX}"
TOP_REGS = int(os.environ.get("TOP_REGS", "25"))
COSINE = os.environ.get("COSINE") == "1"   # D194: per-regulon cosine between the DISEASE and DRUG
                                           # contribution vectors - the pipeline's own per-gene
                                           # edge-weighted quantities, not a reconstruction from z.
NO_TARGETS = os.environ.get("NO_TARGETS") == "1"   # D184: regulator vectors only; the target table at ALL=1 scale
                                                   # would be ~95M rows and the concordance test does not use it.
SHIFT_TOP = int(os.environ.get("SHIFT_TOP", "15"))
CHECK_N = int(os.environ.get("CHECK_N", "200"))
N_CORES = int(os.environ.get("N_CORES", "6"))
MIN_MRS = 2

genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")   # same gene axis as the condition matrix
NG = len(genes)
idx = pd.read_csv(PROC / "drug_condition_profiles_index.tsv", sep="\t")
idx["col"] = np.arange(len(idx))
# D180: an unrecorded time (sentinel -1.0) is "unknown", never early. Only early / late are scored.
_t = idx.pert_time.where(idx.pert_time >= 0)
idx["window"] = np.where(_t.isna(), "unknown", np.where(_t <= 6, "early", np.where(_t >= 24, "late", "other")))
P = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype=np.float32, mode="r", shape=(len(idx), NG))
gpos = pd.Series(np.arange(NG), index=genes.gene_id)
sym = genes.set_index("gene_id").symbol


class Net:
    """Edge arrays of one network's condition-arm regulators, ordered as 79d passed them to viper."""

    def __init__(self, net, edges):
        self.regs = pd.Index((VDIR / f"{net}.regulators.txt").read_text().split())
        e = edges[edges.net == net]
        e = e[e.regulator.isin(set(self.regs)) & e.target.isin(gpos.index)]
        e = e.assign(_k=self.regs.get_indexer(e.regulator)).sort_values("_k", kind="stable")
        assert (e._k >= 0).all()
        self.U = np.sort(gpos[e.target.unique()].values)          # viper eset.filter
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

    def contrib(self, cols):
        """conditions x edges: c_cij, summing over each regulator's edges to its NES."""
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


def net_for(cohort, ovr, available, cohort_of):
    """Network for a cohort, exactly as window_reversal.net_for (D171 production rule)."""
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


G = {}


def run_cohort(args):
    cohort, net_name, sel, full = args
    N, S = G["nets"][net_name], G["S"]
    m = S["mr"][(S["mr"].group_id == cohort) & S["mr"].regulator.isin(N.regs)]
    ci = N.regs.get_indexer(m.regulator)
    nes = np.zeros(len(N.regs)); nes[ci] = m.nes.values
    tested = np.zeros(len(N.regs), bool); tested[ci] = True
    norm = np.linalg.norm(m.nes.values)
    lw = S["lw"]
    wg = lw[(lw.group_id == cohort) & lw.in_key].set_index("ModelID").weight
    sub = idx[idx.ModelID.isin(wg.index) & idx.window.isin(("early", "late"))]
    dz = G["dis"].get(cohort)
    dvec = None
    if COSINE and dz is not None:
        dm = dz.drop_duplicates(["regulator", "target"]).set_index(["regulator", "target"]).disease_contrib_nes
        dvec = np.nan_to_num(dm.reindex(
            pd.MultiIndex.from_arrays([N.regs[N.reg], N.target])).to_numpy(dtype=np.float64))
    R, T = [], []
    for r in sel.itertuples():
        g = sub[(sub.pert_id == r.pert_id) & (sub.window == r.window)]
        if not len(g):
            continue
        ww = wg.reindex(g.ModelID).fillna(0).values
        ww = ww if ww.sum() > 0 else np.ones(len(g))
        ww = ww / ww.sum()
        C = ww @ N.contrib(g.col.values)
        cos_term = None                       # D195: per-gene share of the regulon cosine
        act = np.add.reduceat(C, N.starts)
        rc = -nes * act / norm
        reg = pd.DataFrame({"cohort": cohort, "pert_id": r.pert_id, "cmap_name": r.cmap_name,
                            "window": r.window, "n_conditions": len(g), "n_lines": g.ModelID.nunique(),
                            "regulator": N.regs, "disease_nes": nes, "drug_activity": act,
                            "rev_contrib": rc, "rev": rc[tested].sum(), "rev_stored": r.rev})[tested]
        reg["rank_in_drug"] = reg.rev_contrib.abs().rank(ascending=False, method="first").astype(int)

        # D194 GRANULAR OPPOSITION VIEW. Within each regulon, the cosine between the gene-level
        # DISEASE contributions (50b) and the gene-level DRUG contributions (this script's C).
        # Both already carry each edge's own tfmode x likelihood - that is what a contribution IS -
        # so the per-gene edge weights are maintained end to end rather than rebuilt from z-scores.
        # A NEGATIVE cosine means the drug moves the very genes that account for this regulator's
        # dysregulation, in the opposite direction, gene by gene.
        if COSINE and dvec is not None:
            dot = np.add.reduceat(dvec * C, N.starts)
            nd = np.sqrt(np.add.reduceat(dvec * dvec, N.starts))
            nx = np.sqrt(np.add.reduceat(C * C, N.starts))
            den = nd * nx
            cos = np.divide(dot, den, out=np.full(len(den), np.nan), where=den > 0)
            # reduceat yields one value per regulator group, aligned with N.regs exactly as
            # `act` above is; `reg` is already filtered to the tested regulators.
            assert len(cos) == len(N.regs), (len(cos), len(N.regs))
            reg["cos_contrib"] = cos[tested]
            reg["n_genes_regulon"] = np.diff(np.r_[N.starts, len(N.reg)])[tested]
            # D195: the SAME nesting one level down. Each gene's share of its regulon's cosine,
            #   cos_term_j = d_j * x_j / (||d|| * ||x||),      sum_j in r  cos_term_j = cos_contrib_r
            # so the per-gene scores add back to the regulon score exactly, as the D179 target
            # contributions add back to the NES. NEGATIVE = that gene is pushed opposite to its
            # disease contribution, and |cos_term| says how much of the regulon's opposition it
            # carries. This is the gene-level opposition score, not just the inputs to one.
            cos_term = np.divide(dvec * C, den[N.reg],
                                 out=np.zeros(len(dvec)), where=den[N.reg] > 0)
        R.append(reg)
        if NO_TARGETS:
            continue
        keep = tested.copy() if full else np.isin(
            np.arange(len(N.regs)), N.regs.get_indexer(reg.regulator[reg.rank_in_drug <= TOP_REGS]))
        e = keep[N.reg]
        if e.any():
            t = pd.DataFrame({"cohort": cohort, "pert_id": r.pert_id, "cmap_name": r.cmap_name,
                              "window": r.window, "regulator": N.regs[N.reg[e]], "target": N.target[e],
                              "target_symbol": sym.reindex(N.target[e]).values,
                              "tfmode": N.m[e], "likelihood": N.l[e],
                              "drug_z": ww @ np.asarray(P[g.col.values][:, N.U[N.gene[e]]], dtype=np.float64),
                              "activity_contrib": C[e], "rev_contrib": -nes[N.reg[e]] * C[e] / norm})
            if cos_term is not None:
                # exact decomposition: these sum to the regulon's cos_contrib
                t["cos_term"] = cos_term[e]
                t["opposes"] = cos_term[e] < 0
            if dz is not None:
                t = t.merge(dz, on=["regulator", "target"], how="left")
            T.append(t)
    return (pd.concat(R) if R else pd.DataFrame()), (pd.concat(T) if T else pd.DataFrame())


def check_net(net_name):
    N = G["nets"][net_name]
    cols = np.random.default_rng(0).choice(len(idx), min(CHECK_N, len(idx)), replace=False)
    A = np.add.reduceat(N.contrib(cols), N.starts, axis=1)
    return net_name, float(np.abs(A - N.A[cols]).max())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    w = pd.read_csv(WIN, sep="\t", low_memory=False)
    long = []
    for side in ("early", "late"):
        m = w[f"rev_{side}"].notna()
        long.append(w.loc[m, ["cohort", "network", "pert_id", "cmap_name"]].assign(
            window=side, rev=w.loc[m, f"rev_{side}"].values,
            n_conditions_tbl=w.loc[m, f"n_conditions_{side}"].values))
    L = pd.concat(long, ignore_index=True)

    q = [x.split(":", 1) for x in os.environ.get("QUERY", "").split(",") if x]
    if q:
        sel = pd.concat([L[(L.cohort == c) & ((L.cmap_name == d) | (L.pert_id == d))] for c, d in q])
        assert len(sel), f"QUERY matched nothing: {q}"
    elif os.environ.get("ALL") == "1":
        sel = L
    else:
        fda = pd.read_csv(RES / f"fda_ranking{SFX}.tsv", sep="\t", usecols=["group_id", "cmap_name"], low_memory=False)
        sel = L.merge(fda.drop_duplicates().rename(columns={"group_id": "cohort"}), on=["cohort", "cmap_name"])
    only = [x for x in os.environ.get("ONLY_COHORTS", "").split(",") if x]
    if only:
        sel = sel[sel.cohort.isin(only)]
    assert len(sel), "nothing selected"

    S = {}
    S["mr"] = pd.read_csv(MR_DIR / "group_mr.tsv", sep="\t")
    S["lw"] = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    ovr = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()
    cohort_of = (pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
                 .drop_duplicates("group_id").set_index("group_id").cohort.to_dict())
    avail = {p.stem for p in VDIR.glob("*.f32")}
    sel = sel.assign(net=[net_for(c, ovr, avail, cohort_of) for c in sel.cohort])
    sel = sel[sel.net.notna()]
    G["S"] = S
    edges = pd.read_csv(ROOT / f"data/processed/step5/regulons_aracne_project{SFX}_edges.tsv.gz", sep="\t")
    G["nets"] = {n: Net(n, edges) for n in sorted(sel.net.unique())}
    d = pd.read_csv(MR_DIR / "regulon_contributions.tsv.gz", sep="\t",
                    usecols=["group_id", "regulator", "target", "z", "contrib_nes"])
    d = d[d.group_id.isin(set(sel.cohort))]
    G["dis"] = {u: x.drop(columns="group_id").rename(columns={"z": "disease_z", "contrib_nes": "disease_contrib_nes"})
                for u, x in d.groupby("group_id")}
    print(f"{len(sel):,} (cohort, drug, window) rows over {sel.cohort.nunique()} cohorts, "
          f"{len(G['nets'])} networks; windows {sel.window.value_counts().to_dict()}", flush=True)

    with ProcessPoolExecutor(N_CORES) as ex:
        chk = dict(ex.map(check_net, G["nets"]))
        print(f"check 1 (activity vs the _cond .f32, {CHECK_N} conditions/net): max |diff| {max(chk.values()):.2e}",
              flush=True)
        bad = {k: v for k, v in chk.items() if v > 1e-4}
        assert not bad, f"activity does not reproduce the condition VIPER arm: {bad}"
        jobs = [(c, g.net.iloc[0], g, bool(q)) for c, g in sel.groupby("cohort")]
        out = list(ex.map(run_cohort, jobs))
    R = pd.concat([o[0] for o in out if len(o[0])])
    T = pd.concat([o[1] for o in out if len(o[1])]) if any(len(o[1]) for o in out) else pd.DataFrame()
    u = R.drop_duplicates(["cohort", "pert_id", "window"])
    diff = (u.rev - u.rev_stored).abs()
    print(f"check 2 (regulator terms vs stored rev_<window>, {len(u):,} rows): max |diff| {diff.max():.2e}", flush=True)
    assert diff.max() <= 1e-4, "does not reproduce the window scores"
    R = R.drop(columns="rev_stored")

    if q:
        _rc = ["cohort", "pert_id", "window", "regulator", "disease_nes", "drug_activity",
               "rev_contrib", "rank_in_drug"] + [c for c in ("cos_contrib", "n_genes_regulon") if c in R]
        full = T.merge(R[_rc].rename(columns={"rev_contrib": "regulator_rev_contrib"}),
                       on=["cohort", "pert_id", "window", "regulator"])
        for (c, dn), g in full.groupby(["cohort", "cmap_name"]):
            f = OUT / f"query_window_{c}_{dn}.tsv"
            g.sort_values(["window", "rank_in_drug", "rev_contrib"]).to_csv(f, sep="\t", index=False, float_format="%.6g")
            print(f"-> {f} ({len(g):,} edges)")
        return

    sfx = "_all" if os.environ.get("ALL") == "1" else ""
    R.to_csv(OUT / f"regulator_contributions_by_window{sfx}.tsv.gz", sep="\t", index=False, float_format="%.6g")
    if len(T):
        T.to_csv(OUT / f"target_contributions_by_window{sfx}.tsv.gz", sep="\t", index=False, float_format="%.6g")

    # the point of the exercise: which regulators move between windows, for drugs that have both
    p = R.pivot_table(index=["cohort", "pert_id", "cmap_name", "regulator", "disease_nes"],
                      columns="window", values="rev_contrib", aggfunc="first").reset_index()
    p = p[p.early.notna() & p.late.notna()].copy()
    p["delta"] = p.early - p.late   # not "shift": a column named shift shadows DataFrame.shift
    p["flips_sign"] = np.sign(p.early) != np.sign(p.late)
    p["rank_shift"] = p.groupby(["cohort", "pert_id"])["delta"].transform(
        lambda v: v.abs().rank(ascending=False, method="first"))
    sh = p[p.rank_shift <= SHIFT_TOP].drop(columns="rank_shift").sort_values(
        ["cohort", "cmap_name", "delta"], ascending=[True, True, False])
    sh.to_csv(OUT / "window_shift.tsv", sep="\t", index=False, float_format="%.6g")

    both = p.groupby(["cohort", "pert_id"]).size()
    print(f"\nregulator level {len(R):,} rows, target level {len(T):,} rows (top {TOP_REGS} regulators per row)")
    print(f"window_shift: {len(sh):,} rows over {len(both):,} (cohort, drug) pairs with both windows; "
          f"{int(p.flips_sign.sum()):,} of {len(p):,} regulator terms flip sign between windows "
          f"({100 * p.flips_sign.mean():.1f}%)")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
