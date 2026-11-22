#!/usr/bin/env python
"""Per-target breakdown of drug regulator activity and of the reversal score (user 2026-11-22, D179; runs in `celligner`).

The drug-side half of 50b_regulon_contributions.R. It keeps the nested array drug reversal -> regulator -> regulon target
that 79d / 79e collapse. It re-scores nothing. It rebuilds the exact 79d viper input and splits each number into additive
parts. Read-only on every production file.

1. Drug activity (79d: viper(P, reg, method="none"), no null model, so NES = ES * w exactly, w = sqrt(sum (l/max l)^2)).
   Universe U = LINCS profile genes that are a target of any kept regulator of the network (viper eset.filter); per
   profile rank over U, t2 = qnorm(rank/(|U|+1)), t1 = qnorm(|2 rank/(|U|+1) - 1| + (1 - max)/2). For regulator i and
   target j (mode m_j, likelihood l_j): a_j = m_j l_j t2_j, b_j = (1-|m_j|) l_j t1_j, S1 = sum a, S2 = sum b,
   s = sign(S1) (0 -> +1); activity A_i = w s (|S1| + S2 [S2>0]) / sum l = sum_j c_ij, c_ij = w (a_j + s b_j [S2>0]) / sum l.
   The split is exact, with nothing allocated in proportion.
2. Drug activity per disease unit (79e): weighted mean over the drug's lines (direct: in the tissue key, weighted by the
   unit's Celligner line weights, equal if all 0; borrowed: unweighted over lines outside the key). Linear, so the
   per-target terms carry through: cbar_ij = sum_p ww_p c_pij / sum ww.
3. Reversal (79e SCORE=weighted): rev = -sum_i nes_i a_i / sqrt(sum nes^2) over the regulators tested in the unit.
   rev_contrib_i = -nes_i a_i / norm ; rev_contrib_ij = -nes_i cbar_ij / norm. Both sum to rev exactly.

Checks (stop on failure): (1) A recomputed for CHECK_N random profiles per network vs the production
drug_viper_specific<SFX>/<net>.f32 (|diff| <= 1e-4, float32 storage); (2) sum_i rev_contrib_i vs `rev` in
drug_reversal<SFX>.tsv / stratum_reversal<SFX>.tsv for every row written (|diff| <= 1e-3, rev is stored to 3 dp).

Scope. The full profile x regulator x target array is ~1e9 values per network, so it is not written. Written:
  - every (unit, drug) row of the FDA rankings (fda_ranking<SFX>.tsv + fda_ranking_strata<SFX>.tsv), or ALL=1 for every
    scored drug (regulator level only is then advisable: TOP_REGS=0)
  - regulator level, all tested regulators: regulator_contributions.tsv.gz
      unit_id pert_id cmap_name source n_lines regulator disease_nes drug_activity rev_contrib rev rank_in_drug
  - target level for the TOP_REGS (default 25) regulators per row with the largest |rev_contrib|:
      target_contributions.tsv.gz: unit_id pert_id cmap_name regulator target target_symbol tfmode likelihood
      drug_z (weighted mean MODZ z of the target over the lines used) activity_contrib rev_contrib
      disease_z disease_contrib_nes (the same edge on the disease side, 50b; NaN where the target is outside the LINCS
      or signature universe)
  - QUERY=<unit>:<cmap_name or pert_id>[,...] writes the full array (every regulator, every target) for those rows
    only, to query_<unit>_<drug>.tsv, and skips the bulk tables.
Env: VERSION=v2 REGSET=tf (production) | TOP_REGS | CHECK_N (default 200) | ONLY_UNITS | N_CORES (default 6) | ALL=1 | QUERY
Needs: 50b with EDGES=1 (data/processed/step5/regulons_aracne_project<SFX>_edges.tsv.gz); for disease_z, 50b for groups
and STRATA=1.
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
SFX = ("_v2" if os.environ.get("VERSION", "v2") == "v2" else "") + ("_tf" if os.environ.get("REGSET", "tf") == "tf" else "")
MIN_TARGETS, MIN_FRAC = 25, 0.5                                 # 79d defaults (D139)
TOP_REGS = int(os.environ.get("TOP_REGS", "25"))
CHECK_N = int(os.environ.get("CHECK_N", "200"))
N_CORES = int(os.environ.get("N_CORES", "6"))
OUT = RES / f"regulon_contributions{SFX}"
MR_DIR = ROOT / f"results/step6/aracne_project{SFX}"

genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
idx["col"] = np.arange(len(idx))
NG = len(genes)
P = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype=np.float32, mode="r", shape=(len(idx), NG))
gpos = pd.Series(np.arange(NG), index=genes.gene_id)
sym = genes.set_index("gene_id").symbol


class Net:
    """Edge arrays of one network's scored regulators, pruned and ordered exactly as 79d passes them to viper."""

    def __init__(self, net, edges, regtab):
        e = edges[edges.net == net]
        kept = regtab[(regtab.net == net) & regtab.kept].regulator
        e = e[e.regulator.isin(set(kept)) & e.target.isin(gpos.index)]
        self.U = np.sort(gpos[e.target.unique()].values)          # viper eset.filter (regulator names are symbols: no-op)
        size = e.groupby("regulator").size()
        self.regs = pd.Index((PROC / f"drug_viper_specific{SFX}" / f"{net}.regulators.txt").read_text().split())
        assert set(size[size >= MIN_TARGETS].index) == set(self.regs), f"{net}: scored regulators differ from 79d"
        e = e.assign(_k=self.regs.get_indexer(e.regulator)).sort_values("_k", kind="stable")
        upos = pd.Series(np.arange(len(self.U)), index=gpos.index[self.U])
        self.reg = self.regs.get_indexer(e.regulator)
        self.gene = upos[e.target].values                          # column in the universe
        self.target = e.target.values
        self.m, self.l = e.tfmode.values, e.likelihood.values
        self.starts = np.r_[0, np.flatnonzero(np.diff(self.reg)) + 1]
        L = np.add.reduceat(self.l, self.starts)
        lmax = np.maximum.reduceat(self.l, self.starts)
        self.w = np.sqrt(np.add.reduceat((self.l / lmax[self.reg]) ** 2, self.starts))
        self.k = (self.w / L)[self.reg]                            # per-edge w / sum l
        self.A = np.memmap(PROC / f"drug_viper_specific{SFX}" / f"{net}.f32", dtype=np.float32, mode="r",
                           shape=(len(idx), len(self.regs)))

    def contrib(self, cols):
        """profiles x edges: c_pij, summing over each regulator's edges to its NES."""
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


def setup():
    grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
    netmap = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict() if SFX == "_v2_tf" else {}
    net_of = lambda c: "COADREAD" if c in ("TCGA-COAD", "TCGA-READ") else c.split("-", 1)[1]
    sd = pd.read_csv(RES / "stratum_drugs_v2.tsv", sep="\t", usecols=["group_id", "parent_group"]).drop_duplicates("group_id")
    par = sd.set_index("group_id").parent_group.to_dict()
    unit_net = lambda u: netmap.get(par.get(u, u), net_of(grp.loc[par.get(u, u), "cohort"]))
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    w = pd.concat([pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t"),
                   pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})])
    mr = pd.concat([pd.read_csv(MR_DIR / "group_mr.tsv", sep="\t"), pd.read_csv(f"{MR_DIR}_strata/group_mr.tsv", sep="\t")])
    rv = pd.concat([pd.read_csv(RES / f"drug_reversal{SFX}.tsv", sep="\t", usecols=["group_id", "pert_id", "cmap_name", "source", "rev"]),
                    pd.read_csv(RES / f"stratum_reversal{SFX}.tsv", sep="\t", usecols=["group_id", "pert_id", "cmap_name", "source", "rev"])])
    return dict(unit_net=unit_net, par=par, keys=keys, models=models, w=w, mr=mr, rv=rv)


def key_models(S, u):
    k = S["keys"].loc[S["par"].get(u, u)]
    km = S["models"][S["models"].OncotreeLineage.isin(str(k.lineages).split(";"))]
    if isinstance(k.diseases, str) and k.diseases:
        km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
    return set(km.index)


def line_weights(S, u, pert, source, km):
    """79e: the drug's profiles used for unit u and their weights."""
    p = idx[idx.pert_id == pert]
    wg = S["w"][S["w"].group_id == u].set_index("ModelID").weight
    if source == "direct":
        p = p[p.ModelID.isin(km)]
        ww = wg.reindex(p.ModelID).fillna(0).values
    else:
        p = p[~p.ModelID.isin(km)]
        ww = np.zeros(len(p))
    ww = ww if ww.sum() > 0 else np.ones(len(p))
    return p.col.values, ww / ww.sum()


def run_unit(args):
    u, rows, net_name, full = args
    S, N = G["S"], G["nets"][net_name]
    km = key_models(S, u)
    m = S["mr"][(S["mr"].group_id == u) & S["mr"].regulator.isin(N.regs)]
    ci = N.regs.get_indexer(m.regulator)
    nes = np.zeros(len(N.regs)); nes[ci] = m.nes.values
    tested = np.zeros(len(N.regs), bool); tested[ci] = True
    norm = np.sqrt((m.nes.values ** 2).sum())
    dz = G["dis"].get(u)
    R, T = [], []
    for r in rows.itertuples():
        cols, ww = line_weights(S, u, r.pert_id, r.source, km)
        C = ww @ N.contrib(cols)                                    # edges: cbar_ij (weighted over lines)
        act = np.add.reduceat(C, N.starts)
        rc = -nes * act / norm
        rev = rc[tested].sum()
        reg = pd.DataFrame({"unit_id": u, "pert_id": r.pert_id, "cmap_name": r.cmap_name, "source": r.source,
                            "n_lines": len(cols), "regulator": N.regs, "disease_nes": nes, "drug_activity": act,
                            "rev_contrib": rc, "rev": rev, "rev_stored": r.rev})[tested]
        reg["rank_in_drug"] = reg.rev_contrib.abs().rank(ascending=False, method="first").astype(int)
        R.append(reg)
        keep = tested if full else np.isin(np.arange(len(N.regs)), N.regs.get_indexer(reg.regulator[reg.rank_in_drug <= TOP_REGS]))
        e = keep[N.reg]
        if e.any():
            t = pd.DataFrame({"unit_id": u, "pert_id": r.pert_id, "cmap_name": r.cmap_name, "regulator": N.regs[N.reg[e]],
                              "target": N.target[e], "target_symbol": sym.reindex(N.target[e]).values,
                              "tfmode": N.m[e], "likelihood": N.l[e],
                              "drug_z": ww @ np.asarray(P[cols][:, N.U[N.gene[e]]], dtype=np.float64),
                              "activity_contrib": C[e], "rev_contrib": (-nes[N.reg[e]] * C[e] / norm)})
            if dz is not None:
                t = t.merge(dz, on=["regulator", "target"], how="left")
            T.append(t)
    return pd.concat(R), (pd.concat(T) if T else pd.DataFrame())


def check_net(net_name):
    N = G["nets"][net_name]
    cols = np.random.default_rng(0).choice(len(idx), min(CHECK_N, len(idx)), replace=False)
    A = np.add.reduceat(N.contrib(cols), N.starts, axis=1)
    return net_name, float(np.abs(A - N.A[cols]).max())


G = {}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    S = setup()
    edges = pd.read_csv(ROOT / f"data/processed/step5/regulons_aracne_project{SFX}_edges.tsv.gz", sep="\t")
    regtab = pd.read_csv(RES / f"drug_viper_regulators{SFX}.tsv", sep="\t")
    q = [x.split(":", 1) for x in os.environ.get("QUERY", "").split(",") if x]
    if q:
        sel = pd.concat([S["rv"][(S["rv"].group_id == u) & ((S["rv"].cmap_name == d) | (S["rv"].pert_id == d))] for u, d in q])
        assert len(sel), f"QUERY matched nothing: {q}"
    elif os.environ.get("ALL") == "1":
        sel = S["rv"]
    else:
        fda = pd.concat([pd.read_csv(RES / f"fda_ranking{SFX}.tsv", sep="\t", usecols=["group_id", "cmap_name", "rev"]),
                         pd.read_csv(RES / f"fda_ranking_strata{SFX}.tsv", sep="\t", usecols=["stratum", "cmap_name", "rev"])
                         .rename(columns={"stratum": "group_id"})])
        sel = S["rv"].merge(fda, on=["group_id", "cmap_name", "rev"])
        print(f"FDA ranking rows {len(fda)}, matched to a scored drug {len(sel)}", flush=True)
        assert len(sel) == len(fda), "FDA rows that do not match one drug_reversal row"
    only = [x for x in os.environ.get("ONLY_UNITS", "").split(",") if x]
    if only:
        sel = sel[sel.group_id.isin(only)]
    sel = sel.assign(net=sel.group_id.map(S["unit_net"]))
    G["S"] = S
    G["nets"] = {n: Net(n, edges, regtab) for n in sorted(sel.net.unique())}
    dis = []
    for f in (MR_DIR / "regulon_contributions.tsv.gz", Path(f"{MR_DIR}_strata") / "regulon_contributions.tsv.gz"):
        if f.exists():
            dis.append(pd.read_csv(f, sep="\t", usecols=["group_id", "regulator", "target", "z", "contrib_nes"])
                       .query("group_id in @sel.group_id.unique()", local_dict={"sel": sel}))
    G["dis"] = ({u: d.drop(columns="group_id").rename(columns={"z": "disease_z", "contrib_nes": "disease_contrib_nes"})
                 for u, d in pd.concat(dis).groupby("group_id")} if dis else {})
    print(f"rows {len(sel)} over {sel.group_id.nunique()} units, {len(G['nets'])} networks; disease side "
          f"{'joined' if dis else 'MISSING (run 50b)'}", flush=True)

    with ProcessPoolExecutor(N_CORES) as ex:                       # fork: workers inherit G
        chk = dict(ex.map(check_net, G["nets"]))
        bad = {k: v for k, v in chk.items() if v > 1e-4}
        print(f"check 1 (activity vs 79d .f32, {CHECK_N} profiles/net): max |diff| {max(chk.values()):.2e}", flush=True)
        assert not bad, f"activity does not reproduce 79d: {bad}"
        jobs = [(u, g, g.net.iloc[0], bool(q)) for u, g in sel.groupby("group_id")]
        out = list(ex.map(run_unit, jobs))
    R = pd.concat([o[0] for o in out])
    T = pd.concat([o[1] for o in out])
    d = (R.drop_duplicates(["unit_id", "pert_id"]).rev - R.drop_duplicates(["unit_id", "pert_id"]).rev_stored).abs()
    print(f"check 2 (sum of regulator terms vs stored rev, {len(d)} rows): max |diff| {d.max():.2e}", flush=True)
    assert d.max() <= 1e-3, "reversal does not reproduce 79e"
    R = R.drop(columns="rev_stored")
    for c in ["disease_nes", "drug_activity", "rev_contrib", "rev"]:
        R[c] = R[c].astype(float).map(lambda x: float(f"{x:.6g}"))
    if q:
        full = T.merge(R[["unit_id", "pert_id", "regulator", "disease_nes", "drug_activity", "rev_contrib", "rank_in_drug"]]
                       .rename(columns={"rev_contrib": "regulator_rev_contrib"}), on=["unit_id", "pert_id", "regulator"])
        for (u, d_), g in full.groupby(["unit_id", "cmap_name"]):
            f = OUT / f"query_{u}_{d_}.tsv"
            g.sort_values(["rank_in_drug", "rev_contrib"]).to_csv(f, sep="\t", index=False, float_format="%.6g")
            print(f"-> {f} ({len(g)} edges)")
        return
    R.to_csv(OUT / "regulator_contributions.tsv.gz", sep="\t", index=False, float_format="%.6g")
    T.to_csv(OUT / "target_contributions.tsv.gz", sep="\t", index=False, float_format="%.6g")
    print(f"regulator level {len(R):,} rows, target level {len(T):,} rows (top {TOP_REGS} regulators per drug) -> {OUT}")


if __name__ == "__main__":
    main()
