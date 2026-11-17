#!/usr/bin/env python
"""D133 / D133.1: is the drug-side activity-dimensionality collapse (D131) an artefact of REGULON TARGET OVERLAP?

AUDIT AND QUANTIFICATION ONLY. Reads production inputs, writes only results/lincs/d133/. No score is changed,
nothing is re-ranked or re-tiered, no threshold is tuned, no production file is touched.

Everything measured here was fixed in docs/DECISIONS.md rows D133 and D133.1 before this script was run.

Arms per cohort: consensus (what the shipped score uses), aracne-only, collectri-only; the last two exist only
because the d126_act dump is consensus, so the 79e per-cohort drug aggregation is reproduced here and GATED on
reproducing the shipped consensus dump to float32 (D133.1(2)).
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUT = RES / "d133"
REGDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else None
COHORTS = ["PHEO-kinase-signaling", "PCPG-PHEO", "TCGA-LIHC", "UCEC-EEC", "TGCT-SEM-C1"]
N_BINS = 50          # D133(6): fixed in advance
N_PERM = 10000       # D133(8): fixed in advance
SEED = 0
CONS_SETS = {"aracne": "results/step6/aracne_project_v2_tf", "collectri": "results/step6/collectri_v2"}
t0 = time.time()


def log(*a):
    print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def net_of(cohort):
    return "COADREAD" if cohort in ("TCGA-COAD", "TCGA-READ") else cohort.split("-", 1)[1]


def n_eff(C):
    """D131's participation ratio. Reads as 'this sum's variance behaves as if it were n_eff independent terms',
    NEVER as 'there are only n_eff regulators' worth of information'."""
    p = C.shape[0]
    return float(p ** 2 / np.square(C).sum())


# ---------------------------------------------------------------- inputs (exactly the ones 79e reads)
log("loading production inputs")
idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
idx["col"] = np.arange(len(idx))
grp = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t").set_index("group_id")
keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
NETMAP = pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").set_index("group_id").net.to_dict()

dd_g = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False)
dd_s = pd.read_csv(RES / "stratum_drugs_v2.tsv", sep="\t", low_memory=False)
par = dd_s.drop_duplicates("group_id").set_index("group_id").parent_group.to_dict()
w_g = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
w_s = (pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t")
       .rename(columns={"stratum": "group_id"}))

mr_cache = {}
def mr_table(strata):
    if strata not in mr_cache:
        mr_cache[strata] = pd.concat(
            [pd.read_csv(ROOT / (d + ("_strata" if strata else "")) / "group_mr.tsv", sep="\t")
             .assign(regulator=lambda x, k=k: k + "|" + x.regulator) for k, d in CONS_SETS.items()],
            ignore_index=True)
    return mr_cache[strata]

act_cache = {}
def activity(net):
    """Set-tagged profiles x regulators, exactly as 79e's read()/load() builds it for REGULON=consensus."""
    if net not in act_cache:
        def one(d, n):
            regs = pd.Index((PROC / d / f"{n}.regulators.txt").read_text().split())
            A = np.fromfile(PROC / d / f"{n}.f32", dtype=np.float32).reshape(len(idx), len(regs)).astype(np.float64)
            return regs, A
        ra, Aa = one("drug_viper_specific_v2_tf", net)
        rc, Ac = one("drug_viper_specific_collectri", "ALL")
        act_cache[net] = (pd.Index(["aracne|" + r for r in ra] + ["collectri|" + r for r in rc]), np.hstack([Aa, Ac]))
    return act_cache[net]

# ---------------------------------------------------------------- regulons (79d's own LINCS restriction)
log("loading regulons")
genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
gpos = pd.Index(genes.gene_id)
reg_cache = {}
def regulons(tag):
    if tag not in reg_cache:
        r = pd.read_csv(REGDIR / f"{tag}.tsv.gz", sep="\t")
        r["gi"] = gpos.get_indexer(r.target)
        assert (r.gi >= 0).all()
        reg_cache[tag] = r
    return reg_cache[tag]

# D133.1(1): cross-check the curated set against the independent precursor dump, report rather than resolve
pre = pd.read_csv(PROC / "precursor_regulons/collectri.tsv.gz", sep="\t")
cc = regulons("collectri")
a = set(map(tuple, cc[["regulator", "target"]].values)); b = set(map(tuple, pre[["regulator", "target"]].values))
collectri_check = dict(rds_pairs=len(a), precursor_pairs=len(b), shared=len(a & b),
                       rds_only=len(a - b), precursor_only=len(b - a))
log("collectri rds vs precursor:", collectri_check)

# ---------------------------------------------------------------- OmniPath (symbols on both sides, no mapping)
log("loading omnipath")
op = pd.read_csv(ROOT / "data/raw/omnipath/omnipath_signed_directed.tsv.gz", sep="\t")
log("omnipath columns:", list(op.columns)[:12], "rows", len(op))

rows, zero_tbl, arm_tbl = [], [], []

for g in COHORTS:
    log(f"=== {g} ===")
    strata = g in par and par[g] != g or g in dd_s.group_id.values and g not in dd_g.group_id.values
    strata = g in set(dd_s.group_id) and g not in set(dd_g.group_id)
    P = par.get(g, g)
    net = NETMAP.get(P, net_of(grp.loc[P, "cohort"]))
    dd = dd_s if strata else dd_g
    w = w_s if strata else w_g
    mr = mr_table(strata)
    gd = dd[dd.group_id == g]
    regs0, A = activity(net)
    k = keys.loc[P]
    km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
    if isinstance(k.diseases, str) and k.diseases:
        km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
    wg = w[w.group_id == g].set_index("ModelID").weight
    top_w = wg.reindex(km.index).fillna(0).max()
    by_drug = idx.groupby("pert_id")

    # --- 79e lines 243-266, reproduced
    D0 = np.empty((len(gd), len(regs0)), dtype=np.float64)
    for i, (d, src) in enumerate(zip(gd.pert_id, gd.source)):
        p = by_drug.get_group(d)
        if src == "direct":
            p = p[p.ModelID.isin(km.index)]
            ww = wg.reindex(p.ModelID).fillna(0).values
        else:
            p = p[~p.ModelID.isin(km.index)]
            ww = np.zeros(len(p))
        ww = ww if ww.sum() > 0 else np.ones(len(p))
        D0[i] = (ww[:, None] * A[p.col.values]).sum(0) / ww.sum()

    # --- 79e consensus()
    m = mr[(mr.group_id == g) & mr.regulator.isin(regs0)]
    base = m.regulator.str.split("|", n=1).str[1].values
    sets = m.regulator.str.split("|", n=1).str[0].values
    cols = regs0.get_indexer(m.regulator)
    ub, inv = np.unique(base, return_inverse=True)
    Dc = np.zeros((D0.shape[0], len(ub)))
    cnt = np.bincount(inv, minlength=len(ub))
    np.add.at(Dc.T, inv, D0[:, cols].T)
    Dcons = Dc / cnt
    regs = pd.Index(ub)

    # --- GATE (D133.1(2)): must reproduce the shipped dump exactly
    dump_regs = (PROC / "d126_act" / f"{g}.regulators.txt").read_text().split("\n")
    Ddump = np.fromfile(PROC / "d126_act" / f"{g}.f32", dtype=np.float32).reshape(len(gd), len(dump_regs))
    gate_regs = list(regs) == dump_regs
    gate_diff = float(np.abs(np.asarray(Dcons, dtype=np.float32) - Ddump).max())
    log(f"GATE regulators identical={gate_regs}  max|recon-dump|={gate_diff:.3e}  "
        f"({Dcons.shape[0]} drugs x {Dcons.shape[1]} consensus regulators)")
    assert gate_regs and gate_diff < 1e-4, "D133.1(2) gate failed: arm is void"

    # --- per-set column maps
    set_col = {s: {} for s in ("aracne", "collectri")}
    for s, bn, c in zip(sets, base, cols):
        set_col[s][bn] = c

    # --- regulon vectors per set, restricted to this cohort's consensus regulators
    def vecs(tag, names):
        """Dense genes x regulator vectors: u = tfmode * likelihood (L2-normalised per regulator), b = membership."""
        r = regulons(tag)
        r = r[r.regulator.isin(set(names))]
        pos = pd.Index(names).get_indexer(r.regulator)
        u = np.zeros((len(names), len(gpos)), dtype=np.float32)
        b = np.zeros((len(names), len(gpos)), dtype=np.float32)
        u[pos, r.gi.values] = (r.tfmode.values * r.likelihood.values).astype(np.float32)
        b[pos, r.gi.values] = 1.0
        nrm = np.sqrt((u.astype(np.float64) ** 2).sum(1))
        nrm[nrm == 0] = 1.0
        return (u / nrm[:, None].astype(np.float32)), b

    names = list(regs)
    Ua, Ba = vecs(f"aracne_{net}", [n for n in names if n in set_col["aracne"]])
    Uc, Bc = vecs("collectri", [n for n in names if n in set_col["collectri"]])
    ia = {n: i for i, n in enumerate([n for n in names if n in set_col["aracne"]])}
    ic = {n: i for i, n in enumerate([n for n in names if n in set_col["collectri"]])}
    # consensus weight vector: mean of the present sets' L2-normalised vectors (D133(2d))
    Ucons = np.zeros((len(names), len(gpos)), dtype=np.float32)
    Bcons = np.zeros((len(names), len(gpos)), dtype=np.float32)
    for j, n in enumerate(names):
        nset = 0
        if n in ia:
            Ucons[j] += Ua[ia[n]]; Bcons[j] += Ba[ia[n]]; nset += 1
        if n in ic:
            Ucons[j] += Uc[ic[n]]; Bcons[j] += Bc[ic[n]]; nset += 1
        Ucons[j] /= nset
    Bcons = (Bcons > 0).astype(np.float32)

    # ---------------------------------------------------------------- the three arms
    arms = {"consensus": (Dcons, Ucons, Bcons, names),
            "aracne": (D0[:, [set_col["aracne"][n] for n in names if n in ia]], Ua, Ba, [n for n in names if n in ia]),
            "collectri": (D0[:, [set_col["collectri"][n] for n in names if n in ic]], Uc, Bc,
                          [n for n in names if n in ic])}
    common = [n for n in names if n in ia and n in ic]
    ni = {n: j for j, n in enumerate(names)}
    arms["aracne_common"] = (D0[:, [set_col["aracne"][n] for n in common]],
                             Ua[[ia[n] for n in common]], Ba[[ia[n] for n in common]], common)
    arms["collectri_common"] = (D0[:, [set_col["collectri"][n] for n in common]],
                                Uc[[ic[n] for n in common]], Bc[[ic[n] for n in common]], common)
    arms["consensus_common"] = (Dcons[:, [ni[n] for n in common]],
                                Ucons[[ni[n] for n in common]], Bcons[[ni[n] for n in common]], common)

    rng = np.random.default_rng(SEED)
    for arm, (D, U, B, nm) in arms.items():
        sd = D.std(0, ddof=1)
        ok = sd > 0
        D = D[:, ok]; U = U[ok]; B = B[ok]
        nm = [n for n, o in zip(nm, ok) if o]
        p, nd = D.shape[1], D.shape[0]
        C = np.corrcoef(D - D.mean(0), rowvar=False)
        inter = np.rint(B @ B.T).astype(np.int32)
        sz = np.diag(inter).copy()
        union = sz[:, None] + sz[None, :] - inter
        J = inter / np.maximum(union, 1)
        K = (U @ U.T).astype(np.float64)
        d = np.sqrt(np.clip(np.diag(K), 1e-300, None))
        K = K / d[:, None] / d[None, :]
        iu = np.triu_indices(p, 1)
        c, kk, jj, it = C[iu], K[iu], J[iu], inter[iu]

        ne = n_eff(C)
        ne_iid = float(p / (1 + p / nd))
        rec = dict(cohort=g, arm=arm, net=net, level="stratum" if strata else "group_id", n_drugs=nd, p=p,
                   n_pairs=len(c), n_eff=ne, n_eff_iid=ne_iid, n_eff_ratio=ne / ne_iid,
                   mean_abs_offdiag_corr=float(np.abs(c).mean()), mean_offdiag_corr=float(c.mean()),
                   median_regulon_size=float(np.median(sz)))

        # ---- (5) zero-overlap pairs: a fact about the annotation, not a threshold
        z = it == 0
        rec["n_zero_overlap_pairs"] = int(z.sum())
        rec["frac_zero_overlap_pairs"] = float(z.mean())
        if z.sum():
            cz = c[z]
            rec.update(zero_mean_corr=float(cz.mean()), zero_mean_abs_corr=float(np.abs(cz).mean()),
                       zero_sd_corr=float(cz.std()))
            for q in (1, 5, 25, 50, 75, 95, 99):
                rec[f"zero_p{q}"] = float(np.percentile(cz, q))
        # top overlap decile for contrast
        hi = kk >= np.quantile(kk, 0.9)
        rec.update(top10pct_K_mean_corr=float(c[hi].mean()), top10pct_K_mean_abs_corr=float(np.abs(c[hi]).mean()),
                   top10pct_K_threshold=float(np.quantile(kk, 0.9)))

        # ---- (6) decomposition
        rec["R2_linear_K"] = float(np.corrcoef(c, kk)[0, 1] ** 2)
        rec["R2_linear_J"] = float(np.corrcoef(c, jj)[0, 1] ** 2)
        # nonparametric: 50 equal-count bins of K, fixed in advance
        order = np.argsort(kk, kind="stable")
        binid = np.empty(len(kk), dtype=np.int32)
        binid[order] = np.minimum((np.arange(len(kk)) * N_BINS) // len(kk), N_BINS - 1)
        bm = np.bincount(binid, weights=c, minlength=N_BINS) / np.bincount(binid, minlength=N_BINS)
        r = c - bm[binid]
        rec["R2_binned_K"] = float(1 - r.var() / c.var())
        Cres = np.eye(p)
        Cres[iu] = r
        Cres[(iu[1], iu[0])] = r
        rec["n_eff_residual"] = n_eff(Cres)
        # linear-fit version, reported alongside
        sl, ic_ = np.polyfit(kk, c, 1)
        rl = c - (sl * kk + ic_)
        Cl = np.eye(p); Cl[iu] = rl; Cl[(iu[1], iu[0])] = rl
        rec["n_eff_residual_linear"] = n_eff(Cl)

        # ---- (7) maximal pairwise-zero-overlap subset, greedy ascending degree, seed 0
        deg = (inter > 0).sum(1) - 1
        order2 = np.lexsort((np.array(nm), deg))
        chosen, blocked = [], np.zeros(p, bool)
        for i in order2:
            if not blocked[i]:
                chosen.append(i)
                blocked |= inter[i] > 0
        chosen = np.array(sorted(chosen))
        Cs = C[np.ix_(chosen, chosen)]
        rec.update(n_disjoint_subset=len(chosen), n_eff_disjoint_subset=n_eff(Cs),
                   n_eff_disjoint_iid=float(len(chosen) / (1 + len(chosen) / nd)),
                   mean_abs_corr_disjoint=float(np.abs(Cs[np.triu_indices(len(chosen), 1)]).mean()))

        # ---- (8) cascade check, matched on overlap bin
        sym = pd.Index(nm)
        src, tgt = "source_genesymbol", "target_genesymbol"   # gene symbols, as regulator names are
        e = op[op[src].isin(set(nm)) & op[tgt].isin(set(nm))]
        si, ti = sym.get_indexer(e[src]), sym.get_indexer(e[tgt])
        Adj = np.zeros((p, p), dtype=bool)
        Adj[si, ti] = True
        und = Adj | Adj.T
        two = (und.astype(np.float32) @ und.astype(np.float32)) > 0
        conn1 = und[iu]
        conn2 = (und | two)[iu]
        for lab, conn in (("edge", conn1), ("path2", conn2)):
            nconn = int(conn.sum())
            rec[f"cascade_{lab}_n_pairs"] = nconn
            if nconn < 20:
                continue
            nb = np.bincount(binid, minlength=N_BINS)
            nbc = np.bincount(binid[conn], minlength=N_BINS)
            sc = np.bincount(binid[conn], weights=r[conn], minlength=N_BINS)
            sa = np.bincount(binid, weights=r, minlength=N_BINS)
            use = (nbc > 0) & (nb - nbc > 0)
            mc = sc[use] / nbc[use]
            mu = (sa[use] - sc[use]) / (nb[use] - nbc[use])
            wgt = nbc[use] / nbc[use].sum()
            obs = float((wgt * (mc - mu)).sum())
            sca = np.bincount(binid[conn], weights=np.abs(r[conn]), minlength=N_BINS)
            saa = np.bincount(binid, weights=np.abs(r), minlength=N_BINS)
            obs_abs = float((wgt * (sca[use] / nbc[use] - (saa[use] - sca[use]) / (nb[use] - nbc[use]))).sum())
            # within-bin permutation of the connected label, vectorised per bin
            null = np.zeros(N_PERM)
            for b in np.flatnonzero(use):
                rb = r[binid == b]
                pick = rng.integers(0, len(rb), size=(N_PERM, int(nbc[b])))
                mcb = rb[pick].mean(1)
                tot = rb.sum()
                mub = (tot - mcb * nbc[b]) / (len(rb) - nbc[b])
                null += (nbc[b] / nbc[use].sum()) * (mcb - mub)
            rec[f"cascade_{lab}_resid_diff"] = obs
            rec[f"cascade_{lab}_resid_abs_diff"] = obs_abs
            rec[f"cascade_{lab}_p_perm"] = float((np.abs(null) >= abs(obs)).mean())
            rec[f"cascade_{lab}_mean_resid_connected"] = float((wgt * mc).sum())
            rec[f"cascade_{lab}_mean_resid_unconnected"] = float((wgt * mu).sum())
            rec[f"cascade_{lab}_mean_raw_corr_connected"] = float(c[conn].mean())
        # REPORTED, NOT PRE-REGISTERED (exploratory): the same comparison split by OmniPath edge sign, because a
        # signed directed edge is as often inhibitory as activating and the pooled test cannot see that.
        for sg, lab2 in ((1, "pos"), (-1, "neg")):
            es = e[e["sign"] == sg]
            if len(es) == 0:
                continue
            Ag = np.zeros((p, p), dtype=bool)
            Ag[sym.get_indexer(es[src]), sym.get_indexer(es[tgt])] = True
            cs = (Ag | Ag.T)[iu]
            if cs.sum() < 20:
                continue
            rec[f"xcascade_{lab2}_n_pairs"] = int(cs.sum())
            rec[f"xcascade_{lab2}_mean_raw_corr"] = float(c[cs].mean())
            rec[f"xcascade_{lab2}_mean_resid"] = float(r[cs].mean())
        rows.append(rec)
        log(f"  {arm:<18} p={p:<5} n_eff={ne:6.2f} (iid {ne_iid:6.1f}, ratio {ne/ne_iid:.4f})  "
            f"zero-overlap pairs {rec['n_zero_overlap_pairs']} ({100*rec['frac_zero_overlap_pairs']:.1f}%) "
            f"mean|r|={rec.get('zero_mean_abs_corr', float('nan')):.3f}  "
            f"R2_K={rec['R2_linear_K']:.3f} binned={rec['R2_binned_K']:.3f}  "
            f"n_eff_res={rec['n_eff_residual']:.2f}  disjoint p={rec['n_disjoint_subset']} "
            f"n_eff={rec['n_eff_disjoint_subset']:.2f}")

        if arm == "consensus":
            # the zero-overlap correlation distribution, written out in full for the report
            ov_bins = [(-1, 0), (0, 0.01), (0.01, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.5), (0.5, 1.01)]
            for lo, hi2 in ov_bins:
                s = (jj > lo) & (jj <= hi2) if lo >= 0 else (it == 0)
                if s.sum() == 0:
                    continue
                zero_tbl.append(dict(cohort=g, jaccard_lo=lo, jaccard_hi=hi2, n_pairs=int(s.sum()),
                                     mean_corr=float(c[s].mean()), mean_abs_corr=float(np.abs(c[s]).mean()),
                                     sd_corr=float(c[s].std())))

df = pd.DataFrame(rows)
OUT.mkdir(parents=True, exist_ok=True)
df.to_csv(OUT / "overlap_dimensionality.tsv", sep="\t", index=False)
pd.DataFrame(zero_tbl).to_csv(OUT / "corr_by_overlap_bin.tsv", sep="\t", index=False)
json.dump(dict(cohorts=COHORTS, n_bins=N_BINS, n_perm=N_PERM, seed=SEED, collectri_rds_vs_precursor=collectri_check,
               omnipath_cols=list(op.columns), runtime_s=round(time.time() - t0, 1)),
          open(OUT / "d133_run.json", "w"), indent=1)
log("written", OUT / "overlap_dimensionality.tsv")
