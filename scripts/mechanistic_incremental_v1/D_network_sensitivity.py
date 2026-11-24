#!/usr/bin/env python
"""Experiment D (config block D, amendment 1d): why the same regulator reads differently in different networks.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227 + d3763b7. Identical drug conditions are read through two networks of one tumour family, so the
network is the only thing that changes. Per network, exact per-edge aREA contributions (A_resolvable.Net) are reduced
per regulator into the part carried by targets the partner network shares and the part from targets unique to this
network; their sum is the regulator's activity.
Per (pair, window, regulator shared by both): J (target Jaccard), sign agreement on shared targets, agreement =
Spearman of activity across conditions, mean |delta|, share of |delta| from unique targets.
Test: Spearman across regulators of J vs agreement; bootstrap over regulators; null = J permuted within regulon-size
quintiles. Held-out: drug halves by ingredient. Flags per regulator. ESR1 worked example. Reference: BRCA_LUMA vs LUAD.
Outputs: results/lincs/mechanistic_incremental_v1/D_network_sensitivity/{D_regulators.tsv.gz, D_pairs.tsv,
  D_flags.tsv, D_worked_ESR1.tsv, D_worked_ESR1.md, run.json}
"""
import importlib.util
import json
import os
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("A", HERE / "A_resolvable.py")
A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(A)
CFG, IN, RES, PROC = A.CFG, A.IN, A.RES, A.PROC
OUT = A.WT / CFG["output_root"] / "D_network_sensitivity"
OUT.mkdir(parents=True, exist_ok=True)

FAMILIES = {"TCGA-BRCA": ["BRCA_LUMA", "BRCA_LUMB", "BRCA_BAS", "BRCA"],
            "TCGA-LGG": ["LGG_IDHMUT_CODEL", "LGG_IDHMUT_NONCODEL", "LGG_IDHWT"],
            "TCGA-ESCA": ["ESCA_AD", "ESCA_SCC"], "TCGA-CESC": ["CESC", "CESC_SCC"],
            "TCGA-PCPG": ["PCPG", "PCPG_PHEO"], "TCGA-SARC": ["SARC", "SARC_LMS"], "TCGA-UCEC": ["UCEC", "UCEC_EEC"]}
REFERENCE = ("TCGA-BRCA", "BRCA_LUMA", "LUAD")
ANTI_ESTROGENS = ["fulvestrant", "tamoxifen", "clomifene", "raloxifene"]


def colrank(X):
    return np.apply_along_axis(lambda v: pd.Series(v).rank().to_numpy(), 0, X)


def md_table(df):
    """Markdown table without the optional tabulate dependency."""
    fmt = lambda v: f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v)
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns) + "\n"
    return head + "\n".join("| " + " | ".join(fmt(v) for v in r) + " |" for r in df.itertuples(index=False))


def family_conditions(cohort):
    groups = pd.read_csv(RES / "step1/disease_groups_toil.tsv", sep="\t")
    gids = set(groups[groups.cohort == cohort].group_id)
    lwg = pd.read_csv(IN / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    lws = pd.read_csv(IN / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t")
    lws["in_key"] = lws.in_key.astype(str).str.lower().eq("true")
    lines = set(lwg[lwg.group_id.isin(gids) & lwg.in_key].ModelID) | set(lws[lws.parent.isin(gids) & lws.in_key].ModelID)
    return A.idx[A.idx.ModelID.isin(lines) & A.idx.window.isin(("early", "late"))]


def work(task):
    cohort, nets, pairs = task
    P = np.memmap(PROC / "drug_condition_profiles_specific.f32", dtype=np.float32, mode="r", shape=(len(A.idx), A.NG))
    edges = A.cached("edges", lambda: pd.read_csv(IN / "data/processed/step5/regulons_aracne_project_v2_tf_edges.tsv.gz", sep="\t"))
    cond = family_conditions(cohort)
    cols = cond.col.values
    if len(cols) == 0:      # no LINCS condition in any in-key line of this family: not evaluable, reported in run.json
        return pd.DataFrame(), pd.DataFrame(), 0.0, {"cohort": cohort, "n_conditions": 0}
    N = {n: A.Net(n, edges) for n in nets}
    # regulon target sets and modes, by network
    reg_t = {n: {r: dict(zip(g.target, g.tfmode)) for r, g in edges[edges.net == n].groupby("regulator")} for n in nets}
    partners = {n: [b if a == n else a for a, b in pairs if n in (a, b)] for n in nets}
    SU = {}   # (net, partner) -> (shared, unique) arrays conditions x regs(net)
    for n in nets:
        Nn = N[n]
        masks = {pp: np.array([t in reg_t[pp].get(Nn.regs[k], {}) for k, t in zip(Nn.reg, Nn.target)]) for pp in partners[n]}
        acc = {pp: (np.zeros((len(cols), len(Nn.regs))), np.zeros((len(cols), len(Nn.regs)))) for pp in partners[n]}
        for s in range(0, len(cols), 100):
            C = Nn.contrib(P, cols[s:s + 100])
            for pp, m in masks.items():
                acc[pp][0][s:s + 100] = np.add.reduceat(C * m, Nn.starts, axis=1)
                acc[pp][1][s:s + 100] = np.add.reduceat(C * ~m, Nn.starts, axis=1)
        for pp in partners[n]:
            SU[(n, pp)] = acc[pp]
    # exactness: shared + unique must equal the stored activity (first 50 conditions, every network)
    exact_err = max(float(np.abs((SU[(nn, partners[nn][0])][0] + SU[(nn, partners[nn][0])][1])[:50]
                                 - np.asarray(N[nn].A[cols[:50]], dtype=np.float64)).max()) for nn in nets if partners[nn])
    fm = pd.read_csv(RES / "lincs/fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    ann = pd.read_csv(RES / "lincs/drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id").cmap_name
    ing = cond.pert_id.map(lambda p: fm.get(p, ann.get(p, p))).values
    half = np.array([(zlib.crc32(str(i).encode()) ^ A.SEED) % 2 for i in ing])
    rows, worked = [], []
    for a, b in pairs:
        Na, Nb = N[a], N[b]
        shared = [r for r in Na.regs if r in set(Nb.regs)]
        ia, ib = Na.regs.get_indexer(shared), Nb.regs.get_indexer(shared)
        sa, ua = SU[(a, b)]; sb, ub = SU[(b, a)]
        for w in ("early", "late"):
            m = (cond.window == w).values
            if m.sum() < 20:
                continue
            Aa = (sa + ua)[m][:, ia]; Ab = (sb + ub)[m][:, ib]
            ra, rb = colrank(Aa), colrank(Ab)
            def corr(X, Y):
                X = X - X.mean(0); Y = Y - Y.mean(0)
                d = np.sqrt((X ** 2).sum(0) * (Y ** 2).sum(0))
                return np.divide((X * Y).sum(0), d, out=np.full(X.shape[1], np.nan), where=d > 0)
            agree = corr(ra, rb)
            hh = half[m]
            agree_h = [corr(colrank(Aa[hh == k]), colrank(Ab[hh == k])) for k in (0, 1)]
            dsh = np.abs(sa[m][:, ia] - sb[m][:, ib]).sum(0); dun = np.abs(ua[m][:, ia] - ub[m][:, ib]).sum(0)
            for j, r in enumerate(shared):
                ta, tb = reg_t[a].get(r, {}), reg_t[b].get(r, {})
                inter = set(ta) & set(tb); uni = set(ta) | set(tb)
                rows.append(dict(cohort=cohort, net_a=a, net_b=b, window=w, regulator=r, n_conditions=int(m.sum()),
                                 size_a=len(ta), size_b=len(tb), jaccard=len(inter) / len(uni) if uni else np.nan,
                                 sign_agreement=float(np.mean([np.sign(ta[t]) == np.sign(tb[t]) for t in inter])) if inter else np.nan,
                                 unique_edge_share=1 - 2 * len(inter) / (len(ta) + len(tb)) if (len(ta) + len(tb)) else np.nan,
                                 agreement=agree[j], agreement_h0=agree_h[0][j], agreement_h1=agree_h[1][j],
                                 mean_abs_delta=float(np.abs(Aa[:, j] - Ab[:, j]).mean()),
                                 share_delta_unique=float(dun[j] / (dsh[j] + dun[j])) if (dsh[j] + dun[j]) > 0 else np.nan))
            if {a, b} == {"BRCA_LUMA", "BRCA_BAS"} and "ESR1" in shared:
                j = shared.index("ESR1")
                for dn in ANTI_ESTROGENS:
                    pm = m & (cond.pert_id.map(ann).values == dn)
                    if pm.any():
                        worked.append(dict(drug=dn, window=w, n_conditions=int(pm.sum()),
                                           esr1_activity_a=float((sa + ua)[pm][:, ia[j]].mean()),
                                           esr1_activity_b=float((sb + ub)[pm][:, ib[j]].mean()), net_a=a, net_b=b))
    return pd.DataFrame(rows), pd.DataFrame(worked), exact_err, {"cohort": cohort, "n_conditions": int(len(cols))}


def main():
    t0 = time.time()
    tasks = []
    for cohort, nets in FAMILIES.items():
        pairs = list(combinations(nets, 2))
        if cohort == REFERENCE[0]:
            nets = nets + [REFERENCE[2]]; pairs = pairs + [(REFERENCE[1], REFERENCE[2])]
        tasks.append((cohort, nets, pairs))
    with ProcessPoolExecutor(min(A.NC, len(tasks))) as ex:
        res = list(ex.map(work, tasks))
    R = pd.concat([r[0] for r in res], ignore_index=True); W = pd.concat([r[1] for r in res], ignore_index=True)
    exact_err = max(r[2] for r in res)
    coverage = [r[3] for r in res]
    print("conditions per family:", coverage, flush=True)
    print(f"exactness: max |shared + unique - stored activity| = {exact_err:.2e}", flush=True)
    if exact_err > 1e-3:
        raise SystemExit("decomposition is not exact; stopping")
    R["reference"] = (R.net_a == REFERENCE[1]) & (R.net_b == REFERENCE[2])
    R.to_csv(OUT / "D_regulators.tsv.gz", sep="\t", index=False, float_format="%.5g")
    W.to_csv(OUT / "D_worked_ESR1.tsv", sep="\t", index=False, float_format="%.5g")
    S = []
    for (a, b, w), d in R.groupby(["net_a", "net_b", "window"]):
        d = d.dropna(subset=["jaccard", "agreement"])
        if len(d) < 30:
            continue
        rho = spearmanr(d.jaccard, d.agreement).correlation
        rng = np.random.default_rng(A.rseed("D", a, b, w))
        q = pd.qcut((d.size_a + d.size_b) / 2, 5, labels=False, duplicates="drop").values
        null = []
        for _ in range(A.NP):
            jp = d.jaccard.values.copy()
            for k in np.unique(q):
                ii = np.flatnonzero(q == k); jp[ii] = rng.permutation(jp[ii])
            null.append(spearmanr(jp, d.agreement).correlation)
        nm = float(np.median(null))
        bs = []
        rb = np.random.default_rng(A.rseed("Dboot", a, b, w))
        J, G = d.jaccard.values, d.agreement.values
        for _ in range(A.NB):
            ii = rb.integers(len(d), size=len(d)); bs.append(spearmanr(J[ii], G[ii]).correlation)
        lo, hi = np.nanpercentile(bs, [2.5, 97.5])
        h = d.dropna(subset=["agreement_h0", "agreement_h1"])
        S.append(dict(net_a=a, net_b=b, window=w, reference=bool(d.reference.iloc[0]), n_regulators=len(d),
                      median_jaccard=float(d.jaccard.median()), median_sign_agreement=float(d.sign_agreement.median()),
                      median_agreement=float(d.agreement.median()), rho_jaccard_agreement=rho, null=nm,
                      ci_lo=float(lo), ci_hi=float(hi), resolved=bool(lo > nm or hi < nm),
                      p_perm=float((np.sum(np.array(null) >= rho) + 1) / (A.NP + 1)),
                      rho_sign_agreement=spearmanr(d.sign_agreement, d.agreement, nan_policy="omit").correlation,
                      halves_reliability=spearmanr(h.agreement_h0, h.agreement_h1).correlation,
                      rho_half0=spearmanr(h.jaccard, h.agreement_h0).correlation,
                      rho_half1=spearmanr(h.jaccard, h.agreement_h1).correlation,
                      median_share_delta_unique=float(d.share_delta_unique.median()),
                      median_unique_edge_share=float(d.unique_edge_share.median())))
    SUM = pd.DataFrame(S); SUM.to_csv(OUT / "D_pairs.tsv", sep="\t", index=False, float_format="%.4f")
    # flags (both halves), per pair x window
    F = R[~R.reference].copy()
    med = F.groupby(["net_a", "net_b"]).jaccard.transform("median")
    lowboth = (F.agreement_h0 < 0.3) & (F.agreement_h1 < 0.3)
    F["flag"] = np.select([(F.agreement_h0 >= 0.7) & (F.agreement_h1 >= 0.7), lowboth & (F.jaccard < med), lowboth & (F.jaccard >= med)],
                          ["context_robust", "rewired", "unexplained"], "intermediate")
    FL = F.groupby(["net_a", "net_b", "window", "flag"]).size().unstack(fill_value=0).reset_index()
    FL.to_csv(OUT / "D_flags.tsv", sep="\t", index=False)
    F[["cohort", "net_a", "net_b", "window", "regulator", "jaccard", "agreement", "flag"]].to_csv(OUT / "D_regulator_flags.tsv.gz", sep="\t", index=False, float_format="%.4g")
    esr = R[(R.regulator == "ESR1") & (R.net_a.isin(["BRCA_LUMA", "BRCA_BAS"])) & (R.net_b.isin(["BRCA_LUMA", "BRCA_BAS"]))]
    md = ["# ESR1 through the luminal A and basal breast networks (identical conditions)\n",
          esr[["window", "size_a", "size_b", "jaccard", "sign_agreement", "agreement", "agreement_h0", "agreement_h1",
               "share_delta_unique", "unique_edge_share"]].pipe(md_table), "\n\nAnti-estrogens, mean ESR1 activity:\n",
          md_table(W) if len(W) else "none measured"]
    (OUT / "D_worked_ESR1.md").write_text("\n".join(md) + "\n")
    run = dict(stage="D_network_sensitivity", protocol="5cc2227+d3763b7", regulator_rows=len(R), pairs=int(SUM.shape[0]),
               decomposition_max_error=exact_err, family_conditions=coverage,
               runtime_s=round(time.time() - t0, 1), finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 250)
    print(SUM.to_string(index=False)); print(FL.to_string(index=False)); print("\n".join(md))


if __name__ == "__main__":
    sys.exit(main())
