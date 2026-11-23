#!/usr/bin/env python
"""Experiment C (config block C, amendment 1c): can drug targets be recovered from the drug's regulator response?

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python
Protocol 5cc2227 + cbd2fb5. Per (unit, drug, window), the production window activity over the network's scored
regulators. For each target g of the drug, by role (primary / secondary with measured pChEMBL / annotated):
  exact       percentile of g by |activity| among the drug's regulators (g must be a scored regulator)
  exact_sgn   primary targets with a known orientation: percentile of orient_sign x activity
  down        mean |activity| percentile over H_g = 1-hop OmniPath targets of g among scored regulators (|H| >= 3)
  down_sgn    oriented primary targets: mean percentile of orient_sign x edge_sign x activity over H_g
  down2, down_signor   2-hop and SIGNOR-only H, descriptive
Null: within each (unit, window, role) the drugs' target lists are permuted among the drugs (n_perm), computed as a
matrix M[drug, gene] so each permutation is a row shuffle; cells permuted independently, then pooled.
Uncertainty: cluster bootstrap over target genes. Resolved = interval excludes the permutation median.
Outputs: results/lincs/mechanistic_incremental_v1/C_target_recovery/{C_pairs.tsv.gz, C_summary.tsv, C_potency.tsv,
  C_coverage.tsv, run.json}
"""
import importlib.util
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("A", HERE / "A_resolvable.py")
A = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(A)
CFG, IN, RES = A.CFG, A.IN, A.RES
OUT = A.WT / CFG["output_root"] / "C_target_recovery"
OUT.mkdir(parents=True, exist_ok=True)
MIN_H = 3


def genes_of(s, numeric=False):
    out = {}
    if not isinstance(s, str):
        return out
    for t in re.split(r"[|;,]", s):
        k, _, v = t.strip().partition(":")
        if not k:
            continue
        if numeric:
            m = re.match(r"^([0-9.]+)(\*?)$", v)
            if m:
                out[k] = (float(m.group(1)), m.group(2) == "*")
        else:
            out[k] = None
    return out


def adjacency(path):
    o = pd.read_csv(path, sep="\t", usecols=["source_genesymbol", "target_genesymbol", "sign"])
    adj = defaultdict(dict)
    for s, t, g in zip(o.source_genesymbol, o.target_genesymbol, o.sign):
        if s != t:
            adj[s][t] = int(g)
    return adj


def main():
    t0 = time.time()
    sp = pd.concat([pd.read_csv(RES / f"lincs/{f}", sep="\t", low_memory=False,
                                usecols=lambda c: c in (k, "pert_id", "evidence_targets", "target_source", "secondary_targets"))
                    .rename(columns={k: "unit"}) for f, k in (("drug_reversal_split_v2_tf.tsv", "group_id"),
                                                             ("stratum_reversal_split_v2_tf.tsv", "stratum"))])
    sp = sp.drop_duplicates("pert_id")
    ann = pd.read_csv(RES / "lincs/drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    targets = {}
    for p, e, src, sec in zip(sp.pert_id, sp.evidence_targets, sp.target_source, sp.secondary_targets):
        prim = set(genes_of(e)) if src in ("chembl_mechanism", "confirmed_annotation") else set()
        secd = {g: v for g, v in genes_of(sec, numeric=True).items() if g not in prim}
        hub = set(genes_of(ann.hub_target.get(p))) - prim - set(secd)
        targets[p] = {"primary": {g: None for g in prim}, "secondary": secd, "annotated": {g: None for g in hub}}
    keep = {p for p, d in targets.items() if any(d.values())}
    L = pd.read_csv(RES / "lincs/layered_score_v2_v2_tf.tsv", sep="\t", low_memory=False)
    orient = {(u, p): {"lower": -1, "raise": 1}.get(o, 0) for u, p, o in zip(L.unit, L.pert_id, L.layered_orient)}
    omni = adjacency(IN / "data/raw/omnipath/omnipath_signed_directed.tsv.gz")
    signor = adjacency(IN / "data/raw/omnipath/signor_signed_directed.tsv.gz")
    # median recorded in-key dose per (unit, drug, window)
    idx = A.idx[A.idx.pert_dose > 0]
    lwg = pd.read_csv(IN / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    lws = pd.read_csv(IN / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t").rename(columns={"stratum": "group_id"})
    lws["in_key"] = lws.in_key.astype(str).str.lower().eq("true")
    lw = pd.concat([lwg[lwg.in_key][["group_id", "ModelID"]], lws[lws.in_key][["group_id", "ModelID"]]])
    dose = (idx.merge(lw, on="ModelID").groupby(["group_id", "pert_id", "window"]).pert_dose.median())

    pairs, cells = [], []        # cells: per (level, unit, window, role, metric) numerator/denominator under perms
    rng_master = np.random.default_rng(A.SEED)
    for lv, f in (("group_id", "regulator_contributions_by_window_all.tsv.gz"), ("stratum", "regulator_contributions_by_window_all_strata.tsv.gz")):
        X = pd.concat([ch[ch.pert_id.isin(keep)] for ch in pd.read_csv(
            RES / "lincs/regulon_contributions_v2_tf" / f, sep="\t", chunksize=3_000_000,
            usecols=["cohort", "pert_id", "window", "regulator", "drug_activity"])], ignore_index=True)
        for (u, w), cw in X.groupby(["cohort", "window"], sort=False):
            regs = pd.Index(sorted(cw.regulator.unique())); rpos = {r: i for i, r in enumerate(regs)}
            drugs = sorted(cw.pert_id.unique())
            Amat = cw.pivot_table(index="pert_id", columns="regulator", values="drug_activity").reindex(index=drugs, columns=regs)
            act = Amat.to_numpy(float)
            pabs = pd.DataFrame(act).abs().rank(axis=1, pct=True).to_numpy()
            pup = pd.DataFrame(act).rank(axis=1, pct=True).to_numpy()            # percentile of +activity
            for role in ("primary", "secondary", "annotated"):
                tl = [list(targets[p][role]) for p in drugs]
                pool = sorted({g for t in tl for g in t})
                if not pool:
                    continue
                gi = {g: k for k, g in enumerate(pool)}
                n, G = len(drugs), len(pool)
                Mx = np.full((n, G), np.nan); Md = np.full((n, G), np.nan)
                Md2 = np.full((n, G), np.nan); Mds = np.full((n, G), np.nan)
                Mxs = np.full((n, G), np.nan); Mdsg = np.full((n, G), np.nan)
                ori = np.array([orient.get((u, p), 0) for p in drugs])
                for g, k in gi.items():
                    if g in rpos:
                        Mx[:, k] = pabs[:, rpos[g]]
                        if role == "primary":
                            Mxs[:, k] = np.where(ori > 0, pup[:, rpos[g]], np.where(ori < 0, 1 - pup[:, rpos[g]], np.nan))
                    H = [h for h in omni.get(g, {}) if h in rpos and h != g]
                    if len(H) >= MIN_H:
                        hi = [rpos[h] for h in H]
                        Md[:, k] = np.nanmean(pabs[:, hi], axis=1)
                        if role == "primary":
                            es = np.array([omni[g][h] for h in H])        # +1 stimulation, -1 inhibition
                            sg = ori[:, None] * es[None, :]               # expected direction per drug x regulator
                            v = np.where(sg > 0, pup[:, hi], np.where(sg < 0, 1 - pup[:, hi], np.nan))
                            Mdsg[:, k] = np.where(ori != 0, np.nanmean(v, axis=1), np.nan)
                    H2 = set(H)
                    for h in list(omni.get(g, {})):
                        H2 |= {x for x in omni.get(h, {}) if x in rpos and x != g}
                    H2 = [rpos[h] for h in sorted(H2)[:200]]   # sorted: a set's order is not reproducible
                    if len(H2) >= MIN_H:
                        Md2[:, k] = np.nanmean(pabs[:, H2], axis=1)
                    Hs = [rpos[h] for h in signor.get(g, {}) if h in rpos and h != g]
                    if len(Hs) >= MIN_H:
                        Mds[:, k] = np.nanmean(pabs[:, Hs], axis=1)
                B = np.zeros((n, G))
                for i, t in enumerate(tl):
                    B[i, [gi[g] for g in t]] = 1
                # observed pairs
                for i, p in enumerate(drugs):
                    for g in tl[i]:
                        k = gi[g]
                        pc = targets[p][role].get(g) if role == "secondary" else None
                        pairs.append(dict(level=lv, unit=u, window=w, pert_id=p, role=role, target=g,
                                          is_regulator=g in rpos, exact=Mx[i, k], exact_sgn=Mxs[i, k], down=Md[i, k],
                                          down_sgn=Mdsg[i, k], down2=Md2[i, k], down_signor=Mds[i, k],
                                          pchembl=pc[0] if pc else np.nan, single_measurement=pc[1] if pc else np.nan,
                                          dose_uM=dose.get((u, p, w), np.nan)))
                # permutation null: shuffle which drug carries which target list
                rng = np.random.default_rng(A.rseed("C", lv, u, w, role))
                for name, M in (("exact", Mx), ("exact_sgn", Mxs), ("down", Md), ("down_sgn", Mdsg), ("down2", Md2), ("down_signor", Mds)):
                    ok = ~np.isnan(M)
                    if not (ok & (B > 0)).any():
                        continue
                    Mz = np.where(ok, M, 0.0)
                    num, den = np.empty(A.NP), np.empty(A.NP)
                    for kperm in range(A.NP):
                        Bp = B[rng.permutation(n)]
                        num[kperm] = (Mz * Bp).sum(); den[kperm] = (ok * Bp).sum()
                    cells.append(dict(level=lv, unit=u, window=w, role=role, metric=name, num=num, den=den))
        print(f"{lv}: {time.time() - t0:.0f}s", flush=True)
    P = pd.DataFrame(pairs)
    P["cluster"] = P.target
    P.to_csv(OUT / "C_pairs.tsv.gz", sep="\t", index=False, float_format="%.5g")
    S = []
    for (lv, w, role), d in P.groupby(["level", "window", "role"]):
        for name in ("exact", "exact_sgn", "down", "down_sgn", "down2", "down_signor"):
            x = d.dropna(subset=[name])
            if not len(x):
                continue
            cs = [c for c in cells if c["level"] == lv and c["window"] == w and c["role"] == role and c["metric"] == name]
            num = np.sum([c["num"] for c in cs], axis=0); den = np.sum([c["den"] for c in cs], axis=0)
            nullv = np.divide(num, den, out=np.full(len(num), np.nan), where=den > 0)
            nm = float(np.nanmedian(nullv))
            val = float(x[name].mean())
            lo, hi = A.boot_ci(x[name].values, x.cluster.values, stat=np.mean)
            S.append(dict(level=lv, window=w, role=role, metric=name, n_pairs=len(x), n_targets=x.target.nunique(),
                          n_drugs=x.pert_id.nunique(), value=val, null=nm, lift=val / nm if nm else np.nan,
                          ci_lo=lo, ci_hi=hi, p_perm=float((np.sum(nullv >= val) + 1) / (len(nullv) + 1)),
                          top5pct=float((x[name] >= 0.95).mean()) if name.startswith("exact") else np.nan,
                          resolved=bool(lo > nm) if np.isfinite(lo) else False,
                          feasibility_limited=bool(x.target.nunique() < 30)))
    SUM = pd.DataFrame(S); SUM.to_csv(OUT / "C_summary.tsv", sep="\t", index=False, float_format="%.4f")
    # potency within secondary targets
    sec = P[(P.role == "secondary") & P.pchembl.notna() & P.dose_uM.notna()].copy()
    sec["log_ratio"] = np.log10(sec.dose_uM / 10 ** (6 - sec.pchembl))
    PO = []
    for (lv, w), d in sec.groupby(["level", "window"]):
        for name in ("down", "exact"):
            x = d.dropna(subset=[name])
            if len(x) < 20:
                continue
            r = spearmanr(x.log_ratio, x[name]).correlation
            lo, hi = A.boot_ci(np.arange(len(x)), x.pert_id.values,
                               stat=lambda ii, x=x, name=name: spearmanr(x.log_ratio.values[ii.astype(int)], x[name].values[ii.astype(int)]).correlation)
            PO.append(dict(level=lv, window=w, metric=name, n_pairs=len(x), n_drugs=x.pert_id.nunique(),
                           spearman=r, ci_lo=lo, ci_hi=hi, resolved=bool(lo > 0 or hi < 0)))
    pd.DataFrame(PO).to_csv(OUT / "C_potency.tsv", sep="\t", index=False, float_format="%.4f")
    COV = P.groupby(["level", "window", "role"]).agg(pairs=("target", "size"), targets=("target", "nunique"),
                                                      drugs=("pert_id", "nunique"), regulator_targets=("is_regulator", "sum"),
                                                      with_downstream=("down", lambda s: int(s.notna().sum()))).reset_index()
    COV.to_csv(OUT / "C_coverage.tsv", sep="\t", index=False)
    run = dict(stage="C_target_recovery", protocol="5cc2227+cbd2fb5", pairs=len(P), runtime_s=round(time.time() - t0, 1),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (OUT / "run.json").write_text(json.dumps(run, indent=1))
    pd.set_option("display.width", 250)
    print(COV.to_string(index=False)); print(SUM.to_string(index=False)); print(pd.DataFrame(PO).to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
