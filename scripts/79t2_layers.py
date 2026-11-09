#!/usr/bin/env python
"""D121: three-layer drug-disease score (runs in `celligner`; rule in docs/DECISIONS.md D121, fixed before running).

Layer 1  target role per group x gene: D120 components D / A / O (79t functions).
Layer 2  drug-target engagement: units = ChEMBL mechanism target records (family / complex collapsed) among the drug's
         D105b main targets, else single genes; unit role T_u(w) = weighted mean of the unit's oriented components;
         binding r_u = 10^(p_u - p_best) (unmeasured -> 1, flagged); engaged E_u = 50 + r_u^gamma (T_u - 50);
         chance-corrected S = prod_u G_u(50 + (E_best - 50) / r_u^gamma) (G_u: CDF of random units of the same size and
         direction in the group). Label-free acceptance: median within-group Spearman(S, n_units) in [-0.1, 0.1], else
         fallbacks best E, mean E, D120 T.
Layer 3  drug level = mean S over gate-passing compounds; score = (1 - v) pct(S) + v pct(spec_mech_pct_pregate);
         w (simplex step 0.1 + equal), gamma {0, .5, 1, 1.5, 2}, v {0, .1, ..., 1} fitted by within-group conditional
         logit on Hub (beta profiled 0-20); held out: leave-one-cohort-out on Hub, ChEMBL minus Hub-approved pairs,
         breast controls. Adoption vs baseline (equal w, gamma 1, v 0) by D121's criteria. Reported only.
Outputs results/lincs/target_role/d121/.
"""
import importlib.util
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "target_role/d121"
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402

_s = importlib.util.spec_from_file_location("t79", ROOT / "scripts/79t_target_role.py"); T79 = importlib.util.module_from_spec(_s); _s.loader.exec_module(T79)
KU = T79.KU
SEED, N_RAND, N_BOOT, MARGIN = 0, 2000, 200, 0.01
GAMMAS = [0.0, 0.5, 1.0, 1.5, 2.0]
VS = [round(0.1 * i, 1) for i in range(11)]
BETAS = np.arange(0, 20.0001, 0.25)
EQ = (1 / 3, 1 / 3, 1 / 3)
WS = sorted({(a / 10, b / 10, round(1 - a / 10 - b / 10, 1)) for a in range(11) for b in range(11 - a)}) + [EQ]
POS_CTL, NEG_CTL = T79.POS_CTL, T79.NEG_CTL
ANTI_ER = {"FULVESTRANT", "TAMOXIFEN CITRATE", "RALOXIFENE HYDROCHLORIDE"}


# ------------------------------------------------------------------ layer 2 inputs
def chembl_units(names, main):
    """cmap_name -> list of units {genes, actions, tid}; gene potency per drug {gene: median pChEMBL}."""
    spec = importlib.util.spec_from_file_location("i2", ROOT / "scripts/79i2_target_presence.py")
    i2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i2)
    w = i2.w
    w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
    d2m, _, _, _, _ = w.chembl(sorted(names))
    mech = json.load(open(w.CACHE / "mechanisms.json")) + json.load(open(w.CACHE / "mechanisms_by_parent.json"))
    T = {t["target_chembl_id"]: t for t in json.load(open(w.CACHE / "targets.json")) + json.load(open(w.CACHE / "targets_by_parent.json"))}
    g1 = lambda t: {x["component_synonym"] for c in t.get("target_components") or [] for x in c.get("target_component_synonyms") or [] if x["syn_type"] == "GENE_SYMBOL"}
    one = {k: next(iter(g1(t))) for k, t in T.items() if t.get("target_type") == "SINGLE PROTEIN" and t.get("organism") == "Homo sapiens" and len(g1(t)) == 1}
    act = pd.DataFrame(json.load(open(w.CACHE / "activities.json")) + json.load(open(w.CACHE / "activities_by_parent.json")))
    act = act[act.data_validity_comment.isna() & (act.potential_duplicate.astype(str) != "1") & act.standard_type.isin(["Ki", "Kd", "IC50", "EC50"])
              & act.pchembl_value.notna()].copy()
    act["gene"] = act.target_chembl_id.map(one); act["pchembl_value"] = act.pchembl_value.astype(float)
    act = act.dropna(subset=["gene"])
    units, pot = {}, {}
    for d in names:
        ms = d2m.get(d, set()); M = set(main.get(d, []))
        byt = {}
        for m in mech:
            if m["direct_interaction"] != 1 or not m["target_chembl_id"] or not ({m["molecule_chembl_id"], m["parent_molecule_chembl_id"]} & ms):
                continue
            t = T.get(m["target_chembl_id"])
            if not t:
                continue
            try:
                gs = i2.target_genes(t) & M
            except SystemExit:
                continue
            if gs:
                u = byt.setdefault(m["target_chembl_id"], dict(genes=set(), actions=set(), tid=m["target_chembl_id"], ttype=t.get("target_type")))
                u["genes"] |= gs; u["actions"].add(m["action_type"])
        us = list(byt.values())
        cov = set().union(*[u["genes"] for u in us]) if us else set()
        us += [dict(genes={g}, actions=set(), tid="", ttype="annotation") for g in sorted(M - cov)]
        units[d] = us
        pot[d] = act[act.molecule_chembl_id.isin(ms)].groupby("gene").pchembl_value.median().to_dict()
    return units, pot


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    groups = sorted(sp.group_id.unique())
    _, Pd = T79.comp_dependency(groups)
    _, Pa, _ = T79.comp_alteration(groups)
    _, Po, _ = T79.comp_overactivity(groups)
    print("layer 1 components loaded", flush=True)

    fm = KU.fda_map()
    x = sp[sp.pert_id.isin(fm.index) & sp.passes_gate.map(truthy)].copy()
    x["fda_ingredient"] = x.pert_id.map(fm)
    x = x[x.target_source.isin(T79.PROT_SRC)]
    main_t = {r.cmap_name: [t.split(":")[0] for t in str(r.evidence_targets).split(";") if t and t != "nan"]
              for r in x.drop_duplicates("cmap_name").itertuples()}
    UNITS, POT = chembl_units(set(x.cmap_name), main_t)

    # unit table per drug (binding is drug-level, independent of group)
    ut = []
    for d, us in UNITS.items():
        p = POT[d]
        pu = [max([p[g] for g in u["genes"] if g in p], default=np.nan) for u in us]
        pbest = np.nanmax(pu) if np.isfinite(pu).any() else np.nan
        for u, pv in zip(us, pu):
            o, flag = T79.orient({a for a in u["actions"] if a})
            r = 10 ** (pv - pbest) if np.isfinite(pv) else 1.0
            ut.append(dict(cmap_name=d, tid=u["tid"], ttype=u["ttype"], genes="|".join(sorted(u["genes"])), k=len(u["genes"]),
                           orient=o, action_flag=flag, actions="|".join(sorted(a for a in u["actions"] if a)),
                           p_u=pv, r_u=r, binding_unmeasured=not np.isfinite(pv)))
    UT = pd.DataFrame(ut)

    # gene pools per group x orientation
    pools = {}
    for c in groups:
        gs = sorted(set(Pd.get(c, pd.Series(dtype=float)).index) | set(Pa.get(c, pd.Series(dtype=float)).index) | set(Po.get(c, pd.Series(dtype=float)).index))
        D = Pd[c].reindex(gs).values if c in Pd else np.full(len(gs), np.nan)
        A = Pa[c].reindex(gs).values if c in Pa else np.full(len(gs), np.nan)
        O = Po[c].reindex(gs).values if c in Po else np.full(len(gs), np.nan)
        lo = np.c_[D, A, O]; hi = np.c_[100 - D, np.full(len(gs), np.nan), 100 - O]
        pools[c] = dict(index={g: i for i, g in enumerate(gs)}, lower=lo, raise_=hi)

    # real units in each group: component means (independent of w)
    comp = x[["group_id", "pert_id", "cmap_name", "fda_ingredient"]].drop_duplicates()
    rows = []
    for r in comp.itertuples(index=False):
        P = pools[r.group_id]
        for u in UT[UT.cmap_name == r.cmap_name].itertuples():
            idx = [P["index"][g] for g in u.genes.split("|") if g in P["index"]]
            M = (P["raise_"] if u.orient == "raise" else P["lower"])[idx] if idx else np.full((1, 3), np.nan)
            with np.errstate(all="ignore"):
                cm = np.nanmean(M, axis=0) if len(idx) else np.full(3, np.nan)
            rows.append((r.group_id, r.pert_id, r.cmap_name, r.fda_ingredient, u.Index, u.k, u.orient, u.r_u, *cm))
    U = pd.DataFrame(rows, columns=["group_id", "pert_id", "cmap_name", "fda_ingredient", "unit", "k", "orient", "r_u", "cD", "cA", "cO"])
    U = U.sort_values(["group_id", "pert_id", "unit"]).reset_index(drop=True)
    print(f"units: {len(UT)} drug units ({int(UT.binding_unmeasured.sum())} unmeasured, multi-gene unmeasured "
          f"{int((UT.binding_unmeasured & (UT.k > 1)).sum())}); group x compound x unit rows {len(U)}", flush=True)

    # random units: component means for each group x size x orientation (fixed draws, common to all w)
    RAND = {}
    for (c, k, o), _ in U.groupby(["group_id", "k", "orient"]):
        P = pools[c]
        arr = P["raise_"] if o == "raise" else P["lower"]
        ok = np.where(~np.isnan(arr).all(1))[0]
        idx = ok[rng.integers(0, len(ok), (N_RAND, k))]
        with np.errstate(all="ignore"):
            RAND[(c, k, o)] = np.nanmean(arr[idx], axis=1)   # N_RAND x 3
    cmat = U[["cD", "cA", "cO"]].values
    keyU = list(zip(U.group_id, U.k, U.orient))
    grpU = {}
    for i, kk in enumerate(keyU):
        grpU.setdefault(kk, []).append(i)
    comp_id = (U.group_id + "\t" + U.pert_id).values
    starts = np.r_[0, np.nonzero(comp_id[1:] != comp_id[:-1])[0] + 1]
    CK = U.iloc[starts][["group_id", "pert_id", "cmap_name", "fda_ingredient"]].reset_index(drop=True)

    def wmean(C, w):
        wv = np.asarray(w)
        m = ~np.isnan(C)
        den = (m * wv).sum(1)
        with np.errstate(all="ignore"):
            return np.where(den > 0, np.nansum(np.nan_to_num(C) * wv, 1) / den, np.nan)

    def S_of(w, gamma, mode="chance", r=None):
        """compound-level S for the parameters; returns (S, n_valid_units, best unit row index)."""
        Tu = wmean(cmat, w)
        rr = (U.r_u.values if r is None else r) ** gamma
        E = 50 + rr * (Tu - 50)
        valid = ~np.isnan(E)
        Em = np.where(valid, E, -np.inf)
        best = np.maximum.reduceat(Em, starts)
        nval = np.add.reduceat(valid.astype(int), starts)
        bidx = np.array([s + int(np.argmax(Em[s:e])) for s, e in zip(starts, np.r_[starts[1:], len(U)])])
        if mode == "best":
            S = best
        elif mode == "mean":
            with np.errstate(all="ignore"):
                S = np.add.reduceat(np.where(valid, E, 0), starts) / nval
        else:
            bestu = np.repeat(best, np.diff(np.r_[starts, len(U)]))
            with np.errstate(all="ignore"):
                thr = 50 + (bestu - 50) / np.where(rr > 0, rr, np.nan)
            logG = np.zeros(len(U))
            for key, ii in grpU.items():
                rt = np.sort(np.where(np.isnan(t := wmean(RAND[key], w)), -np.inf, t))
                ii = np.asarray(ii)
                g = np.searchsorted(rt, thr[ii], side="right") / len(rt)
                with np.errstate(divide="ignore"):
                    logG[ii] = np.log(g)
            logG[~valid] = 0
            S = np.exp(np.add.reduceat(logG, starts)) * 100
        S = np.where(nval > 0, S, np.nan)
        return S, nval, bidx

    lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")
    d120 = pd.read_csv(RES / "target_role/compound_target_role.tsv", sep="\t")[["group_id", "pert_id", "T"]]

    def drug_level(S, nval):
        d = CK.assign(S=S, n_units=nval).groupby(["group_id", "fda_ingredient"]).agg(S=("S", "mean"), n_units=("n_units", "mean")).reset_index()
        return d.rename(columns={"group_id": "cohort"})

    def spearman_units(S, nval):
        d = drug_level(S, nval).dropna(subset=["S"])
        return float(np.nanmedian([g.S.corr(g.n_units, method="spearman") for _, g in d.groupby("cohort") if len(g) >= 10]))

    # ---------------- label-free acceptance and fallbacks (equal weights, gamma 1)
    modes = []
    for mode in ("chance", "best", "mean"):
        S, nv, _ = S_of(EQ, 1.0, mode)
        modes.append(dict(mode=mode, spearman_n_units=round(spearman_units(S, nv), 3), S=S, nval=nv))
    tmp = CK.merge(d120, on=["group_id", "pert_id"], how="left")
    _, nv, _ = S_of(EQ, 1.0, "best")
    modes.append(dict(mode="d120", spearman_n_units=round(spearman_units(tmp["T"].values, nv), 3), S=tmp["T"].values, nval=nv))
    chosen = next((m["mode"] for m in modes if abs(m["spearman_n_units"]) <= 0.10), "d120")
    print("acceptance:", [(m["mode"], m["spearman_n_units"]) for m in modes], "->", chosen, flush=True)

    # ---------------- layer 3 grid
    base = drug_level(*S_of(EQ, 1.0, "best")[:2]).merge(
        lab[["cohort", "fda_ingredient", "drug", "anticancer", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl", "eligible", "spec_mech_pct_pregate"]],
        on=["cohort", "fda_ingredient"], how="left")
    scored = base.S.notna() & base.spec_mech_pct_pregate.notna() & base.drug.notna()
    B = base[scored].reset_index(drop=True)
    pct = lambda v: pd.Series(v).groupby(B.cohort.values).rank(pct=True, method="average").values * 100
    spec_p = pct(B.spec_mech_pct_pregate.values)
    combos, SC = [], []
    for w in WS:
        for gm in GAMMAS:
            if chosen == "d120":
                S, nv = modes[3]["S"], modes[3]["nval"]
            else:
                S, nv, _ = S_of(w, gm, chosen)
            d = drug_level(S, nv).set_index(["cohort", "fda_ingredient"]).S.reindex(pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]])).values
            sp_ = np.nan_to_num(pct(d), nan=50.0)   # no weighted component under these weights -> neutral percentile
            for v in VS:
                combos.append((w, gm, v)); SC.append(((1 - v) * sp_ + v * spec_p).astype(np.float32))
    SC = np.vstack(SC)
    ib = combos.index((EQ, 1.0, 0.0))
    print(f"grid: {len(combos)} parameter sets x {SC.shape[1]} drug rows", flush=True)

    coh = B.cohort.values
    cohorts = sorted(set(coh))
    cidx = {c: np.where(coh == c)[0] for c in cohorts}

    def ll_by_cohort(mask, y):
        """(n_combos, n_beta, n_cohorts) conditional log-likelihood contributions."""
        L = np.zeros((len(combos), len(BETAS), len(cohorts)), np.float32)
        for j, c in enumerate(cohorts):
            ii = cidx[c][mask[cidx[c]]]
            if len(ii) == 0 or y[ii].sum() == 0 or (~y[ii]).sum() == 0:
                continue
            X = SC[:, ii].astype(np.float64) / 100
            sp = X[:, y[ii]].sum(1)
            for b, beta in enumerate(BETAS):
                z = beta * X
                m = z.max(1, keepdims=True)
                lse = (m[:, 0] + np.log(np.exp(z - m).sum(1)))
                L[:, b, j] = beta * sp - y[ii].sum() * lse
        return L

    def best_combo(L, wts):
        tot = np.tensordot(L, wts, axes=([2], [0]))
        assert np.isfinite(tot).all(), "non-finite log-likelihood"
        i, b = np.unravel_index(np.argmax(tot), tot.shape)
        return i, BETAS[b]

    def auc_of(score, mask, ycol):
        D = B[mask].assign(_s=score[mask])
        return KU.pooled_auc(D, ycol, "_s")

    pops = {
        "hub": (B.anticancer.astype(bool) & ~B.excl_hub.astype(bool)).values,
        "chembl_minus_hub": (B.anticancer.astype(bool) & ~B.excl_chembl.astype(bool) & ~B.lab_hub.astype(bool)).values,
        "chembl": (B.anticancer.astype(bool) & ~B.excl_chembl.astype(bool)).values,
        "hub_minus_chembl": (B.anticancer.astype(bool) & ~B.excl_hub.astype(bool) & ~B.lab_chembl.astype(bool)).values,
        "pooled": (B.anticancer.astype(bool) & ~(B.excl_hub.astype(bool) | B.excl_chembl.astype(bool))).values}
    Y = {"hub": B.lab_hub.astype(bool).values, "chembl": B.lab_chembl.astype(bool).values}
    Y["pooled"] = Y["hub"] | Y["chembl"]
    Lh = ll_by_cohort(pops["hub"], Y["hub"])
    has = np.array([(Y["hub"][cidx[c]] & pops["hub"][cidx[c]]).any() for c in cohorts])
    i_hub, beta_hub = best_combo(Lh, np.ones(len(cohorts)))
    print("Hub fit:", combos[i_hub], "beta", beta_hub, flush=True)

    # (a) leave-one-cohort-out on Hub
    held = np.full(SC.shape[1], np.nan); held_choice = []
    for j, c in enumerate(cohorts):
        if not has[j]:
            continue
        wts = np.ones(len(cohorts)); wts[j] = 0
        i, _ = best_combo(Lh, wts)
        held[cidx[c]] = SC[i, cidx[c]]; held_choice.append((c, combos[i]))
    m_loco = pops["hub"] & ~np.isnan(held)
    res = []
    def add(test, rows, label, score, mask, ycol):
        a = auc_of(score, mask, ycol); res.append(dict(test=test, rows=rows, score=label, **a)); return a["auc"]
    a_loco = add("(a) Hub leave-one-cohort-out", "hub", "fitted (held-out choice per cohort)", held, m_loco, "lab_hub")
    a_loco_b = add("(a) Hub leave-one-cohort-out", "hub", "baseline", SC[ib].astype(float), m_loco, "lab_hub")
    # (b) ChEMBL minus Hub-approved
    mb = pops["chembl_minus_hub"]
    a_ch = add("(b) ChEMBL minus Hub-approved", "chembl_minus_hub", "Hub-fitted", SC[i_hub].astype(float), mb, "lab_chembl")
    a_ch_b = add("(b) ChEMBL minus Hub-approved", "chembl_minus_hub", "baseline", SC[ib].astype(float), mb, "lab_chembl")
    # (c) breast controls
    def ctl(i):
        s = pd.Series(SC[i], index=pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]]))
        pc = s.groupby(level=0).rank(pct=True) * 100
        pos = [pc.get((g, d), np.nan) for g, d in POS_CTL if d in ANTI_ER]
        neg = [pc.get((g, d), np.nan) for g, d in NEG_CTL]
        return float(np.nanmedian(pos)), float(np.nanmedian(neg))
    c_fit, c_base = ctl(i_hub), ctl(ib)
    crit = dict(a=a_loco >= a_loco_b + MARGIN, b=a_ch >= a_ch_b + MARGIN, c=(c_fit[0] >= c_base[0]) and (c_fit[1] <= c_base[1]))
    adopted = all(crit.values())
    i_final = i_hub if adopted else ib

    # reported: binding effect, reversal effect, D120, reverse and pooled fits, bootstraps
    for g0 in (0.0, 1.0):
        i = combos.index((EQ, g0, 0.0))
        add("binding effect (equal w, v 0)", "hub", f"gamma {g0}", SC[i].astype(float), pops["hub"], "lab_hub")
        add("binding effect (equal w, v 0)", "chembl_minus_hub", f"gamma {g0}", SC[i].astype(float), mb, "lab_chembl")
    w_h, g_h, v_h = combos[i_hub]
    i_v0 = combos.index((w_h, g_h, 0.0))
    add("reversal effect (Hub-fitted w, gamma)", "chembl_minus_hub", "v 0", SC[i_v0].astype(float), mb, "lab_chembl")
    add("reversal effect (Hub-fitted w, gamma)", "chembl_minus_hub", f"v {v_h}", SC[i_hub].astype(float), mb, "lab_chembl")
    add("reversal alone", "chembl_minus_hub", "v 1 (spec only)", spec_p.astype(float), mb, "lab_chembl")
    d120d = drug_level(modes[3]["S"], modes[3]["nval"]).set_index(["cohort", "fda_ingredient"]).S.reindex(pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]])).values
    for p_, yc in (("hub", "lab_hub"), ("chembl_minus_hub", "lab_chembl")):
        add("D120 T on the same rows", p_, "D120 T", d120d, pops[p_], yc)
    for m in modes:
        dl = drug_level(m["S"], m["nval"]).set_index(["cohort", "fda_ingredient"]).S.reindex(pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]])).values
        for p_, yc in (("hub", "lab_hub"), ("chembl_minus_hub", "lab_chembl")):
            add(f"layer-2 mode {m['mode']} (equal w, gamma 1)", p_, f"spearman n_units {m['spearman_n_units']}", dl, pops[p_], yc)
    Lc = ll_by_cohort(pops["chembl"], Y["chembl"]); i_ch, _ = best_combo(Lc, np.ones(len(cohorts)))
    add("reverse fit (ChEMBL -> Hub minus ChEMBL-approved)", "hub_minus_chembl", f"ChEMBL-fitted {combos[i_ch]}", SC[i_ch].astype(float), pops["hub_minus_chembl"], "lab_hub")
    add("reverse fit (ChEMBL -> Hub minus ChEMBL-approved)", "hub_minus_chembl", "baseline", SC[ib].astype(float), pops["hub_minus_chembl"], "lab_hub")
    Lp = ll_by_cohort(pops["pooled"], Y["pooled"]); i_po, _ = best_combo(Lp, np.ones(len(cohorts)))
    B["lab_pooled"] = Y["pooled"]
    add("pooled fit (in-sample)", "pooled", f"pooled-fitted {combos[i_po]}", SC[i_po].astype(float), pops["pooled"], "lab_pooled")
    add("pooled fit (in-sample)", "pooled", "baseline", SC[ib].astype(float), pops["pooled"], "lab_pooled")
    boots = []
    hc = np.where(has)[0]
    for _ in range(N_BOOT):
        cnt = np.bincount(rng.choice(hc, len(hc)), minlength=len(cohorts)).astype(float)
        i, _ = best_combo(Lh, cnt); boots.append(combos[i])
    BT = pd.DataFrame([dict(wD=w[0], wA=w[1], wO=w[2], gamma=g, v=v) for w, g, v in boots])
    bs = BT.describe(percentiles=[0.025, 0.5, 0.975]).T[["2.5%", "50%", "97.5%"]]

    # occupancy at Cmax (reported only): equal w, v 0, r = occupancy where a Cmax exists
    cm = pd.read_csv(RES / "d109/cmax.tsv", sep="\t").dropna(subset=["cmax_uM"])
    cmx = dict(zip(cm.drug.str.lower(), cm.cmax_uM))
    dname = lab.drop_duplicates("fda_ingredient").set_index("fda_ingredient").drug.str.lower()
    U["cmax_uM"] = U.fda_ingredient.map(dname).map(cmx)
    U["p_u"] = UT.p_u.values[U.unit.values]
    occ = np.where(U.cmax_uM.notna() & U.p_u.notna(), (U.cmax_uM * 1e-6) / (U.cmax_uM * 1e-6 + 10 ** (-U.p_u.fillna(0))), np.nan)
    has_c = U.groupby(["group_id", "pert_id"]).cmax_uM.transform(lambda s: s.notna().all()).values
    r_occ = np.where(np.isnan(occ), 1.0, occ)
    S_o, nv_o, _ = S_of(EQ, 1.0, chosen if chosen != "d120" else "best", r=r_occ)
    S_r, nv_r, _ = S_of(EQ, 1.0, chosen if chosen != "d120" else "best")
    sub = CK.assign(hc=np.logical_and.reduceat(has_c, starts)).groupby(["group_id", "fda_ingredient"]).hc.all()
    msub = B.set_index(["cohort", "fda_ingredient"]).index.map(lambda k: bool(sub.get(k, False))).values.astype(bool)
    for nm, S_, nv_ in (("occupancy at Cmax", S_o, nv_o), ("relative binding", S_r, nv_r)):
        dl = drug_level(S_, nv_).set_index(["cohort", "fda_ingredient"]).S.reindex(pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]])).values
        for p_, yc in (("hub", "lab_hub"), ("chembl", "lab_chembl")):
            add("occupancy subset (drugs with Cmax; equal w, gamma 1)", p_, nm, dl, pops[p_] & msub, yc)
    # D123 (reported only): free-drug occupancy (Cmax x fu) in place of relative binding, drugs with it for every measured unit
    ou = RES / "d123/occupancy_units.tsv"
    if ou.exists():
        OU = pd.read_csv(ou, sep="\t")
        key = dict(zip(zip(OU.cmap_name, OU.genes), OU.occ_free))
        ug = UT.genes.values[U.unit.values]
        of = np.array([key.get((c, g), np.nan) for c, g in zip(U.cmap_name, ug)], dtype=float)
        meas = ~UT.binding_unmeasured.values[U.unit.values]
        okc = np.where(meas, ~np.isnan(of), True)
        has_f = np.logical_and.reduceat(okc, starts) & (np.add.reduceat((meas & ~np.isnan(of)).astype(int), starts) > 0)
        r_free = np.where(np.isnan(of), 1.0, of)
        S_f, nv_f, _ = S_of(EQ, 1.0, chosen if chosen != "d120" else "best", r=r_free)
        subf = CK.assign(hf=has_f).groupby(["group_id", "fda_ingredient"]).hf.all()
        mf = B.set_index(["cohort", "fda_ingredient"]).index.map(lambda k: bool(subf.get(k, False))).values.astype(bool)
        for nm, S_, nv_ in (("free occupancy (D123)", S_f, nv_f), ("relative binding", S_r, nv_r)):
            dl = drug_level(S_, nv_).set_index(["cohort", "fda_ingredient"]).S.reindex(pd.MultiIndex.from_frame(B[["cohort", "fda_ingredient"]])).values
            for p_, yc in (("hub", "lab_hub"), ("chembl", "lab_chembl")):
                add("free-occupancy subset (D123; equal w, gamma 1)", p_, nm, dl, pops[p_] & mf, yc)
    R = pd.DataFrame(res)
    R.to_csv(OUT / "eval.tsv", sep="\t", index=False)

    # decision
    fmt = lambda c: f"w(D,A,O)=({c[0][0]:.2f},{c[0][1]:.2f},{c[0][2]:.2f}) gamma {c[1]} v {c[2]}"
    dec = pd.DataFrame([
        dict(item="layer-2 mode (label-free acceptance)", value=chosen),
        dict(item="Hub-fitted parameters", value=f"{fmt(combos[i_hub])} beta {beta_hub}"),
        dict(item="(a) LOCO Hub AUC fitted vs baseline", value=f"{a_loco:.4f} vs {a_loco_b:.4f} -> {'met' if crit['a'] else 'not met'}"),
        dict(item="(b) ChEMBL minus Hub-approved AUC fitted vs baseline", value=f"{a_ch:.4f} vs {a_ch_b:.4f} -> {'met' if crit['b'] else 'not met'}"),
        dict(item="(c) breast controls fitted vs baseline (lum median, basal median)", value=f"({c_fit[0]:.1f}, {c_fit[1]:.1f}) vs ({c_base[0]:.1f}, {c_base[1]:.1f}) -> {'met' if crit['c'] else 'not met'}"),
        dict(item="DECISION", value=("Hub-fitted parameters adopted for layer 3" if adopted else "baseline stays (equal w, gamma 1, v 0)") + f": {fmt(combos[i_final])}"),
        dict(item="ChEMBL-fitted parameters (reverse, reported)", value=fmt(combos[i_ch])),
        dict(item="pooled-fitted parameters (reported)", value=fmt(combos[i_po])),
        dict(item="LOCO choices (cohort: parameters)", value="; ".join(f"{c}: {fmt(p)}" for c, p in held_choice)),
    ] + [dict(item=f"bootstrap {k} (2.5 / 50 / 97.5%)", value=f"{r['2.5%']:.2f} / {r['50%']:.2f} / {r['97.5%']:.2f}") for k, r in bs.iterrows()])
    dec.to_csv(OUT / "decision.tsv", sep="\t", index=False)
    BT.to_csv(OUT / "bootstrap.tsv", sep="\t", index=False)
    pd.DataFrame([{k: v for k, v in m.items() if k in ("mode", "spearman_n_units")} for m in modes]).to_csv(OUT / "acceptance.tsv", sep="\t", index=False)

    # decoded rows at the final parameters
    wF, gF, vF = combos[i_final]
    S, nv, bidx = S_of(wF, gF, chosen if chosen != "d120" else "best")
    Tu = wmean(cmat, wF)
    dec_c = CK.assign(S=S, n_units=nv).join(U.iloc[bidx][["unit", "k", "orient", "r_u", "cD", "cA", "cO"]].reset_index(drop=True))
    dec_c["T_u"] = Tu[bidx]; dec_c["E_u"] = 50 + dec_c.r_u ** gF * (dec_c.T_u - 50)
    dec_c = dec_c.join(UT[["genes", "ttype", "actions", "action_flag", "p_u", "binding_unmeasured"]], on="unit")
    dec_c.to_csv(OUT / "compound_decoded.tsv", sep="\t", index=False, float_format="%.3f")
    top = dec_c.sort_values("S", ascending=False).drop_duplicates(["group_id", "fda_ingredient"])
    F = B[["cohort", "fda_ingredient", "drug", "anticancer", "lab_hub", "lab_chembl", "eligible", "spec_mech_pct_pregate"]].copy()
    F["score"] = SC[i_final]
    F["spec_pct_part"] = vF * spec_p
    F["S_pct_part"] = F.score - F.spec_pct_part
    F["score_pct"] = F.groupby("cohort").score.rank(pct=True) * 100
    F = F.merge(top.rename(columns={"group_id": "cohort"}).drop(columns=["pert_id", "cmap_name"]), on=["cohort", "fda_ingredient"], how="left")
    F.to_csv(OUT / "drug_decoded.tsv", sep="\t", index=False, float_format="%.3f")
    UT.to_csv(OUT / "units.tsv", sep="\t", index=False, float_format="%.3f")
    cov = dict(drug_units=len(UT), measured_units=int((~UT.binding_unmeasured).sum()), unmeasured_units=int(UT.binding_unmeasured.sum()),
               multi_gene_units=int((UT.k > 1).sum()), multi_gene_unmeasured=int((UT.binding_unmeasured & (UT.k > 1)).sum()),
               raise_units=int((UT.orient == "raise").sum()), direction_assumed_units=int(UT.action_flag.eq("direction_assumed").sum()),
               compounds=int(CK.pert_id.nunique()), drug_rows_scored=int(len(B)), occupancy_rows=int(msub.sum()))
    json.dump(dict(seed=SEED, n_rand=N_RAND, n_boot=N_BOOT, margin=MARGIN, n_combos=len(combos), coverage=cov, layer2_mode=chosen,
                   adopted=adopted, final=fmt(combos[i_final])), open(OUT / "run.json", "w"), indent=1)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 140)
    print(json.dumps(cov))
    print(R[["test", "rows", "score", "auc", "n_cohorts", "n_pos", "n_neg"]].to_string(index=False))
    print(dec.to_string(index=False))


if __name__ == "__main__":
    main()
