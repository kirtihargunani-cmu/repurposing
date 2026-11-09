#!/usr/bin/env python
"""D122: the layered score (D121 baseline) as the lead ranking, for the 41 groups and the 14 subgroups (runs in `celligner`).

Layer 1 (D120): D dependency, A alteration, O over-activity per unit x gene (groups: 79t functions; subgroups: the
subgroup's own in-key lines, TCGA tumors, MR tables and step 10 v2 contrasts; A's reference = median over the groups).
Layer 2 (D121): target units, direction, relative binding r = 10^(p_u - p_best) (unmeasured 1), engaged
E_u = 50 + r_u (T_u - 50) with T_u = mean of the unit's available components (equal weights), S = mean E_u over units.
Writes into the production split tables (in place, idempotent): layered_S, layered_pct (rank(max) / n x 100 over the
unit's gate-passing FDA compounds with an S), layered_n_units, layered_genes / _D / _A / _O / _T / _r (highest-E unit),
no_protein_target, minority_alteration, essential_target. Also results/lincs/layered_score_v2.tsv and lead_score_run.json.
"""
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402


def load(name, path):
    s = importlib.util.spec_from_file_location(name, ROOT / path); m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m


T79 = load("t79", "scripts/79t_target_role.py")
T2 = load("t2", "scripts/79t2_layers.py")
N79 = load("n79", "scripts/79n_dependency.py")
COLS = ["layered_S", "layered_pct", "layered_n_units", "layered_genes", "layered_D", "layered_A", "layered_O", "layered_T",
        "layered_r", "no_protein_target", "minority_alteration", "essential_target", "layered_P", "essential_unit", "dna_acting"]
ESSENTIAL = os.environ.get("ESSENTIAL", "0") == "1"   # D124: essential units and DNA-acting drugs scored on proliferation P
SANDBOX = os.environ.get("SANDBOX", "0") == "1"       # D124 test: write layered_score_v2_d124.tsv only, no split tables
PAM50_PROLIF = ["BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C"]   # Nielsen 2010


def proliferation(groups, strata):
    """D124: P per group (percentile among the 41 groups) and per subgroup (placed among the 41 group values)."""
    b = importlib.util.spec_from_file_location("b79", ROOT / "scripts/79b_target_gate.py"); b79 = importlib.util.module_from_spec(b); b.loader.exec_module(b79)
    idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t").query("primary_for_hgnc").drop_duplicates("symbol")
    ids = set(idm[idm.symbol.isin(PAM50_PROLIF)].gene_id)
    gx = b79.group_expression(ids, groups, 1.0)
    sp2 = pd.read_csv(ROOT / "results/step10_v2/stratum_patients.tsv", sep="\t")
    sx = b79.group_expression(ids, None, 1.0, sp2.rename(columns={"group_id": "parent", "stratum": "group_id"})[["group_id", "sample"]])
    ig = lambda x: np.log2(x.pivot_table(index="group_id", columns="gene_id", values="median_tpm") + 1).mean(1)
    G, S = ig(gx), ig(sx[sx.group_id.isin(strata)])
    ref = np.sort(G.values)
    Pg = (G.rank(method="max") / len(G) * 100).to_dict()
    Ps = {u: float((ref <= v).sum() / len(ref) * 100) for u, v in S.items()}
    return Pg, Ps, G.to_dict(), S.to_dict()


def dna_acting(names):
    """D124: cmap names whose ChEMBL mechanism target type is NUCLEIC-ACID."""
    spec = importlib.util.spec_from_file_location("i2b", ROOT / "scripts/79i2_target_presence.py")
    i2 = importlib.util.module_from_spec(spec); spec.loader.exec_module(i2)
    w = i2.w
    w.CACHE = w.RAW / "chembl" / w.RELEASE / "d105_all"
    d2m, _, _, _, _ = w.chembl(sorted(names))
    mech = json.load(open(w.CACHE / "mechanisms.json")) + json.load(open(w.CACHE / "mechanisms_by_parent.json"))
    T = {t["target_chembl_id"]: t for t in json.load(open(w.CACHE / "targets.json")) + json.load(open(w.CACHE / "targets_by_parent.json"))}
    out = set()
    for d in names:
        ms = d2m.get(d, set())
        for m in mech:
            if m["direct_interaction"] == 1 and m["target_chembl_id"] and ({m["molecule_chembl_id"], m["parent_molecule_chembl_id"]} & ms) \
                    and (T.get(m["target_chembl_id"]) or {}).get("target_type") == "NUCLEIC-ACID":
                out.add(d)
    return out
SPLITS = {"group_id": RES / "drug_reversal_split_v2_consensus.tsv", "stratum": RES / "stratum_reversal_split_v2_consensus.tsv"}


def layer1_strata(strata, F_groups):
    """subgroup components: {unit: Series} percentiles for D, A, O and raw alteration shares."""
    G = pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRGeneEffect.csv", index_col=0)
    G.columns = [c.split(" (")[0] for c in G.columns]
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv")
    G = G.loc[G.index.isin(set(models.ModelID[models.OncotreePrimaryDisease != "Non-Cancerous"]))]
    Dz = -G.sub(G.mean(0), axis=1).div(G.std(0).replace(0, np.nan), axis=1)
    W = pd.read_csv(ROOT / "results/celligner/stratum_line_weights_long_v2.tsv", sep="\t")
    W = {u: d.set_index("ModelID").weight for u, d in W[W.in_key & W.stratum.isin(strata)].groupby("stratum")}
    gd = N79.cohort_gene_dep(Dz, W)
    Pd = {u: T79.pct(v.dropna()) for u, v in gd.items()}
    sp = pd.read_csv(ROOT / "results/step10_v2/stratum_patients.tsv", sep="\t")
    alt, prof, _ = T79.alteration_pairs()
    F, _ = T79.alteration_shares(alt, prof, {u: set(d["sample"]) for u, d in sp[sp.stratum.isin(strata)].groupby("stratum")})
    ref = pd.DataFrame(F_groups)
    genes = sorted(set(ref.index) | set(F.index))
    F = F.reindex(genes).fillna(0); ref = ref.reindex(genes).fillna(0)
    med = ref.where(ref >= T79.MIN_F, 0).median(1)
    A = F.where(F >= T79.MIN_F, 0).sub(med, axis=0)
    Pa = {u: T79.pct(A[u]) for u in A}
    _, Po, _ = T79.comp_overactivity(strata, mr_dirs=("aracne_project_v2_tf_strata", "collectri_v2_strata"),
                                     summary="results/step10_v2/signature_summary.tsv", proc="data/processed/step10_v2")
    return Pd, Pa, Po, {u: F[u] for u in F}


def score(x, unit, Pd, Pa, Po, Fr, UT, ess, P=None):
    """x = split rows (FDA, gate-passing, protein target) -> per row S and decoding."""
    gidx = {}
    out = []
    byd = {d: g for d, g in UT.groupby("cmap_name")}
    for r in x[[unit, "pert_id", "cmap_name"]].itertuples(index=False):
        u, d = r[0], r.cmap_name
        if d not in byd:
            out.append(None); continue
        best, Es = None, []
        for t in byd[d].itertuples():
            genes = t.genes.split("|")
            if ESSENTIAL and P is not None and t.ess_frac > 0.5:   # D124: essential unit -> T = P (100 - P for raising)
                pv = P.get(u, np.nan)
                if not np.isfinite(pv):
                    continue
                T = pv if t.orient != "raise" else 100 - pv
                E = 50 + t.r_u * (T - 50)
                Es.append(E)
                if best is None or E > best[0]:
                    best = (E, t.genes, np.nan, np.nan, np.nan, T, t.r_u, np.nan, True)
                continue
            comp = []
            for PC, k in ((Pd, "D"), (Pa, "A"), (Po, "O")):
                if k == "A" and t.orient == "raise":
                    comp.append(np.nan); continue
                s = PC.get(u)
                vals = [s[g] for g in genes if s is not None and g in s.index and np.isfinite(s[g])]
                v = float(np.mean(vals)) if vals else np.nan
                comp.append(100 - v if (t.orient == "raise" and np.isfinite(v)) else v)
            if not np.isfinite(comp).any():
                continue
            T = float(np.nanmean(comp)); E = 50 + t.r_u * (T - 50)
            Es.append(E)
            fmax = max([Fr.get(u, pd.Series(dtype=float)).get(g, np.nan) for g in genes], default=np.nan)
            if best is None or E > best[0]:
                best = (E, t.genes, *comp, T, t.r_u, fmax, False)
        if not Es:
            out.append(None); continue
        E, genes, cD, cA, cO, T, rr, fmax, eu = best
        out.append(dict(layered_S=float(np.mean(Es)), layered_n_units=len(Es), layered_genes=genes, layered_D=cD, layered_A=cA,
                        layered_O=cO, layered_T=T, layered_r=rr, essential_unit=eu, layered_P=P.get(u, np.nan) if P else np.nan,
                        minority_alteration=bool(np.isfinite(cA) and cA >= 90 and np.isfinite(fmax) and fmax < 0.5),
                        essential_target=bool(set(UT[UT.cmap_name == d].genes.str.split("|").explode()) & ess)))
    return out


def main():
    g = pd.read_csv(SPLITS["group_id"], sep="\t", low_memory=False)
    s = pd.read_csv(SPLITS["stratum"], sep="\t", low_memory=False)
    groups, strata = sorted(g.group_id.unique()), sorted(s.stratum.unique())
    _, PdG = T79.comp_dependency(groups)
    FrG, PaG, _ = T79.comp_alteration(groups)
    _, PoG, _ = T79.comp_overactivity(groups)
    PdS, PaS, PoS, FrS = layer1_strata(strata, FrG)
    print(f"layer 1: groups D {len(PdG)} A {len(PaG)} O {len(PoG)}; subgroups D {len(PdS)} A {len(PaS)} O {len(PoS)} of {len(strata)}", flush=True)
    fm = T79.KU.fda_map()
    PROT = T79.PROT_SRC
    names = set(pd.concat([g, s]).query("pert_id in @fm.index").cmap_name)
    main_t = {}
    for df in (g, s):
        for r in df[df.pert_id.isin(fm.index) & df.target_source.isin(PROT)].drop_duplicates("cmap_name").itertuples():
            main_t.setdefault(r.cmap_name, [t.split(":")[0] for t in str(r.evidence_targets).split(";") if t and t != "nan"])
    UNITS, POT = T2.chembl_units(set(main_t), main_t)
    ESS = set(pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv").iloc[:, 0].str.split(" ").str[0])
    ut = []
    for d, us in UNITS.items():
        p = POT[d]
        pu = [max([p[x] for x in u["genes"] if x in p], default=np.nan) for u in us]
        pbest = np.nanmax(pu) if np.isfinite(pu).any() else np.nan
        for u, pv in zip(us, pu):
            o, _ = T79.orient({a for a in u["actions"] if a})
            ut.append(dict(cmap_name=d, genes="|".join(sorted(u["genes"])), orient=o, r_u=10 ** (pv - pbest) if np.isfinite(pv) else 1.0,
                           ess_frac=len(u["genes"] & ESS) / len(u["genes"])))
    UT = pd.DataFrame(ut)
    ess = set(pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv").iloc[:, 0].str.split(" ").str[0])
    run = {}
    allrows = []
    Pg = Ps = None
    DNA = set()
    if ESSENTIAL:
        Pg, Ps, Ig, Is = proliferation(groups, strata)
        run["proliferation_index"] = {**{k: round(v, 3) for k, v in Ig.items()}, **{k: round(v, 3) for k, v in Is.items()}}
        cand = set(pd.concat([g, s]).query("pert_id in @fm.index and target_source == 'non_protein_mechanism'").cmap_name)
        DNA = dna_acting(cand)
        run["dna_acting"] = sorted(DNA)
    for unit, df, (Pd, Pa, Po, Fr), P in (("group_id", g, (PdG, PaG, PoG, FrG), Pg), ("stratum", s, (PdS, PaS, PoS, FrS), Ps)):
        df = df.drop(columns=[c for c in COLS if c in df])
        fda = df.pert_id.isin(fm.index) & df.passes_gate.map(truthy)
        prot = fda & df.target_source.isin(PROT)
        res = score(df[prot], unit, Pd, Pa, Po, Fr, UT, ess, P)
        R = pd.DataFrame([r if r else {} for r in res], index=df.index[prot])
        for c in COLS:
            if c in R:
                df[c] = R[c]
        if ESSENTIAL:   # D124: DNA-acting drugs without a protein target: S = P, flag dna_acting
            dn = fda & df.target_source.eq("non_protein_mechanism") & df.cmap_name.isin(DNA)
            pv = df.loc[dn, unit].map(P)
            df.loc[dn, "layered_S"] = pv; df.loc[dn, "layered_T"] = pv; df.loc[dn, "layered_P"] = pv
            df.loc[dn, "layered_r"] = 1.0; df.loc[dn, "layered_n_units"] = 1; df.loc[dn, "layered_genes"] = "DNA"
            df["dna_acting"] = dn
        df["no_protein_target"] = fda & df["layered_S"].isna()
        for c in COLS:
            if c not in df:
                df[c] = np.nan
        for c in ("minority_alteration", "essential_target", "essential_unit", "dna_acting"):
            if c not in df:
                df[c] = False
            df[c] = df[c].fillna(False).astype(bool)
        has = fda & df.layered_S.notna()
        df["layered_pct"] = np.nan
        df.loc[has, "layered_pct"] = df[has].groupby(unit).layered_S.transform(lambda v: v.rank(method="max") / len(v) * 100)
        if not SANDBOX:
            df.to_csv(SPLITS[unit], sep="\t", index=False)
        run[unit] = dict(rows=len(df), fda_gate_passing=int(fda.sum()), scored=int(has.sum()), no_protein_target=int(df.no_protein_target.sum()),
                         minority_alteration=int(df.minority_alteration.sum()), essential_target=int((df.essential_target & has).sum()))
        allrows.append(df[fda][[unit, "pert_id", "cmap_name"] + COLS].rename(columns={unit: "unit"}).assign(level=unit))
        print(unit, run[unit], flush=True)
    sfx = "_d124" if SANDBOX else ""
    pd.concat(allrows).to_csv(RES / f"layered_score_v2{sfx}.tsv", sep="\t", index=False, float_format="%.3f")
    json.dump(run, open(RES / f"lead_score_run{sfx}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
