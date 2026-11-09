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
        "layered_r", "no_protein_target", "minority_alteration", "essential_target"]
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


def score(x, unit, Pd, Pa, Po, Fr, UT, ess):
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
            comp = []
            for P, k in ((Pd, "D"), (Pa, "A"), (Po, "O")):
                if k == "A" and t.orient == "raise":
                    comp.append(np.nan); continue
                s = P.get(u)
                vals = [s[g] for g in genes if s is not None and g in s.index and np.isfinite(s[g])]
                v = float(np.mean(vals)) if vals else np.nan
                comp.append(100 - v if (t.orient == "raise" and np.isfinite(v)) else v)
            if not np.isfinite(comp).any():
                continue
            T = float(np.nanmean(comp)); E = 50 + t.r_u * (T - 50)
            Es.append(E)
            fmax = max([Fr.get(u, pd.Series(dtype=float)).get(g, np.nan) for g in genes], default=np.nan)
            if best is None or E > best[0]:
                best = (E, t.genes, *comp, T, t.r_u, fmax)
        if not Es:
            out.append(None); continue
        E, genes, cD, cA, cO, T, rr, fmax = best
        out.append(dict(layered_S=float(np.mean(Es)), layered_n_units=len(Es), layered_genes=genes, layered_D=cD, layered_A=cA,
                        layered_O=cO, layered_T=T, layered_r=rr,
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
    ut = []
    for d, us in UNITS.items():
        p = POT[d]
        pu = [max([p[x] for x in u["genes"] if x in p], default=np.nan) for u in us]
        pbest = np.nanmax(pu) if np.isfinite(pu).any() else np.nan
        for u, pv in zip(us, pu):
            o, _ = T79.orient({a for a in u["actions"] if a})
            ut.append(dict(cmap_name=d, genes="|".join(sorted(u["genes"])), orient=o, r_u=10 ** (pv - pbest) if np.isfinite(pv) else 1.0))
    UT = pd.DataFrame(ut)
    ess = set(pd.read_csv(ROOT / "data/raw/depmap_24q4/CRISPRInferredCommonEssentials.csv").iloc[:, 0].str.split(" ").str[0])
    run = {}
    allrows = []
    for unit, df, (Pd, Pa, Po, Fr) in (("group_id", g, (PdG, PaG, PoG, FrG)), ("stratum", s, (PdS, PaS, PoS, FrS))):
        df = df.drop(columns=[c for c in COLS if c in df])
        fda = df.pert_id.isin(fm.index) & df.passes_gate.map(truthy)
        prot = fda & df.target_source.isin(PROT)
        res = score(df[prot], unit, Pd, Pa, Po, Fr, UT, ess)
        R = pd.DataFrame([r if r else {} for r in res], index=df.index[prot])
        for c in COLS:
            if c in R:
                df[c] = R[c]
        df["no_protein_target"] = fda & df["layered_S"].isna()
        for c in ("minority_alteration", "essential_target"):
            df[c] = df[c].fillna(False).astype(bool)
        has = fda & df.layered_S.notna()
        df["layered_pct"] = np.nan
        df.loc[has, "layered_pct"] = df[has].groupby(unit).layered_S.transform(lambda v: v.rank(method="max") / len(v) * 100)
        df.to_csv(SPLITS[unit], sep="\t", index=False)
        run[unit] = dict(rows=len(df), fda_gate_passing=int(fda.sum()), scored=int(has.sum()), no_protein_target=int(df.no_protein_target.sum()),
                         minority_alteration=int(df.minority_alteration.sum()), essential_target=int((df.essential_target & has).sum()))
        allrows.append(df[fda][[unit, "pert_id", "cmap_name"] + COLS].rename(columns={unit: "unit"}).assign(level=unit))
        print(unit, run[unit], flush=True)
    pd.concat(allrows).to_csv(RES / "layered_score_v2.tsv", sep="\t", index=False, float_format="%.3f")
    json.dump(run, open(RES / "lead_score_run.json", "w"), indent=1)


if __name__ == "__main__":
    main()
