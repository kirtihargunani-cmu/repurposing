#!/usr/bin/env python
"""D124 test (rule: docs/DECISIONS.md D124). Compares the D122 lead (results/lincs/layered_score_v2.tsv) with the D124 sandbox
(layered_score_v2_d124.tsv, from 79t3 ESSENTIAL=1 SANDBOX=1) on the 41 groups. Runs in `celligner`.
Drug level = mean layered_pct over the ingredient's gate-passing compounds; labels / anticancer / eligible (D122 gates) /
exclusions from known_use/labels.tsv. Essential-unit drugs = ingredients whose best unit is essential in the sandbox in
any group (implementation of 'drugs with an essential unit'). Outputs results/lincs/d124/{criteria.tsv, between.tsv, top10.tsv}.
"""
import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/lincs"
NEW = os.environ.get("NEW", "_d124")
OUT = RES / NEW.lstrip("_")
TIE_MAX = 3   # D124a criterion (4)
sys.path.insert(0, str(ROOT / "scripts"))
_s = importlib.util.spec_from_file_location("ku", ROOT / "scripts/79l_known_use.py"); KU = importlib.util.module_from_spec(_s); _s.loader.exec_module(KU)
_c = importlib.util.spec_from_file_location("ec", ROOT / "scripts/checks/essential_check.py"); EC = importlib.util.module_from_spec(_c); _c.loader.exec_module(EC)
ANTI_ER = ["FULVESTRANT", "TAMOXIFEN CITRATE", "RALOXIFENE HYDROCHLORIDE"]
rng = np.random.default_rng(0)


def drug_level(f):
    x = pd.read_csv(f, sep="\t")
    x = x[x.level == "group_id"].rename(columns={"unit": "cohort"})
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    x["fda_ingredient"] = x.pert_id.map(fm)
    agg = dict(pct=("layered_pct", "mean"))
    for c in ("essential_unit", "dna_acting", "essential_target"):
        if c in x:
            agg[c] = (c, lambda v: bool(pd.Series(v).astype(str).eq("True").any()))
    agg["statin_cyclo"] = ("layered_genes", lambda v: bool(pd.Series(v).astype(str).str.contains(r"\bHMGCR\b|\bPPIA\b").any()))
    return x.groupby(["cohort", "fda_ingredient"]).agg(**agg).reset_index()


def between(D, score, drugs):
    B = D[D.fda_ingredient.isin(drugs)].dropna(subset=[score])
    return EC.pooled(B.fda_ingredient.values, B.y.values.astype(bool), B[score].astype(float).values)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    old, new = drug_level(RES / "layered_score_v2.tsv"), drug_level(RES / f"layered_score_v2{NEW}.tsv")
    old = old.rename(columns={"statin_cyclo": "statin_cyclo_old"})
    lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")[["cohort", "fda_ingredient", "drug", "anticancer", "eligible", "lab_hub",
                                                                  "lab_chembl", "excl_hub", "excl_chembl"]]
    D = lab.merge(old[["cohort", "fda_ingredient", "pct", "statin_cyclo_old"]].rename(columns={"pct": "old"}), how="left").merge(
        new.rename(columns={"pct": "new"}), how="left")
    ess_drugs = set(D[D.essential_unit.fillna(False).astype(bool)].fda_ingredient)
    dna_drugs = set(D[D.dna_acting.fillna(False).astype(bool)].fda_ingredient)
    crit, btw = [], []
    for src in ("hub", "chembl"):
        y = f"lab_{src}"
        A = D[D.anticancer.astype(bool) & ~D[f"excl_{src}"].astype(bool)].assign(y=lambda d: d[y].astype(bool))
        E = A[A.eligible.astype(bool)]
        ko = KU.pooled_auc(E.dropna(subset=["old"]), y, "old", rng, 1000, 1000)
        kn = KU.pooled_auc(E.dropna(subset=["new"]), y, "new", rng, 1000, 1000)
        crit.append(dict(source=src, test="(2) known use, primary rows", old=ko["auc"], new=kn["auc"], n_pos_old=ko["n_pos"], n_pos_new=kn["n_pos"],
                         ci_new=f"{kn.get('ci_lo')}-{kn.get('ci_hi')}"))
        both = A.dropna(subset=["old", "new"])
        bo, pb = between(both, "old", ess_drugs); bn, _ = between(both, "new", ess_drugs)
        crit.append(dict(source=src, test="(1) between cancers, essential-unit drugs", old=round(bo, 4), new=round(bn, 4), pairs=int(pb)))
        dn, pdn = between(A.dropna(subset=["new"]), "new", dna_drugs)
        btw.append(dict(source=src, drugs="DNA-acting (new only)", auc=round(dn, 4), pairs=int(pdn), n_drugs=len(dna_drugs)))
        mk = pd.read_csv(RES / "essential_check/markers.tsv", sep="\t").set_index("group_id")
        A2 = A.join(mk[["SLFN11", "proliferation"]], on="cohort")
        for m in ("SLFN11", "proliferation"):
            v, pp = between(A2, m, dna_drugs)
            btw.append(dict(source=src, drugs="DNA-acting", score=m, auc=round(v, 4), pairs=int(pp)))
    # (3) breast controls
    def ctl(col):
        X = D.dropna(subset=[col]).copy()   # all scored rows, as the D122 breast check (anti-estrogens are excluded in Basal by the target-absent gate)
        X["p"] = X.groupby("cohort")[col].rank(pct=True) * 100
        lum = X[X.cohort.isin(["BRCA-LumA", "BRCA-LumB"]) & X.fda_ingredient.isin(ANTI_ER)].p.median()
        bas = X[(X.cohort == "BRCA-Basal") & X.fda_ingredient.isin(ANTI_ER)].p.median()
        return round(lum, 1), round(bas, 1)
    co, cn = ctl("old"), ctl("new")
    crit.append(dict(source="both", test="(3) anti-estrogens (luminal median, basal median)", old=str(co), new=str(cn)))
    C = pd.DataFrame(crit)
    k = {r.source: r for r in C.itertuples() if r.test.startswith("(1)")}
    q = {r.source: r for r in C.itertuples() if r.test.startswith("(2)")}
    ok1 = k["chembl"].new >= k["chembl"].old + 0.03 and k["hub"].new >= k["hub"].old
    ok2 = all(q[s].new >= q[s].old - 0.01 for s in ("hub", "chembl"))
    ok3 = cn[0] >= co[0] and cn[1] <= co[1]
    # (4) ties among non-DNA-acting drugs in each group's top 10 (eligible, drug level)
    ties = []
    for coh, d in D[D.eligible.astype(bool) & D.new.notna()].groupby("cohort"):
        t10 = d.nlargest(10, "new", keep="all")
        nd = t10[~t10.dna_acting.fillna(False).astype(bool)]
        ties.append(dict(cohort=coh, largest_tie=int(nd.new.round(6).value_counts().max()) if len(nd) else 0,
                         statin_cyclo_top10_new=int(t10.statin_cyclo.fillna(False).astype(bool).sum())))
    for coh, d in D[D.eligible.astype(bool) & D.old.notna()].groupby("cohort"):
        t10o = d.nlargest(10, "old")
        for r in ties:
            if r["cohort"] == coh:
                r["statin_cyclo_top10_old"] = int(t10o.statin_cyclo_old.fillna(False).astype(bool).sum())
    TI = pd.DataFrame(ties); TI.to_csv(OUT / "ties.tsv", sep="\t", index=False)
    ok4 = bool((TI.largest_tie <= TIE_MAX).all())
    crit.append(dict(source="both", test=f"(4) largest non-DNA tie in a top 10 (max over groups; <= {TIE_MAX})", old="", new=str(int(TI.largest_tie.max()))))
    crit.append(dict(source="both", test="statin / cyclosporine rows in top 10s (reported)", old=str(int(TI.get("statin_cyclo_top10_old", pd.Series(dtype=int)).sum())),
                     new=str(int(TI.statin_cyclo_top10_new.sum()))))
    C = pd.DataFrame(crit)
    use4 = NEW != "_d124"
    allok = ok1 and ok2 and ok3 and (ok4 or not use4)
    C.loc[len(C)] = dict(source="DECISION", test=("adopt " + NEW.lstrip("_").upper()) if allok else "current objective stays",
                         old=f"(1) {'met' if ok1 else 'not met'}; (2) {'met' if ok2 else 'not met'}; (3) {'met' if ok3 else 'not met'}"
                             + (f"; (4) {'met' if ok4 else 'not met'}" if use4 else ""))
    C.to_csv(OUT / "criteria.tsv", sep="\t", index=False)
    pd.DataFrame(btw).to_csv(OUT / "between.tsv", sep="\t", index=False)
    # top 10 (eligible drugs)
    T = []
    for coh, d in D[D.eligible.astype(bool)].groupby("cohort"):
        a = set(d.nlargest(10, "old").fda_ingredient); b = d.nlargest(10, "new")
        T.append(dict(cohort=coh, kept=len(a & set(b.fda_ingredient)), new_ess=int(b.essential_unit.fillna(False).astype(bool).sum()),
                      new_dna=int(b.dna_acting.fillna(False).astype(bool).sum()), top3_new=", ".join(b.drug.head(3))))
    T = pd.DataFrame(T); T.to_csv(OUT / "top10.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 90)
    print(C.to_string(index=False)); print(pd.DataFrame(btw).to_string(index=False))
    print("top-10 kept median", T.kept.median(), "; essential-unit rows in new top 10s", T.new_ess.sum(), "; DNA-acting", T.new_dna.sum())
    print(T.to_string(index=False))


if __name__ == "__main__":
    main()
