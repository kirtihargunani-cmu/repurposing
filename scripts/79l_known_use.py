#!/usr/bin/env python
"""Known-use recovery: are drugs already used for a cancer ranked high for it? (C28; redesigned D118, 2026-11-08). Reported only.

Unit (D118): cohort (disease) x FDA ingredient. Score = mean post-gate spec_mech_pct over the ingredient's LINCS compounds that pass
the target gate in the cohort (fda_matches; was: the best compound). Eligible = any compound passes the D115 test (fda_rows.eligible).
Labels (config/known_use.yaml), per drug: a drug is a known use of a cohort when one of its indication terms (Repurposing Hub
`indication` pipe entries, Launched; ChEMBL 37 drug_indication mesh_heading / efo_term separately, max_phase_for_ind >= 3) matches the
cohort's include and not its exclude regex AND matches cancer_term. ChEMBL ids of a drug = every id of every compound (annotation, 79w
name lookup, D113 resolved ids, ChEMBL parents); Hub = any compound's normalised name. Subtype table: a labelled drug of a listed class
keeps the label only in its allowed cohorts; its rows in the other siblings are left out of the tests.
Anticancer = any compound with ATC L01 / L02, or Hub disease area oncology.
Tests:
  PRIMARY  eligible rows of anticancer drugs; pooled within-cohort AUC (every labelled x unlabelled pair in a cohort once; cohorts with
           >= 1 of each); p from N_PERM within-cohort permutations; 95% CI from N_BOOT cohort bootstraps; all / matched / borrowed.
  Shown with it: class baseline (anticancer 0/1 as the score, all eligible rows), lead score on all eligible rows, the primary
           metric on pre-gate percentiles (all gate-passing rows).
  A1 / A2  (C28, as before, on the new labels and drug-level scores; A1 = median of per-cohort AUCs, >= MIN_LAB / MIN_UNLAB).
Outputs (results/lincs/known_use/): labels.tsv (drug level), label_changes.tsv (vs labels_pre_d118.tsv), primary.tsv,
         a1_per_cohort.tsv, a1_summary.tsv, a2_per_drug.tsv, a2_summary.tsv, coverage.tsv, run.json
`drug_table(split_df)` and `pooled_auc(...)` are imported by d100_decide.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chembl_lib as C  # noqa: E402
from fda_rows import eligible, truthy, drug_name  # noqa: E402
from route_flag import apply_route  # noqa: E402
from evidence_flags import apply_evidence  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "known_use"
CH = ROOT / "data/raw/chembl/ChEMBL_37"
CFG = yaml.safe_load(open(ROOT / "config/known_use.yaml"))
CANCER = re.compile(CFG["cancer_term"], re.I)
SCORES = ["spec_mech_pct", "spec_pct", "rev_pct", "target_frac", "mean_pct"]
LEAD = "spec_mech_pct"
MIN_LAB, MIN_UNLAB, MIN_COH, N_PERM, N_BOOT, SEED = 3, 20, 10, 1000, 1000, 0


def matches(term, pat):
    if not CANCER.search(term) or not re.search(pat["include"], term, re.I):
        return False
    return not (pat.get("exclude") and re.search(pat["exclude"], term, re.I))


def auc(pos, neg):
    """P(a labelled drug scores above an unlabelled one), ties at 0.5."""
    v = np.concatenate([pos, neg])
    r = pd.Series(v).rank().values
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


# ---------------------------------------------------------------- drug-level scores (shared with d100_decide)
def fda_map():
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id")
    return fm.set_index("pert_id").fda_ingredient


def ing_compounds():
    """FDA ingredient -> (all its LINCS pert_ids, all their cmap_names), from fda_matches (drug level, any cohort)."""
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id")
    return {g: (set(d.pert_id), set(d.cmap_name.astype(str))) for g, d in fm.groupby("fda_ingredient")}


def drug_table(sp, unit="group_id"):
    """split table -> one row per unit x FDA ingredient: mean scores over gate-passing compounds, eligible (any compound), evidence."""
    fm = fda_map()
    x = sp[sp.pert_id.isin(fm.index) & sp.passes_gate.map(truthy)].copy()
    x["fda_ingredient"] = x.pert_id.map(fm)
    x = apply_evidence(apply_route(x, unit), unit)
    x["elig"] = eligible(x, unit)
    x["direct"] = x.source.eq("direct")
    agg = dict(spec_mech_pct=("spec_mech_pct", "mean"), n_compounds=("pert_id", "nunique"), eligible=("elig", "any"),
               any_direct_elig=("direct", lambda s: bool((s & x.loc[s.index, "elig"]).any())), any_direct=("direct", "any"),
               cmap_names=("cmap_name", lambda s: "|".join(sorted(set(s)))), pert_ids=("pert_id", lambda s: "|".join(sorted(set(s)))))
    for c in ("spec_mech_pct_pregate", "spec_pct", "rev_pct", "target_frac"):
        if c in x:
            agg[c] = (c, "mean")
    D = x.groupby([unit, "fda_ingredient"]).agg(**agg).reset_index().rename(columns={unit: "cohort"})
    D["drug"] = D.fda_ingredient.map(drug_name)
    D["evidence"] = np.where(np.where(D.eligible, D.any_direct_elig, D.any_direct), "matched", "borrowed")
    return D


def pooled_auc(D, lab, score, rng=None, n_perm=0, n_boot=0):
    """pooled within-cohort AUC of `score` (labelled vs unlabelled pairs inside each cohort); returns dict."""
    per = []
    for _, d in D.dropna(subset=[score]).groupby("cohort"):
        y = d[lab].values.astype(bool); v = d[score].values.astype(float)
        if y.sum() and (~y).sum():
            per.append((v, y))
    if not per:
        return dict(auc=np.nan, n_cohorts=0, n_pos=0, n_neg=0, pairs=0)

    def u(v, y):
        r = pd.Series(v).rank().values
        return r[y].sum() - y.sum() * (y.sum() + 1) / 2, y.sum() * (~y).sum()
    UP = np.array([u(v, y) for v, y in per])
    obs = UP[:, 0].sum() / UP[:, 1].sum()
    out = dict(auc=round(obs, 4), n_cohorts=len(per), n_pos=int(sum(y.sum() for _, y in per)),
               n_neg=int(sum((~y).sum() for _, y in per)), pairs=int(UP[:, 1].sum()))
    if rng is not None and n_perm:
        null = []
        for _ in range(n_perm):
            P = np.array([u(v, rng.permutation(y)) for v, y in per])
            null.append(P[:, 0].sum() / P[:, 1].sum())
        out["p_perm"] = round((1 + np.sum(np.array(null) >= obs)) / (1 + n_perm), 4)
    if rng is not None and n_boot:
        bs = []
        for _ in range(n_boot):
            i = rng.integers(0, len(UP), len(UP))
            bs.append(UP[i, 0].sum() / UP[i, 1].sum())
        out["ci_lo"], out["ci_hi"] = round(np.percentile(bs, 2.5), 4), round(np.percentile(bs, 97.5), 4)
    return out


# ---------------------------------------------------------------- labels
def drug_ids(D, ann):
    """FDA ingredient -> every ChEMBL id of its compounds (annotation, name lookup, D113, parents)."""
    nl = json.load(open(CH / "d105_all/name_lookup.json"))
    d113 = json.load(open(CH / "d105_all/d113_ids.json")) if (CH / "d105_all/d113_ids.json").exists() else {}
    par = {f["molecule_chembl_id"]: f["parent_chembl_id"] for f in json.load(open(CH / "d105_all/forms.json"))}
    a = ann.set_index("pert_id").chembl_id
    out = {}
    IC = ing_compounds()
    for ing in D.fda_ingredient.unique():
        perts, names = IC[ing]
        ids = {a.get(p) for p in perts} | {nl.get(n) for n in names}
        for n in names:
            ids |= set(d113.get(n, {}).get("ids", []))
        ids = {i for i in ids if isinstance(i, str)}
        out[ing] = ids | {par.get(i, i) for i in ids}
    return out


def indications(ids):
    """ChEMBL drug_indication for these molecules (cache extended with ids not fetched before; ChEMBL 37 checked live)."""
    cache = CH / "indications"; f = cache / "drug_indication.json"; fi = cache / "fetched_ids.json"
    ind = json.load(open(f)) if f.exists() else []
    done = set(json.load(open(fi))) if fi.exists() else {i["molecule_chembl_id"] for i in ind}
    need = sorted(set(ids) - done)
    if need:
        live = C.get(f"{C.API}/status.json")["chembl_db_version"]
        assert live == C.RELEASE, f"ChEMBL API is {live}, config pins {C.RELEASE}"
        cache.mkdir(parents=True, exist_ok=True)
        ind += sum([C.paged("drug_indication", "drug_indications", {"molecule_chembl_id__in": ",".join(c),
                    "only": "molecule_chembl_id,mesh_heading,efo_term,max_phase_for_ind"}) for c in C.chunks(need, 30)], [])
        json.dump(ind, open(f, "w")); json.dump(sorted(done | set(need)), open(fi, "w"))
    return ind, len(need)


def labels(D, ann):
    cohorts = sorted(D.cohort.unique())
    # Hub: name -> cohorts, per indication term
    h = pd.read_csv(ROOT / "data/raw/repurposing_hub/repurposing_drugs_20200324.txt", sep="\t", comment="!")
    hda = h.set_index(h.pert_iname.map(lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower()))).disease_area
    h = h[h.clinical_phase.isin(CFG["sources"]["hub"]["phases"]) & h.indication.notna()]
    hub = {}
    for r in h.itertuples():
        k = re.sub(r"[^a-z0-9]", "", r.pert_iname.lower())
        for term in r.indication.split("|"):
            for coh in cohorts:
                if coh in CFG["cohorts"] and matches(term, CFG["cohorts"][coh]):
                    hub.setdefault(k, set()).add(coh)
    ids = drug_ids(D, ann)
    ind, n_new = indications({i for s in ids.values() for i in s})
    lo = CFG["sources"]["chembl"]["min_phase"]
    chem = {}
    for i in ind:
        try:
            if float(i.get("max_phase_for_ind") or 0) < lo:
                continue
        except (TypeError, ValueError):
            continue
        for term in (i.get("mesh_heading"), i.get("efo_term")):
            if not term:
                continue
            for coh in cohorts:
                if coh in CFG["cohorts"] and matches(str(term), CFG["cohorts"][coh]):
                    chem.setdefault(i["molecule_chembl_id"], set()).add(coh)
    IC = ing_compounds()
    keys = lambda names: {re.sub(r"[^a-z0-9]", "", n.lower()) for n in names}
    hub_d = {g: set().union(*[hub.get(k, set()) for k in keys(IC[g][1])]) for g in D.fda_ingredient.unique()}
    chem_d = {g: set().union(*[chem.get(i, set()) for i in s]) if s else set() for g, s in ids.items()}
    D["lab_hub"] = [c in hub_d.get(g, ()) for g, c in zip(D.fda_ingredient, D.cohort)]
    D["lab_chembl"] = [c in chem_d.get(g, ()) for g, c in zip(D.fda_ingredient, D.cohort)]
    # anticancer: any compound ATC L01 / L02, or Hub disease area oncology
    a = ann.drop_duplicates("pert_id").set_index("pert_id")
    def onc(g):
        ps, ns = IC[g]
        atc = any(re.search(r"\bL0[12]", str(a.atc.get(p, ""))) for p in ps)
        hd = any("oncology" in str(hda.get(k, "")) for k in keys(ns)) or any("oncology" in str(a.hub_disease_area.get(p, "")) for p in ps)
        return atc or hd
    on = {g: onc(g) for g in D.fda_ingredient.unique()}
    D["anticancer"] = D.fda_ingredient.map(on)
    # subtype table (main targets = D105b evidence targets of any compound)
    sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False, usecols=["pert_id", "evidence_targets"])
    mt = sp.dropna().drop_duplicates("pert_id").set_index("pert_id").evidence_targets.map(lambda s: {t.split(":")[0] for t in s.split(";") if t})
    tg = {g: set().union(*[mt.get(p, set()) for p in IC[g][0]]) for g in D.fda_ingredient.unique()}
    for src in ("hub", "chembl"):
        D[f"excl_{src}"] = False
        D[f"subtype_{src}"] = ""
        for rule in CFG.get("subtype", []):
            cls = {g for g, t in tg.items() if t & set(rule.get("main_targets", []))} | \
                  {g for g in tg if drug_name(g) in set(rule.get("drugs", []))}
            sib = D.cohort.isin(rule["siblings"]) & D.fda_ingredient.isin(cls)
            labd = set(D.fda_ingredient[sib & D[f"lab_{src}"]])
            hit = sib & D.fda_ingredient.isin(labd)
            off = hit & ~D.cohort.isin(rule["allowed"])
            D.loc[off, f"excl_{src}"] = True
            D.loc[off, f"lab_{src}"] = False
            D.loc[hit, f"subtype_{src}"] = rule["name"]
    return D, dict(chembl_indication_records=len(ind), chembl_ids_newly_fetched=n_new)


def label_changes(D):
    f = OUT / "labels_pre_d118.tsv"
    if not f.exists():
        return None
    o = pd.read_csv(f, sep="\t")
    fm = fda_map()
    o["fda_ingredient"] = o.pert_id.map(fm)
    o = o.dropna(subset=["fda_ingredient"]).groupby(["cohort", "fda_ingredient"])[["lab_hub", "lab_chembl"]].any().reset_index()
    m = D[["cohort", "fda_ingredient", "drug", "lab_hub", "lab_chembl", "excl_hub", "excl_chembl"]].merge(
        o, on=["cohort", "fda_ingredient"], how="outer", suffixes=("", "_old"))
    rows = []
    for src in ("hub", "chembl"):
        new, old = m[f"lab_{src}"].fillna(False).astype(bool), m[f"lab_{src}_old"].fillna(False).astype(bool)
        exc = m[f"excl_{src}"].fillna(False).astype(bool)
        for kind, sel in (("gained", new & ~old), ("lost", ~new & old & ~exc), ("subtype_removed", old & exc)):
            for r in m[sel].itertuples():
                rows.append(dict(source=src, change=kind, cohort=r.cohort, fda_ingredient=r.fda_ingredient,
                                 drug=r.drug if isinstance(r.drug, str) else drug_name(r.fda_ingredient)))
    ch = pd.DataFrame(rows)
    ch.to_csv(OUT / "label_changes.tsv", sep="\t", index=False)
    return ch


# ---------------------------------------------------------------- C28 A1 / A2 (as before)
def tests(D, src, rng):
    a1, a1s, a2, a2s = [], [], [], []
    for sub, name in ((D, "all"), (D[D.evidence == "matched"], "matched"), (D[D.evidence == "borrowed"], "borrowed")):
        per = []
        for coh, d in sub.groupby("cohort"):
            if d.label.sum() < MIN_LAB or (~d.label).sum() < MIN_UNLAB:
                continue
            row = dict(source=src, subset=name, cohort=coh, n_labelled=int(d.label.sum()), n_unlabelled=int((~d.label).sum()))
            for c in SCORES:
                v = d[c].values.astype(float); ok = ~np.isnan(v)
                if ok.sum() < MIN_LAB + MIN_UNLAB or d.label[ok].sum() < MIN_LAB or (~d.label[ok]).sum() < MIN_UNLAB:
                    row[f"auc_{c}"] = np.nan; continue
                row[f"auc_{c}"] = auc(v[ok & d.label.values], v[ok & ~d.label.values])
            per.append((row, d))
            if name == "all":
                a1.append(row)
        if not per:
            continue
        R = pd.DataFrame([r for r, _ in per])
        for c in SCORES:
            a = R[f"auc_{c}"].dropna()
            if len(a) < 3:
                continue
            obs = a.median(); null = []
            for _ in range(N_PERM):
                vals = []
                for row, d in per:
                    v = d[c].values.astype(float); ok = ~np.isnan(v)
                    lab = rng.permutation(d.label.values[ok])
                    if lab.sum() < MIN_LAB or (~lab).sum() < MIN_UNLAB:
                        continue
                    vals.append(auc(v[ok][lab], v[ok][~lab]))
                null.append(np.median(vals) if vals else np.nan)
            null = np.array(null, float)
            bs = [np.median(rng.choice(a.values, len(a))) for _ in range(N_BOOT)]
            a1s.append(dict(source=src, subset=name, score=c, n_cohorts=len(a), median_auc=round(obs, 4),
                            ci_lo=round(np.percentile(bs, 2.5), 4), ci_hi=round(np.percentile(bs, 97.5), 4),
                            frac_cohorts_above_half=round((a > 0.5).mean(), 3),
                            p_perm=(1 + np.nansum(null >= obs)) / (1 + np.isfinite(null).sum())))
    for drug, d in D.groupby("drug"):
        if len(d) < MIN_COH or d.label.sum() == 0 or (~d.label).sum() == 0:
            continue
        here, there = d[LEAD][d.label].median(), d[LEAD][~d.label].median()
        a2.append(dict(source=src, drug=drug, n_cohorts=len(d), n_labelled=int(d.label.sum()),
                       pct_labelled=round(here, 2), pct_other=round(there, 2), diff=round(here - there, 2)))
    if a2:
        A = pd.DataFrame(a2)
        k, n = int((A["diff"] > 0).sum()), int((A["diff"] != 0).sum())
        from scipy.stats import binomtest
        a2s.append(dict(source=src, n_drugs=len(A), median_diff=round(A["diff"].median(), 3),
                        frac_higher_in_own_cohort=round(k / max(n, 1), 3),
                        p_sign=binomtest(k, n, 0.5, alternative="greater").pvalue if n else np.nan))
    return a1, a1s, a2, a2s


def main():
    rng = np.random.default_rng(SEED)
    OUT.mkdir(parents=True, exist_ok=True)
    sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_targets", "targets_expressed"]).drop_duplicates(["group_id", "pert_id"])
    nexp = dd.targets_expressed.fillna("").map(lambda x: len([t for t in str(x).split("|") if t]))
    dd["target_frac"] = np.where(dd.n_targets > 0, nexp / dd.n_targets.where(dd.n_targets > 0), np.nan)
    sp = sp.merge(dd[["group_id", "pert_id", "target_frac"]], on=["group_id", "pert_id"], how="left")
    D = drug_table(sp)
    D["mean_pct"] = D.groupby("drug")[LEAD].transform("mean")
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False)
    D, info = labels(D, ann)
    ch = label_changes(D)
    if ch is not None and len(ch):
        print("label changes vs the pre-D118 labels (cohort x drug):")
        print(ch.groupby(["source", "change"]).size().to_string())

    # PRIMARY and companions
    prim = []
    E = D[D.eligible]
    for src in ("hub", "chembl"):
        lab = f"lab_{src}"
        keep = lambda X: X[~X[f"excl_{src}"]]
        for subset in ("all", "matched", "borrowed"):
            sel = lambda X: X if subset == "all" else X[X.evidence == subset]
            runs = [("primary: anticancer, eligible, lead score", sel(keep(E[E.anticancer])), LEAD, True),
                    ("baseline: anticancer 0/1, all eligible", sel(keep(E.assign(onc=E.anticancer.astype(float)))), "onc", subset == "all"),
                    ("lead score, all eligible", sel(keep(E)), LEAD, subset == "all"),
                    ("primary on pre-gate percentiles (anticancer, gate-passing)", sel(keep(D[D.anticancer])), "spec_mech_pct_pregate", subset == "all")]
            for name, X, score, full in runs:
                if score not in X:
                    continue
                r = pooled_auc(X, lab, score, rng, N_PERM if full else 0, N_BOOT if full else 0)
                prim.append(dict(source=src, subset=subset, test=name, score=score, **r))
    P = pd.DataFrame(prim)
    P.to_csv(OUT / "primary.tsv", sep="\t", index=False)

    runs, a1, a1s, a2, a2s = {}, [], [], [], []
    for src in ("hub", "chembl"):
        d = D[~D[f"excl_{src}"]].assign(label=lambda x: x[f"lab_{src}"])
        runs[src] = dict(labelled_rows=int(d.label.sum()), labelled_drugs=int(d.drug[d.label].nunique()),
                         cohorts_with_labels=int(d.groupby("cohort").label.sum().gt(0).sum()),
                         subtype_removed_rows=int(D[f"excl_{src}"].sum()))
        r = tests(d, src, rng)
        a1 += r[0]; a1s += r[1]; a2 += r[2]; a2s += r[3]
    D.drop(columns=["any_direct_elig", "any_direct"]).to_csv(OUT / "labels.tsv", sep="\t", index=False)
    pd.DataFrame(a1).to_csv(OUT / "a1_per_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(a1s).to_csv(OUT / "a1_summary.tsv", sep="\t", index=False)
    pd.DataFrame(a2).to_csv(OUT / "a2_per_drug.tsv", sep="\t", index=False)
    pd.DataFrame(a2s).to_csv(OUT / "a2_summary.tsv", sep="\t", index=False)
    cov = D.groupby("cohort").agg(n_drugs=("drug", "nunique"), n_eligible=("eligible", "sum"), n_anticancer_eligible=(
        "anticancer", lambda s: int((s & D.loc[s.index, "eligible"]).sum())), hub=("lab_hub", "sum"), chembl=("lab_chembl", "sum"))
    cov.to_csv(OUT / "coverage.tsv", sep="\t")
    json.dump(dict(rule="D118", rows=len(D), drugs=int(D.drug.nunique()), cohorts=int(D.cohort.nunique()),
                   eligible_rows=int(D.eligible.sum()), anticancer_drugs=int(D.drop_duplicates("drug").anticancer.sum()),
                   sources=runs, min_lab=MIN_LAB, min_unlab=MIN_UNLAB, n_perm=N_PERM, n_boot=N_BOOT, seed=SEED, **info),
              open(OUT / "run.json", "w"), indent=1)
    pd.set_option("display.width", 250)
    print(P.to_string(index=False))
    print(pd.DataFrame(a1s)[lambda x: x.subset == "all"].to_string(index=False))
    print(pd.DataFrame(a2s).to_string(index=False))


if __name__ == "__main__":
    main()
