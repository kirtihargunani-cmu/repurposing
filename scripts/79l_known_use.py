#!/usr/bin/env python
"""Known-use recovery: are drugs already used for a cancer ranked high for it? (protocol C28, runs in `celligner`).

Labels (config/known_use.yaml, written before any result): a drug is a known use of a cohort when its recorded
indication text matches the cohort's include regex and not its exclude regex. Two independent label sets:
  hub     Repurposing Hub `indication`, Launched drugs (data/raw/repurposing_hub, 2020-03-24 snapshot)
  chembl  ChEMBL 37 drug_indication, max_phase_for_ind >= 3 (cached in data/raw/chembl/ChEMBL_37/indications/)
Rows: the published FDA rows (fda_rows.py), diseases only (subgroups inherit the parent's labels, so they cannot
separate siblings). Matching is on the LINCS cmap_name / FDA ingredient for the Hub and on ChEMBL ids for ChEMBL.
Tests (C28):
  A1  per cohort with >= MIN_LAB labelled and >= MIN_UNLAB unlabelled drugs: AUC of the score, labelled vs unlabelled;
      median over cohorts, permutation of labels within cohort, bootstrap CI over cohorts
  A2  per labelled drug in >= MIN_COH cohorts: percentile in its labelled cohorts vs its median elsewhere; sign test
  Baselines: target_frac, rev, and the drug's mean percentile (high everywhere); split matched / borrowed.
Outputs (results/lincs/known_use/): labels.tsv, a1_per_cohort.tsv, a1_summary.tsv, a2_per_drug.tsv, a2_summary.tsv,
         coverage.tsv, run.json
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import chembl_lib as C
from fda_rows import fda_rows

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
OUT = RES / "known_use"
CFG = yaml.safe_load(open(ROOT / "config/known_use.yaml"))
SCORES = ["spec_mech_pct", "spec_pct", "rev_pct", "target_frac", "mean_pct"]
LEAD = "spec_mech_pct"
MIN_LAB, MIN_UNLAB, MIN_COH, N_PERM, N_BOOT, SEED = 3, 20, 10, 1000, 1000, 0


def matches(text, pat):
    inc, exc = pat["include"], pat.get("exclude")
    if not re.search(inc, text, re.I):
        return False
    return not (exc and re.search(exc, text, re.I))


def auc(pos, neg):
    """P(a labelled drug scores above an unlabelled one), ties at 0.5."""
    v = np.concatenate([pos, neg])
    r = pd.Series(v).rank().values
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def hub_labels():
    h = pd.read_csv(ROOT / "data/raw/repurposing_hub/repurposing_drugs_20200324.txt", sep="\t", comment="!")
    h = h[h.clinical_phase.isin(CFG["sources"]["hub"]["phases"]) & h.indication.notna()]
    out = {}
    for r in h.itertuples():
        for coh, pat in CFG["cohorts"].items():
            if matches(r.indication, pat):
                out.setdefault(re.sub(r"[^a-z0-9]", "", r.pert_iname.lower()), set()).add(coh)
    return out


def chembl_labels(chembl_ids):
    """ChEMBL drug_indication for the FDA rows' molecules -> {chembl_id: {cohort}} (cached JSON)."""
    cache = ROOT / "data/raw/chembl/ChEMBL_37/indications"
    f = cache / "drug_indication.json"
    if f.exists():
        ind = json.load(open(f))
    else:
        live = C.get(f"{C.API}/status.json")["chembl_db_version"]
        assert live == C.RELEASE, f"ChEMBL API is {live}, config pins {C.RELEASE}"
        cache.mkdir(parents=True, exist_ok=True)
        ind = sum([C.paged("drug_indication", "drug_indications", {
            "molecule_chembl_id__in": ",".join(c),
            "only": "molecule_chembl_id,mesh_heading,efo_term,max_phase_for_ind"})
            for c in C.chunks(sorted(chembl_ids), 30)], [])
        json.dump(ind, open(f, "w"))
    lo = CFG["sources"]["chembl"]["min_phase"]
    out = {}
    for i in ind:
        try:
            if float(i.get("max_phase_for_ind") or 0) < lo:
                continue
        except (TypeError, ValueError):
            continue
        txt = " ".join(str(i.get(k) or "") for k in ("mesh_heading", "efo_term"))
        for coh, pat in CFG["cohorts"].items():
            if matches(txt, pat):
                out.setdefault(i["molecule_chembl_id"], set()).add(coh)
    return out, len(ind)


def tests(D, src, rng):
    """A1 and A2 for one label set over the score columns; D has columns cohort, drug, label, evidence + SCORES."""
    a1, a1s, a2, a2s = [], [], [], []
    for sub, name in ((D, "all"), (D[D.evidence == "matched"], "matched"), (D[D.evidence == "borrowed"], "borrowed")):
        per = []
        for coh, d in sub.groupby("cohort"):
            if d.label.sum() < MIN_LAB or (~d.label).sum() < MIN_UNLAB:
                continue
            row = dict(source=src, subset=name, cohort=coh, n_labelled=int(d.label.sum()), n_unlabelled=int((~d.label).sum()))
            for c in SCORES:
                v = d[c].values.astype(float)
                ok = ~np.isnan(v)
                if ok.sum() < MIN_LAB + MIN_UNLAB or d.label[ok].sum() < MIN_LAB or (~d.label[ok]).sum() < MIN_UNLAB:
                    row[f"auc_{c}"] = np.nan
                    continue
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
            obs = a.median()
            null = []
            for _ in range(N_PERM):
                vals = []
                for row, d in per:
                    v = d[c].values.astype(float)
                    ok = ~np.isnan(v)
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
    # A2: a labelled drug's percentile in its labelled cohorts vs elsewhere (lead score)
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
    D, _ = fda_rows()
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t")
    chem = ann.dropna(subset=["chembl_id"]).groupby("pert_id").chembl_id.apply(lambda s: set(s)).to_dict()
    dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False,
                     usecols=["group_id", "pert_id", "n_targets", "targets_expressed"]).drop_duplicates(["group_id", "pert_id"])
    nexp = dd.targets_expressed.fillna("").map(lambda x: len([t for t in str(x).split("|") if t]))
    dd["target_frac"] = np.where(dd.n_targets > 0, nexp / dd.n_targets.where(dd.n_targets > 0), np.nan)
    D = D.merge(dd[["group_id", "pert_id", "target_frac"]], on=["group_id", "pert_id"], how="left")
    D["cohort"] = D.group_id
    D["evidence"] = np.where(D.source == "direct", "matched", "borrowed")
    D["mean_pct"] = D.groupby("drug")[LEAD].transform("mean")
    D["key"] = D.cmap_name.map(lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower()))

    hub = hub_labels()
    cl, n_ind = chembl_labels({i for p in D.pert_id.unique() for i in chem.get(p, ())})
    D["lab_hub"] = [bool(hub.get(k, set()) & {c}) for k, c in zip(D.key, D.cohort)]
    D["lab_chembl"] = [any(c in cl.get(i, ()) for i in chem.get(p, ())) for p, c in zip(D.pert_id, D.cohort)]

    runs, a1, a1s, a2, a2s = {}, [], [], [], []
    for src in ("hub", "chembl"):
        d = D.assign(label=D[f"lab_{src}"])
        runs[src] = dict(labelled_rows=int(d.label.sum()), labelled_drugs=int(d.drug[d.label].nunique()),
                         cohorts_with_labels=int(d.groupby("cohort").label.sum().gt(0).sum()))
        r = tests(d, src, rng)
        a1 += r[0]; a1s += r[1]; a2 += r[2]; a2s += r[3]
    D[["cohort", "drug", "pert_id", "cmap_name", "evidence", "lab_hub", "lab_chembl"] + SCORES].to_csv(
        OUT / "labels.tsv", sep="\t", index=False)
    pd.DataFrame(a1).to_csv(OUT / "a1_per_cohort.tsv", sep="\t", index=False)
    pd.DataFrame(a1s).to_csv(OUT / "a1_summary.tsv", sep="\t", index=False)
    pd.DataFrame(a2).to_csv(OUT / "a2_per_drug.tsv", sep="\t", index=False)
    pd.DataFrame(a2s).to_csv(OUT / "a2_summary.tsv", sep="\t", index=False)
    cov = D.groupby("cohort").agg(n_drugs=("drug", "nunique"), hub=("lab_hub", "sum"), chembl=("lab_chembl", "sum"))
    cov.to_csv(OUT / "coverage.tsv", sep="\t")
    json.dump(dict(rows=len(D), drugs=int(D.drug.nunique()), cohorts=int(D.cohort.nunique()),
                   chembl_indication_records=n_ind, sources=runs, min_lab=MIN_LAB, min_unlab=MIN_UNLAB,
                   n_perm=N_PERM, n_boot=N_BOOT, seed=SEED), open(OUT / "run.json", "w"), indent=1)
    print(pd.DataFrame(a1s).to_string(index=False))
    print(pd.DataFrame(a2s).to_string(index=False))
    print(cov.to_string())


if __name__ == "__main__":
    main()
