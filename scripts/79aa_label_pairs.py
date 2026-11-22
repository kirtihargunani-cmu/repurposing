#!/usr/bin/env python
"""Labelled (unit, drug) pairs over the WHOLE scored universe, graded by clinical phase (D181, 2026-11-22).

WHY. `79l_known_use.py` builds labels for the validation metric only: FDA ingredient x cohort, 41 cohorts,
391 FDA-matched drugs, binary at phase >= 3. A learning layer needs the labels over the rows that actually
carry scores - 41 groups + 14 subgroups, every scored compound - and needs the PHASE, not a binary cut,
because phase 1 / 2 / 3 / 4 is an ordinal outcome and throwing it away discards most of the signal.
This script only LABELS. No score is read, nothing is fitted, no row is filtered.

UNIVERSE = every (unit, pert_id) row of drug_reversal<SFX>.tsv + stratum_reversal<SFX>.tsv.

SOURCES (config/known_use.yaml supplies the per-cohort regex and the subtype table; nothing here is derived
from our scores):
  ids      pert_id -> ChEMBL ids: annotation, D105 name lookup, D113 ids, the D181b name fetch
           (`scripts/79ab_fetch_missing_ids.py` -> d105_all/d181_ids.json), and parent forms.
  chembl   data/raw/chembl_indications/drug_indication.tsv (D178 fetch, ChEMBL 37: 60,055 molecule x
           indication rows with max_phase_for_ind and the NCT ids). A pair is labelled when any indication
           term of the drug (mesh_heading or efo_term) matches the cohort's include, not its exclude, and
           the shared cancer_term. `chembl_phase` = the HIGHEST max_phase_for_ind over the matching terms.
  hub      Repurposing Hub 2020-03-24 `indication` field, Launched compounds only -> `hub_label`.
  trials   n_nct summed over the matching ChEMBL indication rows -> `n_trials` (volume, not outcome: D176
           showed our ranking tracks trial volume and NOT trial success, so this is a covariate to control
           for, never a target).

DRUG IDENTITY. pert_id -> ChEMBL ids via the annotation `chembl_id`, the D105 name lookup, the D113 resolved
ids and ChEMBL parent forms (same resolution as 79l), so salt and parent forms are pooled.

SUBGROUPS inherit their parent cohort's regex (stated in known_use.yaml). The subtype table still applies:
a labelled drug of a listed class (endocrine, her2, lung_driver) keeps its label only in the allowed
siblings; in the other siblings the row is marked `subtype_excluded` and its label is cleared.

HOW TO READ THE OUTPUT, and the one thing that must not be done with it:
  * `label_phase` 4 / 3 / 2 / 1 / 0.5 = the highest phase this drug reached FOR THIS DISEASE (ChEMBL writes
    0.5 for early phase 1); 0 = no record.
  * **0 IS NOT A NEGATIVE.** It means "no trial record matched", which for a repurposing candidate is the
    expected state of a true positive that nobody has tested. Treating 0 as a negative class trains a model
    to reproduce what has already been tried, which is the opposite of the goal. Use positives-only /
    PU-style framing, or restrict negatives to drugs trialled for OTHER cancers and not this one
    (`tried_elsewhere`), which at least conditions on having been tried at all.
  * `tried_elsewhere` = the drug has a phase >= 1 cancer indication for some OTHER cohort of a different
    tissue key. That is the only defensible negative pool here, and it is still weak.

WINDOW (D182, user 2026-11-22: "the label table should carry the window columns"). A label is a
(disease, drug) fact and does NOT depend on timepoint - it is orthogonal to the window, the target gate and
target dependency, so nothing here is filtered by them. But the EVIDENCE behind a labelled pair may or may
not be timepoint-resolved, and that must not be invisible when the labels are used. Carried per row from the
split tables (79j / D180):
  window_coverage        both | early_only | late_only | unknown
  window_unknown_reason  no_window_run | no_condition_profile | no_condition_in_lines (unknown rows)
  has_window             window_coverage in {both, early_only, late_only}
Why so many labelled pairs have none: the window arm scores IN-KEY (matched) lines only, while 79e also
scores borrowed drugs from non-matched lines, and no subgroup has window data at all. Counting rows with
`has_window` therefore gives a much smaller - and honest - training set than the label count alone.

Env: SFX (default _v2_tf). Output: results/lincs/labels/pair_labels<SFX>.tsv.gz, label_counts<SFX>.tsv, run.json
"""
import json
import os
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
CH = ROOT / "data/raw/chembl/ChEMBL_37"
IND = ROOT / "data/raw/chembl_indications/drug_indication.tsv"
OUT = RES / "labels"
SFX = os.environ.get("SFX", "_v2_tf")
CFG = yaml.safe_load(open(ROOT / "config/known_use.yaml"))
CANCER = re.compile(CFG["cancer_term"], re.I)
norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def matches(term, pat):
    if not CANCER.search(term) or not re.search(pat["include"], term, re.I):
        return False
    return not (pat.get("exclude") and re.search(pat["exclude"], term, re.I))


def universe():
    """Every scored (unit, drug) row; `parent` is the cohort whose label regex a subgroup inherits."""
    g = pd.read_csv(RES / f"drug_reversal{SFX}.tsv", sep="\t", low_memory=False,
                    usecols=["group_id", "pert_id", "cmap_name", "source", "passes_gate"])
    g["unit_kind"], g["parent"] = "group", g.group_id
    s = pd.read_csv(RES / f"stratum_reversal{SFX}.tsv", sep="\t", low_memory=False,
                    usecols=["group_id", "pert_id", "cmap_name", "source", "passes_gate"])
    par = (pd.read_csv(RES / "stratum_drugs_v2.tsv", sep="\t", usecols=["group_id", "parent_group"], low_memory=False)
           .drop_duplicates("group_id").set_index("group_id").parent_group)
    s["unit_kind"], s["parent"] = "subgroup", s.group_id.map(par)
    d = pd.concat([g, s], ignore_index=True).rename(columns={"group_id": "unit_id"})
    assert d.parent.notna().all(), "subgroup with no parent_group"
    return d


def chembl_ids(perts, ann):
    """pert_id -> {ChEMBL ids}, POOLED BY DRUG NAME.

    D181 fix: the same drug appears under several BRD ids and usually only ONE of them carries the
    annotation (doxorubicin has 4 pert_ids, 1 with CHEMBL53463). Resolving per pert_id therefore lost
    69 of 555 approved drugs. Ids are collected per normalised name - as 79l pools an FDA ingredient's
    compounds - and every pert_id of that name gets the union.
    """
    # D181b: these caches are keyed by the RAW cmap_name ("CPCCOEt", "MLN-9708"), so a normalised lookup
    # silently misses. Key each one both ways and try both.
    def bothkeys(d):
        return {**{k: v for k, v in d.items()}, **{norm(k): v for k, v in d.items()}}
    nl = bothkeys(json.load(open(CH / "d105_all/name_lookup.json")))
    d113 = bothkeys(json.load(open(CH / "d105_all/d113_ids.json")) if (CH / "d105_all/d113_ids.json").exists() else {})
    d181 = bothkeys(json.load(open(CH / "d105_all/d181_ids.json")) if (CH / "d105_all/d181_ids.json").exists() else {})
    par = {f["molecule_chembl_id"]: f["parent_chembl_id"] for f in json.load(open(CH / "d105_all/forms.json"))}
    a = ann.drop_duplicates("pert_id").set_index("pert_id")
    names_of = {p: {n for n in (a.cmap_name.get(p), a.hub_name.get(p), a.chembl_pref_name.get(p))
                    if isinstance(n, str)} for p in perts}
    by_name = {}
    for p in perts:
        ids = {a.chembl_id.get(p)}
        for nm in names_of[p]:
            for k in (nm, norm(nm)):
                ids |= {nl.get(k)} | set(d113.get(k, {}).get("ids", [])) | set(d181.get(k, {}).get("ids", []))
        ids = {i for i in ids if isinstance(i, str)}
        for nm in names_of[p]:
            by_name.setdefault(norm(nm), set()).update(ids)
    out = {}
    for p in perts:
        ids = set().union(*[by_name.get(norm(nm), set()) for nm in names_of[p]]) if names_of[p] else set()
        out[p] = ids | {par.get(i, i) for i in ids}
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    D = universe()
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False)
    ids = chembl_ids(sorted(D.pert_id.unique()), ann)   # includes the D181b name fetch (79ab)
    # FDA-approved subsets (user 2026-11-22: "we want to only look at fda approved drugs first").
    # fda_matched is the strict one - a Drugs@FDA ingredient match, which is what the FDA rankings are built on.
    # approved_anywhere is the looser ChEMBL phase 4 / Hub Launched flag 79e uses. Both are reported.
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t", low_memory=False)
    fda_name = {norm(n) for n in fm.cmap_name.dropna()}
    a1 = ann.drop_duplicates("pert_id").set_index("pert_id")
    appr_name = {norm(n) for p, n in a1.cmap_name.items() if isinstance(n, str)
                 and (a1.chembl_max_phase.get(p) == 4 or a1.hub_clinical_phase.get(p) == "Launched")}
    D["fda_matched"] = D.cmap_name.map(lambda n: norm(n) in fda_name)
    D["approved_anywhere"] = D.cmap_name.map(lambda n: norm(n) in appr_name) | D.fda_matched
    D["fda_ingredient"] = D.cmap_name.map(
        fm.assign(k=fm.cmap_name.map(norm)).drop_duplicates("k").set_index("k").fda_ingredient.to_dict().get)
    cohorts = {c: p for c, p in CFG["cohorts"].items()}
    keys = sorted(set(D.parent.unique()) & set(cohorts))
    missing = sorted(set(D.parent.unique()) - set(cohorts))
    tkey = (pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#")
            .set_index("group_id").key.to_dict())

    # ---- ChEMBL: (molecule, cohort) -> highest phase, and the trial count behind it
    ind = pd.read_csv(IND, sep="\t", low_memory=False)
    ind["ph"] = pd.to_numeric(ind.max_phase_for_ind, errors="coerce").fillna(0).clip(lower=0)
    term_cols = ["mesh_heading", "efo_term"]
    terms = sorted({t for c in term_cols for t in ind[c].dropna().unique() if CANCER.search(str(t))})
    term_coh = {t: {c for c in keys if matches(str(t), cohorts[c])} for t in terms}
    term_coh = {t: c for t, c in term_coh.items() if c}
    hit = ind[ind[term_cols].isin(set(term_coh)).any(axis=1)]   # set(): DataFrame.isin(dict) keys on COLUMNS
    phase, ntr = {}, {}
    for r in hit.itertuples():
        cs = set().union(*[term_coh.get(getattr(r, c), set()) for c in term_cols])
        for mol in {r.molecule_chembl_id, r.parent_molecule_chembl_id}:
            if not isinstance(mol, str):
                continue
            for c in cs:
                k = (mol, c)
                phase[k] = max(phase.get(k, 0), float(r.ph))
                ntr[k] = ntr.get(k, 0) + int(r.n_nct or 0)

    # ---- Repurposing Hub: normalised name -> cohorts (Launched only)
    h = pd.read_csv(ROOT / "data/raw/repurposing_hub/repurposing_drugs_20200324.txt", sep="\t", comment="!")
    h = h[h.clinical_phase.isin(CFG["sources"]["hub"]["phases"]) & h.indication.notna()]
    hub = {}
    for r in h.itertuples():
        for term in str(r.indication).split("|"):
            cs = {c for c in keys if matches(term, cohorts[c])}
            if cs:
                hub.setdefault(norm(r.pert_iname), set()).update(cs)
    name_of = ann.drop_duplicates("pert_id").set_index("pert_id")

    # ---- per row
    def row_phase(p, parent):
        return max([phase.get((i, parent), 0) for i in ids[p]] or [0])

    def row_trials(p, parent):
        return sum(ntr.get((i, parent), 0) for i in ids[p])

    def row_hub(p, parent):
        ns = {norm(n) for n in (name_of.cmap_name.get(p), name_of.hub_name.get(p)) if isinstance(n, str)}
        return any(parent in hub.get(n, set()) for n in ns)

    D["chembl_phase"] = [row_phase(p, c) for p, c in zip(D.pert_id, D.parent)]
    D["n_trials"] = [row_trials(p, c) for p, c in zip(D.pert_id, D.parent)]
    D["hub_label"] = [row_hub(p, c) for p, c in zip(D.pert_id, D.parent)]

    # tried_elsewhere: phase >= 1 cancer indication in a cohort of a DIFFERENT tissue key
    any_coh = {}
    for (mol, c), ph in phase.items():
        if ph >= 1:
            any_coh.setdefault(mol, set()).add(c)
    def elsewhere(p, parent):
        cs = set().union(*[any_coh.get(i, set()) for i in ids[p]]) if ids[p] else set()
        return any(tkey.get(c) != tkey.get(parent) for c in cs)
    D["tried_elsewhere"] = [elsewhere(p, c) for p, c in zip(D.pert_id, D.parent)]

    # ---- subtype restriction (endocrine / her2 / lung_driver): clear the label outside the allowed siblings
    D["subtype_excluded"] = ""
    sub = ann.drop_duplicates("pert_id").set_index("pert_id")
    for st in CFG.get("subtype", []):
        tg, drugs = set(st.get("main_targets", [])), {norm(x) for x in st.get("drugs", [])}
        sibs, allow = set(st["siblings"]), set(st["allowed"])
        def is_class(p):
            t = str(sub.hub_target.get(p, "") or "") + "|" + str(sub.lincs_target.get(p, "") or "")
            return bool(tg & set(re.split(r"[|,;\s]+", t))) or norm(sub.cmap_name.get(p, "")) in drugs
        m = D.parent.isin(sibs - allow) & D.pert_id.map(is_class) & (D.chembl_phase.gt(0) | D.hub_label)
        D.loc[m, ["chembl_phase", "hub_label"]] = [0.0, False]
        D.loc[m, "subtype_excluded"] = st["name"]

    # ---- window evidence (D182): carried, never used to filter or to define a label
    wcols = ["window_coverage", "window_unknown_reason"]
    wp = []
    for f, unit in ((RES / f"drug_reversal_split{SFX}.tsv", "group_id"),
                    (RES / f"stratum_reversal_split{SFX}.tsv", "stratum")):
        if not f.exists():
            continue
        h = pd.read_csv(f, sep="\t", nrows=0).columns
        use = [unit, "pert_id"] + [c for c in wcols if c in h]
        x = pd.read_csv(f, sep="\t", usecols=use, low_memory=False).rename(columns={unit: "unit_id"})
        wp.append(x.drop_duplicates(["unit_id", "pert_id"]))
    if wp:
        D = D.merge(pd.concat(wp, ignore_index=True), on=["unit_id", "pert_id"], how="left")
    for c in wcols:
        if c not in D.columns:
            D[c] = pd.NA
    D["window_coverage"] = D.window_coverage.fillna("unknown")
    D["has_window"] = D.window_coverage.isin(["both", "early_only", "late_only"])

    D["label_phase"] = D[["chembl_phase"]].max(axis=1)
    D.loc[D.hub_label & D.label_phase.lt(4), "label_phase"] = 4.0   # Hub Launched = approved for it
    D["labelled"] = D.label_phase > 0
    D["label_source"] = ["both" if h and c else "hub" if h else "chembl" if c else ""
                         for h, c in zip(D.hub_label, D.chembl_phase > 0)]
    D["has_chembl_id"] = [bool(ids[p]) for p in D.pert_id]

    # ---- counts
    f = OUT / f"pair_labels{SFX}.tsv.gz"
    D.to_csv(f, sep="\t", index=False)
    rows = []
    gp = D[D.passes_gate.astype(str).isin(["True", "TRUE", "1"])]
    for kind, g in [("all", D), ("group", D[D.unit_kind == "group"]), ("subgroup", D[D.unit_kind == "subgroup"]),
                    ("gate-passing", gp),
                    ("fda_matched", D[D.fda_matched]),
                    ("fda_matched+gate", gp[gp.fda_matched]),
                    ("fda_matched group", D[D.fda_matched & D.unit_kind.eq("group")]),
                    ("fda_matched subgroup", D[D.fda_matched & D.unit_kind.eq("subgroup")]),
                    ("approved_anywhere", D[D.approved_anywhere])]:
        r = dict(subset=kind, pairs=len(g), drugs=g.pert_id.nunique(), units=g.unit_id.nunique(),
                 labelled=int(g.labelled.sum()), labelled_drugs=int(g[g.labelled].pert_id.nunique()),
                 units_with_label=int(g[g.labelled].unit_id.nunique()),
                 tried_elsewhere=int(g.tried_elsewhere.sum()),
                 has_window=int(g.has_window.sum()),
                 labelled_with_window=int((g.labelled & g.has_window).sum()),
                 labelled_both_windows=int((g.labelled & g.window_coverage.eq("both")).sum()))
        for p in (0.5, 1, 2, 3, 4):   # ChEMBL uses 0.5 for early phase 1
            r[f"phase{p}"] = int((g.label_phase == p).sum())
        rows.append(r)
    cnt = pd.DataFrame(rows)
    cnt.to_csv(OUT / f"label_counts{SFX}.tsv", sep="\t", index=False)
    print(cnt.to_string(index=False))
    print(f"\nby source: {D[D.labelled].label_source.value_counts().to_dict()}")
    L = D[D.labelled]
    print(f"window evidence behind the labels: {int(L.has_window.sum()):,} of {len(L):,} labelled rows "
          f"({100 * L.has_window.mean():.1f}%); both windows {int(L.window_coverage.eq('both').sum()):,}")
    print(f"  labelled rows with NO window, by reason: {L[~L.has_window].window_unknown_reason.value_counts().to_dict()}")
    uq = D.drop_duplicates("pert_id")
    print(f"drugs with any ChEMBL id: {uq.has_chembl_id.sum():,} of {len(uq):,}  "
          f"(FDA-matched {int(uq[uq.fda_matched].has_chembl_id.sum()):,} of {int(uq.fda_matched.sum()):,}; "
          f"approved-anywhere {int(uq[uq.approved_anywhere].has_chembl_id.sum()):,} of {int(uq.approved_anywhere.sum()):,})")
    if missing:
        print(f"units whose parent has NO label regex in known_use.yaml ({len(missing)}): {missing}")
    per_unit = (D.groupby(["unit_kind", "unit_id"])
                .agg(pairs=("labelled", "size"), labelled=("labelled", "sum"),
                     phase4=("label_phase", lambda v: int((v == 4).sum())),
                     phase3plus=("label_phase", lambda v: int((v >= 3).sum())),
                     tried_elsewhere=("tried_elsewhere", "sum"))
                .reset_index().sort_values("labelled", ascending=False))
    per_unit.to_csv(OUT / f"label_counts_per_unit{SFX}.tsv", sep="\t", index=False)
    top = (D[D.labelled].groupby("unit_id").labelled.size().sort_values(ascending=False))
    print(f"\nlabelled pairs per unit: median {int(top.median())}, min {int(top.min())}, max {int(top.max())}; "
          f"{D.unit_id.nunique() - len(top)} of {D.unit_id.nunique()} units have none")
    json.dump(dict(sfx=SFX, pairs=len(D), drugs=int(D.pert_id.nunique()), units=int(D.unit_id.nunique()),
                   labelled=int(D.labelled.sum()), chembl_indication_rows=len(ind),
                   cancer_terms_matched=len(term_coh), units_without_regex=missing),
              open(OUT / f"run{SFX}.json", "w"), indent=1)
    print(f"\nwrote {f.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
