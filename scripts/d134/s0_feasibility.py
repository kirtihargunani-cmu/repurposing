#!/usr/bin/env python
"""D134 stage S0 -- FEASIBILITY AUDIT for docs/SPEC_target_recoverability.md. Counts only.

NO RETRIEVAL IS COMPUTED HERE. No similarity between two different compounds' profiles, no ranking, no AUROC,
no permutation null, no outcome label (CTRP / PRISM / GDSC) read, no production score touched, no threshold
chosen. Every quantity is either an annotation count, a presence count, or a PER-PROFILE statistic. The point
is that this output cannot be fitted to a result, because no result exists yet.

All definitions were fixed in docs/DECISIONS.md row D134-S0, committed and pushed before this ran.

Inputs
  results/lincs/target_role/d121/units.tsv          target units (rows with non-empty tid)
  results/lincs/drug_annotations.tsv                inchi_key (skeleton hash = first 14 chars)
  data/processed/lincs/drug_line_profiles_index.tsv row order of the two consensus matrices
  data/processed/lincs/drug_line_profiles_bing.f32          raw arm      (9,678 x 10,112 float32, memory-mapped)
  data/processed/lincs/drug_line_profiles_bing_specific.f32 specific arm (same shape, memory-mapped)
  data/processed/lincs/level5_bing_hiq_sigs.tsv.gz  pert_id -> cmap_name, per-signature TAS

Outputs (results/lincs/d134/)
  strength_profiles.tsv      per (cmap_name, ModelID): per-profile strength measures
  feasibility.tsv            per target unit: annotation, shared-line, analog-collapsed and strength counts
  chance_by_line_target.tsv  per (line, unit) analytic per-query chance, both candidate pools
  attrition.tsv              testable-unit counts at every filter and every candidate strength floor
  d134_s0_run.json           headline numbers, provenance and every definition actually used
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "8"

import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
PROC = ROOT / "data/processed/lincs"
RES = ROOT / "results/lincs"
OUT = RES / "d134"
OUT.mkdir(parents=True, exist_ok=True)

N_GENES = 10112          # drug_line_profiles_genes.tsv, checked below
CHUNK = 512
NZ2_FLOORS = [0, 50, 100, 200, 400, 800]        # candidate floors, COSTED NOT CHOSEN (D134-S0(6))
TAS_FLOORS = [0.4, 0.5, 0.6, 0.7]               # 0.4 is the floor the universe already carries
CHANCE_FLAG = 0.05                              # reporting flag only, gates nothing (D134-S0(7))


def qstat(x, name):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0, "note": f"{name}: no finite values"}
    q = np.percentile(x, [0, 5, 25, 50, 75, 95, 100])
    return {"n": int(x.size), "mean": round(float(x.mean()), 4),
            "min": round(float(q[0]), 4), "p05": round(float(q[1]), 4), "p25": round(float(q[2]), 4),
            "median": round(float(q[3]), 4), "p75": round(float(q[4]), 4), "p95": round(float(q[5]), 4),
            "max": round(float(q[6]), 4)}


def unit_stats(members, present, skel):
    """members: {tid: [cmap_name]}; present: {cmap_name: set(ModelID)}.

    Returns one row per unit. A PAIR {A,B} is shared-line testable iff A and B both have a surviving profile in
    at least one common line. A CROSS-SKELETON pair additionally has different InChIKey skeletons: two analogs
    are not two independent pieces of evidence (D134-S0(5)).
    """
    rows = []
    for tid, comps in members.items():
        comps = sorted(comps)
        prof = [c for c in comps if present.get(c)]
        n_pairs = n_sl = n_sl_cross = 0
        in_sl, in_sl_cross, lines_sl, lines_sl_cross = set(), set(), set(), set()
        for a, b in combinations(prof, 2):
            n_pairs += 1
            sh = present[a] & present[b]
            if not sh:
                continue
            n_sl += 1
            in_sl.update((a, b))
            lines_sl |= sh
            if skel[a] != skel[b]:
                n_sl_cross += 1
                in_sl_cross.update((a, b))
                lines_sl_cross |= sh
        n_pairs_annot = len(comps) * (len(comps) - 1) // 2
        rows.append(dict(
            tid=tid,
            n_comp_annot=len(comps),
            n_comp_profiled=len(prof),
            n_eff_annot=len({skel[c] for c in comps}),
            n_eff_profiled=len({skel[c] for c in prof}),
            n_pairs_annot=n_pairs_annot,
            n_pairs_profiled=n_pairs,
            n_pairs_shared_line=n_sl,
            n_pairs_shared_line_cross_skeleton=n_sl_cross,
            n_comp_shared_line=len(in_sl),
            n_eff_shared_line=len({skel[c] for c in in_sl}),
            n_comp_shared_line_cross=len(in_sl_cross),
            n_eff_shared_line_cross=len({skel[c] for c in in_sl_cross}),
            n_lines_any=len(set().union(*[present[c] for c in prof])) if prof else 0,
            n_lines_shared_line=len(lines_sl),
            n_lines_shared_line_cross=len(lines_sl_cross),
            testable_annot=len(comps) >= 2,
            testable_shared_line=len(in_sl) >= 2,
            testable_shared_line_eff=len({skel[c] for c in in_sl_cross}) >= 2,
        ))
    return pd.DataFrame(rows)


def main():
    t0 = time.time()
    prov, notes = {}, []

    # ---------------------------------------------------------------- 1. annotation level (spec Section 2)
    u_all = pd.read_csv(RES / "target_role/d121/units.tsv", sep="\t")
    u = u_all[u_all.tid.notna() & (u_all.tid.astype(str).str.strip() != "")].copy()
    prov["units_rows_total"] = int(len(u_all))
    prov["units_rows_dropped_empty_tid"] = int(len(u_all) - len(u))
    members = {t: sorted(set(g)) for t, g in u.groupby("tid").cmap_name}
    comps_annot = sorted(set(u.cmap_name))
    sec2 = dict(compounds=len(comps_annot), units=len(members),
                genes=len({g for s in u.genes.dropna() for g in str(s).split("|")}),
                units_ge2=sum(len(v) >= 2 for v in members.values()),
                units_ge3=sum(len(v) >= 3 for v in members.values()),
                units_ge5=sum(len(v) >= 5 for v in members.values()),
                singletons=sum(len(v) == 1 for v in members.values()))
    sec2["compounds_in_units_ge2"] = len({c for v in members.values() if len(v) >= 2 for c in v})
    prov["section2_reproduced"] = sec2

    # ---------------------------------------------------------------- 2. profile universe
    idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t")
    genes = pd.read_csv(PROC / "drug_line_profiles_genes.tsv", sep="\t")
    assert len(genes) == N_GENES, f"gene count {len(genes)} != {N_GENES}"
    sig = pd.read_csv(PROC / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
    p2c = sig.drop_duplicates("pert_id").set_index("pert_id").cmap_name
    idx["cmap_name"] = idx.pert_id.map(p2c)
    assert idx.cmap_name.notna().all()
    prov["profiles"] = dict(rows=int(len(idx)), pert_ids=int(idx.pert_id.nunique()),
                            compounds=int(idx.cmap_name.nunique()), lines=int(idx.ModelID.nunique()),
                            genes=N_GENES)
    # what the TAS >= 0.4 floor the universe ALREADY carries has cost (D134-S0(2))
    pre = sig.groupby(["pert_id", "ModelID"]).size()
    post = sig[sig.tas >= 0.4].groupby(["pert_id", "ModelID"]).size()
    prov["already_applied_tas_floor"] = dict(
        definition="79c keeps is_hiq signatures with per-signature tas >= 0.4; a (pert_id, line) cell survives "
                   "if >= 1 of its signatures does",
        cells_hiq_only=int(len(pre)), cells_after_tas_0p4=int(len(post)),
        cells_lost=int(len(pre) - len(post)),
        frac_lost=round(1 - len(post) / len(pre), 4))

    # ---------------------------------------------------------------- 3. per-profile strength (NO comparison)
    raw = np.memmap(PROC / "drug_line_profiles_bing.f32", dtype=np.float32, mode="r", shape=(len(idx), N_GENES))
    spc = np.memmap(PROC / "drug_line_profiles_bing_specific.f32", dtype=np.float32, mode="r",
                    shape=(len(idx), N_GENES))
    nz2 = np.empty(len(idx), np.int32); nz4 = np.empty(len(idx), np.int32)
    mab = np.empty(len(idx), np.float32); nz2s = np.empty(len(idx), np.int32)
    mabs = np.empty(len(idx), np.float32)
    for i in range(0, len(idx), CHUNK):
        j = min(i + CHUNK, len(idx))
        a = np.abs(np.asarray(raw[i:j], dtype=np.float32))
        nz2[i:j] = (a >= 2).sum(1); nz4[i:j] = (a >= 4).sum(1); mab[i:j] = a.mean(1)
        b = np.abs(np.asarray(spc[i:j], dtype=np.float32))
        nz2s[i:j] = (b >= 2).sum(1); mabs[i:j] = b.mean(1)
    idx["nz2"], idx["nz4"], idx["mean_abs"] = nz2, nz4, mab
    idx["nz2_specific"], idx["mean_abs_specific"] = nz2s, mabs
    tas = sig[sig.tas >= 0.4].groupby(["pert_id", "ModelID"]).tas.agg(["max", "mean"])
    idx = idx.merge(tas.rename(columns={"max": "tas_max", "mean": "tas_mean"}),
                    left_on=["pert_id", "ModelID"], right_index=True, how="left")
    prov["tas_available"] = dict(source="level5_bing_hiq_sigs.tsv.gz (from raw/lincs2020/siginfo_beta.txt)",
                                 profiles_without_tas=int(idx.tas_max.isna().sum()))

    # a compound's value in a line is the MAX over its pert_ids there (D134-S0(2))
    cl = idx.groupby(["cmap_name", "ModelID"]).agg(
        n_pert_ids=("pert_id", "nunique"), n_sig=("n_sig", "sum"),
        nz2=("nz2", "max"), nz4=("nz4", "max"), mean_abs=("mean_abs", "max"),
        nz2_specific=("nz2_specific", "max"), mean_abs_specific=("mean_abs_specific", "max"),
        tas_max=("tas_max", "max"), tas_mean=("tas_mean", "max")).reset_index()
    cl["target_annotated"] = cl.cmap_name.isin(comps_annot)
    cl.to_csv(OUT / "strength_profiles.tsv", sep="\t", index=False)

    ann_cl = cl[cl.target_annotated]
    prov["strength_distribution"] = {
        "all_profiles": {k: qstat(cl[k], k) for k in
                         ["nz2", "nz4", "mean_abs", "nz2_specific", "tas_max"]},
        "target_annotated_profiles": {k: qstat(ann_cl[k], k) for k in
                                      ["nz2", "nz4", "mean_abs", "nz2_specific", "tas_max"]},
        "note": "per-profile only; no profile is compared with any other profile",
    }

    # ---------------------------------------------------------------- 4. InChIKey skeletons
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False)
    ann["skel"] = ann.inchi_key.astype(str).str.slice(0, 14).where(ann.inchi_key.notna())
    prof_pids = set(idx.pert_id)
    sub = ann[ann.pert_id.isin(prof_pids) & ann.cmap_name.isin(comps_annot)]
    skel, no_key, disagree = {}, [], []
    for c, g in sub.groupby("cmap_name"):
        s = g.skel.dropna()
        if s.empty:
            skel[c] = f"NOKEY::{c}"; no_key.append(c); continue
        vc = s.value_counts()
        skel[c] = vc.index[0]
        if len(vc) > 1:
            disagree.append(c)
    for c in comps_annot:                      # compounds with no annotation row at all
        if c not in skel:
            skel[c] = f"NOKEY::{c}"; no_key.append(c)
    prov["skeletons"] = dict(
        definition="first 14 characters of inchi_key in results/lincs/drug_annotations.tsv",
        compounds=len(skel), distinct_skeletons=len(set(skel.values())),
        compounds_without_inchikey=sorted(set(no_key)),
        n_compounds_without_inchikey=len(set(no_key)),
        compounds_with_disagreeing_pert_ids=sorted(disagree),
        handling="a compound with no InChIKey gets its own singleton skeleton and is flagged; it is never "
                 "imputed into another skeleton (D134-S0(5)/(9))")
    coll = {}
    for c, k in skel.items():
        coll.setdefault(k, []).append(c)
    prov["skeletons"]["collapsed_groups"] = {k: sorted(v) for k, v in coll.items() if len(v) > 1}
    prov["skeletons"]["n_collapsed_groups"] = sum(len(v) > 1 for v in coll.values())
    prov["skeletons"]["caveat"] = (
        "MEASURED LIMIT OF THE SPEC'S CONTROL 5.1 AS WRITTEN. The first 14 characters of an InChIKey are the "
        "hash of the molecular SKELETON (connectivity). Two compounds share it only when they are the same "
        "connectivity differing in salt form, stereochemistry, isotope or charge. Structural ANALOGS in the "
        "ordinary medicinal-chemistry sense - two statins, two tinibs, two taxanes - have DIFFERENT first "
        "blocks and are NOT collapsed by this rule. The compound set is already keyed by cmap_name, so "
        "near-duplicates are mostly gone before the rule is applied. The effective-compound count below is "
        "therefore a near-duplicate check, not an analog check, and it should be read as a LOWER BOUND on "
        "analog redundancy. A real analog control needs a structural fingerprint (e.g. Morgan / Tanimoto) and "
        "a similarity cut-off: a new dependency and a threshold, both of which belong to S1, not to S0.")

    # ---------------------------------------------------------------- 5. presence and the filters
    def presence(mask=None):
        d = cl if mask is None else cl[mask]
        d = d[d.target_annotated]
        return {c: set(g) for c, g in d.groupby("cmap_name").ModelID}

    pres_all = presence()
    base = unit_stats(members, pres_all, skel)
    base = base.merge(u.drop_duplicates("tid")[["tid", "ttype", "genes"]], on="tid", how="left")

    # ---------------------------------------------------------------- 6. candidate strength floors (COSTED)
    attr, floor_unit = [], {}

    def testable_sets(st, pres, col="n_eff_shared_line_cross"):
        """distinct compounds and distinct lines carried by the units still testable under `col` >= 2.

        Counted as SETS, never as a sum of per-unit counts: a compound or a line that serves two target units
        must be counted once (the earlier draft summed and over-counted).
        """
        comps, lines = set(), set()
        for tid, r in st.set_index("tid").iterrows():
            if r[col] >= 2:
                cs = [c for c in members[tid] if c in pres]
                comps |= set(cs)
                for a, b in combinations(sorted(cs), 2):
                    if skel[a] != skel[b]:
                        lines |= pres[a] & pres[b]
        return comps, lines

    attr.append(dict(step="1 annotation (spec Section 2)", filter="target unit has >= 2 annotated compounds",
                     profiles_kept=None, compounds_kept=len(comps_annot),
                     units_testable=sec2["units_ge2"],
                     units_testable_eff=int((base.n_eff_annot >= 2).sum()),
                     pairs=int(base.n_pairs_annot.sum()),
                     compounds_in_testable=sec2["compounds_in_units_ge2"],
                     lines_in_testable=None, unit_line_cells_testable=None))
    attr.append(dict(step="2 profiled", filter="compound has >= 1 consensus profile",
                     profiles_kept=int(len(cl[cl.target_annotated])), compounds_kept=int(len(pres_all)),
                     units_testable=int((base.n_comp_profiled >= 2).sum()),
                     units_testable_eff=int((base.n_eff_profiled >= 2).sum()),
                     pairs=int(base.n_pairs_profiled.sum()),
                     compounds_in_testable=int(len({c for t in members if base.set_index("tid").at[t, "n_comp_profiled"] >= 2
                                                    for c in members[t] if c in pres_all})),
                     lines_in_testable=None, unit_line_cells_testable=None))
    attr.append(dict(step="3 shared line", filter="pair shares >= 1 cell line",
                     profiles_kept=int(len(cl[cl.target_annotated])), compounds_kept=int(len(pres_all)),
                     units_testable=int(base.testable_shared_line.sum()),
                     units_testable_eff=int((base.n_eff_shared_line >= 2).sum()),
                     pairs=int(base.n_pairs_shared_line.sum()),
                     compounds_in_testable=len(testable_sets(base, pres_all, "n_comp_shared_line")[0]),
                     lines_in_testable=len(testable_sets(base, pres_all, "n_comp_shared_line")[1]),
                     unit_line_cells_testable=int(base.loc[base.testable_shared_line, "n_lines_shared_line"].sum())))
    attr.append(dict(step="4 shared line + analog collapse",
                     filter="pair shares a line AND the two compounds have different InChIKey skeletons",
                     profiles_kept=int(len(cl[cl.target_annotated])), compounds_kept=int(len(pres_all)),
                     units_testable=int((base.n_comp_shared_line_cross >= 2).sum()),
                     units_testable_eff=int(base.testable_shared_line_eff.sum()),
                     pairs=int(base.n_pairs_shared_line_cross_skeleton.sum()),
                     compounds_in_testable=len(testable_sets(base, pres_all)[0]),
                     lines_in_testable=len(testable_sets(base, pres_all)[1]),
                     unit_line_cells_testable=int(base.loc[base.testable_shared_line_eff,
                                                    "n_lines_shared_line_cross"].sum())))

    for col, floors in (("nz2", NZ2_FLOORS), ("tas_max", TAS_FLOORS)):
        for f in floors:
            m = cl[col] >= f
            p = presence(m)
            st = unit_stats(members, p, skel)
            key = f"{col}>={f}"
            floor_unit[key] = st.set_index("tid")
            ma = cl.target_annotated & m
            attr.append(dict(
                step="5 shared line + analog collapse + candidate strength floor", filter=key,
                profiles_kept=int(ma.sum()), compounds_kept=int(len(p)),
                units_testable=int((st.n_comp_shared_line_cross >= 2).sum()),
                units_testable_eff=int(st.testable_shared_line_eff.sum()),
                pairs=int(st.n_pairs_shared_line_cross_skeleton.sum()),
                compounds_in_testable=len(testable_sets(st, p)[0]),
                lines_in_testable=len(testable_sets(st, p)[1]),
                unit_line_cells_testable=int(st.loc[st.testable_shared_line_eff, "n_lines_shared_line_cross"].sum()),
                profiles_dropped=int((cl.target_annotated & ~m).sum())))
    pd.DataFrame(attr).to_csv(OUT / "attrition.tsv", sep="\t", index=False)

    # bottom line: the shape of the surviving set, with no floor and at illustrative floor points
    bl = {}
    for key, st in [("no strength floor", base.copy())] + [(k, v.reset_index()) for k, v in floor_unit.items()]:
        t = st[st.n_eff_shared_line_cross >= 2]
        pres = pres_all if key == "no strength floor" else presence(cl[key.split(">=")[0]] >=
                                                                    float(key.split(">=")[1]))
        comps, lines = testable_sets(st, pres)
        bl[key] = dict(units_testable=int(len(t)),
                       units_eff_ge3=int((t.n_eff_shared_line_cross >= 3).sum()),
                       units_eff_ge5=int((t.n_eff_shared_line_cross >= 5).sum()),
                       eff_compounds_per_unit=qstat(t.n_eff_shared_line_cross, "eff"),
                       lines_per_unit=qstat(t.n_lines_shared_line_cross, "lines"),
                       distinct_compounds=len(comps), distinct_lines=len(lines))
    prov["bottom_line"] = bl
    prov["bottom_line_note"] = ("floor points are ILLUSTRATIVE and costed only; no floor is chosen at S0 "
                                "(D134-S0(6)). 'eff' = distinct InChIKey skeletons among the unit's compounds "
                                "that appear in at least one cross-skeleton shared-line pair.")

    # per-unit strength summary + per-floor survival, into feasibility.tsv
    cl_a = cl[cl.target_annotated]
    c2t = u.groupby("cmap_name").tid.apply(list)
    rec = []
    for _, r in cl_a.iterrows():
        for t in c2t.get(r.cmap_name, []):
            rec.append((t, r.nz2, r.tas_max))
    sdf = pd.DataFrame(rec, columns=["tid", "nz2", "tas_max"])
    agg = sdf.groupby("tid").agg(nz2_median=("nz2", "median"), nz2_min=("nz2", "min"),
                                 nz2_max=("nz2", "max"), tas_max_median=("tas_max", "median")).reset_index()
    base = base.merge(agg, on="tid", how="left")
    for key, st in floor_unit.items():
        base[f"testable_{key.replace('>=', '_ge_')}"] = base.tid.map(st.testable_shared_line_eff)
        base[f"n_eff_{key.replace('>=', '_ge_')}"] = base.tid.map(st.n_eff_shared_line_cross)
    base = base.sort_values(["testable_shared_line_eff", "n_eff_shared_line_cross", "n_comp_annot"],
                            ascending=False)
    base.to_csv(OUT / "feasibility.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- 7. chance baselines (spec Section 4)
    pool_all = cl.groupby("ModelID").cmap_name.nunique()          # every profiled compound in the line
    pool_ann = cl_a.groupby("ModelID").cmap_name.nunique()        # annotated-only sensitivity pool
    ch = []
    for tid, comps in members.items():
        byline = {}
        for c in comps:
            for L in pres_all.get(c, ()):
                byline.setdefault(L, []).append(c)
        for L, cs in byline.items():
            if len(cs) < 2:
                continue
            na, nb = int(pool_all[L]), int(pool_ann[L])
            ch.append(dict(ModelID=L, tid=tid, n_comp_in_line=len(cs),
                           n_eff_in_line=len({skel[c] for c in cs}),
                           n_candidates_all=na, n_candidates_annotated=nb,
                           chance_all=(len(cs) - 1) / (na - 1) if na > 1 else np.nan,
                           chance_annotated=(len(cs) - 1) / (nb - 1) if nb > 1 else np.nan))
    chd = pd.DataFrame(ch)
    if len(chd):
        chd["flag_uninformative_all"] = chd.chance_all >= CHANCE_FLAG
        chd["flag_uninformative_annotated"] = chd.chance_annotated >= CHANCE_FLAG
        chd = chd.sort_values("chance_all", ascending=False)
    chd.to_csv(OUT / "chance_by_line_target.tsv", sep="\t", index=False)
    prov["chance"] = dict(
        definition="(compounds in line L carrying T, excluding self) / (candidate compounds in L - 1); "
                   "cells with >= 2 of T's compounds in L only",
        cells=int(len(chd)),
        distinct_lines=int(chd.ModelID.nunique()) if len(chd) else 0,
        distinct_units=int(chd.tid.nunique()) if len(chd) else 0,
        candidate_pool_all=qstat(pool_all.values, "pool_all"),
        candidate_pool_annotated=qstat(pool_ann.values, "pool_ann"),
        chance_all=qstat(chd.chance_all, "chance_all") if len(chd) else {},
        chance_annotated=qstat(chd.chance_annotated, "chance_annotated") if len(chd) else {},
        n_flag_all=int(chd.flag_uninformative_all.sum()) if len(chd) else 0,
        n_flag_annotated=int(chd.flag_uninformative_annotated.sum()) if len(chd) else 0,
        flag_threshold=CHANCE_FLAG,
        flag_note="reporting flag only; gates nothing, chooses nothing")

    # ---------------------------------------------------------------- 8. Section 2 coverage cross-check
    lines_per_comp = pd.Series({c: len(v) for c, v in pres_all.items()})
    testable_comps = {c for v in members.values() if len(v) >= 2 for c in v}
    lp_t = lines_per_comp[lines_per_comp.index.isin(testable_comps)]
    annt = ann[ann.cmap_name.isin(testable_comps)]
    annsum = annt.groupby("cmap_name").n_active_lines.sum()
    prov["coverage_crosscheck"] = dict(
        spec_claim="median 4 active cell lines (IQR 1-8); 206 testable compounds active in >= 3 lines, 163 in >= 5",
        distinct_lines_per_compound=dict(median=float(lp_t.median()), p25=float(lp_t.quantile(.25)),
                                         p75=float(lp_t.quantile(.75)), ge3=int((lp_t >= 3).sum()),
                                         ge5=int((lp_t >= 5).sum()), n=int(len(lp_t))),
        sum_n_active_lines_over_pert_ids=dict(median=float(annsum.median()), p25=float(annsum.quantile(.25)),
                                              p75=float(annsum.quantile(.75)), ge3=int((annsum >= 3).sum()),
                                              ge5=int((annsum >= 5).sum()), n=int(len(annsum))),
        n_active_lines_per_pert_id=dict(median=float(annt.n_active_lines.median()),
                                        p25=float(annt.n_active_lines.quantile(.25)),
                                        p75=float(annt.n_active_lines.quantile(.75)),
                                        ge3=int((annt.n_active_lines >= 3).sum()),
                                        ge5=int((annt.n_active_lines >= 5).sum()), n=int(len(annt))),
        verdict="NOT REPRODUCED",
        note="The spec Section 2 coverage line matches NONE of the three readings computed here: distinct lines "
             "per cmap_name, n_active_lines summed over a compound's pert_ids, or n_active_lines per pert_id. "
             "Reported as unreproduced rather than guessed (D134-S0(9)); whoever drafted the spec should say "
             "which object was counted. The figure retrieval can actually use is DISTINCT LINES PER cmap_name, "
             "because a compound is one cmap_name and the same line profiled under two pert_ids of the same "
             "drug is one line, not two.")

    prov["constraints_held"] = ["no retrieval computed", "no between-compound similarity computed",
                                "no ranking", "no AUROC", "no permutation null",
                                "no CTRP / PRISM / GDSC or any outcome label read",
                                "no production score changed, nothing re-ranked",
                                "no threshold chosen; candidate floors costed only"]
    prov["threads"] = 8
    prov["runtime_s"] = round(time.time() - t0, 1)
    json.dump(prov, open(OUT / "d134_s0_run.json", "w"), indent=1)
    print(json.dumps({k: prov[k] for k in ["section2_reproduced", "profiles", "already_applied_tas_floor",
                                           "skeletons", "coverage_crosscheck"]}, indent=1)[:4000])
    print(pd.DataFrame(attr).to_string(index=False))
    print(f"done in {prov['runtime_s']}s")


if __name__ == "__main__":
    main()
