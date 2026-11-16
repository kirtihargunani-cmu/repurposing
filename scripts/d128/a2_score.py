#!/usr/bin/env python
"""D128 A2: score M and R per (cohort, pert_id, source line, target unit) from the A1 footprints (plan Section 4 A2).

Loads the frozen inputs once and scores whatever subset the flags select, so the same code path serves the smoke
run, the ESR1 control (A4) and later the full population. Production is never written; everything lands under
results/lincs/d128/.

  --cohorts / --compounds / --lines   restrict the run (comma-separated); default: everything eligible
  --arm specific|raw                  P5 primary is `specific`; `raw` is the mandatory second arm of the control
  --nulls N                           degree-matched, sign-shuffled and gene-label nulls per scored row
  --out NAME                          subdirectory under results/lincs/d128/

Outputs: scores_by_line.parquet, gene_contributions.parquet, abstentions.tsv, nulls.parquet, score_run.json
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import d128io as IO  # noqa: E402
import nulls as NU  # noqa: E402
import score as SC  # noqa: E402

T0 = time.time()
LOG = []


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def load_inputs(arm="specific"):
    """Everything A2 needs, loaded once. Read-only."""
    idx = IO.profile_index()
    genes = IO.profile_genes()
    ann = pd.read_csv(C.DRUG_ANNOT, sep="\t", usecols=["pert_id", "cmap_name"]).drop_duplicates("pert_id")
    idx = idx.merge(ann, on="pert_id", how="left").reset_index(drop=True)
    path = C.PROFILES_SPECIFIC if arm == "specific" else C.PROFILES_RAW
    M = IO.profiles(path, len(idx), len(genes))
    col = {g: i for i, g in enumerate(genes.gene_id)}
    land = set(genes.gene_id[genes.feature_space == "landmark"])
    F = pd.read_parquet(C.OUT / "a1/target_gene_footprints.parquet")
    amb = pd.read_csv(C.OUT / "a1/path_ambiguity.tsv", sep="\t")
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str}).fillna({"action_flag": ""})
    U["genes_str"] = U.genes
    U["unit_id"] = U.cmap_name + "::" + U.tid.fillna("annot:" + U.genes_str)
    dup = U.unit_id.duplicated(keep=False)
    if dup.any():
        U.loc[dup, "unit_id"] = U.unit_id[dup] + "#" + U[dup].groupby("unit_id").cumcount().add(1).astype(str)
    # which TF carries each gene, for the hub diagnostic
    P = pd.read_parquet(C.OUT / "a1/target_tf_paths.parquet")
    inp = dict(idx=idx, genes=genes, M=M, col=col, land=land, F=F, amb=amb, U=U, P=P, arm=arm, path=path)

    # --- D128-A3 indexes: the same values, looked up instead of re-scanned ------------------------------------------
    # A3 scores ~114k rows. The first version filtered the 368k-row footprint frame and the unit table once per row,
    # which is O(rows x frame) and dominates everything else. Nothing below changes a score: it is the same gene list,
    # the same expected signs and the same unit record, prepared once.
    FP = {}
    for uid, d in F.groupby("unit_id", sort=False):
        res = d[d.status == "resolved"]
        g = res.gene_id.to_numpy()
        FP[uid] = dict(
            genes=g.tolist(), cols=np.array([col[x] for x in g], dtype=np.int64),
            expected=res.expected_sign_pre_action.to_numpy().astype(float),
            lm=np.array([x in land for x in g], dtype=bool), n_unresolved=int((d.status != "resolved").sum()),
            all_genes=d.gene_id.to_numpy().tolist(),
            # D128-A3: identity of the resolved footprint, so a degree-matched null does not draw the same target
            # unit several times under different drug names (the 497-unit pool is full of such repeats).
            sig=hash(tuple(sorted(zip(g.tolist(), res.expected_sign_pre_action.to_numpy().tolist())))))
        # D128-A3 landmark sensitivity (post-freeze, reported only): the same footprint restricted to the 978
        # physically measured landmark genes. P4 (all BING) stays the primary arm.
        k = FP[uid]["lm"]
        FP[uid].update(lm_genes=[x for x, t in zip(FP[uid]["genes"], k) if t], lm_cols=FP[uid]["cols"][k],
                       lm_expected=FP[uid]["expected"][k], lm_mask=np.ones(int(k.sum()), bool),
                       all_genes_landmark=[x for x in d.gene_id.to_numpy().tolist() if x in land])
    inp["FP"] = FP
    inp["U_by_id"] = {r.unit_id: r for r in U.itertuples(index=False)}
    inp["prow"] = {(p_, m_): i for i, (p_, m_) in enumerate(zip(idx.pert_id, idx.ModelID))}
    inp["cand_cache"] = {}
    return inp


def disease_genes(cohort, level):
    """P1/P2/P8: signed effect and the significant gene set for one cohort."""
    proc, summ = (C.DE_GROUPS if level == "group_id" else C.DE_STRATA)
    S = pd.read_csv(summ, sep="\t")
    key = "group_id" if "group_id" in S.columns else S.columns[0]
    ok = S[(S[key] == cohort) & (S.status == "ok")].contrast.tolist()
    con = next((x for x in C.CONTRAST_ORDER if x in ok), None)
    if con is None:
        return None, None
    d = IO.de_table(proc, cohort, con)
    if d is None:
        return None, None
    sig = d[d[C.DE_P_COL] < C.R_P_MAX]
    return con, dict(zip(sig.gene_id, np.sign(sig[C.DE_SIGN_COL])))


def score_one(inp, pert_id, model_id, unit_id, de_sign=None, want_nulls=0, cov=None, gene_contrib=True,
              gene_arm="all_bing"):
    """M and R for one (drug, line, unit). Returns (Score_M, Score_R, per-gene rows, null summaries).

    `gene_arm` is the D128-A3 landmark sensitivity switch: `all_bing` is the frozen primary (P4), `landmark_only`
    runs the identical scoring on the subset of footprint genes that L1000 physically measures. The floors apply to
    whichever gene set the arm uses, so a unit can clear in one arm and abstain in the other.
    """
    i_row = inp["prow"].get((pert_id, model_id))
    if i_row is None:
        return (SC.Score(np.nan, "abstained", "no_profile_for_line"),
                SC.Score(np.nan, "abstained", "no_profile_for_line"), [], {})
    zrow = np.asarray(inp["M"][i_row], dtype=float)

    fp = inp["FP"].get(unit_id)
    if fp is None or not len(fp["genes"]):
        return (SC.Score(np.nan, "abstained", "no_footprint"), SC.Score(np.nan, "abstained", "no_footprint"), [], {})
    if gene_arm == "landmark_only":
        gids, z, lm, expected = fp["lm_genes"], zrow[fp["lm_cols"]], fp["lm_mask"], fp["lm_expected"]
        r_universe = fp["all_genes_landmark"]
    else:
        gids, z, lm, expected = fp["genes"], zrow[fp["cols"]], fp["lm"], fp["expected"]
        r_universe = fp["all_genes"]
    if not len(gids):
        return (SC.Score(np.nan, "abstained", "no_landmark_gene_in_footprint"),
                SC.Score(np.nan, "abstained", "no_landmark_gene_in_footprint"), [], {})
    u = inp["U_by_id"][unit_id]
    action = 1 if u.orient == "raise" else -1

    m = SC.score_mechanism(gids, expected, action, z, landmark_mask=lm,
                           n_unresolved=fp["n_unresolved"], action_flag=u.action_flag, tf_of_gene=None)
    # R: the same footprint genes (resolved or not: sign ambiguity blocks M, not direction), restricted to the
    # cohort's significant genes (P8)
    if de_sign:
        rg = [g for g in r_universe if g in de_sign]
        rz = np.array([zrow[inp["col"][g]] if g in inp["col"] else np.nan for g in rg])
        rlm = np.array([g in inp["land"] for g in rg])
        r = SC.score_rescue(rg, [de_sign[g] for g in rg], rz, landmark_mask=rlm,
                            n_unresolved=len(r_universe) - len(rg))
    else:
        r = SC.Score(np.nan, "abstained", "no_de_contrast")

    grows = []
    if gene_contrib:
        e_signed = expected * action
        pos = {g: k for k, g in enumerate(gids)}
        for g, c in m.contributions.items():
            k = pos[g]
            grows.append(dict(unit_id=unit_id, pert_id=pert_id, ModelID=model_id, gene_id=g, score="M",
                              expected_sign=float(e_signed[k]), z=float(z[k]), f_z=float(SC.bounded(z[k])),
                              contribution=float(c), feature_space="landmark" if lm[k] else "best inferred"))

    nsum = {}
    if want_nulls and m.scored:
        e = expected * action
        if NU.exact_shortcut_ok(e, z):                       # the fast path is exact here (see nulls.py)
            fz = SC.bounded(z)
            nsum["sign_shuffle"] = NU.summarise(NU.m_null_sign_shuffle(e, fz, want_nulls), m.value)
            nsum["gene_label"] = NU.summarise(NU.m_null_gene_label(e, fz, lm, want_nulls), m.value)
        else:
            nsum["sign_shuffle"] = NU.summarise(
                [SC.score_mechanism(gids, NU.shuffled_edge_signs(e, i), 1, z, landmark_mask=lm).value
                 for i in range(want_nulls)], m.value)
            nsum["gene_label"] = NU.summarise(
                [SC.score_mechanism(gids, e, 1, NU.gene_label_permutation(z, lm, i), landmark_mask=lm).value
                 for i in range(want_nulls)], m.value)
        if cov is not None:
            cand = inp["cand_cache"].get(unit_id)
            if cand is None:
                ntf = int(cov.loc[cov.unit_id == unit_id, "n_resolved_tfs"].iloc[0]) if (cov.unit_id == unit_id).any() else 0
                raw = list(NU.degree_matched_units(unit_id, cov, int(m.n_genes), ntf, k=min(want_nulls, 20),
                                                   exclude_genes=str(u.genes_str).split("|")).unit_id)
                seen, cand = set(), []                      # D128-A3: one draw per distinct footprint
                for cu in raw:
                    cf = inp["FP"].get(cu)
                    if cf is None or not len(cf["genes"]) or cf["sig"] in seen:
                        continue
                    seen.add(cf["sig"]); cand.append(cu)
                inp["cand_cache"][unit_id] = cand
            gk, ck, ek = (("lm_genes", "lm_cols", "lm_expected") if gene_arm == "landmark_only"
                          else ("genes", "cols", "expected"))
            dvals = []
            for cu in cand:
                cf = inp["FP"].get(cu)
                if cf is None or not len(cf[gk]):
                    continue
                cz = zrow[cf[ck]]
                ce = cf[ek] * action
                if NU.exact_shortcut_ok(ce, cz):
                    dvals.append(float((ce * SC.bounded(cz)).mean()))
                else:
                    dvals.append(SC.score_mechanism(cf[gk], cf[ek], action, cz).value)
            d_summ = NU.summarise(dvals, m.value)
            # D128-A3 reporting rule: fewer than 3 distinct alternative footprints is a degenerate null, not a number
            d_summ["degenerate"] = bool(d_summ["n_null"] < 3)
            if d_summ["degenerate"]:
                d_summ["z"] = np.nan
            nsum["degree_matched_target"] = d_summ
    return m, r, grows, nsum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--cohorts", default="")
    ap.add_argument("--compounds", default="")
    ap.add_argument("--lines", default="")
    ap.add_argument("--arm", default=C.PROFILE_ARM, choices=["specific", "raw"])
    ap.add_argument("--nulls", type=int, default=0)
    ap.add_argument("--gene-contrib", default="none", choices=["none", "all", "sample"],
                    help="write gene_contributions.parquet: none, all (small runs only), or a seeded sample")
    ap.add_argument("--gene-contrib-n", type=int, default=2000)
    ap.add_argument("--status", default="", help="restrict to eligibility rows with these a0 status values")
    ap.add_argument("--gene-arms", default="all_bing",
                    help="comma list of gene-space arms: all_bing (P4 primary) and/or landmark_only (D128-A3 "
                         "landmark sensitivity amendment, reported only)")
    ap.add_argument("--out", default="a2")
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen: A2 must not run (config.FROZEN is False)")
    out = C.OUT / a.out
    out.mkdir(parents=True, exist_ok=True)

    inp = load_inputs(a.arm)
    cov = pd.read_csv(C.OUT / "a0/coverage.tsv", sep="\t")
    E = pd.read_csv(C.OUT / "a0/eligibility.tsv.gz", sep="\t", low_memory=False)
    E = E[(E.unit_id != "") & E.unit_id.notna()]
    if a.cohorts:
        E = E[E.cohort.isin(a.cohorts.split(","))]
    if a.compounds:
        E = E[E.cmap_name.isin(a.compounds.split(","))]
    if a.lines:
        E = E[E.ModelID.isin(a.lines.split(","))]
    E = E[E.ModelID.notna() & (E.ModelID != "")]
    if a.status:
        E = E[E.status.isin(a.status.split(","))]
    E = E.sort_values(["pert_id", "ModelID"], kind="stable").reset_index(drop=True)
    gc_rows = set()
    if a.gene_contrib == "all":
        gc_rows = set(range(len(E)))
    elif a.gene_contrib == "sample":
        n_gc = min(a.gene_contrib_n, len(E))
        gc_rows = set(np.random.default_rng(C.SEED).choice(len(E), n_gc, replace=False).tolist())
    gene_arms = [x for x in a.gene_arms.split(",") if x]
    assert set(gene_arms) <= {"all_bing", "landmark_only"}, f"unknown gene arm in {gene_arms}"
    log(f"profile arm {a.arm}; gene arms {gene_arms}; scoring {len(E)} (cohort, pert_id, line, unit) rows over "
        f"{E.cohort.nunique()} cohorts, {E.pert_id.nunique()} compounds, {E.ModelID.nunique()} lines")

    de_cache = {}
    rows, grows, nrows = [], [], []
    for n, r in enumerate(E.itertuples(index=False), 1):
        if r.cohort not in de_cache:
            de_cache[r.cohort] = disease_genes(r.cohort, r.level)
        con, dsign = de_cache[r.cohort]
        u = inp["U_by_id"][r.unit_id]
        for ga in gene_arms:
            m, rr, gr, ns = score_one(inp, r.pert_id, r.ModelID, r.unit_id, dsign, a.nulls, cov,
                                      gene_contrib=((n - 1) in gc_rows and ga == "all_bing"), gene_arm=ga)
            rows.append(dict(cohort=r.cohort, level=r.level, pert_id=r.pert_id, cmap_name=r.cmap_name,
                             ModelID=r.ModelID, line_weight=r.line_weight, source=r.source_a0, unit_id=r.unit_id,
                             unit_genes=u.genes_str, unit_kind=r.unit_kind, action=u.orient,
                             action_flag=u.action_flag, arm=a.arm, gene_arm=ga, contrast=con,
                             M=m.value, M_status=m.status, M_abstain_reason=m.abstain_reason, M_n_genes=m.n_genes,
                             M_n_landmark=m.n_landmark, M_n_inferred=m.n_inferred,
                             M_n_excluded_unresolved=m.n_excluded_unresolved, M_n_missing=m.n_missing_response,
                             M_hub_fraction=m.hub_fraction,
                             R=rr.value, R_status=rr.status, R_abstain_reason=rr.abstain_reason,
                             R_n_genes=rr.n_genes, R_n_landmark=rr.n_landmark, R_n_inferred=rr.n_inferred))
            grows.extend(gr)
            for k, v in ns.items():
                nrows.append(dict(cohort=r.cohort, pert_id=r.pert_id, ModelID=r.ModelID, unit_id=r.unit_id,
                                  gene_arm=ga, null=k, observed=m.value, **v))
        if n % 5000 == 0 or n == len(E):
            log(f"scored {n}/{len(E)}")

    S = pd.DataFrame(rows)
    # E17 hub diagnostic: A1 already measured the largest TF's share of each unit's resolved footprint
    S = S.merge(inp["amb"][["unit_id", "hub_tf", "hub_fraction"]].rename(
        columns={"hub_fraction": "hub_fraction_structural"}), on="unit_id", how="left")
    S.to_parquet(out / "scores_by_line.parquet", index=False)
    if grows:
        pd.DataFrame(grows).to_parquet(out / "gene_contributions.parquet", index=False)
    if nrows:
        pd.DataFrame(nrows).to_parquet(out / "nulls.parquet", index=False)
    ab = S[(S.M_status == "abstained") | (S.R_status == "abstained")]
    ab.groupby(["gene_arm", "M_abstain_reason", "R_abstain_reason"]).size().rename("n").reset_index() \
        .to_csv(out / "abstentions.tsv", sep="\t", index=False)
    defined = S.M[S.M_status == "scored"]
    assert defined.between(-1, 1).all(), "an M value fell outside [-1, 1]"
    log(f"M scored {int((S.M_status == 'scored').sum())} / {len(S)}; R scored {int((S.R_status == 'scored').sum())}")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (out / "score_run.json").write_text(json.dumps(dict(
        stage="D128-A2", git_sha=sha, frozen=C.FROZEN, arm=a.arm, gene_arms=gene_arms,
        profile_file=str(inp["path"].name),
        filters=dict(cohorts=a.cohorts, compounds=a.compounds, lines=a.lines, status=a.status), nulls=a.nulls,
        counts=dict(rows=len(S), M_scored=int((S.M_status == "scored").sum()),
                    R_scored=int((S.R_status == "scored").sum()),
                    by_gene_arm=S.groupby("gene_arm").apply(
                        lambda d: dict(rows=len(d), M_scored=int((d.M_status == "scored").sum()),
                                       R_scored=int((d.R_status == "scored").sum()))).to_dict()),
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
