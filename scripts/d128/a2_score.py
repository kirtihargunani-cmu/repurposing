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
    return dict(idx=idx, genes=genes, M=M, col=col, land=land, F=F, amb=amb, U=U, P=P, arm=arm, path=path)


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


def score_one(inp, pert_id, model_id, unit_id, de_sign=None, want_nulls=0, cov=None):
    """M and R for one (drug, line, unit). Returns (Score_M, Score_R, per-gene rows, null summaries)."""
    row = inp["idx"].index[(inp["idx"].pert_id == pert_id) & (inp["idx"].ModelID == model_id)]
    if not len(row):
        return (SC.Score(np.nan, "abstained", "no_profile_for_line"),
                SC.Score(np.nan, "abstained", "no_profile_for_line"), [], {})
    zrow = np.asarray(inp["M"][int(row[0])], dtype=float)

    f = inp["F"][(inp["F"].unit_id == unit_id)]
    if not len(f):
        return (SC.Score(np.nan, "abstained", "no_footprint"), SC.Score(np.nan, "abstained", "no_footprint"), [], {})
    res = f[f.status == "resolved"]
    n_unres = int((f.status != "resolved").sum())
    gids = res.gene_id.tolist()
    z = np.array([zrow[inp["col"][g]] if g in inp["col"] else np.nan for g in gids])
    lm = np.array([g in inp["land"] for g in gids])
    u = inp["U"][inp["U"].unit_id == unit_id].iloc[0]
    action = 1 if u.orient == "raise" else -1

    m = SC.score_mechanism(gids, res.expected_sign_pre_action.to_numpy(), action, z, landmark_mask=lm,
                           n_unresolved=n_unres, action_flag=u.action_flag)
    # R: the same footprint genes, restricted to the cohort's significant genes (P8)
    if de_sign:
        rg = [g for g in f.gene_id.tolist() if g in de_sign]
        rz = np.array([zrow[inp["col"][g]] if g in inp["col"] else np.nan for g in rg])
        rlm = np.array([g in inp["land"] for g in rg])
        r = SC.score_rescue(rg, [de_sign[g] for g in rg], rz, landmark_mask=rlm,
                            n_unresolved=len(f) - len(rg))
    else:
        r = SC.Score(np.nan, "abstained", "no_de_contrast")

    grows = []
    for g, c in m.contributions.items():
        grows.append(dict(unit_id=unit_id, pert_id=pert_id, ModelID=model_id, gene_id=g, score="M",
                          expected_sign=float(res.expected_sign_pre_action.to_numpy()[gids.index(g)] * action),
                          z=float(zrow[inp["col"][g]]), f_z=float(SC.bounded(zrow[inp["col"][g]])),
                          contribution=float(c), feature_space="landmark" if g in inp["land"] else "best inferred"))

    nsum = {}
    if want_nulls and m.scored:
        e = res.expected_sign_pre_action.to_numpy() * action
        vals = [SC.score_mechanism(gids, NU.shuffled_edge_signs(e, i), 1, z, landmark_mask=lm).value
                for i in range(want_nulls)]
        nsum["sign_shuffle"] = NU.summarise(vals, m.value)
        vals = [SC.score_mechanism(gids, e, 1, NU.gene_label_permutation(z, lm, i), landmark_mask=lm).value
                for i in range(want_nulls)]
        nsum["gene_label"] = NU.summarise(vals, m.value)
        if cov is not None:
            cand = NU.degree_matched_units(unit_id, cov, int(m.n_genes),
                                           int(cov.loc[cov.unit_id == unit_id, "n_resolved_tfs"].iloc[0])
                                           if (cov.unit_id == unit_id).any() else 0,
                                           k=min(want_nulls, 20), exclude_genes=str(u.genes_str).split("|"))
            dvals = []
            for cu in cand.unit_id:
                cf = inp["F"][(inp["F"].unit_id == cu) & (inp["F"].status == "resolved")]
                cg = cf.gene_id.tolist()
                cz = np.array([zrow[inp["col"][g]] if g in inp["col"] else np.nan for g in cg])
                dvals.append(SC.score_mechanism(cg, cf.expected_sign_pre_action.to_numpy(), action, cz).value)
            nsum["degree_matched_target"] = NU.summarise(dvals, m.value)
    return m, r, grows, nsum


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--cohorts", default="")
    ap.add_argument("--compounds", default="")
    ap.add_argument("--lines", default="")
    ap.add_argument("--arm", default=C.PROFILE_ARM, choices=["specific", "raw"])
    ap.add_argument("--nulls", type=int, default=0)
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
    log(f"arm {a.arm}; scoring {len(E)} (cohort, pert_id, line, unit) rows over {E.cohort.nunique()} cohorts, "
        f"{E.pert_id.nunique()} compounds, {E.ModelID.nunique()} lines")

    de_cache = {}
    rows, grows, nrows = [], [], []
    for n, r in enumerate(E.itertuples(index=False), 1):
        if r.cohort not in de_cache:
            de_cache[r.cohort] = disease_genes(r.cohort, r.level)
        con, dsign = de_cache[r.cohort]
        m, rr, gr, ns = score_one(inp, r.pert_id, r.ModelID, r.unit_id, dsign, a.nulls, cov)
        u = inp["U"][inp["U"].unit_id == r.unit_id].iloc[0]
        rows.append(dict(cohort=r.cohort, level=r.level, pert_id=r.pert_id, cmap_name=r.cmap_name,
                         ModelID=r.ModelID, line_weight=r.line_weight, source=r.source_a0, unit_id=r.unit_id,
                         unit_genes=u.genes_str, unit_kind=r.unit_kind, action=u.orient, action_flag=u.action_flag,
                         arm=a.arm, contrast=con,
                         M=m.value, M_status=m.status, M_abstain_reason=m.abstain_reason, M_n_genes=m.n_genes,
                         M_n_landmark=m.n_landmark, M_n_inferred=m.n_inferred,
                         M_n_excluded_unresolved=m.n_excluded_unresolved, M_n_missing=m.n_missing_response,
                         M_hub_fraction=m.hub_fraction,
                         R=rr.value, R_status=rr.status, R_abstain_reason=rr.abstain_reason, R_n_genes=rr.n_genes,
                         R_n_landmark=rr.n_landmark, R_n_inferred=rr.n_inferred))
        grows.extend(gr)
        for k, v in ns.items():
            nrows.append(dict(cohort=r.cohort, pert_id=r.pert_id, ModelID=r.ModelID, unit_id=r.unit_id,
                              null=k, observed=m.value, **v))
        if n % 200 == 0 or n == len(E):
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
    ab.groupby(["M_abstain_reason", "R_abstain_reason"]).size().rename("n").reset_index() \
        .to_csv(out / "abstentions.tsv", sep="\t", index=False)
    defined = S.M[S.M_status == "scored"]
    assert defined.between(-1, 1).all(), "an M value fell outside [-1, 1]"
    log(f"M scored {int((S.M_status == 'scored').sum())} / {len(S)}; R scored {int((S.R_status == 'scored').sum())}")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (out / "score_run.json").write_text(json.dumps(dict(
        stage="D128-A2", git_sha=sha, frozen=C.FROZEN, arm=a.arm, profile_file=str(inp["path"].name),
        filters=dict(cohorts=a.cohorts, compounds=a.compounds, lines=a.lines), nulls=a.nulls,
        counts=dict(rows=len(S), M_scored=int((S.M_status == "scored").sum()),
                    R_scored=int((S.R_status == "scored").sum())),
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
