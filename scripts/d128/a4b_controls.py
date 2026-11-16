#!/usr/bin/env python
"""D128 A4 actions 2, 3 and 5: the remaining mechanism controls, the genotype negatives, and the failure analysis.

Action 1 (the mandatory MCF7 ESR1 control) is `a4_controls.py` and already ran; it is not repeated here. This script
is a second file rather than an extension of it so that action 1's committed output stays byte-reproducible.

  action 2  HER2, BRAF and EGFR compounds, in the lines whose genotype A0 confirmed. EGFR wave A4 is **not
            computable** (no HCC827 / PC9 / HCC4006 in the LINCS lines; NCIH1975 is T790M-resistant) and is recorded
            as a stated gap - no substitute line is quietly used in its place. The BRAF footprint's MYC
            hub-dominance is quantified: what fraction of the signed numerator comes from MYC-regulated genes.
  action 3  mutation / line-context negatives: the RAF-inhibitor paradox (RAF inhibitor in a BRAF-WT, RAS-active
            line) and EGFR inhibitors in EGFR-wild-type, KRAS-active lines, plus HER2 drugs in ERBB2-non-amplified
            lines. The question is whether a negative becomes a *confident* positive on a generic growth-arrest
            response.
  action 5  failures analysed by MECHANISM, not by rank: `mechanism_case_reports/`.

Outputs results/lincs/d128/a4/: control_results.tsv, control_summary.tsv, braf_hub_dominance.tsv,
braf_hub_genes.tsv, mechanism_case_reports/.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a2_score as A2  # noqa: E402
import config as C  # noqa: E402
import score as SC  # noqa: E402

A4 = C.OUT / "a4"
T0 = time.time()
LOG = []

# Wave definitions. `positive` and `negative` are genotype RULES evaluated against a0/source_line_context.tsv, never
# hand-picked line lists, so the membership of each set is reproducible and auditable.
WAVES = [
    dict(wave="A2-HER2", anchor="ERBB2",
         compounds=["lapatinib", "neratinib", "afatinib", "canertinib", "tucatinib"],
         positive=("ERBB2_amp_candidate", True), negative=("ERBB2_amp_candidate", False),
         note="A0: executable but THIN - 12 signatures on SKBR3. Do not over-read it."),
    dict(wave="A3-BRAF", anchor="BRAF",
         compounds=["vemurafenib", "dabrafenib", "encorafenib", "PLX-4720", "SB-590885"],
         positive=("BRAF_V600", True), negative=("KRAS_NRAS_hotspot", True),
         note="A0: the richest wave (104 signatures). Negative = RAF-inhibitor paradox context "
              "(BRAF-WT, RAS-active): a RAF inhibitor must not read as a confident positive there."),
    dict(wave="A4-EGFR", anchor="EGFR",
         compounds=["gefitinib", "erlotinib", "afatinib", "osimertinib", "dacomitinib"],
         positive=("EGFR_sensitizing", True), negative=("KRAS_NRAS_hotspot", True),
         note="A0: wave A4 is NOT COMPUTABLE as designed - no HCC827 / PC9 / HCC4006 among the LINCS lines. "
              "NCIH1975 carries a sensitizing mutation AND T790M (resistant), so it is reported as a "
              "resistance-genotype line, not as the missing sensitive line. STATED GAP."),
]
CONFIDENT = dict(min_pct=95.0)    # a 'confident positive' = M > 0 and BOTH permutation nulls at pct >= 95


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def units_for(inp, cmap_name, anchor):
    """Every target unit of this compound whose annotated gene set contains the wave's anchor gene."""
    out = []
    for u in inp["U"].itertuples(index=False):
        if u.cmap_name == cmap_name and anchor in str(u.genes_str).split("|"):
            out.append(u.unit_id)
    return out


def run_controls(inp, ctx, nulls, cov, gene_arms=("all_bing",)):
    rows = []
    prof = inp["idx"]
    for w in WAVES:
        for ctxname, (flag, want) in (("positive", w["positive"]), ("negative", w["negative"])):
            lines = set(ctx.ModelID[ctx[flag] == want])
            if ctxname == "negative" and w["anchor"] == "BRAF":
                lines &= set(ctx.ModelID[~ctx.BRAF_V600.astype(bool)])          # paradox needs BRAF-WT
            if ctxname == "negative" and w["anchor"] == "EGFR":
                lines &= set(ctx.ModelID[~ctx.EGFR_sensitizing.astype(bool)])   # EGFR-WT
            for cmp_name in w["compounds"]:
                uids = units_for(inp, cmp_name, w["anchor"])
                p = prof[(prof.cmap_name == cmp_name) & (prof.ModelID.isin(lines))]
                if not uids or not len(p):
                    rows.append(dict(wave=w["wave"], context=ctxname, anchor=w["anchor"], gene_arm="", cmap_name=cmp_name,
                                     pert_id="", ModelID="", unit_id=uids[0] if uids else "",
                                     M=np.nan, M_status="not_run",
                                     M_abstain_reason="no_target_unit" if not uids else "no_profile_in_context",
                                     note=w["note"]))
                    continue
                for uid in uids:
                  for ga in gene_arms:
                    for r in p.itertuples(index=False):
                        m, _, _, ns = A2.score_one(inp, r.pert_id, r.ModelID, uid, None, nulls, cov,
                                                   gene_contrib=False, gene_arm=ga)
                        cx = ctx[ctx.ModelID == r.ModelID].iloc[0]
                        d = dict(wave=w["wave"], context=ctxname, anchor=w["anchor"], gene_arm=ga, cmap_name=cmp_name,
                                 pert_id=r.pert_id, ModelID=r.ModelID, line=cx.StrippedCellLineName,
                                 lineage=cx.OncotreeLineage, n_sig=r.n_sig, unit_id=uid,
                                 unit_genes=inp["U_by_id"][uid].genes_str,
                                 action=inp["U_by_id"][uid].orient, action_flag=inp["U_by_id"][uid].action_flag,
                                 BRAF_V600=bool(cx.BRAF_V600), EGFR_sensitizing=bool(cx.EGFR_sensitizing),
                                 EGFR_T790M=bool(cx.EGFR_T790M), KRAS_NRAS=bool(cx.KRAS_NRAS_hotspot),
                                 ERBB2_amp=bool(cx.ERBB2_amp_candidate),
                                 M=m.value, M_status=m.status, M_abstain_reason=m.abstain_reason,
                                 M_n_genes=m.n_genes, M_n_landmark=m.n_landmark, M_n_inferred=m.n_inferred,
                                 M_n_excluded_unresolved=m.n_excluded_unresolved, M_n_missing=m.n_missing_response,
                                 note=w["note"])
                        for k, v in ns.items():
                            d[f"{k}_pct"] = round(v["percentile"], 1) if np.isfinite(v["percentile"]) else np.nan
                            d[f"{k}_z"] = round(v["z"], 2) if np.isfinite(v["z"]) else np.nan
                            d[f"{k}_n"] = v["n_null"]
                        rows.append(d)
    D = pd.DataFrame(rows)
    sc = D.M_status == "scored"
    D["confident_positive"] = (sc & (D.M > 0)
                               & (D.get("sign_shuffle_pct", pd.Series(np.nan, index=D.index)) >= CONFIDENT["min_pct"])
                               & (D.get("gene_label_pct", pd.Series(np.nan, index=D.index)) >= CONFIDENT["min_pct"]))
    return D


def braf_hub(inp, ctx, out_dir, gene_arm="all_bing"):
    """The check the user asked for: is BRAF's M measuring BRAF mechanism or generic MYC-driven proliferation?"""
    E = pd.read_parquet(C.OUT / "a1/regulon_edges.parquet")
    E = E[E.in_lincs_gene_space]
    reg = {t: set(d.gene_id) for t, d in E.groupby("regulator")}
    P = pd.read_parquet(C.OUT / "a1/target_tf_paths.parquet")
    anch = P[(P.target_gene == "BRAF") & (P.network == "omnipath") & (P.hops <= C.MAX_HOPS)
             & (P.path_status == "resolved") & P.tf_has_regulon]
    tfs = sorted(set(anch.tf))
    rows, grows = [], []
    for uid, fp in inp["FP"].items():
        if not uid.startswith(tuple(c + "::" for c in
                                    ["vemurafenib", "dabrafenib", "PLX-4720", "SB-590885", "encorafenib"])) \
                and "BRAF" not in str(inp["U_by_id"][uid].genes_str).split("|"):
            continue
        if "BRAF" not in str(inp["U_by_id"][uid].genes_str).split("|"):
            continue
        gk, ck, ek = (("lm_genes", "lm_cols", "lm_expected") if gene_arm == "landmark_only"
                      else ("genes", "cols", "expected"))
        g = np.array(fp[gk])
        if not len(g):
            continue
        memb = {t: np.array([x in reg.get(t, ()) for x in g]) for t in tfs}
        share_genes = {t: float(memb[t].mean()) for t in tfs}
        top = sorted(share_genes, key=share_genes.get, reverse=True)[:5]
        cm = inp["U_by_id"][uid].cmap_name
        for r in inp["idx"][(inp["idx"].cmap_name == cm)].itertuples(index=False):
            cx = ctx[ctx.ModelID == r.ModelID]
            if not len(cx):
                continue
            m, _, _, _ = A2.score_one(inp, r.pert_id, r.ModelID, uid, None, 0, None, gene_contrib=False,
                                      gene_arm=gene_arm)
            if not m.scored:
                continue
            e = fp[ek] * (1 if inp["U_by_id"][uid].orient == "raise" else -1)
            contrib = e * SC.bounded(np.asarray(inp["M"][inp["prow"][(r.pert_id, r.ModelID)]], float)[fp[ck]])
            tot_abs, tot_sig = float(np.abs(contrib).sum()), float(contrib.sum())
            d = dict(unit_id=uid, gene_arm=gene_arm, cmap_name=cm, pert_id=r.pert_id, ModelID=r.ModelID,
                     line=cx.StrippedCellLineName.iloc[0], BRAF_V600=bool(cx.BRAF_V600.iloc[0]),
                     n_genes=len(g), n_anchored_tfs=len(tfs), M=m.value)
            for t in ["MYC"] + [x for x in top if x != "MYC"][:3]:
                k = memb.get(t)
                if k is None or not k.any():
                    continue
                d[f"share_genes_{t}"] = round(float(k.mean()), 4)
                d[f"share_abs_numerator_{t}"] = round(float(np.abs(contrib[k]).sum() / tot_abs), 4) if tot_abs else np.nan
                d[f"share_signed_numerator_{t}"] = round(float(contrib[k].sum() / tot_sig), 4) if tot_sig else np.nan
                d[f"M_{t}_only"] = round(float(contrib[k].mean()), 5)
                d[f"M_without_{t}"] = round(float(contrib[~k].mean()), 5) if (~k).any() else np.nan
                d[f"n_genes_{t}"] = int(k.sum())
            # genes carried by MYC alone (no other anchored, resolved TF of BRAF regulates them)
            if "MYC" in memb:
                others = np.zeros(len(g), bool)
                for t in tfs:
                    if t != "MYC":
                        others |= memb[t]
                only = memb["MYC"] & ~others
                d["n_genes_MYC_only_source"] = int(only.sum())
                d["share_signed_numerator_MYC_only_source"] = round(float(contrib[only].sum() / tot_sig), 4) \
                    if tot_sig and only.any() else np.nan
                d["M_excluding_MYC_only_source"] = round(float(contrib[~only].mean()), 5) if (~only).any() else np.nan
            rows.append(d)
        for t in tfs:
            grows.append(dict(unit_id=uid, gene_arm=gene_arm, tf=t, n_regulon_genes_in_lincs=len(reg.get(t, ())),
                              n_in_footprint=int(memb[t].sum()), share_of_footprint=round(share_genes[t], 4)))
    H = pd.DataFrame(rows)
    pd.DataFrame(grows).drop_duplicates().sort_values("share_of_footprint", ascending=False) \
        .to_csv(out_dir / f"braf_hub_genes_{gene_arm}.tsv", sep="\t", index=False)
    return H


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--nulls", type=int, default=200)
    ap.add_argument("--gene-arms", default="all_bing,landmark_only",
                    help="the D128-A3 landmark amendment failed its criterion, so the controls are reported on "
                         "both gene spaces rather than on the frozen primary alone")
    a = ap.parse_args()
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen")
    A4.mkdir(parents=True, exist_ok=True)

    inp = A2.load_inputs("specific")
    cov = pd.read_csv(C.OUT / "a0/coverage.tsv", sep="\t")
    ctx = pd.read_csv(C.OUT / "a0/source_line_context.tsv", sep="\t", low_memory=False)
    for c in ["BRAF_V600", "EGFR_sensitizing", "EGFR_T790M", "KRAS_NRAS_hotspot", "ERBB2_amp_candidate",
              "RB1_null_candidate"]:
        ctx[c] = ctx[c].astype(bool)

    D = run_controls(inp, ctx, a.nulls, cov, gene_arms=tuple(x for x in a.gene_arms.split(",") if x))
    D.to_csv(A4 / "control_results.tsv", sep="\t", index=False)
    sc = D[D.M_status == "scored"]
    SU = sc.groupby(["wave", "context", "gene_arm"]).agg(
        n_rows=("M", "size"), n_compounds=("cmap_name", "nunique"), n_lines=("ModelID", "nunique"),
        n_signatures=("n_sig", "sum"), M_median=("M", "median"), M_min=("M", "min"), M_max=("M", "max"),
        share_M_positive=("M", lambda x: round(float((x > 0).mean()), 3)),
        median_sign_shuffle_pct=("sign_shuffle_pct", "median"),
        median_gene_label_pct=("gene_label_pct", "median"),
        n_confident_positive=("confident_positive", "sum"),
        share_confident_positive=("confident_positive", lambda x: round(float(x.mean()), 3)),
        median_n_genes=("M_n_genes", "median"), median_n_landmark=("M_n_landmark", "median")).reset_index()
    SU.to_csv(A4 / "control_summary.tsv", sep="\t", index=False)
    log("control summary:\n" + SU.to_string())
    log("not run / abstained rows:\n" + D[D.M_status != "scored"][
        ["wave", "context", "cmap_name", "ModelID", "M_status", "M_abstain_reason"]].to_string())

    H = pd.concat([braf_hub(inp, ctx, A4, ga) for ga in (x for x in a.gene_arms.split(",") if x)],
                  ignore_index=True)
    H.to_csv(A4 / "braf_hub_dominance.tsv", sep="\t", index=False)
    cols = [c for c in H.columns if "MYC" in c or c in ("cmap_name", "line", "BRAF_V600", "M", "n_genes",
                                                        "n_anchored_tfs", "pert_id", "gene_arm")]
    log("BRAF MYC hub dominance (per profile):\n" + H[cols].to_string())

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A4 / "a4b_run.json").write_text(json.dumps(dict(
        stage="D128-A4 actions 2,3,5", git_sha=sha, frozen=C.FROZEN, nulls=a.nulls,
        confident_positive_rule=CONFIDENT, waves=[{k: v for k, v in w.items()} for w in WAVES],
        summary=SU.to_dict("records"), log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
