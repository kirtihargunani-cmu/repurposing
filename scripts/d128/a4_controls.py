#!/usr/bin/env python
"""D128 A4 action 1: the mandatory MCF7 anti-estrogen control (plan Section 4 A4; D0.2; run before any other control).

Gated on M alone (D0.1: M carries no disease term; breast R is reported against a prediction recorded in advance).
Fulvestrant is the primary positive set, each of its pert_ids reported separately and as a mean. The SERM set is
scored and reported SEPARATELY, never pooled, and only in breast-lineage lines (D0.2). Both profile arms are run:
a pass on raw alone is "generic removal removes the target signal", not a pass (P5).

Outputs results/lincs/d128/a4/:
  esr1_control.tsv           one row per (set, compound, pert_id, arm) with M, coverage and nulls
  esr1_control_summary.tsv   per set and arm, with the gate verdict
  esr1_breast_R.tsv          breast R beside the D0.1 prediction, per cohort and arm
  a4_run.json
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
import nulls as NU  # noqa: E402

A4 = C.OUT / "a4"
T0 = time.time()
LOG = []
MCF7 = "ACH-000019"
BREAST_COHORTS = ["BRCA-LumA", "BRCA-LumB", "BRCA-Her2", "BRCA-Basal"]


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--nulls", type=int, default=200)
    a = ap.parse_args()
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen: the control must not run")
    A4.mkdir(parents=True, exist_ok=True)

    cov = pd.read_csv(C.OUT / "a0/coverage.tsv", sep="\t")
    lines = pd.read_csv(C.OUT / "a0/source_line_context.tsv", sep="\t")
    breast_lines = set(lines.ModelID[lines.OncotreeLineage == C.SERM_LINEAGE])
    rows, brows = [], []

    for arm in ("specific", "raw"):
        inp = A2.load_inputs(arm)
        idx = inp["idx"]
        # the ESR1 unit of each compound; tamoxifen's ESRRG unit is excluded by D0.2
        esr1_units = {}
        for u in inp["U"].itertuples(index=False):
            gs = str(u.genes_str).split("|")
            if "ESR1" in gs and not set(gs) & set(C.ESR1_EXCLUDED_UNITS):
                esr1_units.setdefault(u.cmap_name, []).append(u.unit_id)

        for setname, compounds in (("primary_fulvestrant", C.ESR1_PRIMARY), ("serm_set", C.ESR1_SERM_SET)):
            for cmp_name in compounds:
                for unit_id in esr1_units.get(cmp_name, []):
                    prof = idx[(idx.cmap_name == cmp_name) & (idx.ModelID == MCF7)]
                    if not len(prof):
                        rows.append(dict(set=setname, cmap_name=cmp_name, pert_id="", arm=arm, unit_id=unit_id,
                                         M=np.nan, M_status="abstained", M_abstain_reason="no_MCF7_profile"))
                        continue
                    for p in prof.itertuples(index=False):
                        m, r, _, ns = A2.score_one(inp, p.pert_id, MCF7, unit_id, None, a.nulls, cov)
                        u = inp["U"][inp["U"].unit_id == unit_id].iloc[0]
                        # D0.2: SERMs count as ESR1 antagonists in breast-lineage lines only; MCF7 is one
                        serm_rule = setname == "serm_set" and MCF7 in breast_lines
                        if m.status == "abstained" and "ambiguous_action" in m.abstain_reason and serm_rule:
                            m, r, _, ns = A2.score_one(
                                inp, p.pert_id, MCF7, unit_id, None, a.nulls, cov)
                            # re-score with the tissue-context antagonist rule applied
                            import score as SC
                            f = inp["F"][(inp["F"].unit_id == unit_id) & (inp["F"].status == "resolved")]
                            row_i = int(idx.index[(idx.pert_id == p.pert_id) & (idx.ModelID == MCF7)][0])
                            z = np.array([np.asarray(inp["M"][row_i], float)[inp["col"][g]]
                                          if g in inp["col"] else np.nan for g in f.gene_id])
                            lm = np.array([g in inp["land"] for g in f.gene_id])
                            m = SC.score_mechanism(f.gene_id.tolist(), f.expected_sign_pre_action.to_numpy(), -1, z,
                                                   landmark_mask=lm)
                            e = f.expected_sign_pre_action.to_numpy() * -1
                            ns = {"sign_shuffle": NU.summarise(
                                [SC.score_mechanism(f.gene_id.tolist(), NU.shuffled_edge_signs(e, i), 1, z).value
                                 for i in range(a.nulls)], m.value),
                                "gene_label": NU.summarise(
                                [SC.score_mechanism(f.gene_id.tolist(), e, 1, NU.gene_label_permutation(z, lm, i)).value
                                 for i in range(a.nulls)], m.value)}
                        d = dict(set=setname, cmap_name=cmp_name, pert_id=p.pert_id, n_sig=p.n_sig, arm=arm,
                                 unit_id=unit_id, unit_genes=u.genes_str, action=u.orient,
                                 action_flag=u.action_flag, serm_rule_applied=bool(serm_rule and u.action_flag),
                                 M=m.value, M_status=m.status, M_abstain_reason=m.abstain_reason,
                                 M_n_genes=m.n_genes, M_n_landmark=m.n_landmark, M_n_inferred=m.n_inferred,
                                 M_n_excluded_unresolved=m.n_excluded_unresolved, M_n_missing=m.n_missing_response)
                        for k, v in ns.items():
                            d[f"{k}_null_median"] = round(v["null_median"], 5)
                            d[f"{k}_pct"] = round(v["percentile"], 1)
                            d[f"{k}_z"] = round(v["z"], 2) if np.isfinite(v["z"]) else np.nan
                        rows.append(d)
                        log(f"{arm} {setname} {cmp_name} {p.pert_id}: M={m.value if m.scored else 'NA'} "
                            f"({m.status}{':' + m.abstain_reason if m.abstain_reason else ''})")

        # breast R, against the D0.1 prediction (reported, never a gate)
        for cohort in BREAST_COHORTS:
            con, dsign = A2.disease_genes(cohort, "group_id")
            for unit_id in esr1_units.get("fulvestrant", []):
                for p in idx[(idx.cmap_name == "fulvestrant") & (idx.ModelID == MCF7)].itertuples(index=False):
                    m, r, _, _ = A2.score_one(inp, p.pert_id, MCF7, unit_id, dsign, 0, None)
                    brows.append(dict(cohort=cohort, contrast=con, arm=arm, pert_id=p.pert_id, network=C.NETWORK,
                                      R=r.value, R_status=r.status, R_n_genes=r.n_genes,
                                      R_n_landmark=r.n_landmark, M=m.value,
                                      d01_prediction="uninformative on CollecTRI (inverted predicted on ARACNe)"))

    D = pd.DataFrame(rows)
    D.to_csv(A4 / "esr1_control.tsv", sep="\t", index=False)
    B = pd.DataFrame(brows)
    B.to_csv(A4 / "esr1_breast_R.tsv", sep="\t", index=False)

    sc = D[D.M_status == "scored"]
    summ = []
    for (setname, arm), g in sc.groupby(["set", "arm"]):
        summ.append(dict(set=setname, arm=arm, n_pert_ids=len(g), n_compounds=g.cmap_name.nunique(),
                         M_mean=round(g.M.mean(), 4), M_min=round(g.M.min(), 4), M_max=round(g.M.max(), 4),
                         all_positive=bool((g.M > 0).all()),
                         median_sign_shuffle_pct=round(g.sign_shuffle_pct.median(), 1) if "sign_shuffle_pct" in g else np.nan,
                         median_gene_label_pct=round(g.gene_label_pct.median(), 1) if "gene_label_pct" in g else np.nan,
                         n_genes=int(g.M_n_genes.max()), n_landmark=int(g.M_n_landmark.max())))
    SU = pd.DataFrame(summ)
    SU.to_csv(A4 / "esr1_control_summary.tsv", sep="\t", index=False)

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A4 / "a4_run.json").write_text(json.dumps(dict(
        stage="D128-A4-esr1", git_sha=sha, frozen=C.FROZEN, nulls=a.nulls, line=MCF7,
        primary=list(C.ESR1_PRIMARY), serm_set=list(C.ESR1_SERM_SET), excluded_units=list(C.ESR1_EXCLUDED_UNITS),
        gap="no agonist sign-flip control: estradiol has 0 MCF7 signatures",
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
