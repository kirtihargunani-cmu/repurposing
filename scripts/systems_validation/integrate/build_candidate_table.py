#!/usr/bin/env python
"""Side-by-side candidate table: the LAYERED SCORE and the MECHANISM avenues, per candidate pair.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

Purpose (user, 2026-11-18): "output the results from the layered scoring function and then the work
we are completing now calculating F, O, S, the networks and traces ... to have a better
understanding of what limitations each presents."

This does NOT merge the avenues into a score. It places them beside each other on the SPEC section 0
entry set - existing drug, different established indication, matched source line - and records for
every pair whether each avenue can say anything at all, and why not when it cannot. The
not-evaluable reasons are the point of the table.

  LAYERED  production lead ranking (layered_pct and its D/A/O/T/r components, D122)
  REV      transcriptional reversal, ARACNe-only arm at minsize 25 (rev_pct, spec_mech_pct)
  F        causal fault isolation: is the target in the frozen signed PKN with outgoing edges
  O        observability: D134 target recovery / D135 CRISPR pairing
  S        state transitions: >=3 timepoints for this compound

Output: results/systems_validation/integrated/{candidate_table.tsv, coverage_summary.tsv}
"""
import os, sys
from pathlib import Path
import pandas as pd, numpy as np

ROOT = Path("/home/kirti/projects/repurposing")
RES = ROOT / "results/lincs"
OUT = ROOT / "results/systems_validation/integrated"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "4"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lab = pd.read_csv(RES / "known_use/labels.tsv", sep="\t")
    lab["known_use"] = lab.lab_hub.fillna(False) | lab.lab_chembl.fillna(False)
    cand = lab[(~lab.known_use) & lab.eligible & (lab.evidence == "matched")].copy()
    cand = cand.assign(cmap_name=cand.cmap_names.astype(str).str.split("|")).explode("cmap_name")
    cand["cmap_name"] = cand.cmap_name.str.strip()
    cand = cand.rename(columns={"cohort": "group_id"})[
        ["group_id", "fda_ingredient", "cmap_name", "drug", "anticancer"]]

    sp = pd.read_csv(RES / "drug_reversal_split_v2_tf_ms25.tsv", sep="\t", low_memory=False)
    keep = [c for c in ["group_id", "pert_id", "cmap_name", "source", "passes_gate", "rev_pct",
                        "spec_mech_pct", "layered_S", "layered_pct", "layered_D", "layered_A",
                        "layered_O", "layered_T", "layered_r", "n_mrs", "line_match"]
            if c in sp.columns]
    sp = sp[keep].drop_duplicates(["group_id", "cmap_name"])
    df = cand.merge(sp, on=["group_id", "cmap_name"], how="left")
    has_layered = df["layered_S"].notna() if "layered_S" in df.columns else pd.Series(False, index=df.index)
    df["layered_status"] = np.where(has_layered, "ok", "not_evaluable:no_layered_score")
    df["rev_status"] = np.where(df.rev_pct.notna(), "ok", "not_evaluable:no_reversal_score")

    pkn = pd.read_csv(ROOT / "results/systems_validation/pkn/signed_pkn.tsv", sep="\t")
    pkn_sources, pkn_nodes = set(pkn.source), set(pkn.source) | set(pkn.target)
    units = pd.read_csv(RES / "target_role/d121/units.tsv", sep="\t",
                        dtype={"tid": str, "action_flag": str})
    g = units[units.k == 1].groupby("cmap_name").genes.apply(lambda s: sorted(set(s))).to_dict()
    df["targets"] = df.cmap_name.map(g)
    df["n_targets"] = df.targets.map(lambda t: len(t) if isinstance(t, list) else 0)
    df["targets_in_pkn"] = df.targets.map(
        lambda t: sum(x in pkn_nodes for x in t) if isinstance(t, list) else 0)
    df["targets_with_outgoing"] = df.targets.map(
        lambda t: sum(x in pkn_sources for x in t) if isinstance(t, list) else 0)
    df["F_status"] = np.select(
        [df.n_targets == 0, df.targets_in_pkn == 0, df.targets_with_outgoing == 0],
        ["not_evaluable:no_single_gene_target", "not_evaluable:target_absent_from_PKN",
         "not_evaluable:target_has_no_outgoing_edge"], default="ok")

    testable = set()
    f134 = RES / "d134/feasibility.tsv"
    if f134.exists():
        fe = pd.read_csv(f134, sep="\t")
        col = next((c for c in fe.columns if c.lower() in ("genes", "gene", "target", "unit_genes")), None)
        if col:
            testable = set(fe[col].astype(str))
    crispr_drugs = set()
    f135 = RES / "d135/paired_cells.tsv"
    if f135.exists():
        pc = pd.read_csv(f135, sep="\t")
        dcol = next((c for c in pc.columns if c.lower() in ("cmap_name", "drug", "pert_iname")), None)
        if dcol:
            crispr_drugs = set(pc[dcol].astype(str))
    df["O_target_recoverable_unit"] = df.targets.map(
        lambda t: bool(set(t) & testable) if isinstance(t, list) else False)
    df["O_crispr_paired"] = df.cmap_name.isin(crispr_drugs)
    df["O_status"] = np.where(df.O_target_recoverable_unit | df.O_crispr_paired, "ok",
                              "not_evaluable:no_independent_perturbation_evidence")

    sig = pd.read_csv(ROOT / "data/raw/lincs2020/siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["pert_id", "pert_type", "pert_time", "is_hiq", "tas"])
    sig = sig[(sig.pert_type == "trt_cp") & (sig.is_hiq == 1) & (sig.tas >= 0.4)]
    tp = sig.groupby("pert_id").pert_time.nunique()
    df["S_max_timepoints"] = df.pert_id.map(tp).fillna(0).astype(int)
    df["S_status"] = np.where(df.S_max_timepoints >= 3, "ok",
                              "not_evaluable:<3_timepoints_and_no_washout_data_exists")

    cols = ["group_id", "fda_ingredient", "drug", "cmap_name", "pert_id", "anticancer", "source",
            "line_match", "n_mrs", "layered_pct", "layered_S", "layered_D", "layered_A",
            "layered_O", "layered_T", "layered_r", "layered_status", "rev_pct", "spec_mech_pct",
            "rev_status", "n_targets", "targets_in_pkn", "targets_with_outgoing", "F_status",
            "O_target_recoverable_unit", "O_crispr_paired", "O_status", "S_max_timepoints",
            "S_status"]
    df = df[[c for c in cols if c in df.columns]]
    sortcol = "layered_pct" if "layered_pct" in df.columns else "rev_pct"
    df = df.sort_values(["group_id", sortcol], ascending=[True, False])
    df.to_csv(OUT / "candidate_table.tsv", sep="\t", index=False)

    rows = []
    for a in ["layered", "rev", "F", "O", "S"]:
        st = df[f"{a}_status"]
        rows.append(dict(avenue=a, n_pairs=len(df), n_ok=int((st == "ok").sum()),
                         pct_ok=round(100 * (st == "ok").mean(), 1),
                         main_reason_not_evaluable=(st[st != "ok"].value_counts().index[0]
                                                    if (st != "ok").any() else "-")))
    cov = pd.DataFrame(rows)
    cov.to_csv(OUT / "coverage_summary.tsv", sep="\t", index=False)

    print(f"SPEC section 0 entry set: {len(df):,} pairs, {df.group_id.nunique()} cohorts, "
          f"{df.cmap_name.nunique():,} compounds\n")
    print("WHAT EACH AVENUE CAN SAY ANYTHING ABOUT")
    print(cov.to_string(index=False))
    allthree = ((df.layered_status == "ok") & (df.F_status == "ok") & (df.O_status == "ok"))
    print(f"\nlayered AND F AND O all evaluable : {int(allthree.sum()):,} "
          f"({100*allthree.mean():.1f}%)")
    print(f"plus S                            : "
          f"{int((allthree & (df.S_status=='ok')).sum()):,}")
    print("\nwrote", OUT / "candidate_table.tsv")


if __name__ == "__main__":
    sys.exit(main())
