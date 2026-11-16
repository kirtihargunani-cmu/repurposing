#!/usr/bin/env python
"""D128 A3: cohort aggregation, uncertainty and baseline joins (plan Section 4 A3; decision row D128-A3).

Reads the full A2 run (`results/lincs/d128/a2/scores_by_line.parquet`) and produces one row per (cohort, drug) with
M and R, their percentiles within the declared comparison population, their spread across source lines and target
units, and the production baseline columns read (never recomputed) beside them.

Outputs results/lincs/d128/a3/:
  scores_by_cohort.tsv / .parquet   one row per (cohort, level, pert_id)
  line_heterogeneity.tsv            one row per (cohort, level, pert_id, unit): per-line spread and line genotype
  baseline_join_audit.tsv           rows in / matched / unmatched / duplicated, plus the value-identity checks
  abstentions.tsv                   100% accounting of the A0 candidate rows, engine verdict vs A0 prediction
  a3_run.json
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

import aggregate as AG  # noqa: E402
import baselines as BL  # noqa: E402
import config as C  # noqa: E402

A3 = C.OUT / "a3"
T0 = time.time()
LOG = []


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--a2", default="a2")
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen")
    A3.mkdir(parents=True, exist_ok=True)

    S = pd.read_parquet(C.OUT / a.a2 / "scores_by_line.parquet")
    E = pd.read_csv(C.OUT / "a0/eligibility.tsv.gz", sep="\t", low_memory=False)
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str})
    ctx = pd.read_csv(C.OUT / "a0/source_line_context.tsv", sep="\t", low_memory=False)
    log(f"A2 rows {len(S)}; A0 candidate rows {len(E)}")

    # --- P14 re-assertion: the direct/borrowed call must still reproduce production on every pair ------------------
    pairs = E.drop_duplicates(["cohort", "pert_id"])
    agree = int((pairs.source_prod == pairs.source_a0).sum())
    log(f"P14 check: direct/borrowed agrees with production on {agree} / {len(pairs)} pairs")
    if agree != len(pairs):
        sys.exit(f"P14 broken: {len(pairs) - agree} pairs disagree with the production split table")

    # --- step 1: source lines -> target unit -----------------------------------------------------------------------
    UN = AG.lines_to_unit(S)
    log(f"unit level: {len(UN)} (cohort, pert_id, unit) rows; M {int(UN.M.notna().sum())}, R {int(UN.R.notna().sum())}")

    # line heterogeneity, with the source line's genotype so spread from genotype mismatch is visible
    gcols = ["ModelID", "StrippedCellLineName", "OncotreeLineage", "BRAF_V600", "EGFR_sensitizing", "EGFR_T790M",
             "KRAS_NRAS_hotspot", "ERBB2_amp_candidate", "RB1_null_candidate", "ESR1_tpm_log1p"]
    sc = S[S.M_status == "scored"].merge(ctx[gcols], on="ModelID", how="left")
    geno = sc.groupby(["cohort", "level", "pert_id", "unit_id"]).agg(
        lines_used="|".join if False else ("StrippedCellLineName", lambda x: "|".join(map(str, x))),
        lineages=("OncotreeLineage", lambda x: "|".join(sorted(set(map(str, x))))),
        n_BRAF_V600=("BRAF_V600", "sum"), n_EGFR_sens=("EGFR_sensitizing", "sum"),
        n_KRAS_NRAS=("KRAS_NRAS_hotspot", "sum"), n_ERBB2_amp=("ERBB2_amp_candidate", "sum"),
        n_RB1_null=("RB1_null_candidate", "sum")).reset_index()
    H = UN.merge(geno, on=["cohort", "level", "pert_id", "unit_id"], how="left")
    H["borrowed_fraction"] = (H.source == "borrowed").astype(float)
    H = H.merge(E.drop_duplicates(["cohort", "pert_id"])[["cohort", "pert_id", "line_match", "transfer_conf",
                                                          "passes_gate"]], on=["cohort", "pert_id"], how="left")
    H.to_csv(A3 / "line_heterogeneity.tsv", sep="\t", index=False)
    log(f"line_heterogeneity: {len(H)} rows; M sd defined on {int(H.M_sd.notna().sum())}")

    # --- step 2: target units -> drug ------------------------------------------------------------------------------
    n_units_ann = U.groupby("cmap_name").size().to_dict()
    D = AG.units_to_drug(UN, n_units_ann)
    log(f"drug level: {len(D)} (cohort, pert_id) rows; M {int(D.M.notna().sum())}, R {int(D.R.notna().sum())}")

    # --- baselines, read not recomputed ----------------------------------------------------------------------------
    B = BL.load()
    D = BL.join(D, B)
    # production excludes gate-failing rows from the lead ranking; D128 flags them and excludes them from percentiles
    D = D.merge(E.drop_duplicates(["cohort", "pert_id"])[["cohort", "pert_id", "passes_gate", "line_match",
                                                          "transfer_conf", "status"]].rename(
        columns={"status": "a0_status"}), on=["cohort", "pert_id"], how="left", validate="one_to_one")
    D["in_comparison_population"] = D.passes_gate.fillna(False).astype(bool)
    D = AG.add_percentiles(D, mask=D.in_comparison_population)
    AU = pd.DataFrame(BL.AUDIT)
    ID = BL.identity_checks(D)
    AU = pd.concat([AU, ID.assign(source="value identity")], ignore_index=True)
    AU.to_csv(A3 / "baseline_join_audit.tsv", sep="\t", index=False)
    log("baseline join audit:\n" + AU.to_string())
    bad = ID[(ID.n_disagree > 0)]
    if len(bad):
        log("WARNING: a baseline value-identity check disagrees:\n" + bad.to_string())

    cols = [c for c in ["cohort", "level", "pert_id", "cmap_name", "fda_ingredient", "source", "split_source",
                        "contrast", "in_comparison_population", "passes_gate", "a0_status",
                        "M", "M_pct", "R", "R_pct", "n_M_ranked", "n_R_ranked",
                        "n_units_annotated", "n_units_in_run", "n_units_M", "n_units_R", "unit_scored_fraction",
                        "M_unit_sd", "M_unit_min", "M_unit_max", "R_unit_sd", "M_line_sd_max",
                        "n_lines_M", "n_lines_R", "M_n_genes", "M_n_landmark", "R_n_genes",
                        "hub_fraction_structural", "hub_tf", "units", "unit_kinds",
                        "M_abstain_reason", "R_abstain_reason",
                        "layered_pct", "layered_S", "layered_T", "layered_P", "layered_r", "layered_rho",
                        "minority_alteration", "essential_target", "dna_acting", "no_protein_target",
                        "split_rev", "split_rev_pct", "split_spec_mech", "split_spec_mech_pct",
                        "split_spec_mech_pct_pregate", "split_prolif_effect_pct", "split_prolif_driven",
                        "split_tier_lead", "split_target_absent", "split_receptor_status",
                        "target_role_pct", "T", "T_unsigned", "D", "A", "O", "dep_z", "alt_share",
                        "n_main_targets", "target_source", "line_match", "transfer_conf"] if c in D.columns]
    D[cols].to_csv(A3 / "scores_by_cohort.tsv", sep="\t", index=False)
    D[cols].to_parquet(A3 / "scores_by_cohort.parquet", index=False)

    # --- abstentions: 100% of the A0 candidate rows ----------------------------------------------------------------
    key = ["cohort", "pert_id", "ModelID", "unit_id"]
    Ek = E.copy()
    Ek["unit_id"] = Ek.unit_id.fillna("")
    Ek["ModelID"] = Ek.ModelID.fillna("")
    Sk = S[key].assign(entered_a2=True)
    Ek = Ek.merge(Sk, on=key, how="left")
    Ek["entered_a2"] = Ek.entered_a2.fillna(False)
    Sm = S[key + ["M_status", "M_abstain_reason", "R_status", "R_abstain_reason"]]
    Ek = Ek.merge(Sm, on=key, how="left", suffixes=("_a0", "_engine"))
    Ek["engine_M"] = Ek.M_status.fillna("not_run")
    Ek["engine_M_reason"] = Ek.M_abstain_reason.fillna("")
    AB = (Ek.groupby(["level", "status", "M_status_a0" if "M_status_a0" in Ek.columns else "M_status",
                      "engine_M", "engine_M_reason", "reason"], dropna=False).size().rename("n_rows")
          .reset_index().sort_values("n_rows", ascending=False))
    AB.to_csv(A3 / "abstentions.tsv", sep="\t", index=False)
    assert int(AB.n_rows.sum()) == len(E), "abstentions.tsv does not account for every A0 candidate row"
    log(f"abstentions: {len(AB)} categories covering all {len(E)} A0 candidate rows; "
        f"entered A2 {int(Ek.entered_a2.sum())}")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A3 / "a3_run.json").write_text(json.dumps(dict(
        stage="D128-A3", git_sha=sha, frozen=C.FROZEN, a2_run=a.a2,
        p14_pairs_agreeing=agree, p14_pairs=len(pairs),
        counts=dict(a2_rows=len(S), unit_rows=len(UN), drug_rows=len(D),
                    drug_rows_in_population=int(D.in_comparison_population.sum()),
                    M_scored=int(D.M.notna().sum()), R_scored=int(D.R.notna().sum()),
                    M_ranked=int(D.M_pct.notna().sum()), R_ranked=int(D.R_pct.notna().sum())),
        identity_checks=ID.to_dict("records"), log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
