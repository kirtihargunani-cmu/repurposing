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
import landmark as LM  # noqa: E402

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
    ap.add_argument("--nulls-run", default="a2_nulls",
                    help="directory of the gated null pass (D128-A3 landmark amendment); used if present")
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
    if "gene_arm" not in S.columns:
        S["gene_arm"] = "all_bing"
    n_units_ann = U.groupby("cmap_name").size().to_dict()
    arms = {}
    for ga, d in S.groupby("gene_arm"):
        un = AG.lines_to_unit(d)
        dr = AG.units_to_drug(un, n_units_ann)
        arms[ga] = (un, dr)
        log(f"gene arm {ga}: {len(un)} unit rows (M {int(un.M.notna().sum())}, R {int(un.R.notna().sum())}); "
            f"{len(dr)} drug rows (M {int(dr.M.notna().sum())}, R {int(dr.R.notna().sum())})")
    S = S[S.gene_arm == "all_bing"]           # everything downstream is the frozen primary arm (P4)
    UN = arms["all_bing"][0]
    # Engine-vs-frozen-floor audit. P10 as frozen by D128-A0.1 floors each SUBSET of a unit separately (a direct-TF
    # subset on its zero-hop regulon genes, a signalling subset on resolved genes AND >= 2 anchored TFs). score.py
    # applies only the >= 10 usable genes part, to the merged footprint. A unit that fails the frozen per-subset
    # floor can therefore still clear the engine's. These rows are FLAGGED, not silently kept or silently dropped.
    a0u = E.dropna(subset=["unit_id"]).drop_duplicates(["cohort", "pert_id", "unit_id"])[
        ["cohort", "pert_id", "unit_id", "M_status"]].rename(columns={"M_status": "a0_M_unit_status"})
    UN = UN.merge(a0u, on=["cohort", "pert_id", "unit_id"], how="left")
    UN["scored_below_a0_floor"] = UN.M.notna() & (UN.a0_M_unit_status == "unscorable")
    log(f"engine vs frozen P10 floor: {int(UN.scored_below_a0_floor.sum())} unit rows scored by the engine although "
        f"the frozen per-subset floor marks the unit unscorable "
        f"({UN.loc[UN.scored_below_a0_floor, 'cmap_name'].nunique()} compounds)")
    arms["all_bing"] = (UN, arms["all_bing"][1])
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
    D = arms["all_bing"][1].merge(
        UN.groupby(["cohort", "level", "pert_id"]).scored_below_a0_floor.sum()
          .rename("n_units_scored_below_a0_floor").reset_index(),
        on=["cohort", "level", "pert_id"], how="left", validate="one_to_one")
    D["n_units_scored_below_a0_floor"] = D.n_units_scored_below_a0_floor.fillna(0).astype(int)
    log(f"drug-level rows affected by the floor mismatch: "
        f"{int((D.n_units_scored_below_a0_floor > 0).sum())} of {len(D)} "
        f"(M scored on at least one unit the frozen floor would not admit)")
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
    bad = ID[ID.fails]
    if len(bad):
        sys.exit("a baseline value-identity check FAILED:\n" + bad.to_string())
    noted = ID[(ID.n_outside_tol > 0) & ~ID.fails]
    if len(noted):
        log("baseline copies differ but the difference is explained (recorded, not a failure):\n"
            + noted[["check", "n", "n_exact", "n_outside_tol", "max_abs_diff", "why"]].to_string())

    # --- D128-A3 landmark sensitivity (amendment 2026-11-16): pre-stated criterion, evaluated here ----------------
    if "landmark_only" in arms:
        # the landmark arm is ranked in exactly the same population as the primary arm, so a percentile gap is a
        # gap in the score and not in the denominator
        dr_lm = arms["landmark_only"][1].merge(
            D[["cohort", "level", "pert_id", "in_comparison_population"]], on=["cohort", "level", "pert_id"],
            how="left", validate="one_to_one")
        dr_lm["in_comparison_population"] = dr_lm.in_comparison_population.fillna(False).astype(bool)
        dr_lm = AG.add_percentiles(dr_lm, mask=dr_lm.in_comparison_population)
        T, P = LM.compare(arms["all_bing"][0], arms["landmark_only"][0], D, dr_lm)
        T.to_csv(A3 / "landmark_sensitivity.tsv", sep="\t", index=False)
        SUM = LM.summary(T, P)
        SUM.to_csv(A3 / "landmark_sensitivity_summary.tsv", sep="\t", index=False)
        log("landmark sensitivity (pre-stated criterion: rho >= 0.90, flips <= 5%, median |dpct| <= 5):\n"
            + SUM.to_string())
        crit_met = bool(SUM.loc[(SUM.subset == "all units") & (SUM.score == "M"), "interchangeable"].iloc[0])
        log(f"CRITERION for M over all units: {'MET' if crit_met else 'NOT MET'} "
            f"-> population-wide nulls are {'run' if crit_met else 'NOT run'} (decision row, fixed in advance)")
        D = D.merge(arms["landmark_only"][1][["cohort", "level", "pert_id", "M", "R", "n_units_M"]].rename(
            columns={"M": "M_landmark", "R": "R_landmark", "n_units_M": "n_units_M_landmark"}),
            on=["cohort", "level", "pert_id"], how="left", validate="one_to_one")
        D = D.merge(P[["cohort", "level", "pert_id", "M_pct_landmark", "R_pct_landmark"]].drop_duplicates(
            ["cohort", "level", "pert_id"]), on=["cohort", "level", "pert_id"], how="left", validate="one_to_one")
    else:
        crit_met = None

    # --- null power: how many degree-matched candidates actually survive (user question, A3 amendment) ------------
    npath = C.OUT / a.a2 / "nulls.parquet"
    nullp = C.OUT / a.nulls_run / "nulls.parquet"
    src = nullp if nullp.exists() else npath
    NP = None
    if src.exists():
        N = pd.read_parquet(src)
        NP = N.groupby(["gene_arm", "null"] if "gene_arm" in N.columns else ["null"]).agg(
            n_rows=("n_null", "size"), n_null_min=("n_null", "min"), n_null_median=("n_null", "median"),
            n_null_max=("n_null", "max"), frac_lt_10=("n_null", lambda x: round(float((x < 10).mean()), 4)),
            frac_lt_3=("n_null", lambda x: round(float((x < 3).mean()), 4)),
            median_pct=("percentile", "median"), median_z=("z", "median"),
            frac_pct_ge_95=("percentile", lambda x: round(float((x >= 95).mean()), 4))).reset_index()
        NP.to_csv(A3 / "null_power.tsv", sep="\t", index=False)
        log("null power:\n" + NP.to_string())
        nn = N[N.null == "degree_matched_target"].drop_duplicates(["unit_id"])
        log(f"degree-matched null: {len(nn)} distinct units; "
            f"{int((nn.n_null < 10).sum())} with < 10 surviving candidates, "
            f"{int((nn.n_null < 3).sum())} degenerate (< 3)")
        for k in ("sign_shuffle", "gene_label", "degree_matched_target"):
            sub = N[N.null == k]
            if len(sub):
                D = D.merge(sub.groupby(["cohort", "pert_id"]).agg(
                    **{f"{k}_pct": ("percentile", "mean"), f"{k}_z": ("z", "mean"),
                       f"{k}_n": ("n_null", "median")}).reset_index(), on=["cohort", "pert_id"], how="left")

    cols = [c for c in ["cohort", "level", "pert_id", "cmap_name", "fda_ingredient", "source", "split_source",
                        "contrast", "in_comparison_population", "passes_gate", "a0_status",
                        "M", "M_pct", "R", "R_pct", "n_M_ranked", "n_R_ranked",
                        "M_landmark", "M_pct_landmark", "R_landmark", "R_pct_landmark", "n_units_M_landmark",
                        "sign_shuffle_pct", "sign_shuffle_z", "sign_shuffle_n", "gene_label_pct", "gene_label_z",
                        "gene_label_n", "degree_matched_target_pct", "degree_matched_target_z",
                        "degree_matched_target_n",
                        "n_units_annotated", "n_units_in_run", "n_units_M", "n_units_R", "unit_scored_fraction",
                        "n_units_scored_below_a0_floor",
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
    Sm = S[key + ["M_status", "M_abstain_reason", "R_status", "R_abstain_reason"]].rename(
        columns={"M_status": "engine_M", "M_abstain_reason": "engine_M_reason",
                 "R_status": "engine_R", "R_abstain_reason": "engine_R_reason"})
    Ek = Ek.merge(Sm, on=key, how="left")
    for c in ("engine_M", "engine_R"):
        Ek[c] = Ek[c].fillna("not_run")
    for c in ("engine_M_reason", "engine_R_reason"):
        Ek[c] = Ek[c].fillna("")
    AB = (Ek.rename(columns={"M_status": "a0_predicted_M", "R_status": "a0_predicted_R"})
          .groupby(["level", "status", "a0_predicted_M", "engine_M", "engine_M_reason", "a0_predicted_R",
                    "engine_R", "reason"], dropna=False).size().rename("n_rows")
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
        gene_arms=sorted(arms), landmark_criterion_met=crit_met,
        landmark_summary=SUM.to_dict("records") if "landmark_only" in arms else None,
        null_power=NP.to_dict("records") if NP is not None else None,
        identity_checks=ID.to_dict("records"), log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
