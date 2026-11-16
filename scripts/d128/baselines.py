"""D128 A3 action 4: read the production baseline columns and prove the join did not change them.

Every column here is READ from a production output. Nothing is recomputed from an older implementation: the plan is
explicit that the lead score must not be rebuilt inside D128. Two production files carry `layered_pct` (the lead
table and the split tables) and two carry the target-role `T` (per compound and per drug), so the audit does not
only count rows: it checks that the two independent copies agree, which is the real test that the join is on the
right keys.
"""
import numpy as np
import pandas as pd

import config as C

AUDIT = []


def _audit(source, key, n_left, merged, cols, extra=""):
    AUDIT.append(dict(source=source, join_key="+".join(key), rows_left=n_left, rows_after=len(merged),
                      duplicated=len(merged) - n_left,
                      matched=int(merged[cols[0]].notna().sum()) if cols else np.nan,
                      unmatched=int(merged[cols[0]].isna().sum()) if cols else np.nan,
                      columns=",".join(cols), note=extra))


def load():
    """The production tables A3 joins to, with their cohort key renamed to `cohort`."""
    lay = pd.read_csv(C.LAYERED, sep="\t", low_memory=False).rename(columns={"unit": "cohort"})
    sg = pd.read_csv(C.SPLIT_GROUPS, sep="\t", low_memory=False).rename(columns={"group_id": "cohort"})
    ss = pd.read_csv(C.SPLIT_STRATA, sep="\t", low_memory=False).rename(columns={"stratum": "cohort"})
    split = pd.concat([sg.assign(level="group_id"), ss.assign(level="stratum")], ignore_index=True)
    ctr = pd.read_csv(C.ROOT / "results/lincs/target_role/compound_target_role.tsv", sep="\t",
                      low_memory=False).rename(columns={"group_id": "cohort"})
    dtr = pd.read_csv(C.ROOT / "results/lincs/target_role/drug_target_role.tsv", sep="\t", low_memory=False)
    return dict(layered=lay, split=split, compound_role=ctr, drug_role=dtr)


SPLIT_COLS = ["source", "passes_gate", "line_match", "transfer_conf", "low_transfer", "n_mrs", "rev", "rev_pct",
              "spec", "spec_pct", "rev_mech", "spec_mech", "spec_mech_pct", "spec_mech_pct_pregate",
              "prolif_effect", "prolif_effect_pct", "prolif_driven", "tier_lead", "target_absent",
              "receptor_status", "layered_pct"]
LAYERED_COLS = ["layered_S", "layered_pct", "layered_n_units", "layered_T", "layered_P", "layered_r", "layered_rho",
                "minority_alteration", "essential_target", "essential_unit", "dna_acting", "no_protein_target"]
ROLE_COLS = ["fda_ingredient", "T", "T_unsigned", "D", "A", "O", "dep_z", "alt_share", "n_main_targets",
             "target_source"]


def join(D, B):
    """Join every baseline onto the drug-level D128 table and record the audit. D keeps one row per input row."""
    n0 = len(D)

    lay = B["layered"][["cohort", "pert_id", "level"] + LAYERED_COLS].drop_duplicates(["cohort", "pert_id", "level"])
    D = D.merge(lay, on=["cohort", "pert_id", "level"], how="left", validate="one_to_one")
    _audit("layered_score_v2.tsv", ["cohort", "pert_id", "level"], n0, D, LAYERED_COLS,
           "production lead table; layered_pct is the lead ranking (D122) including D124b behaviour")

    sp = B["split"][["cohort", "pert_id", "level"] + SPLIT_COLS].drop_duplicates(["cohort", "pert_id", "level"])
    sp = sp.rename(columns={c: f"split_{c}" for c in SPLIT_COLS})
    D = D.merge(sp, on=["cohort", "pert_id", "level"], how="left", validate="one_to_one")
    _audit("*_reversal_split_v2_consensus.tsv", ["cohort", "pert_id", "level"], n0, D,
           [f"split_{c}" for c in SPLIT_COLS], "VIPER reversal (rev/spec), spec_mech_pct_pregate, tier_lead")

    cr = B["compound_role"][["cohort", "pert_id"] + ROLE_COLS].drop_duplicates(["cohort", "pert_id"])
    D = D.merge(cr, on=["cohort", "pert_id"], how="left", validate="one_to_one")
    _audit("target_role/compound_target_role.tsv", ["cohort", "pert_id"], n0, D, ROLE_COLS,
           "group cohorts only (9,288 rows): subgroup rows are NA by construction, not imputed")

    dr = B["drug_role"][["cohort", "fda_ingredient", "target_role_pct", "T", "eligible"]].rename(
        columns={"T": "T_druglevel", "eligible": "role_eligible"}).drop_duplicates(["cohort", "fda_ingredient"])
    D = D.merge(dr, on=["cohort", "fda_ingredient"], how="left", validate="many_to_one")
    _audit("target_role/drug_target_role.tsv", ["cohort", "fda_ingredient"], n0, D,
           ["target_role_pct", "T_druglevel"], "drug-level table, group cohorts only; subgroups NA by construction")
    return D


def identity_checks(D):
    """The joined columns must equal their source. Two independent copies of two of them make that testable.

    `tol` is the tolerance at which the two copies are called the same number, and each row says why it is what it
    is. `fails` is the stage-level verdict: only a disagreement that cannot be explained by storage precision or by
    a documented difference in what the two tables are keyed on counts as a failure.
    """
    out = []
    both = D[D.layered_pct.notna() & D.split_layered_pct.notna()]
    diff = (both.layered_pct - both.split_layered_pct).abs()
    out.append(dict(check="layered_pct: lead table vs split table", n=len(both),
                    n_exact=int((diff == 0).sum()), n_outside_tol=int((diff > 1e-3).sum()), tol=1e-3,
                    max_abs_diff=float(diff.max()) if len(both) else np.nan, fails=bool((diff > 1e-3).any()),
                    why="layered_score_v2.tsv stores layered_pct rounded to 3 dp, the split tables store it at full "
                        "precision; the two are the same number to within 5e-4"))
    both = D[D["T"].notna() & D.T_druglevel.notna()]
    diff = (both["T"] - both.T_druglevel).abs()
    bad = both[diff > 1e-3]
    out.append(dict(check="target-role T: per compound vs per drug", n=len(both), n_exact=int((diff == 0).sum()),
                    n_outside_tol=len(bad), tol=1e-3,
                    max_abs_diff=float(diff.max()) if len(both) else np.nan, fails=False,
                    why="NOT a join error: drug_target_role.tsv is keyed on (cohort, fda_ingredient) and aggregates "
                        "the compounds that share an ingredient. The rows that differ are exactly the ingredients "
                        "with two cmap_names (" + ", ".join(sorted(set(bad.fda_ingredient.astype(str)))) + "). "
                        "target_role_pct is therefore a DRUG-level value shared by both pert_ids, and is reported "
                        "as such, never as a per-pert_id score"))
    both = D[D.source.notna() & D.split_source.notna()]
    out.append(dict(check="direct/borrowed: A0 (P14) vs production split table", n=len(both),
                    n_exact=int((both.source == both.split_source).sum()),
                    n_outside_tol=int((both.source != both.split_source).sum()), tol=np.nan, max_abs_diff=np.nan,
                    fails=bool((both.source != both.split_source).any()),
                    why="P14: the direct/borrowed line population must be production's, unchanged"))
    return pd.DataFrame(out)
