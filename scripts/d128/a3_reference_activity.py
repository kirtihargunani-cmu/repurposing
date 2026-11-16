#!/usr/bin/env python
"""D128-A3.2: the reference-activity flag - is R even readable for this (target unit, cohort)?

R asks "does the drug push the tumor back toward normal?", with the cohort's P2 reference (GTEx normal tissue,
or adjacent normal) standing for "normal". When the target's own pathway is ACTIVE in that reference, R inverts by
construction: moving away from a pathway the normal tissue runs looks like moving away from normal, so R penalises
exactly the drug that works. D0.1 predicted this in advance for ESR1 in breast and A4 confirmed it. This script
computes the same thing as a general, falsifiable quantity for every (target unit, cohort) pair.

The quantity reuses the frozen scoring machinery, fed the DISEASE contrast instead of the drug profile:

    M_disease(unit, cohort, contrast) = sum_g w_g e_g f(z_disease_g) / sum_g w_g

same footprint, same expected signs `e` (A1's pre-action sign x the unit's frozen action sign), same weights
(wM = 1, P7), same bounded transform f = tanh(z / 2) (P9), and the frozen P1 columns of the step4_v2 / step10_v2
tables (`logFC` for the sign, `z` for the magnitude; the two agree in sign on every gene of every cohort, checked).
The frozen per-subset P10 floor (D128-A0.1, engine fix D128-A2.1) applies, so a unit that cannot carry M cannot
carry a flag either.

    M_disease > 0  the tumor ALREADY looks target-inhibited against the reference -> an inhibitor pushes further
                   the way the disease went, i.e. AWAY from normal -> R is misleading          reference_inverted
    M_disease < 0  the tumor looks target-ACTIVE against the reference -> the drug reverses it  reference_aligned
    M_disease ~ 0  the reference carries no target-pathway signal                               reference_neutral

Significance is decided against the two nulls A3 showed are well powered (sign-shuffle and gene-label, 200 seeded
draws each); the degree-matched null is degenerate for 95% of units and is not used. Nothing here gates, filters or
re-ranks anything: M, R and every A3/A4 output keep their definitions and values. Columns and a report only.

Outputs results/lincs/d128/a3/:
  reference_activity.tsv              one row per (cohort, level, unit_id, contrast)
  reference_activity_rows.tsv.gz      the scored (cohort, drug, unit) rows with the flag and the centred R joined
  reference_activity_cohort.tsv       per-cohort R baseline (median R over all drugs) and the tau summary
  reference_activity_summary.tsv      level distribution and the four pre-stated validation items
  reference_activity_run.json
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

A3 = C.OUT / "a3"
T0 = time.time()
LOG = []

# the significance rule, fixed in the D128-A3.2 decision row BEFORE the run
PCT_HI, PCT_LO = 97.5, 2.5          # two-sided 5%, and BOTH nulls must agree
N_NULL = 200                        # the frozen A3 null budget, same seeded generator (SEED = 128)


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def footprints():
    """A1's resolved footprints, in the form the engine takes (same preparation as a2_score.load_inputs)."""
    F = pd.read_parquet(C.OUT / "a1/target_gene_footprints.parquet")
    amb = pd.read_csv(C.OUT / "a1/path_ambiguity.tsv", sep="\t")
    if "n_sig_subset_tfs" not in amb.columns:
        raise SystemExit("a1/path_ambiguity.tsv predates D128-A2.1: re-run a1_build_paths.py")
    n_sig_tfs = dict(zip(amb.unit_id, amb.n_sig_subset_tfs))
    genes = IO.profile_genes()
    land = set(genes.gene_id[genes.feature_space == "landmark"])
    FP = {}
    for uid, d in F.groupby("unit_id", sort=False):
        res = d[d.status == "resolved"]
        g = res.gene_id.to_numpy()
        FP[uid] = dict(genes=g, expected=res.expected_sign_pre_action.to_numpy().astype(float),
                       subset=res.subset.to_numpy().astype(str),
                       lm=np.array([x in land for x in g], dtype=bool),
                       n_unresolved=int((d.status != "resolved").sum()),
                       n_sig_tfs=int(n_sig_tfs.get(uid, 0)))
    return FP, amb


def de_tables():
    """Every (cohort, contrast) disease contrast on disk, with its signature_summary status (P1 columns)."""
    out = {}
    for level, (proc, summ) in (("group_id", C.DE_GROUPS), ("stratum", C.DE_STRATA)):
        S = pd.read_csv(summ, sep="\t")
        key = "group_id" if "group_id" in S.columns else S.columns[0]
        for r in S.itertuples(index=False):
            cohort, con, status = getattr(r, key), r.contrast, r.status
            d = IO.de_table(proc, cohort, con)
            if d is None:
                continue
            z = d[C.DE_WEIGHT_COL].to_numpy(dtype=float)
            fc = d[C.DE_SIGN_COL].to_numpy(dtype=float)
            # P1: sign from logFC, magnitude from z. They agree on every gene of every cohort; the line below makes
            # the frozen convention explicit rather than relying on that.
            zz = np.sign(fc) * np.abs(z)
            out[(cohort, con)] = dict(level=level, status=status, z=dict(zip(d.gene_id, zz)),
                                      n_genes=len(d), n_sign_disagree=int((np.sign(fc) != np.sign(z)).sum()))
    return out


def gtex_tau():
    """tau(gene, cohort) = mean log2-CPM in the cohort's GTEx reference tissue - median over all GTEx tissues."""
    p = A3 / "gtex_tissue_logcpm.tsv.gz"
    if not p.exists():
        log("WARNING: gtex_tissue_logcpm.tsv.gz missing (run scripts/d128/gtex_tissue_expression.R); tau is NA")
        return None, {}
    E = pd.read_csv(p, sep="\t")
    E = E.drop_duplicates("gene_id").set_index("gene_id")
    med = E.median(axis=1)
    TAU = E.sub(med, axis=0)
    ref = pd.read_csv(C.ROOT / "results/step4_v2/gtex_reference.tsv", sep="\t")
    ref = ref[ref.use.astype(str).str.upper().isin(("TRUE", "1"))]
    cohort_tissues = {}
    for r in ref.itertuples(index=False):
        for grp in str(r.groups).split(";"):
            cohort_tissues.setdefault(grp.strip(), set()).add(r.tissue)
    return TAU, cohort_tissues


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    ap.add_argument("--a3", default="a3")
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    if not C.FROZEN:
        sys.exit("Section 0 is not frozen")
    A3.mkdir(parents=True, exist_ok=True)

    FP, amb = footprints()
    DE = de_tables()
    TAU, cohort_tissues = gtex_tau()
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str}).fillna({"action_flag": ""})
    U["genes_str"] = U.genes
    U["unit_id"] = U.cmap_name + "::" + U.tid.fillna("annot:" + U.genes_str)
    dup = U.unit_id.duplicated(keep=False)
    if dup.any():
        U.loc[dup, "unit_id"] = U.unit_id[dup] + "#" + U[dup].groupby("unit_id").cumcount().add(1).astype(str)
    U_by_id = {r.unit_id: r for r in U.itertuples(index=False)}

    # the cohorts and the contrast each one's R actually uses (P2, as a2_score selects it)
    SC_COH = pd.read_csv(C.OUT / a.a3 / "scores_by_cohort.tsv", sep="\t", low_memory=False)
    cohorts = SC_COH[["cohort", "level", "contrast"]].drop_duplicates().dropna(subset=["cohort"])
    p2 = {r.cohort: r.contrast for r in cohorts.itertuples(index=False) if isinstance(r.contrast, str)}
    levels = dict(zip(cohorts.cohort, cohorts.level))
    log(f"{len(FP)} units with a resolved footprint; {len(p2)} cohorts with a P2 contrast; "
        f"{len(DE)} (cohort, contrast) disease tables on disk")

    rows = []
    for (cohort, con), de in sorted(DE.items()):
        if cohort not in levels:
            continue
        zmap = de["z"]
        tiss = sorted(cohort_tissues.get(cohort, cohort_tissues.get(str(cohort).split(":")[0], set())))
        tau_col = None
        if TAU is not None and tiss:
            cols = [t for t in tiss if t in TAU.columns]
            tau_col = TAU[cols].mean(axis=1) if cols else None
        for uid, fp in FP.items():
            u = U_by_id[uid]
            action = 1 if u.orient == "raise" else -1
            keep = np.array([g in zmap for g in fp["genes"]], dtype=bool)
            gids = fp["genes"][keep]
            z = np.array([zmap[g] for g in gids], dtype=float)
            m = SC.score_mechanism(gids, fp["expected"][keep], action, z, landmark_mask=fp["lm"][keep],
                                   action_flag=u.action_flag, subset=fp["subset"][keep],
                                   n_signaling_tfs=fp["n_sig_tfs"])
            row = dict(cohort=cohort, level=levels[cohort], contrast=con, is_p2_contrast=(p2.get(cohort) == con),
                       de_status=de["status"], unit_id=uid, cmap_name=u.cmap_name, unit_genes=u.genes_str,
                       action=u.orient, action_flag=u.action_flag,
                       M_disease=m.value, status=m.status, abstain_reason=m.abstain_reason,
                       n_genes=m.n_genes, n_landmark=m.n_landmark, n_footprint=len(fp["genes"]),
                       n_missing_de=int((~keep).sum()), gtex_tissues="|".join(tiss))
            if m.scored:
                e = fp["expected"][keep] * action
                fz = SC.bounded(z)
                ss = NU.summarise(NU.m_null_sign_shuffle(e, fz, N_NULL), m.value)
                gl = NU.summarise(NU.m_null_gene_label(e, fz, fp["lm"][keep], N_NULL), m.value)
                row.update(sign_shuffle_pct=ss["percentile"], sign_shuffle_z=ss["z"],
                           gene_label_pct=gl["percentile"], gene_label_z=gl["z"],
                           null_median=ss["null_median"], null_sd=ss["null_sd"],
                           # how much of f() saturates on the disease scale: the disease z is far larger than a
                           # LINCS z, so M_disease is close to a signed sign-concordance. Reported, not corrected.
                           frac_saturated=float((np.abs(fz) > 0.9).mean()),
                           frac_small_z=float((np.abs(z) < 1.0).mean()))
                if m.value > 0 and ss["percentile"] >= PCT_HI and gl["percentile"] >= PCT_HI:
                    flag = "reference_inverted"
                elif m.value < 0 and ss["percentile"] <= PCT_LO and gl["percentile"] <= PCT_LO:
                    flag = "reference_aligned"
                else:
                    flag = "reference_neutral"
                if tau_col is not None:
                    t = tau_col.reindex(list(m.contributions)).to_numpy(dtype=float)
                    row.update(tau_weighted=float(np.nanmean(t)) if np.isfinite(t).any() else np.nan,
                               tau_n_genes=int(np.isfinite(t).sum()))
            else:
                flag = "reference_unavailable"
            row["r_reference_flag"] = flag
            rows.append(row)
        log(f"{cohort} / {con} ({de['status']}): {len(rows)} (cohort, unit, contrast) rows so far")
    REF = pd.DataFrame(rows)
    REF.to_csv(A3 / "reference_activity.tsv", sep="\t", index=False)
    log(f"reference_activity: {len(REF)} rows; "
        + ", ".join(f"{k} {v}" for k, v in REF.r_reference_flag.value_counts().items()))

    # --- companion diagnostic 1: the cohort R baseline ---------------------------------------------------------
    H = pd.read_csv(C.OUT / a.a3 / "line_heterogeneity.tsv", sep="\t", low_memory=False)
    base = SC_COH.groupby(["cohort", "level"]).agg(
        n_drugs=("pert_id", "nunique"), n_R=("R", "count"), R_median=("R", "median"),
        R_q25=("R", lambda x: x.quantile(0.25)), R_q75=("R", lambda x: x.quantile(0.75)),
        M_median=("M", "median"), n_M=("M", "count")).reset_index()
    base_unit = H.groupby("cohort").R.median().rename("R_median_unit_level").reset_index()
    base = base.merge(base_unit, on="cohort", how="left")
    # the flag mix per cohort, on the P2 contrast
    mix = (REF[REF.is_p2_contrast].groupby(["cohort", "r_reference_flag"]).size().unstack(fill_value=0)
           .reset_index())
    base = base.merge(mix, on="cohort", how="left")
    if TAU is not None:
        base = base.merge(REF[REF.is_p2_contrast].groupby("cohort").tau_weighted.median()
                          .rename("tau_weighted_median").reset_index(), on="cohort", how="left")
    base.sort_values("R_median", ascending=False).to_csv(A3 / "reference_activity_cohort.tsv", sep="\t", index=False)
    log("cohort R baseline (top and bottom 6 by median R):\n"
        + base.sort_values("R_median", ascending=False)[
            ["cohort", "level", "n_R", "R_median", "R_q25", "R_q75", "M_median"]].head(6).to_string()
        + "\n...\n"
        + base.sort_values("R_median", ascending=False)[
            ["cohort", "level", "n_R", "R_median", "R_q25", "R_q75", "M_median"]].tail(6).to_string())

    # --- validation 3: does the flag predict R? ------------------------------------------------------------------
    key = ["cohort", "unit_id"]
    J = H.merge(REF[REF.is_p2_contrast][key + ["r_reference_flag", "M_disease", "sign_shuffle_pct",
                                               "gene_label_pct", "tau_weighted", "n_genes"]]
                .rename(columns={"n_genes": "n_genes_disease"}), on=key, how="left", validate="many_to_one")
    J = J.merge(base[["cohort", "R_median"]].rename(columns={"R_median": "cohort_R_baseline"}),
                on="cohort", how="left")
    J["R_centred"] = J.R - J.cohort_R_baseline
    J.to_csv(A3 / "reference_activity_rows.tsv.gz", sep="\t", index=False)

    def effect(x, y):
        """Cliff's delta for x against y, from the Mann-Whitney U, with a tie-corrected normal-approximation p.

        |delta| 0.11 / 0.28 / 0.43 = small / medium / large by the usual convention. delta > 0 means x tends to be
        larger than y. scipy is not in the environment, so U is computed from average ranks (which is what
        scipy.stats.mannwhitneyu does) rather than imported.
        """
        x = np.asarray([v for v in x if np.isfinite(v)], float)
        y = np.asarray([v for v in y if np.isfinite(v)], float)
        nx, ny = len(x), len(y)
        if nx < 2 or ny < 2:
            return dict(n_x=nx, n_y=ny, delta=np.nan, p=np.nan)
        allv = np.concatenate([x, y])
        order = np.argsort(allv, kind="mergesort")
        ranks = np.empty(len(allv), float)
        sv = allv[order]
        i = 0
        tie_term = 0.0
        while i < len(sv):                                  # average ranks within each tie group
            j = i
            while j + 1 < len(sv) and sv[j + 1] == sv[i]:
                j += 1
            ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            t = j - i + 1
            tie_term += t ** 3 - t
            i = j + 1
        u = float(ranks[:nx].sum() - nx * (nx + 1) / 2.0)
        n = nx + ny
        mu = nx * ny / 2.0
        sd = np.sqrt(nx * ny / 12.0 * ((n + 1) - tie_term / (n * (n - 1))))
        from math import erfc, sqrt
        p = float(erfc(abs(u - mu) / (sd * sqrt(2)))) if sd > 0 else np.nan
        return dict(n_x=nx, n_y=ny, delta=float(2 * u / (nx * ny) - 1), p=p)

    def spearman(x, y):
        """Rank correlation without scipy (which is not in the environment): Pearson on average ranks."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 3:
            return np.nan
        rx = pd.Series(x[ok]).rank().to_numpy()
        ry = pd.Series(y[ok]).rank().to_numpy()
        return float(np.corrcoef(rx, ry)[0, 1])

    V = []
    sc = J[J.M.notna() | J.R.notna()]
    for col in ("R", "R_centred", "M"):
        for lv in ("reference_inverted", "reference_aligned", "reference_neutral"):
            d = sc.loc[sc.r_reference_flag == lv, col]
            V.append(dict(check=f"{col} by flag level", level=lv, n=int(d.notna().sum()),
                          median=float(d.median()) if d.notna().any() else np.nan,
                          mean=float(d.mean()) if d.notna().any() else np.nan,
                          q25=float(d.quantile(0.25)) if d.notna().any() else np.nan,
                          q75=float(d.quantile(0.75)) if d.notna().any() else np.nan))
        e = effect(sc.loc[sc.r_reference_flag == "reference_inverted", col],
                   sc.loc[sc.r_reference_flag == "reference_aligned", col])
        V.append(dict(check=f"{col}: inverted vs aligned", level="effect size (Cliff's delta)",
                      n=e["n_x"] + e["n_y"], median=e["delta"], mean=e["p"], q25=e["n_x"], q75=e["n_y"]))
    # supporting evidence, pre-stated: within inverted rows M and R should oppose, within aligned rows agree
    for lv in ("reference_inverted", "reference_aligned", "reference_neutral"):
        d = sc[(sc.r_reference_flag == lv) & sc.M.notna() & sc.R.notna()]
        V.append(dict(check="Spearman rho(M, R) within flag level", level=lv, n=len(d),
                      median=spearman(d.M.to_numpy(), d.R.to_numpy()) if len(d) > 10 else np.nan,
                      mean=np.nan, q25=np.nan, q75=np.nan))
    VAL = pd.DataFrame(V)

    # --- the pre-stated named cases ------------------------------------------------------------------------------
    named = []
    cases = [("ESR1 x BRCA-LumA (validation 1: must NOT be aligned)", "ESR1", ("BRCA-LumA",)),
             ("ESR1 x BRCA-LumB (validation 1: must NOT be aligned)", "ESR1", ("BRCA-LumB",)),
             ("ESR1 x BRCA-Her2 (reported)", "ESR1", ("BRCA-Her2",)),
             ("ESR1 x BRCA-Basal (reported)", "ESR1", ("BRCA-Basal",)),
             ("ERBB2 x BRCA-Her2 (validation 4: expected aligned/neutral)", "ERBB2", ("BRCA-Her2",)),
             ("BRAF x SKCM (validation 4: expected aligned/neutral)", "BRAF", ("TCGA-SKCM",))]
    for label, gene, cohs in cases:
        d = REF[REF.is_p2_contrast & REF.cohort.isin(cohs)
                & REF.unit_genes.fillna("").apply(lambda g: gene in str(g).split("|"))]
        for r in d.itertuples(index=False):
            named.append(dict(case=label, cohort=r.cohort, unit=r.unit_id, genes=r.unit_genes, action=r.action,
                              M_disease=r.M_disease, flag=r.r_reference_flag,
                              sign_shuffle_pct=getattr(r, "sign_shuffle_pct", np.nan),
                              gene_label_pct=getattr(r, "gene_label_pct", np.nan),
                              n_genes=r.n_genes, tau=getattr(r, "tau_weighted", np.nan),
                              abstain_reason=r.abstain_reason))
    NAMED = pd.DataFrame(named)

    dist = (REF.groupby(["contrast", "r_reference_flag"]).size().rename("n").reset_index())
    dist["share_of_contrast"] = dist.n / dist.groupby("contrast").n.transform("sum")
    elig = REF[REF.status == "scored"]
    dist_el = (elig.groupby(["contrast", "r_reference_flag"]).size().rename("n").reset_index())
    dist_el["share_of_eligible"] = dist_el.n / dist_el.groupby("contrast").n.transform("sum")
    SUM = pd.concat([
        dist.assign(table="all (unit, cohort, contrast) pairs"),
        dist_el.assign(table="eligible pairs only (M_disease scored)")], ignore_index=True)
    SUM.to_csv(A3 / "reference_activity_summary.tsv", sep="\t", index=False)
    with open(A3 / "reference_activity_validation.tsv", "w") as fh:
        VAL.to_csv(fh, sep="\t", index=False)
    NAMED.to_csv(A3 / "reference_activity_named_cases.tsv", sep="\t", index=False)
    log("level distribution:\n" + SUM.to_string())
    log("validation:\n" + VAL.to_string())
    log("named cases:\n" + (NAMED.to_string() if len(NAMED) else "none found"))

    max_share = float(dist_el.share_of_eligible.max()) if len(dist_el) else np.nan
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    (A3 / "reference_activity_run.json").write_text(json.dumps(dict(
        stage="D128-A3.2", git_sha=sha, frozen=C.FROZEN, n_null=N_NULL, pct_hi=PCT_HI, pct_lo=PCT_LO,
        rows=len(REF), units=int(REF.unit_id.nunique()), cohorts=int(REF.cohort.nunique()),
        contrasts=sorted(REF.contrast.unique()),
        flag_counts=REF.r_reference_flag.value_counts().to_dict(),
        flag_counts_eligible=elig.r_reference_flag.value_counts().to_dict(),
        max_level_share_of_eligible=max_share,
        discrimination_ok=bool(max_share <= 0.90) if np.isfinite(max_share) else None,
        de_sign_disagreements={f"{k[0]}|{k[1]}": v["n_sign_disagree"] for k, v in DE.items() if v["n_sign_disagree"]},
        validation=VAL.to_dict("records"), named_cases=NAMED.to_dict("records"),
        log=LOG, seconds=round(time.time() - T0, 1)), indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
