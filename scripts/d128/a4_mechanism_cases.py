#!/usr/bin/env python
"""D128 A4 action 5: analyse failures by MECHANISM, not by rank.

The plan is explicit that "it ranked badly" is not a diagnosis. Each case below is a *reason a mechanism score
cannot be trusted or cannot be computed*, defined from data that does not involve the drug's rank:

  missing_network_path        the annotated target reaches no regulator in the frozen prior, so no footprint exists
  below_floor                 a footprint exists but is smaller than the P10 floor
  ambiguous_action            ChEMBL gives no unambiguous direction (P11 forbids assuming one for M)
  sign_contradiction          the footprint's own routes disagree about a large share of its genes
  inferred_only_genes         the score rests entirely on L1000-inferred genes; nothing in it was measured
  absent_target_expression    the annotated target is not expressed in the source line the profile came from
  cytotoxic_signature         production already reads the drug's response as proliferation-driven
  line_mismatch               the score is borrowed from lines outside the disease's tissue key
  unknown_genotype            the control's genotype is not established for the source line
  polypharmacology            the drug's own target units disagree with each other

Writes one markdown report per case under results/lincs/d128/a4/mechanism_case_reports/, plus an index with the
counts. Reading, counting and writing only: nothing here changes a score or a parameter.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import d128io as IO  # noqa: E402

OUTDIR = C.OUT / "a4/mechanism_case_reports"
TPM_FLOOR = 1.0     # log2(TPM+1) below this counts as "not expressed" (declared in the D128-A4 decision row)


def md_table(d, n=12):
    d = d.head(n)
    if not len(d):
        return "_(no rows)_\n"
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    body = "".join("| " + " | ".join("" if pd.isna(v) else (f"{v:.4g}" if isinstance(v, float) else str(v))
                                     for v in r) + " |\n" for r in d.itertuples(index=False))
    return head + body


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a3", default="a3")
    a = ap.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    S = pd.read_csv(C.OUT / a.a3 / "scores_by_cohort.tsv", sep="\t", low_memory=False)
    H = pd.read_csv(C.OUT / a.a3 / "line_heterogeneity.tsv", sep="\t", low_memory=False)
    E = pd.read_csv(C.OUT / "a0/eligibility.tsv.gz", sep="\t", low_memory=False)
    cov = pd.read_csv(C.OUT / "a0/coverage.tsv", sep="\t")
    ctx = pd.read_csv(C.OUT / "a0/source_line_context.tsv", sep="\t", low_memory=False)
    ctl = pd.read_csv(C.OUT / "a4/control_results.tsv", sep="\t", low_memory=False)

    N = len(S)
    cases, index = {}, []

    def case(name, d, what, means, cols):
        cases[name] = (d, what, means, cols)
        index.append(dict(case=name, n_rows=len(d), share_of_drug_rows=round(len(d) / N, 4) if N else np.nan,
                          what_it_is=what))

    # --- cases that stop a score from existing at all --------------------------------------------------------------
    ab = S[S.M.isna()]
    case("missing_network_path", ab[ab.M_abstain_reason.fillna("").str.contains("no_footprint|no_anchored")],
         "the annotated target reaches no regulator in the frozen OmniPath prior at 0-2 hops, so there is no "
         "footprint to score",
         "not a scoring failure but a coverage limit of the prior; M must abstain and does",
         ["cohort", "pert_id", "cmap_name", "n_units_annotated", "M_abstain_reason", "layered_pct"])
    case("below_floor", ab[ab.M_abstain_reason.fillna("").str.contains("below_floor")],
         f"a footprint exists but holds fewer than the P10 floor of {C.FLOOR_GENES} sign-resolved genes",
         "the abstention is correct; raising the floor would lose more rows, lowering it would score noise",
         ["cohort", "pert_id", "cmap_name", "M_abstain_reason", "layered_pct"])
    case("ambiguous_action", ab[ab.M_abstain_reason.fillna("").str.contains("ambiguous_action")],
         "ChEMBL gives no unambiguous direction for the drug on the target (MODULATOR, mixed), and P11 forbids "
         "assuming one for M",
         "a real annotation gap, not a model failure. R is still computable for these rows: it carries no action sign",
         ["cohort", "pert_id", "cmap_name", "M_abstain_reason", "R", "layered_pct"])

    # --- cases where a score exists but should not be trusted -------------------------------------------------------
    sc = S[S.M.notna()].copy()
    cf = cov.set_index("unit_id")
    sc["conflicted_share"] = np.nan
    case("inferred_only_genes", sc[sc.M_n_landmark.fillna(0) == 0],
         "every gene behind the score is inferred by the L1000 model; none of it was physically measured",
         "the score may be real, but it cannot be separated from the imputation model. Reported, never promoted "
         "on its own",
         ["cohort", "pert_id", "cmap_name", "M", "M_pct", "M_n_genes", "M_n_landmark", "M_landmark"])

    amb = pd.read_csv(C.OUT / "a1/path_ambiguity.tsv", sep="\t")
    amb_hi = set(amb.unit_id[amb.gene_ambiguity_fraction >= 0.5])
    hu = H[H.unit_id.isin(amb_hi) & H.M.notna()]
    case("sign_contradiction", hu.merge(amb[["unit_id", "gene_ambiguity_fraction", "n_conflicted"]], on="unit_id"),
         "at least half of the genes the target reaches are excluded from the signed numerator because the routes "
         "to them disagree",
         "M is computed on the minority of genes whose sign survives; the coverage denominator records it, and a "
         "high value here means the prior is ambiguous about this target, not that the drug failed",
         ["cohort", "pert_id", "cmap_name", "unit_genes", "M", "M_n_genes", "gene_ambiguity_fraction",
          "n_conflicted"])

    if "split_prolif_driven" in S.columns:
        cyt = sc[(sc.split_prolif_driven.astype(str).str.lower().isin(["true", "1"]))
                 | (sc.get("split_prolif_effect_pct", pd.Series(np.nan, index=sc.index)) >= 90)]
        case("cytotoxic_signature", cyt[cyt.M > 0],
             "production already reads the drug's transcriptional response as proliferation-driven, and M is "
             "positive anyway",
             "the most dangerous confound for M: a growth-arrest response can agree with almost any route that "
             "runs through a proliferation hub. These rows must not be promoted on M alone",
             ["cohort", "pert_id", "cmap_name", "M", "M_pct", "hub_tf", "hub_fraction_structural",
              "split_prolif_effect_pct", "layered_pct"])

    case("line_mismatch", sc[(sc.source == "borrowed") | (sc.line_match.fillna(1) < 0.1)],
         "the drug has no profiled line inside the disease's tissue key, so its response is borrowed from other "
         "tissues (or the best matched line is a poor match)",
         "the mechanism may hold in the profiled line and not in the disease's tissue; `source` and `line_match` "
         "travel with every row so this is never invisible",
         ["cohort", "pert_id", "cmap_name", "source", "line_match", "transfer_conf", "M", "M_pct", "layered_pct"])

    poly = sc[(sc.n_units_M.fillna(0) > 1)]
    case("polypharmacology", poly[(poly.M_unit_min.fillna(0) < 0) & (poly.M_unit_max.fillna(0) > 0)],
         "the drug's own annotated target units disagree: at least one gives a positive M and another a negative "
         "one, and P13 averages them",
         "the averaged M hides a real disagreement about which target the response belongs to; the per-unit values "
         "are kept in `units` so the average is never the only number",
         ["cohort", "pert_id", "cmap_name", "n_units_M", "M", "M_unit_min", "M_unit_max", "units"])

    # --- absent target expression, from DepMap, per source line -----------------------------------------------------
    genes = sorted({g for s in cov.genes.dropna() for g in str(s).split("|")})
    ex = IO.depmap_matrix(C.DEPMAP / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", set(genes))
    rows = []
    for r in H[H.M.notna()].itertuples(index=False):
        gs = [g for g in str(r.unit_genes).split("|") if g in ex.columns]
        used = [m.split(":")[0] for m in str(r.lines).split("|") if m]
        if not gs or not used:
            continue
        v = ex.reindex(used)[gs]
        if v.isna().all().all():
            continue
        rows.append(dict(cohort=r.cohort, pert_id=r.pert_id, cmap_name=r.cmap_name, unit_genes=r.unit_genes,
                         M=r.M, n_lines=len(used),
                         max_target_tpm_log1p=round(float(np.nanmax(v.max(axis=1).values)), 3),
                         mean_target_tpm_log1p=round(float(np.nanmean(v.values)), 3)))
    EXP = pd.DataFrame(rows)
    case("absent_target_expression",
         EXP[EXP.max_target_tpm_log1p < TPM_FLOOR] if len(EXP) else EXP,
         f"the annotated target is not expressed (log2(TPM+1) < {TPM_FLOOR} in every source line the profile came "
         f"from), so the drug cannot have engaged it there",
         "a positive M on such a row is a false mechanism by construction; this is the cleanest internal negative "
         "control D128 has, and it needs no outcome label",
         ["cohort", "pert_id", "cmap_name", "unit_genes", "M", "max_target_tpm_log1p", "n_lines"])

    unk = ctl[(ctl.M_status == "scored") & (~ctl.BRAF_V600) & (~ctl.EGFR_sensitizing) & (~ctl.ERBB2_amp)
              & (~ctl.KRAS_NRAS)]
    case("unknown_genotype", unk,
         "a control row whose source line carries none of the audited alterations for its wave: neither the "
         "sensitising genotype nor the resistance one",
         "it belongs to neither the positive nor the negative set and is reported separately rather than counted "
         "as either",
         ["wave", "context", "cmap_name", "line", "M", "sign_shuffle_pct"])

    IDX = pd.DataFrame(index)
    IDX.to_csv(OUTDIR / "index.tsv", sep="\t", index=False)
    for name, (d, what, means, cols) in cases.items():
        cols = [c for c in cols if c in d.columns]
        txt = [f"# {name.replace('_', ' ')}", "",
               f"**What this case is.** {what}.", "",
               f"**Why it matters.** {means}.", "",
               f"**Rows:** {len(d)}" + (f" of {N} drug-level rows ({len(d)/N:.1%})" if N and name not in
                                        ("unknown_genotype", "sign_contradiction") else ""), ""]
        if len(d):
            srt = "M" if "M" in cols else cols[0]
            txt += ["**Largest cases** (by M where defined):", "",
                    md_table(d[cols].sort_values(srt, ascending=False, na_position="last"))]
        else:
            txt += ["_No row falls in this case._", ""]
        (OUTDIR / f"{name}.md").write_text("\n".join(txt))
        d[cols].to_csv(OUTDIR / f"{name}.tsv", sep="\t", index=False)
    print(IDX.to_string())


if __name__ == "__main__":
    main()
