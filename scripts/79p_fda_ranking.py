#!/usr/bin/env python
"""Ranking of FDA-approved drugs per disease and subgroup (user 2026-10-29; runs in `celligner`).

FDA approval from Drugs@FDA (data/raw/drugsatfda/extracted, FDA download https://www.fda.gov/media/89850/download,
fetched 2026-10-29; small molecules and CDER biologics). The `approved` column of 79e means approved anywhere
(ChEMBL max phase 4 or Repurposing Hub "Launched"); this file is FDA only.
Matching: a LINCS drug's names (cmap_name, ChEMBL preferred name, Repurposing Hub name) against each Drugs@FDA active
ingredient (combination products split on ";"), both normalized: upper case, salt / hydrate words removed (SALTS),
non-alphanumerics removed. fda_status = "marketed" if any product of the ingredient is prescription or OTC now, else
"discontinued" (only discontinued products; still FDA-approved historically). fda_routes = all routes of the
ingredient's products; local_route = none of them systemic (D107). first_fda_approval = earliest original
(ORIG) approval of an NDA / BLA containing the ingredient.
Ranking (D80: specificity first): per disease, FDA-approved drugs that pass the target gate, are not flagged
(unmatched_lines / low_transfer / local_route, D107: no systemic formulation, scripts/route_flag.py; D109: dose_inconsistent / inactive_matched_line / low_exposure, scripts/evidence_flags.py), not demoted by D102c / D104c (target_absent, receptor_mismatch) and have rev > 0, ordered by spec_mech (then spec, rev as context); subgroups by
spec_sib_mech (vs the sibling), with spec_mech. A drug listed under several LINCS ids keeps its best row. Literature /
trial counts are joined from 79o where that pair was searched.
Outputs: results/lincs/fda_ranking_v2_consensus.tsv (41 groups), fda_ranking_strata_v2_consensus.tsv (14 subgroups),
         fda_matches.tsv (LINCS drug -> FDA ingredient, status, first approval, application numbers)
"""
import re
from pathlib import Path

import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from route_flag import apply_route, fda_routes, is_local  # noqa: E402
from evidence_flags import apply_evidence  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
FDA = ROOT / "data/raw/drugsatfda/extracted"
from fda_rows import SALTS, FDA_ALIAS, FDA_EXCLUDE, drug_name  # noqa: E402  (D112: one list for 79p, fda_rows, 79o2)


def norm(s):
    if not isinstance(s, str) or not s.strip():
        return ""
    w = [x for x in re.split(r"[\s\-,()]+", s.upper()) if x and x not in SALTS]
    return re.sub(r"[^A-Z0-9]", "", "".join(w))


def fda_table():
    prod = pd.read_csv(FDA / "Products.txt", sep="\t", dtype=str, encoding="latin-1", on_bad_lines="skip")
    ms = pd.read_csv(FDA / "MarketingStatus.txt", sep="\t", dtype=str)
    app = pd.read_csv(FDA / "Applications.txt", sep="\t", dtype=str, encoding="latin-1")
    sub = pd.read_csv(FDA / "Submissions.txt", sep="\t", dtype=str, encoding="latin-1")
    p = prod.merge(ms, on=["ApplNo", "ProductNo"], how="left").merge(app[["ApplNo", "ApplType"]], on="ApplNo", how="left")
    p["ingredient"] = p.ActiveIngredient.str.split(";")
    p = p.explode("ingredient")
    p["key"] = p.ingredient.map(norm)
    p = p[p.key != ""]
    orig = sub[(sub.SubmissionType == "ORIG") & (sub.SubmissionStatus == "AP")][["ApplNo", "SubmissionStatusDate"]]
    orig = orig.groupby("ApplNo").SubmissionStatusDate.min()
    p["orig_date"] = p.ApplNo.map(orig)
    nb = p[p.ApplType.isin(["NDA", "BLA"])]
    t = p.groupby("key").agg(fda_ingredient=("ingredient", lambda v: v.str.strip().value_counts().index[0]),
                             marketed=("MarketingStatusID", lambda v: bool(v.isin(["1", "2"]).any())))
    t["first_fda_approval"] = nb.dropna(subset=["orig_date"]).groupby("key").orig_date.min().reindex(t.index).str[:10]
    t["fda_appl_nda_bla"] = nb.groupby("key").ApplNo.apply(lambda v: ";".join(sorted(set(v))[:5])).reindex(t.index)
    t["fda_status"] = t.marketed.map({True: "marketed", False: "discontinued"})
    rt = fda_routes(norm).reindex(t.index)   # D107: routes over all products of the ingredient
    t["fda_routes"] = rt.map(lambda r: "|".join(sorted(r)) if isinstance(r, set) else "")
    t["local_route"] = rt.map(lambda r: is_local(r) if isinstance(r, set) else None)
    return t.drop(columns="marketed")


def main():
    t = fda_table()
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id")
    rows = []
    for r in ann.itertuples(index=False):
        if r.cmap_name in FDA_EXCLUDE:   # D112
            continue
        for col in ("fda_alias", "cmap_name", "chembl_pref_name", "hub_name"):
            k = norm(FDA_ALIAS.get(r.cmap_name)) if col == "fda_alias" else norm(getattr(r, col))
            if k in t.index:
                rows.append(dict(pert_id=r.pert_id, cmap_name=r.cmap_name, matched_on=col, **t.loc[k].to_dict()))
                break
    m = pd.DataFrame(rows)
    m.to_csv(RES / "fda_matches.tsv", sep="\t", index=False)
    print(f"LINCS drugs: {ann.pert_id.nunique()}; matched to Drugs@FDA: {len(m)} ({(m.fda_status == 'marketed').sum()} marketed); "
          f"79e 'approved' (anywhere) but not in Drugs@FDA: "
          f"{int(((ann.chembl_max_phase == 4) | (ann.hub_clinical_phase == 'Launched')).sum() - ann.pert_id.isin(m.pert_id[m.pert_id.isin(ann.pert_id[(ann.chembl_max_phase == 4) | (ann.hub_clinical_phase == 'Launched')])]).sum())}")
    keep_ann = ann.set_index("pert_id")[["approved_indications", "hub_indication"]]
    lit = RES / "literature_v2_consensus.tsv"
    lit = pd.read_csv(lit, sep="\t") if lit.exists() else None
    summ = pd.read_csv(RES / "drug_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")

    def rank(df, unit, order, extra):
        ok = df[df.passes_gate & ~df.unmatched_lines & ~df.low_transfer & (df.rev > 0) & (~df.get("target_absent", pd.Series(False, index=df.index)).astype(str).str.lower().eq("true")) & ~df.get("receptor_status", pd.Series("", index=df.index)).eq("receptor_mismatch")]   # D102c / D104c demotions excluded
        ok = apply_route(ok, unit)
        ok = ok[~ok.local_route]   # D107: no systemic formulation
        ok = apply_evidence(ok, unit)
        ok = ok[~ok.weak_d109]   # D109: dose_inconsistent / inactive_matched_line / low_exposure
        x = ok.merge(m.drop(columns=["cmap_name", "local_route"]), on="pert_id").join(keep_ann, on="pert_id")
        x = x.sort_values(order, ascending=False).drop_duplicates([unit, "fda_ingredient"])
        x["fda_rank"] = x.groupby(unit).cumcount() + 1
        x["drug"] = x.fda_ingredient.map(drug_name)
        x["n_fda_eligible"] = x.groupby(unit)[unit].transform("size")
        x["drugsatfda_url"] = x.fda_appl_nda_bla.fillna("").str.split(";").str[0].map(
            lambda a: f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo={a}" if a else "")
        cols = [unit] + extra + ["fda_rank", "n_fda_eligible", "drug", "cmap_name", "fda_ingredient", "fda_status", "first_fda_approval",
                "hub_moa", "approved_indications", "source", "transfer_conf"] + order + ["rev", "rev_pct", "spec_mech_pct",
                "drugsatfda_url"]
        return x[[c for c in dict.fromkeys(cols) if c in x]].sort_values([unit, "fda_rank"])

    g = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    gr = rank(g, "group_id", ["spec_mech", "spec"], [])
    gr.insert(1, "evidence_status", gr.group_id.map(summ.evidence_status))
    if lit is not None:
        gr = gr.merge(lit[["group_id", "cmap_name", "pubmed_count", "pubmed_search_url", "trial_count", "trials_search_url"]],
                      on=["group_id", "cmap_name"], how="left")
    gr.to_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t", index=False)
    s = pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    sr = rank(s, "stratum", ["spec_sib_mech", "spec_mech"], ["parent_group"])
    sr.to_csv(RES / "fda_ranking_strata_v2_consensus.tsv", sep="\t", index=False)
    print(gr.groupby("group_id").fda_rank.size().describe().round(0).to_dict())
    for gid, d in gr.groupby("group_id"):
        print(f"{gid:<20} {d.evidence_status.iloc[0]:<14} " + ", ".join(f"{a} ({b})" for a, b in zip(d.drug.head(6), d.fda_status.head(6).str[0])))


if __name__ == "__main__":
    main()
