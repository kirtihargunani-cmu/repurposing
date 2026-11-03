#!/usr/bin/env python
"""Browsable page: LINCS compounds tested in each cohort's own lines and their closest FDA drugs (user 2026-11-02).

Per disease (41) and subgroup (14): compounds that pass the target gate and were tested in the cohort's matched lines
(direct, no unmatched_lines flag), ranked by the cohort's lead score (diseases spec_mech, subgroups spec_sib_mech),
each with up to 3 closest FDA drugs from 79y (shared annotated target, similar regulator-activity change over shared
profiled lines). For each FDA drug: similarity, shared lines and targets, the drug's own score and data source in the
cohort, whether it needs a stand-in (borrowed / weak evidence), and its literature counts (79o MODE=fda_all).
Cohorts without matched compounds say so. The 79y validation is shown on the page.
Inputs: results/lincs/standin_pairs_v2.tsv, standin_run.json, drug_reversal_split_v2_consensus.tsv,
        stratum_reversal_split_v2_consensus.tsv, drug_annotations.tsv, drug_reversal_summary_v2_consensus.tsv,
        stratum_reversal_summary_v2_consensus.tsv, literature_fda_all_v2.tsv (if present), standin_vs_borrowed_summary_v2.tsv
        (79y2, if present), results/step10_v2/strata.tsv;
        template scripts/templates/standin_compounds.html (placeholder __DATA__)
Output: results/lincs/standin_compounds_v2.html
"""
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import truthy  # noqa: E402


def clean(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else v


def rnd(v, k=1):
    v = clean(v)
    return None if v is None else round(float(v), k)


def main():
    sp = pd.read_csv(RES / "standin_pairs_v2.tsv", sep="\t")
    run = json.load(open(RES / "standin_run.json"))
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id")
    fm = set(pd.read_csv(RES / "fda_matches.tsv", sep="\t").pert_id)
    fa = RES / "literature_fda_all_v2.tsv"
    lit = {}
    if fa.exists():
        for r in pd.read_csv(fa, sep="\t").itertuples(index=False):
            lit[(r.unit, r.pert_id)] = [None if pd.isna(r.pubmed_count) else int(r.pubmed_count),
                                        None if pd.isna(r.trial_count) else int(r.trial_count),
                                        r.pubmed_search_url, r.trials_search_url]
    full = {"group": pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False),
            "subgroup": pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
            .rename(columns={"stratum": "group_id"})}
    gs = pd.read_csv(RES / "drug_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    ss = pd.read_csv(RES / "stratum_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    st = pd.read_csv(ROOT / "results/step10_v2/strata.tsv", sep="\t").set_index("stratum")
    close = sp[sp.closest_rank <= 3].sort_values("closest_rank")
    C = {}
    for r in close.itertuples(index=False):
        own = r.fda_spec_sib_mech_pct if r.level == "subgroup" else r.fda_spec_mech_pct
        C.setdefault((r.group_id, r.cmp_pert), []).append(
            [r.fda_drug, rnd(r.similarity, 2), int(r.n_shared_lines), r.shared_targets.replace("|", ", "), rnd(own),
             r.fda_source[0], int(bool(r.fda_unmatched)), int(bool(r.fda_low_transfer)), int(bool(r.fda_needs_standin)),
             lit.get((r.group_id, r.fda_pert))])
    rows, meta = {}, {}
    for kind, df in full.items():
        df = df[df.passes_gate.map(truthy) & (df.source == "direct") & ~df.unmatched_lines.map(truthy)].copy()
        lead = "spec_sib_mech" if kind == "subgroup" else "spec_mech"
        units = sorted(gs.index if kind == "group" else ss.index)
        for u in units:
            x = df[df.group_id == u].sort_values(lead, ascending=False).reset_index(drop=True)
            out = []
            for k, r in enumerate(x.itertuples(index=False)):
                cl = C.get((u, r.pert_id))
                if not cl:
                    continue
                a = ann.loc[r.pert_id] if r.pert_id in ann.index else None
                tg = sorted({t for c in ("lincs_target", "hub_target") for t in str(a[c]).split("|") if t not in ("nan", "")}) if a is not None else []
                out.append([k + 1, r.cmap_name, (clean(a.hub_moa) or clean(a.lincs_moa) or "") if a is not None else "",
                            ", ".join(tg[:8]) + (" …" if len(tg) > 8 else ""), int(r.pert_id in fm),
                            rnd(getattr(r, "spec_sib_mech_pct", None)) if kind == "subgroup" else None,
                            rnd(r.spec_mech_pct), rnd(r.spec_pct), rnd(r.rev_pct), rnd(r.line_match, 2), cl])
            rows[u] = out
            sm = gs if kind == "group" else ss
            meta[u] = dict(kind=kind, parent=str(st.loc[u, "group_id"]) if kind == "subgroup" else "",
                           status=sm.loc[u, "evidence_status"], note=sm.loc[u, "evidence_note"],
                           bl=clean(sm.loc[u, "best_line"]) or "", n_matched=len(x),
                           desc=(clean(st.loc[u, "description"]) or "") if kind == "subgroup" else "")
    val = {v["set"]: v for v in run["validation"]}
    data = dict(rows=rows, meta=meta, v=dict(rho=val["accepted stand-ins"]["spearman"], lo=val["accepted stand-ins"]["ci_lo"],
                hi=val["accepted stand-ins"]["ci_hi"], rand=val["random shared-target pairs (median of 200)"]["spearman"],
                n=val["accepted stand-ins"]["n"], r=run["config"]["min_similarity"], k=run["config"]["min_shared_lines"],
                ge3=val["shared lines >= 3"]["spearman"], two=val["shared lines = 2"]["spearman"]))
    cmpf = RES / "standin_vs_borrowed_summary_v2.tsv"
    if cmpf.exists():   # 79y2: stand-in vs the drug's own borrowed score, same drugs
        c = pd.read_csv(cmpf, sep="\t").set_index("estimate")
        data["v"].update(b_n=int(c.loc["standin", "n"]), b_si=float(c.loc["standin", "spearman"]),
                         b_bo=float(c.loc["borrowed_pct", "spearman"]), b_av=float(c.loc["average", "spearman"]),
                         b_lo=float(c.loc["standin", "diff_ci_lo"]), b_hi=float(c.loc["standin", "diff_ci_hi"]))
    s = json.dumps(data, separators=(",", ":"), ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    page = (ROOT / "scripts/templates/standin_compounds.html").read_text().replace("__DATA__", s, 1)
    (RES / "standin_compounds_v2.html").write_text(page)
    n = sum(map(len, rows.values()))
    print(f"{len(rows)} cohorts ({sum(1 for v in rows.values() if v)} with compounds), {n} compound rows; "
          f"{len(page) / 1e6:.2f} MB -> results/lincs/standin_compounds_v2.html")


if __name__ == "__main__":
    main()
