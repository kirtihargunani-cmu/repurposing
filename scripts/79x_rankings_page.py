#!/usr/bin/env python
"""Browsable page of the full disease x drug rankings (user 2026-11-02; runs in `celligner`).

One self-contained HTML file (data embedded): 41 diseases and 14 subgroups; per cohort every FDA-approved drug
(LINCS compound matched to Drugs@FDA in 79p; one per ingredient, the best on the lead score, as 79p) that passed the
target gate (79b) with reversal (rev), specificity (spec), mechanism specificity (spec_mech, the lead score, D80) and,
for subgroups, sibling specificity (spec_sib, spec_sib_mech): raw values and within-cohort percentiles; data source,
transfer confidence and the flags (unmatched_lines, low_transfer, multi_target, prolif_driven, low_mrs); FDA status;
the receptor evidence gate (79w) per drug; and literature where a search was run: the LINCS-name search of each
disease's top 10 by spec_mech and by rev (79o) and the FDA-name search of each cohort's top 3 FDA drugs (79o MODE=fda).
Pairs never searched link to a PubMed search instead.
Inputs: results/lincs/drug_reversal_split_v2_consensus.tsv, stratum_reversal_split_v2_consensus.tsv,
        drug_reversal_summary_v2_consensus.tsv, stratum_reversal_summary_v2_consensus.tsv, drug_annotations.tsv,
        fda_matches.tsv, literature_v2_consensus.tsv, literature_fda_top_v2.tsv, receptor_gate_v2.tsv,
        results/step10_v2/strata.tsv; template scripts/templates/disease_drug_rankings.html (placeholder __DATA__)
Output: results/lincs/disease_drug_rankings_v2.html
"""
import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
FLAGS = ["unmatched_lines", "low_transfer", "multi_target", "prolif_driven", "low_mrs"]   # bit order in the page
GATE = {"pass": "p", "binds_but_absent_in_matched_lines": "a", "assembly_incomplete_in_tumors": "s"}


def clean(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else v


def rnd(v, k=2):
    v = clean(v)
    return None if v is None else round(float(v), k)


def truthy(v):
    return str(v) in ("True", "true", "1", "1.0")


def links(lr):
    pa = [[lr[f"paper{j}"], clean(lr[f"paper{j}_title"]) or ""] for j in (1, 2, 3) if isinstance(lr[f"paper{j}"], str) and lr[f"paper{j}"]]
    tr = [[lr[f"trial{j}"], clean(lr[f"trial{j}_info"]) or ""] for j in (1, 2, 3) if isinstance(lr[f"trial{j}"], str) and lr[f"trial{j}"]]
    pm, ct = clean(lr["pubmed_count"]), clean(lr["trial_count"])
    return [None if pm is None else int(pm), None if ct is None else int(ct), lr["pubmed_search_url"], lr["trials_search_url"], pa, tr]


def main():
    g = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    s = pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False).rename(columns={"stratum": "group_id"})
    g, s = g[g.passes_gate.map(truthy)], s[s.passes_gate.map(truthy)]
    # FDA-approved drugs only (user 2026-11-02): LINCS compounds matched to Drugs@FDA (79p; marketed or discontinued);
    # several LINCS IDs of one ingredient -> keep the best on the lead score, as 79p does
    fdm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    n_all = (len(g), len(s))
    g = (g[g.pert_id.isin(fdm.index)].assign(ing=lambda x: x.pert_id.map(fdm))
         .sort_values(["spec_mech", "spec", "rev"], ascending=False).drop_duplicates(["group_id", "ing"]))
    s = (s[s.pert_id.isin(fdm.index)].assign(ing=lambda x: x.pert_id.map(fdm))
         .sort_values(["spec_sib_mech", "spec_mech", "rev"], ascending=False).drop_duplicates(["group_id", "ing"]))
    print(f"FDA filter: diseases {n_all[0]} -> {len(g)} rows, subgroups {n_all[1]} -> {len(s)} rows")
    # flags not carried into the split tables come from the reversal tables
    for df, f in ((g, "drug_reversal_v2_consensus.tsv"), (s, "stratum_reversal_v2_consensus.tsv")):
        r = pd.read_csv(RES / f, sep="\t", low_memory=False, usecols=["group_id", "pert_id", "n_lines_used"])
        df["n_lines_used"] = pd.MultiIndex.from_frame(df[["group_id", "pert_id"]]).map(
            r.drop_duplicates(["group_id", "pert_id"]).set_index(["group_id", "pert_id"]).n_lines_used)
    gs = pd.read_csv(RES / "drug_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    ss = pd.read_csv(RES / "stratum_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    st = pd.read_csv(ROOT / "results/step10_v2/strata.tsv", sep="\t").set_index("stratum")
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id")
    fda = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id")
    lit = pd.read_csv(RES / "literature_v2_consensus.tsv", sep="\t")
    top = pd.read_csv(RES / "literature_fda_top_v2.tsv", sep="\t")
    rg = pd.read_csv(RES / "receptor_gate_v2.tsv", sep="\t", low_memory=False)

    # ---- drug table
    pids = sorted(set(g.pert_id) | set(s.pert_id))
    D, di = [], {}
    for p in pids:
        a = ann.loc[p] if p in ann.index else None
        tg = sorted({t for c in ("lincs_target", "hub_target") for t in str(a[c]).split("|") if a is not None and t not in ("nan", "")})
        f = fda.loc[p] if p in fda.index else None
        di[p] = len(D)
        D.append(dict(c=(clean(a.cmap_name) if a is not None else None) or p, p=p,
                      moa=(clean(a.hub_moa) or clean(a.lincs_moa) or "") if a is not None else "",
                      ph=(clean(a.hub_clinical_phase) or "") if a is not None else "",
                      tg=tg, ind=((clean(a.approved_indications) or "")[:200]) if a is not None else "",
                      fda=(f.fda_ingredient.title(), f.fda_status[0], str(clean(f.first_fda_approval) or "")[:4]) if f is not None else None))

    # ---- literature: (cohort, pert_id or name) -> link record; FDA-name search wins (most specific)
    L, li = [], {}
    for df, unit, by in ((lit, "group_id", "l"), (top, "unit", "f")):
        name2p = {}
        for p in pids:
            name2p.setdefault(D[di[p]]["c"], []).append(p)
        for r in df.itertuples(index=False):
            r = r._asdict()
            ps = [r["cmap_name"]] if r["cmap_name"] in di else name2p.get(r["cmap_name"], [])
            for p in ps:
                li[(r[unit], p)] = len(L)
            L.append(links(r) + [by])

    # ---- receptor gate per cohort x drug
    R = {}
    for r in rg.itertuples(index=False):
        R.setdefault((r.group_id, r.pert_id), []).append(
            [r.target, r.receptor_class, r.binding_evidence, rnd(r.pchembl_target), rnd(r.window_log10),
             rnd(r.frac_matched_lines_expr), int(r.n_matched_lines), GATE.get(r.receptor_gate, "n"), int(bool(r.single_line)),
             r.receptor_gate.split(":")[1] if r.receptor_gate.startswith("no_binding") else ""])

    def rows(df, sub):
        out = {}
        for u, d in df.groupby("group_id"):
            lst = []
            for r in d.itertuples(index=False):
                fl = sum(1 << i for i, f in enumerate(FLAGS) if truthy(getattr(r, f, False)))
                row = [di[r.pert_id], r.source[0], rnd(r.rev), rnd(r.rev_pct, 1), rnd(r.spec), rnd(r.spec_pct, 1),
                       rnd(r.spec_mech), rnd(r.spec_mech_pct, 1), rnd(r.prolif_share), rnd(r.transfer_conf),
                       rnd(r.line_match), fl, li.get((u, r.pert_id), -1), R.get((u, r.pert_id), []),
                       int(clean(r.n_lines_used) or 0)]
                if sub:
                    row += [rnd(r.spec_sib), rnd(r.spec_sib_pct, 1), rnd(r.spec_sib_mech), rnd(r.spec_sib_mech_pct, 1)]
                lst.append(row)
            out[u] = lst
        return out

    G, S = rows(g, False), rows(s, True)
    meta = {u: dict(kind="group", parent="", status=gs.loc[u, "evidence_status"], note=gs.loc[u, "evidence_note"],
                    bl=clean(gs.loc[u, "best_line"]) or "", desc="") for u in G}
    meta.update({u: dict(kind="subgroup", parent=st.loc[u, "group_id"], status=ss.loc[u, "evidence_status"],
                         note=ss.loc[u, "evidence_note"], bl=clean(ss.loc[u, "best_line"]) or "",
                         desc=clean(st.loc[u, "description"]) or "") for u in S})
    data = json.dumps(dict(drugs=D, lit=L, rows={**G, **S}, meta=meta, flags=FLAGS), separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    page = (ROOT / "scripts/templates/disease_drug_rankings.html").read_text().replace("__DATA__", data, 1)
    (RES / "disease_drug_rankings_v2.html").write_text(page)
    n = sum(map(len, G.values())) + sum(map(len, S.values()))
    nl = sum(1 for v in [*G.values(), *S.values()] for r in v if r[12] >= 0)
    print(f"{len(G)} diseases, {len(S)} subgroups, {n} rows ({nl} with literature), {len(D)} drugs; "
          f"{len(page) / 1e6:.2f} MB -> results/lincs/disease_drug_rankings_v2.html")


if __name__ == "__main__":
    main()
