#!/usr/bin/env python
"""Browsable page of the FDA-approved drug rankings (user 2026-10-29; runs in `celligner`).

One self-contained HTML file (data embedded, no server needed): cohort list (41 diseases, 14 subgroups) with evidence
status, and per cohort the ranked FDA drugs (79p) with scores, data source, FDA status, Drugs@FDA link, PubMed / trial
counts and links, and a PubMed reviews link per drug. Literature: the top 3 FDA drugs per cohort from the FDA-name
search (79o MODE=fda, literature_fda_top_v2.tsv); other rows from the LINCS-name search (79o, literature_v2_consensus.tsv)
where that pair was searched.
Inputs: results/lincs/fda_ranking_v2_consensus.tsv, fda_ranking_strata_v2_consensus.tsv, literature_fda_top_v2.tsv,
        literature_v2_consensus.tsv, drug_reversal_summary_v2_consensus.tsv, stratum_reversal_summary_v2_consensus.tsv,
        results/step10_v2/strata.tsv; template scripts/templates/fda_rankings.html (placeholder __DATA__)
Output: results/lincs/fda_rankings_v2.html (published copy: https://claude.ai/artifact/4UVF5F6SxD639VKKJEZenw)
"""
import json
import math
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"


def clean(v):
    return None if v is None or (isinstance(v, float) and math.isnan(v)) else v


def decode(r):
    """D122: the layered score's parts for the best target unit."""
    if clean(getattr(r, "layered_genes", None)) is None:
        return ""
    b = "" if clean(r.layered_r) is None else f"{r.layered_r:.2f}"
    return f"{r.layered_genes}: D {r1(r.layered_D)} A {r1(r.layered_A)} O {r1(r.layered_O)}, binding {b}"


def r1(v):
    v = clean(v)
    return None if v is None else round(float(v), 1)


def links(lr):
    pa = [[lr[f"paper{j}"], clean(lr[f"paper{j}_title"]) or ""] for j in (1, 2, 3) if isinstance(lr[f"paper{j}"], str) and lr[f"paper{j}"]]
    tr = [[lr[f"trial{j}"], clean(lr[f"trial{j}_info"]) or ""] for j in (1, 2, 3) if isinstance(lr[f"trial{j}"], str) and lr[f"trial{j}"]]
    pm, ct = clean(lr.pubmed_count), clean(lr.trial_count)
    return dict(pm=None if pm is None else int(pm), ct=None if ct is None else int(ct), pmu=lr.pubmed_search_url,
                ctu=lr.trials_search_url, pa=pa, tr=tr)


def main():
    g = pd.read_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t")
    s = pd.read_csv(RES / "fda_ranking_strata_v2_consensus.tsv", sep="\t")
    lit = pd.read_csv(RES / "literature_v2_consensus.tsv", sep="\t").set_index(["group_id", "cmap_name"])
    top = pd.read_csv(RES / "literature_fda_top_v2.tsv", sep="\t").set_index(["unit", "drug"])
    gs = pd.read_csv(RES / "drug_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    ss = pd.read_csv(RES / "stratum_reversal_summary_v2_consensus.tsv", sep="\t").set_index("group_id")
    st = pd.read_csv(ROOT / "results/step10_v2/strata.tsv", sep="\t").set_index("stratum")

    def rows(df, unit, scorecol):
        out = {}
        for u, d in df.groupby(unit):
            lst = []
            for r in d.itertuples(index=False):
                if (u, r.drug) in top.index:
                    L = links(top.loc[(u, r.drug)]); L["by"] = "fda"
                elif unit == "group_id" and (u, r.cmap_name) in lit.index:
                    L = links(lit.loc[(u, r.cmap_name)]); L["by"] = "lincs"
                else:
                    L = dict(pm=None, ct=None, pmu="", ctu="", pa=[], tr=[], by="")
                ind = clean(r.approved_indications) or ""
                appr = clean(r.first_fda_approval)
                lst.append(dict(k=int(r.fda_rank), d=r.drug, c=r.cmap_name, st=r.fda_status[0], y=str(appr)[:4] if appr else "",
                                moa=clean(r.hub_moa) or "", ind=ind[:160], src=r.source[0], tc=r1(r.transfer_conf),
                                s1=r1(r.layered_pct), sp=r1(getattr(r, scorecol)), rv=r1(r.rev_pct),
                                dec=decode(r),
                                fl=", ".join(x.replace("_", " ") for x in [*(str(r.lincs_flags).split("|") if clean(r.lincs_flags) else []),
                                     *(["minority alteration"] if str(r.minority_alteration) == "True" else []),
                                     *(["essential target"] if str(r.essential_target) == "True" else [])]),
                                fda=clean(r.drugsatfda_url) or "", **L))
            out[u] = lst
        return out

    G, S = rows(g, "group_id", "spec_mech_pct"), rows(s, "stratum", "spec_sib_mech_pct")   # D122: lead = layered_pct (s1); sp = the former lead
    meta = {u: dict(kind="group", parent="", status=gs.loc[u, "evidence_status"], note=gs.loc[u, "evidence_note"],
                    bl=clean(gs.loc[u, "best_line"]) or "", desc="") for u in G}
    meta.update({u: dict(kind="subgroup", parent=st.loc[u, "group_id"], status=ss.loc[u, "evidence_status"],
                         note=ss.loc[u, "evidence_note"], bl=clean(ss.loc[u, "best_line"]) or "", desc=st.loc[u, "description"])
                 for u in S})
    data = json.dumps(dict(rows={**G, **S}, meta=meta), separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")
    page = (ROOT / "scripts/templates/fda_rankings.html").read_text().replace("__DATA__", data, 1)
    (RES / "fda_rankings_v2.html").write_text(page)
    print(f"{len(G)} diseases, {len(S)} subgroups, {sum(map(len, G.values())) + sum(map(len, S.values()))} rows; "
          f"{len(page) / 1e6:.2f} MB -> results/lincs/fda_rankings_v2.html")


if __name__ == "__main__":
    main()
