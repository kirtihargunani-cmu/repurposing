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
        fda_matches.tsv, literature_v2_consensus.tsv, literature_fda_top_v2.tsv, literature_fda_all_v2.tsv (79o
        MODE=fda_all, every row; preferred where present), receptor_gate_v2.tsv, standin_pairs_v2.tsv + standin_run.json
        (79y: up to 3 stand-ins for borrowed / weak-evidence FDA drugs; where one exists the page also shows the
        average of the drug's own lead percentile and the best stand-in's, and can rank by it: user 2026-11-02),
        results/step10_v2/strata.tsv, results/lincs/d100/timing_scores.tsv (D100 timing class, optional);
        template scripts/templates/disease_drug_rankings.html (placeholder __DATA__)
Output: results/lincs/disease_drug_rankings_v2.html (published copy with the D97 networks and D100 timing column:
        https://claude.ai/artifact/DkoreTiuHXRkqaEYMG5hfk ; earlier copy https://claude.ai/artifact/LnfCejaAaBThJGpgM91naA)
"""
import json
import os
import math
from pathlib import Path

import sys

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from fda_rows import drug_name, fda_rows  # noqa: E402
RES = ROOT / "results/lincs"
PRESENT_ONLY = os.environ.get("PRESENT_ONLY", "0") == "1"
OUT_HTML = "disease_drug_rankings_v2_target_present.html" if PRESENT_ONLY else "disease_drug_rankings_v2.html"
FLAGS = ["unmatched_lines", "low_transfer", "multi_target", "prolif_driven", "low_mrs", "target_absent", "minority_target", "receptor_mismatch", "tumor_environment", "unverified_binding"]   # bit order in the page
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


def standin_vs_borrowed():
    """79y2 comparison (stand-in vs the drug's own borrowed score), if run"""
    f = RES / "standin_vs_borrowed_summary_v2.tsv"
    if not f.exists():
        return {}
    c = pd.read_csv(f, sep="\t").set_index("estimate")
    return dict(b_si=float(c.loc["standin", "spearman"]), b_bo=float(c.loc["borrowed_pct", "spearman"]),
                b_av=float(c.loc["average", "spearman"]), b_n=int(c.loc["standin", "n"]))


def main():
    g = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    s = pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False).rename(columns={"stratum": "group_id"})
    # FDA-approved drugs only (user 2026-11-02): scripts/fda_rows.py (Drugs@FDA match, one LINCS id per ingredient)
    n_all = (int(g.passes_gate.map(truthy).sum()), int(s.passes_gate.map(truthy).sum()))
    g, s = fda_rows()
    for df in (g, s):   # D102c (79i2): target absent in this cohort's typical tumor; minority = >= 3 tumors carry it
        df["minority_target"] = df.target_absent.map(truthy) & (df.get("target_carriers", 0).fillna(0) >= 3)
        rs = df.get("receptor_status", pd.Series("", index=df.index)).fillna("")   # D104c (79i3)
        df["receptor_mismatch"], df["tumor_environment"] = rs.eq("receptor_mismatch"), rs.eq("acts_on_tumor_environment")
    print(f"FDA filter: diseases {n_all[0]} -> {len(g)} rows, subgroups {n_all[1]} -> {len(s)} rows")
    if PRESENT_ONLY:   # user 2026-11-04: only drugs whose main target (D105b) is present in the cancer's tumors
        keep = lambda df: df[df.get("evidence_targets", pd.Series("", index=df.index)).fillna("").ne("") & ~df.target_absent.map(truthy)
                             & ~df.receptor_mismatch & df.get("target_source", pd.Series("", index=df.index)).isin(["chembl_mechanism", "confirmed_annotation", "annotation"])]
        n0 = (len(g), len(s)); g, s = keep(g), keep(s)
        # user 2026-11-04: presence must rest on data: >= 1 main target with a tumor value that passes the rule
        # (receptors: median >= 1 TPM and >= 10% of the typical high R; other targets: median >= 1 TPM)
        import importlib.util
        sp_ = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w79 = importlib.util.module_from_spec(sp_); sp_.loader.exec_module(w79)
        REC = set(w79.classes(w79.hgnc())[0])
        def by_data(ev):
            for item in str(ev).split(";"):
                gene, _, v = item.partition(":"); m, _, R = v.partition("/")
                try:
                    m, R = float(m), float(R)
                except ValueError:
                    continue
                if m == m and m >= 1 and (gene not in REC or (R == R and m >= 0.1 * R)):
                    return True
            return False
        n1 = (len(g), len(s)); g, s = g[g.evidence_targets.map(by_data)], s[s.evidence_targets.map(by_data)]
        print(f"target present only: diseases {n0[0]} -> {n1[0]} -> {len(g)} rows (presence backed by data), subgroups {n0[1]} -> {n1[1]} -> {len(s)}")
        # user 2026-11-04: percentiles re-ranked among the FDA drugs on this page, per cohort / subgroup
        for df in (g, s):
            for c in ["rev", "spec", "spec_mech", "spec_sib", "spec_sib_mech"]:
                if c in df:
                    df[c + "_pct"] = (df.groupby("group_id")[c].rank(method="max", pct=True) * 100).round(1)
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
    fa = RES / "literature_fda_all_v2.tsv"
    lall = pd.read_csv(fa, sep="\t") if fa.exists() else None
    sp = pd.read_csv(RES / "standin_pairs_v2.tsv", sep="\t")
    sp = sp[sp.fda_needs_standin & (sp.standin_rank <= 3)].sort_values("standin_rank")
    SI = {}
    for r in sp.itertuples(index=False):
        SI.setdefault((r.group_id, r.fda_pert), []).append(
            [r.cmp_name, rnd(r.similarity), int(r.n_shared_lines), r.shared_targets, rnd(r.cmp_spec_mech_pct, 1),
             rnd(r.cmp_spec_pct, 1), rnd(r.cmp_rev_pct, 1), rnd(r.cmp_spec_sib_mech_pct, 1), int(bool(r.cmp_is_fda))])
    vrun = json.load(open(RES / "standin_run.json"))
    rg = pd.read_csv(RES / "receptor_gate_v2.tsv", sep="\t", low_memory=False)
    # D100 (user 2026-11-03): early (6 h) / late (24 h) scores and timing class per disease x drug, where available
    tf = RES / "d100/timing_scores.tsv"
    TC = {"both": "b", "early only": "e", "late only": "l", "neither": "n"}
    TM = {}
    if tf.exists():
        for t in pd.read_csv(tf, sep="\t").itertuples(index=False):
            if t.timing_class in TC:
                TM[(t.group_id, t.pert_id)] = [TC[t.timing_class], rnd(t.spec_mech_early_pct, 1), rnd(t.spec_mech_late_pct, 1)]

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
                      fda=(drug_name(f.fda_ingredient), f.fda_status[0], str(clean(f.first_fda_approval) or "")[:4]) if f is not None else None))

    # ---- literature: (cohort, pert_id or name) -> link record; FDA-name search wins (most specific)
    L, li = [], {}
    for df, unit, by in ((lit, "group_id", "l"), (top, "unit", "f")) + (((lall, "unit", "a"),) if lall is not None else ()):
        name2p = {}
        for p in pids:
            name2p.setdefault(D[di[p]]["c"], []).append(p)
        for r in df.itertuples(index=False):
            r = r._asdict()
            ps = ([r["pert_id"]] if r.get("pert_id") in di else []) or ([r["cmap_name"]] if r["cmap_name"] in di else name2p.get(r["cmap_name"], []))
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
                row += ([rnd(r.spec_sib), rnd(r.spec_sib_pct, 1), rnd(r.spec_sib_mech), rnd(r.spec_sib_mech_pct, 1)]
                        if sub else [None] * 4)
                si = SI.get((u, r.pert_id), [])
                row.append(si)
                # user 2026-11-02 (C27): where a stand-in exists, show the average of the drug's own lead percentile and
                # the best stand-in's (79y2: average 0.68 vs own borrowed 0.56 Spearman with the matched-line truth)
                own = rnd(r.spec_sib_mech_pct, 1) if sub else rnd(r.spec_mech_pct, 1)
                sip = (si[0][7] if sub else si[0][4]) if si else None
                avg = round((own + sip) / 2, 1) if si and own is not None and sip is not None and not PRESENT_ONLY else None   # stand-in percentiles are on the all-compound scale
                row += [avg, avg if avg is not None else own]
                row += TM.get((u, r.pert_id), [None, None, None])   # 22-24: timing class, early / late percentile (D100)
                row += [clean(getattr(r, "target_source", "")) or "", clean(getattr(r, "evidence_targets", "")) or "",
                        clean(getattr(r, "secondary_targets", "")) or ""]   # 25-27: main-target source, main targets, secondary (D105b)
                lst.append(row)
            out[u] = lst
        return out

    G, S = rows(g, False), rows(s, True)
    meta = {u: dict(kind="group", parent="", status=gs.loc[u, "evidence_status"], note=gs.loc[u, "evidence_note"],
                    bl=clean(gs.loc[u, "best_line"]) or "", desc="") for u in G}
    meta.update({u: dict(kind="subgroup", parent=st.loc[u, "group_id"], status=ss.loc[u, "evidence_status"],
                         note=ss.loc[u, "evidence_note"], bl=clean(ss.loc[u, "best_line"]) or "",
                         desc=clean(st.loc[u, "description"]) or "") for u in S})
    val = {v["set"]: v for v in vrun["validation"]}
    data = json.dumps(dict(drugs=D, lit=L, rows={**G, **S}, meta=meta, flags=FLAGS,
                           standin=dict(adopted=vrun["adopted"], rho=val["accepted stand-ins"]["spearman"],
                                        lo=val["accepted stand-ins"]["ci_lo"], hi=val["accepted stand-ins"]["ci_hi"],
                                        rand=val["random shared-target pairs (median of 200)"]["spearman"],
                                        n=val["accepted stand-ins"]["n"], r=vrun["config"]["min_similarity"],
                                        k=vrun["config"]["min_shared_lines"],
                                        **standin_vs_borrowed())), separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    page = (ROOT / "scripts/templates/disease_drug_rankings.html").read_text().replace("__DATA__", data, 1)
    if PRESENT_ONLY:
        page = page.replace("<title>Disease Drug Rankings</title>", "<title>Target-Present Drug Rankings</title>", 1).replace(
            "<h1>Disease Drug Rankings</h1>", "<h1>Target-Present Drug Rankings</h1>\n    <p class=\"sub\">Only FDA drugs with at least one main target (ChEMBL curated mechanism, binding-confirmed or annotated target) that tumor data show is made by this cancer's typical tumor, and, for receptor drugs, by the cell lines the score came from (or by the tumor's vessels / immune cells). Drugs acting on DNA or pathogens and drugs without target data are left out. Percentiles are among the FDA drugs on this page, per cohort. Hidden by default (untick 'Hide weak-evidence drugs'): unmatched or low-transfer lines and unverified binding.</p>", 1).replace(
            'const HIDE = ["unmatched_lines", "low_transfer", "target_absent", "receptor_mismatch"]',
            'const HIDE = ["unmatched_lines", "low_transfer", "target_absent", "receptor_mismatch", "unverified_binding"]', 1)
        assert '"receptor_mismatch", "unverified_binding"]' in page
    (RES / OUT_HTML).write_text(page)
    n = sum(map(len, G.values())) + sum(map(len, S.values()))
    nl = sum(1 for v in [*G.values(), *S.values()] for r in v if r[12] >= 0)
    ns = sum(1 for v in [*G.values(), *S.values()] for r in v if r[19])
    print(f"rows with a stand-in: {ns}")
    print(f"{len(G)} diseases, {len(S)} subgroups, {n} rows ({nl} with literature), {len(D)} drugs; "
          f"{len(page) / 1e6:.2f} MB -> results/lincs/{OUT_HTML}")


if __name__ == "__main__":
    main()
