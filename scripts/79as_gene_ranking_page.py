#!/usr/bin/env python
"""Data for the production FDA ranking page (D211): filter on D/A/O, rank by gene-level reversal per timepoint.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-23: "create a new html of FDA rankings with this as production. show the breakdown of the D/A/O
scores, timepoint, etc. Show the target annotations as well."

One record per (unit, FDA ingredient), carrying BOTH timepoints, the D/A/O breakdown of the scored target unit,
the raw per-gene evidence behind it (79t gene_components: dependency z, alteration share, overactivity and its
source; groups only), and the target annotations from three independent sources:
  ChEMBL   evidence_targets (D105b: curated mechanism or binding; value = tumour expression / threshold)
  Hub      Drug Repurposing Hub target list and MOA
  LINCS    LINCS target annotation
plus FDA status, approved indications, the known-use label for this unit (79aa label_phase) and the LINCS flags.

Inputs: gene_ranking[_chemo]<SFX>.tsv (79ar), drug/stratum split tables, drug_annotations.tsv,
        fda_ranking[_strata]<SFX>.tsv, labels/pair_labels<SFX>.tsv.gz, target_role/gene_components.tsv.gz,
        step10_v2/strata.tsv, gene_ranking_run<SFX>.json
Output: results/lincs/gene_rankings_page/{index.html (copied template), data.json}
"""
import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
RES = ROOT / "results/lincs"
SFX = os.environ.get("SFX", "_v2_tf")
if "consensus" in SFX:
    raise SystemExit("D208: the consensus arm has been removed")
OUT = RES / "gene_rankings_page"
TEMPLATE = ROOT / "scripts/templates/gene_rankings.html"


def r(x, d=1):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) or pd.isna(x) else round(float(x), d)


def s(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) or pd.isna(x) or str(x) in ("", "nan") else str(x)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    G = pd.read_csv(RES / f"gene_ranking{SFX}.tsv", sep="\t", low_memory=False).assign(table="main")
    C = pd.read_csv(RES / f"gene_ranking_chemo{SFX}.tsv", sep="\t", low_memory=False).assign(table="chemo")
    X = pd.concat([G, C], ignore_index=True)
    X["ing"] = X.fda_ingredient.fillna(X.cmap_name)

    # target annotations
    sp = []
    for f, u in ((f"drug_reversal_split{SFX}.tsv", "group_id"), (f"stratum_reversal_split{SFX}.tsv", "stratum")):
        d = pd.read_csv(RES / f, sep="\t", low_memory=False,
                        usecols=lambda c: c in (u, "pert_id", "evidence_targets", "target_source", "secondary_targets",
                                                "hub_moa", "target_carrier_share", "unverified_binding", "target_absent"))
        sp.append(d.rename(columns={u: "unit"}))
    SP = pd.concat(sp, ignore_index=True).drop_duplicates(["unit", "pert_id"]).set_index(["unit", "pert_id"])
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    fr = []
    for f, u in ((f"fda_ranking{SFX}.tsv", "group_id"), (f"fda_ranking_strata{SFX}.tsv", "stratum")):
        d = pd.read_csv(RES / f, sep="\t", low_memory=False)
        fr.append(d.rename(columns={u: "unit"}))
    FR = pd.concat(fr, ignore_index=True).drop_duplicates(["unit", "cmap_name"]).set_index(["unit", "cmap_name"])
    lab = pd.read_csv(RES / f"labels/pair_labels{SFX}.tsv.gz", sep="\t", low_memory=False,
                      usecols=["unit_id", "pert_id", "label_phase", "label_source"])
    LAB = lab.drop_duplicates(["unit_id", "pert_id"]).set_index(["unit_id", "pert_id"])
    gc = pd.read_csv(RES / "target_role/gene_components.tsv.gz", sep="\t", low_memory=False).set_index(["group_id", "gene"])
    strata = pd.read_csv(ROOT / "results/step10_v2/strata.tsv", sep="\t").drop_duplicates("stratum").set_index("stratum")

    recs = []
    for (unit, level, ing), g in X.groupby(["unit", "level", "ing"], sort=False):
        b = g.iloc[0]
        win = {}
        for _, w in g.iterrows():
            if w.window in ("early", "late"):
                win[w.window] = dict(p=r(w.gene_pct), k=None if pd.isna(w.gene_rank) else int(w.gene_rank),
                                     rg=r(w.rev_gene, 4), rv=r(w.rev, 3), nc=None if pd.isna(w.n_conditions) else int(w.n_conditions),
                                     nl=None if pd.isna(w.n_lines) else int(w.n_lines), pid=w.pert_id)
        unknown = s(g[g.window == "unknown"].window_unknown_reason.iloc[0]) if (g.window == "unknown").any() else None
        pid = (win.get("late") or win.get("early") or {}).get("pid", b.pert_id)
        a = ann.loc[pid] if pid in ann.index else None
        spr = SP.loc[(unit, pid)] if (unit, pid) in SP.index else None
        f = FR.loc[(unit, b.cmap_name)] if (unit, b.cmap_name) in FR.index else None
        lb = LAB.loc[(unit, pid)] if (unit, pid) in LAB.index else None
        genes = [x for x in str(b.layered_genes).split("|") if x and x != "nan"]
        graw = []
        for gn in genes:
            if (unit, gn) in gc.index:
                q = gc.loc[(unit, gn)]
                graw.append(dict(g=gn, dz=r(q.dep_z, 2), D=r(q.D), af=r(q.alt_share, 3), A=r(q.A), o=r(q.o_raw, 2),
                                 O=r(q.O), os=s(q.o_source)))
            else:
                graw.append(dict(g=gn))
        recs.append(dict(
            u=unit, lv="subgroup" if level == "stratum" else "disease", t=b.table, drug=s(b.fda_ingredient) or s(b.cmap_name),
            cm=s(b.cmap_name), path=s(b.path), el=bool(b.eligible) if b.table == "main" else None,
            T=r(b.layered_T), D=r(b.layered_D), A=r(b.layered_A), O=r(b.layered_O),
            nax=None if pd.isna(b.layered_n_axes) else int(b.layered_n_axes), ori=s(b.layered_orient), genes=genes,
            graw=graw, lead=r(b.layered_pct), w=win, unk=unknown,
            ev=s(spr.evidence_targets) if spr is not None else None,
            tsrc=s(spr.target_source) if spr is not None else None,
            sec=s(spr.secondary_targets) if spr is not None else None,
            car=r(spr.target_carrier_share, 3) if spr is not None and "target_carrier_share" in spr else None,
            ub=bool(spr.unverified_binding) if spr is not None and s(spr.get("unverified_binding")) in ("True", "true", "1") else False,
            moa=s(a.hub_moa) if a is not None else None, hubt=s(a.hub_target) if a is not None else None,
            lt=s(a.lincs_target) if a is not None else None,
            ind=(s(f.approved_indications) if f is not None else None) or (s(a.approved_indications) if a is not None else None),
            fs=s(f.fda_status) if f is not None else None,
            fy=s(f.first_fda_approval) if f is not None else (s(a.chembl_first_approval) if a is not None else None),
            fl=s(f.lincs_flags) if f is not None else None,
            ph=r(lb.label_phase, 1) if lb is not None and pd.notna(lb.label_phase) else None,
            phs=s(lb.label_source) if lb is not None else None,
            url=s(f.drugsatfda_url) if f is not None else None,
        ))

    units = []
    for (unit, level), g in X.groupby(["unit", "level"]):
        m = g[(g.table == "main")]
        st = strata.loc[unit] if unit in strata.index else None
        units.append(dict(
            id=unit, lv="subgroup" if level == "stratum" else "disease",
            parent=s(st.group_id) if st is not None else None, desc=s(st.description) if st is not None else None,
            n_el=int(m[m.eligible].ing.nunique()),
            n_e=int(m[(m.window == "early") & m.gene_pct.notna()].ing.nunique()),
            n_l=int(m[(m.window == "late") & m.gene_pct.notna()].ing.nunique())))
    run = json.loads((RES / f"gene_ranking_run{SFX}.json").read_text())
    d127 = RES / "d127/d211_gene_summary.tsv"
    meta = dict(run=run, built=pd.Timestamp.now().strftime("%Y-%m-%d"),
                d127=(lambda t: t.astype(object).where(t.notna(), None).to_dict("records"))(pd.read_csv(d127, sep="\t"))
                if d127.exists() else [])
    # one small index plus one file per unit, fetched when the unit is opened (the whole set is ~10 MB)
    (OUT / "data").mkdir(exist_ok=True)
    for old in (OUT / "data").glob("*.json"):
        old.unlink()
    (OUT / "data/index.json").write_text(json.dumps(dict(meta=meta, units=units), separators=(",", ":"), allow_nan=False))
    by = {}
    for rec in recs:
        by.setdefault(rec["u"], []).append(rec)
    for u, rr in by.items():
        (OUT / f"data/{u}.json").write_text(json.dumps(rr, separators=(",", ":"), allow_nan=False))
    if (OUT / "data.json").exists():
        (OUT / "data.json").unlink()
    shutil.copy(TEMPLATE, OUT / "index.html")
    sizes = [f.stat().st_size for f in (OUT / "data").glob("*.json")]
    print(f"{len(recs):,} (unit, ingredient) records over {len(units)} units -> {OUT}/data "
          f"({len(sizes)} files, largest {max(sizes) / 1e3:.0f} kB, total {sum(sizes) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
