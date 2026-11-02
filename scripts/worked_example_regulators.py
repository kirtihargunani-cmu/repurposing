#!/usr/bin/env python
"""Per-regulator decomposition of the published consensus rev for two drug-disease pairs (read-only; mirrors 79e / 79i).

Not a pipeline step (figure helper, 2026-11-01). Reads existing matrices only; writes results/lincs/worked_example_regulators.tsv.
Columns: disease_nes = ARACNe TF NES; weight_used = consensus NES (the w_i in rev); drug_activity = consensus a_i;
contribution = -w_i a_i / sqrt(sum w^2) (sums to rev); prolif_linked = ARACNe TF P1 list; prolif_linked_79i = P1 or P2
in any set (the 79i definition behind rev_prolif / rev_nonprolif).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC, RES = ROOT / "data/processed/lincs", ROOT / "results/lincs"
PAIRS = [("TCGA-LUAD", "LUAD", "BRD-K70401845", "erlotinib", ["MYBL2", "FOXM1", "FOXA1", "ETV4", "E2F1"]),
         ("BRCA-LumB", "BRCA", "BRD-K51313569", "palbociclib", ["FOXM1", "MYBL2", "E2F1", "ESR1"])]
SETS = {"aracne": ("results/step6/aracne_project_v2_tf", "_v2_tf"), "collectri": ("results/step6/collectri_v2", "_collectri")}

idx = pd.read_csv(PROC / "drug_line_profiles_index.tsv", sep="\t"); idx["col"] = np.arange(len(idx))
keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
dd = pd.read_csv(RES / "disease_drugs.tsv", sep="\t", low_memory=False)
pub = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
pub79e = pd.read_csv(RES / "drug_reversal_v2_consensus.tsv", sep="\t", low_memory=False)
p1 = pd.read_csv(RES / "prolif_regulators_v2_tf.tsv", sep="\t")
p2 = pd.read_csv(RES / "prolif_regulators_p2_v2_consensus.tsv", sep="\t")
mr = {k: pd.read_csv(ROOT / d / "group_mr.tsv", sep="\t") for k, (d, _) in SETS.items()}


def read(tag, net):
    d = PROC / f"drug_viper_specific{tag}"
    regs = pd.Index((d / f"{net}.regulators.txt").read_text().split())
    raw = np.fromfile(d / f"{net}.f32", dtype=np.float32)
    assert raw.size == len(idx) * len(regs), f"{tag}/{net}: size {raw.size} != {len(idx)} x {len(regs)}"
    return regs, raw.reshape(len(idx), len(regs)).astype(np.float64)


out, report = [], []
for g, net, pid, name, named in PAIRS:
    src = dd[(dd.group_id == g) & (dd.pert_id == pid)].source
    assert len(src) == 1, (g, pid, len(src)); src = src.iloc[0]
    k = keys.loc[g]
    km = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
    if isinstance(k.diseases, str) and k.diseases:
        km = km[km.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
    wg = w[w.group_id == g].set_index("ModelID").weight
    p = idx[idx.pert_id == pid]
    p = p[p.ModelID.isin(km.index)] if src == "direct" else p[~p.ModelID.isin(km.index)]
    ww = wg.reindex(p.ModelID).fillna(0).values if src == "direct" else np.zeros(len(p))
    equal = ww.sum() == 0
    ww = np.ones(len(p)) if equal else ww
    lines = ", ".join(f"{models.CellLineName.get(m, m)} w={x:.4f}" for m, x in zip(p.ModelID, wg.reindex(p.ModelID).fillna(0)))

    per_set = {}
    for s, (_, tag) in SETS.items():
        regs, A = read(tag, net if s == "aracne" else "ALL")
        a = (ww[:, None] * A[p.col.values]).sum(0) / ww.sum()              # Celligner-weighted mean over lines
        m = mr[s][(mr[s].group_id == g) & mr[s].regulator.isin(regs)]       # tested in disease AND scored on drug side
        per_set[s] = pd.DataFrame({"regulator": m.regulator.values, f"nes_{s}": m.nes.values,
                                   f"act_{s}": a[regs.get_indexer(m.regulator)]})
    t = per_set["aracne"].merge(per_set["collectri"], on="regulator", how="outer")
    t["in_aracne"] = t.nes_aracne.notna(); t["in_collectri"] = t.nes_collectri.notna()
    t["n_sets"] = t.in_aracne.astype(int) + t.in_collectri.astype(int)
    t["weight"] = t[["nes_aracne", "nes_collectri"]].mean(axis=1)           # consensus disease NES = w_i
    t["drug_activity"] = t[["act_aracne", "act_collectri"]].mean(axis=1)    # consensus drug activity = a_i
    norm = np.sqrt((t.weight ** 2).sum())
    t["contribution"] = -t.weight * t.drug_activity / norm
    rev = t.contribution.sum()

    p1s = set(p1[(p1.net == net) & p1.prolif_linked].regulator)
    t["prolif_linked"] = t.regulator.isin(p1s)                               # as requested: ARACNe TF P1 list only
    la = set(p2[(p2.set == "aracne") & (p2.net == net) & p2.prolif_linked].regulator)
    lc = set(p2[(p2.set == "collectri") & (p2.net == "ALL") & p2.prolif_linked].regulator)
    t["prolif_linked_79i"] = (t.in_aracne & t.regulator.isin(la)) | (t.in_collectri & t.regulator.isin(lc))

    pr = pub[(pub.group_id == g) & (pub.pert_id == pid)].iloc[0]
    pe = pub79e[(pub79e.group_id == g) & (pub79e.pert_id == pid)].iloc[0]
    split = {}
    for col in ["prolif_linked", "prolif_linked_79i"]:
        sp, sn = t[t[col]], t[~t[col]]
        own = lambda x: -(x.weight * x.drug_activity).sum() / np.sqrt((x.weight ** 2).sum())
        split[col] = dict(n_p=len(sp), n_np=len(sn), contrib_p=sp.contribution.sum(), contrib_np=sn.contribution.sum(),
                          rev_prolif=own(sp), rev_nonprolif=own(sn))
    report.append(dict(g=g, name=name, src=src, equal=equal, lines=lines, n=len(t), rev=rev, rev_pub=pr.rev, rev_79e=pe.rev,
                       n_mrs_pub=pr.n_mrs, spec_mech=pr.spec_mech, rp_pub=pr.rev_prolif, rn_pub=pr.rev_nonprolif,
                       np_pub=pr.n_mrs_prolif, nn_pub=pr.n_mrs_nonprolif, split=split, t=t, named=named, norm=norm,
                       pos_ok=bool(((t.weight > 0) & (t.drug_activity < 0) <= (t.contribution > 0)).all()
                                   and ((t.weight < 0) & (t.drug_activity > 0) <= (t.contribution > 0)).all())))
    t.insert(0, "cmap_name", name); t.insert(0, "group_id", g)
    t = t.rename(columns={"nes_aracne": "disease_nes", "nes_collectri": "disease_nes_collectri",
                          "act_aracne": "drug_activity_aracne", "act_collectri": "drug_activity_collectri", "weight": "weight_used"})
    out.append(t[["group_id", "cmap_name", "regulator", "disease_nes", "disease_nes_collectri", "weight_used",
                  "drug_activity", "drug_activity_aracne", "drug_activity_collectri", "contribution",
                  "prolif_linked", "prolif_linked_79i", "in_aracne", "in_collectri", "n_sets"]]
               .sort_values("contribution", ascending=False))

res = pd.concat(out, ignore_index=True)
res.round(6).to_csv(RES / "worked_example_regulators.tsv", sep="\t", index=False)

pd.set_option("display.width", 200)
cols = ["regulator", "disease_nes", "weight_used", "drug_activity", "contribution", "n_sets", "prolif_linked_79i"]
for r in report:
    t = res[res.group_id == r["g"]]
    print(f"\n===== {r['g']} / {r['name']} ({r['src']}; lines: {r['lines']}{'; equal weights' if r['equal'] else ''})")
    print(f"regulators {r['n']} (published n_mrs {r['n_mrs_pub']}); norm sqrt(sum w^2) = {r['norm']:.3f}")
    print(f"CHECK 1  sum(contribution) = {r['rev']:.4f}   published rev = {r['rev_pub']:.3f} (79e file {r['rev_79e']:.3f})   "
          f"diff {r['rev'] - r['rev_pub']:+.4f}")
    print(f"CHECK 2  positive contribution <=> drug moves regulator opposite to its disease direction: {r['pos_ok']}")
    print("top 10 by contribution:"); print(t.head(10)[cols].round(3).to_string(index=False))
    print("bottom 5:"); print(t.tail(5)[cols].round(3).to_string(index=False))
    print("named regulators:"); print(t[t.regulator.isin(r["named"])].set_index("regulator").reindex(r["named"]).reset_index()[cols].round(3).to_string(index=False))
    for col, s in r["split"].items():
        print(f"split by {col:<18} prolif {s['n_p']:>4} / other {s['n_np']:>4} regulators | contribution sums "
              f"{s['contrib_p']:.3f} + {s['contrib_np']:.3f} = {s['contrib_p'] + s['contrib_np']:.3f} | each on its own norm "
              f"(79i): rev_prolif {s['rev_prolif']:.3f}, rev_nonprolif {s['rev_nonprolif']:.3f}")
    print(f"published: rev_prolif {r['rp_pub']}, rev_nonprolif {r['rn_pub']} (n {r['np_pub']} / {r['nn_pub']}); spec_mech {r['spec_mech']}")
