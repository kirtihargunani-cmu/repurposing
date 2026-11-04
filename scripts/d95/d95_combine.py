#!/usr/bin/env python
"""D95 (rule fixed 2026-11-03, docs/DECISIONS.md): metaVIPER-style combination of the breast networks BRCA (pooled),
BRCA_LUM and BRCA_BAS, as viper() combines a list of interactomes (mvws = 1): per regulator,
NES = sum_k NES_k |NES_k| / sum_k |NES_k| over the networks that score it.
Disease side: group_mr of the four breast groups on each network (same seeds: _d92bpooled, _d95lum, _d95bas) ->
  results/step6/aracne_project_v2_tf_d95/{group_mr,regulator_filter}.tsv (filter: kept on any network).
Drug side: activity on each network (BRCA from _d85; BRCA_LUM / BRCA_BAS from 79d runs _d95lum / _d95bas) ->
  data/processed/lincs/drug_viper_specific_v2_tf_d95/BRCA_META.{f32,regulators.txt}; other networks linked from _d85.
P1: a regulator is proliferation-linked on BRCA_META if linked on any of the three -> results/lincs/prolif_regulators_v2_tf_d95extra.tsv
Reported: weight share per network -> results/lincs/d95/weight_share.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
S6 = ROOT / "results/step6"; P = ROOT / "data/processed/lincs"; OUT = ROOT / "results/lincs/d95"; OUT.mkdir(parents=True, exist_ok=True)
G = ["BRCA-LumA", "BRCA-LumB", "BRCA-Basal", "BRCA-Her2"]
SRC = {"BRCA": "aracne_project_v2_tf_d92bpooled", "BRCA_LUM": "aracne_project_v2_tf_d95lum", "BRCA_BAS": "aracne_project_v2_tf_d95bas"}


def combine(mats):
    """mats: dict net -> DataFrame (rows x regulators); returns combined NES and weight shares."""
    regs = sorted(set().union(*(m.columns for m in mats.values())))
    num = 0; den = 0; w = {}
    for k, m in mats.items():
        x = m.reindex(columns=regs).fillna(0.0).values
        num = num + x * np.abs(x); den = den + np.abs(x); w[k] = np.abs(x)
    nes = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
    share = {k: np.where(den > 0, v / np.where(den > 0, den, 1), np.nan) for k, v in w.items()}
    return regs, nes, share


# disease side
mr = {k: pd.read_csv(S6 / d / "group_mr.tsv", sep="\t") for k, d in SRC.items()}
rows, ws = [], []
for g in G:
    mats = {k: m[m.group_id == g].set_index("regulator").nes.to_frame().T for k, m in mr.items()}
    regs, nes, share = combine(mats)
    fdr = pd.concat([m[m.group_id == g].set_index("regulator").fdr for m in mr.values()], axis=1).min(axis=1).reindex(regs)
    p = pd.concat([m[m.group_id == g].set_index("regulator").p for m in mr.values()], axis=1).min(axis=1).reindex(regs)
    rows.append(pd.DataFrame({"group_id": g, "regulator": regs, "nes": nes[0], "p": p.values, "fdr": fdr.values}))
    ws.append(pd.DataFrame({"side": "disease", "unit": g, **{k: np.nanmean(v) for k, v in share.items()}}, index=[0]))
od = S6 / "aracne_project_v2_tf_d95"; od.mkdir(exist_ok=True)
pd.concat(rows).to_csv(od / "group_mr.tsv", sep="\t", index=False)
flt = pd.concat([pd.read_csv(S6 / d / "regulator_filter.tsv", sep="\t") for d in SRC.values()])
flt = flt[flt.group_id.isin(G)].sort_values("keep", ascending=False).drop_duplicates(["group_id", "regulator"])
flt.to_csv(od / "regulator_filter.tsv", sep="\t", index=False)

# drug side
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t")
def act(d, net):
    regs = (P / d / f"{net}.regulators.txt").read_text().split()
    return pd.DataFrame(np.fromfile(P / d / f"{net}.f32", dtype=np.float32).reshape(len(idx), len(regs)), columns=regs)
A = {"BRCA": act("drug_viper_specific_v2_tf_d85", "BRCA"), "BRCA_LUM": act("drug_viper_specific_v2_tf_d95lum", "BRCA_LUM"),
     "BRCA_BAS": act("drug_viper_specific_v2_tf_d95bas", "BRCA_BAS")}
regs, nes, share = combine(A)
dd = P / "drug_viper_specific_v2_tf_d95"; dd.mkdir(exist_ok=True)
nes.astype(np.float32).tofile(dd / "BRCA_META.f32"); (dd / "BRCA_META.regulators.txt").write_text("\n".join(regs) + "\n")
for f in (P / "drug_viper_specific_v2_tf_d85").iterdir():
    t = dd / f.name
    if not t.exists():
        t.symlink_to(f.resolve())
cd = P / "drug_viper_specific_collectri_d95"; cd.mkdir(exist_ok=True)
for f in (P / "drug_viper_specific_collectri_d85").iterdir():
    if not (cd / f.name).exists():
        (cd / f.name).symlink_to(f.resolve())
ws.append(pd.DataFrame({"side": "drug", "unit": "all profiles", **{k: np.nanmean(v) for k, v in share.items()}}, index=[0]))
# MCF7 anti-estrogens, per regulator: combined NES and where it comes from
cmp = pd.read_csv(ROOT / "results/lincs/compounds.tsv", sep="\t")[["pert_id", "cmap_name"]]
r = idx.reset_index().merge(cmp, on="pert_id")
r = r[(r.ModelID == "ACH-000019") & r.cmap_name.isin(["tamoxifen", "fulvestrant", "raloxifene", "toremifene", "bazedoxifene", "clomifene"])]["index"].values
mm = []
for reg in ["ESR1", "FOXA1", "GATA3", "XBP1", "PGR"]:
    if reg in regs:
        j = regs.index(reg)
        mm.append(dict(regulator=reg, combined=round(float(np.median(nes[r, j])), 2), share_negative=round(float((nes[r, j] < 0).mean()), 2),
                       **{f"{k}_nes": round(float(np.median(A[k][reg].values[r])), 2) if reg in A[k] else np.nan for k in A},
                       **{f"{k}_weight": round(float(np.nanmean(share[k][r, j])), 2) for k in A}))
pd.DataFrame(mm).to_csv(OUT / "mcf7_antier.tsv", sep="\t", index=False)
pd.concat(ws).to_csv(OUT / "weight_share.tsv", sep="\t", index=False, float_format="%.3f")
# P1
pl = pd.concat([pd.read_csv(ROOT / "results/lincs/prolif_regulators_v2_tf.tsv", sep="\t").query("net in ['BRCA', 'BRCA_LUM']"),
                pd.read_csv(ROOT / "results/lincs/prolif_regulators_v2_tf_d92bextra.tsv", sep="\t").query("net == 'BRCA_BAS'")])
pl = pl.groupby("regulator").prolif_linked.any().reset_index().assign(net="BRCA_META")
pl.to_csv(ROOT / "results/lincs/prolif_regulators_v2_tf_d95extra.tsv", sep="\t", index=False)
pd.DataFrame({"group_id": G, "net": "BRCA_META"}).to_csv(ROOT / "results/step5/d95_group_net.tsv", sep="\t", index=False)
print(pd.concat(ws).round(3).to_string(index=False)); print(pd.DataFrame(mm).to_string(index=False))
