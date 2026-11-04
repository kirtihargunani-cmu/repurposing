#!/usr/bin/env python
"""D104 sizing (descriptive, no ranking change): for each cohort x drug x receptor the drug binds (79w binding evidence
mechanism / potent), tumor median TPM vs the weighted mean DepMap TPM of the lines its score came from (79e / 79w line
rule: direct = in-key lines with Celligner weights, borrowed = out-of-key lines, equal weights). Provisional classes:
line_only (line >= 1 TPM, tumor < 10% of line), tumor_only (tumor >= 1, line < 1), both, neither.
Poorly matched = unmatched_lines (direct, line_match < 0.1) or borrowed. Groups (41) only.
Outputs results/lincs/d104/size_rows.tsv, size_summary.tsv"""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; RES = ROOT / "results/lincs"; OUT = RES / "d104"; RAW = ROOT / "data/raw"
spec = importlib.util.spec_from_file_location("w79", ROOT / "scripts/79w_receptor_gate.py"); w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)
r = pd.read_csv(RES / "receptor_gate_v2.tsv", sep="\t", low_memory=False)
r = r[(r.level == "disease") & r.binding_evidence.isin(["mechanism", "potent"])].copy()
sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False,
                 usecols=["group_id", "pert_id", "cmap_name", "source", "line_match", "unmatched_lines", "spec_mech_pct", "spec_mech_pct_pregate", "target_absent", "tier_spec_mech", "approved"])
sp = sp.drop_duplicates(["group_id", "pert_id"])
r = r.merge(sp.drop(columns=["cmap_name", "source"]), on=["group_id", "pert_id"], how="inner")
idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
by_drug = {p: d.ModelID.values for p, d in idx.groupby("pert_id")}
keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
models = pd.read_csv(RAW / "depmap_24q4/Model.csv", index_col=0); models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
W = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
wmap = {g: d.set_index("ModelID").weight for g, d in W.groupby("group_id")}
f = RAW / "depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv"; hdr = pd.read_csv(f, nrows=0).columns
genes = set(r.target)
X = pd.read_csv(f, usecols=[hdr[0]] + [c for c in hdr[1:] if c.split(" (")[0] in genes], index_col=0)
X.columns = [c.split(" (")[0] for c in X.columns]; X = X.loc[:, ~X.columns.duplicated()]; X = np.power(2.0, X) - 1
ks = {}
def km(g):
    if g not in ks:
        k = keys.loc[g]; m = models[models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        ks[g] = set(m.index)
    return ks[g]
lt, nl = [], []
for x in r.itertuples(index=False):
    L = by_drug.get(x.pert_id, np.array([]))
    if x.source == "direct":
        L = [m for m in L if m in km(x.group_id)]; ww = wmap.get(x.group_id, pd.Series(dtype=float)).reindex(L).fillna(0).values
    else:
        L = [m for m in L if m not in km(x.group_id)]; ww = np.zeros(len(L))
    ww = ww if ww.sum() > 0 else np.ones(len(L))
    s = pd.Series(ww, index=L); s = s[s.index.isin(X.index)]; s = s[~s.index.duplicated()]
    if len(s) and x.target in X.columns and s.sum() > 0:
        lt.append(float((s * X.loc[s.index, x.target]).sum() / s.sum())); nl.append(len(s))
    else:
        lt.append(np.nan); nl.append(len(s))
r["line_tpm"] = np.round(lt, 2); r["n_lines"] = nl
r["poor_match"] = (r.source == "borrowed") | r.unmatched_lines.astype(str).str.lower().eq("true")
t, l = r.median_tpm, r.line_tpm
r["cls"] = np.select([l.isna(), (l >= 1) & (t < 0.1 * l), (t >= 1) & (l < 1), (t >= 1) & (l >= 1)], ["no_line_data", "line_only", "tumor_only", "both"], "neither")
r["log2_tumor_over_line"] = np.round(np.log2((t + 1) / (l + 1)), 2)
r.to_csv(OUT / "size_rows.tsv", sep="\t", index=False)
dp = r.groupby(["group_id", "pert_id"]).agg(poor=("poor_match", "first"), cls=("cls", lambda c: "line_only_all" if (c == "line_only").all() else ("tumor_only_any" if (c == "tumor_only").any() else "other")),
                                             absent=("target_absent", "first"), pct=("spec_mech_pct", "first"), tier=("tier_spec_mech", "first"))
S = []
for poor, d in r.groupby("poor_match"):
    for c, n in d.cls.value_counts().items():
        S.append(dict(level="receptor_row", poor_match=poor, cls=c, n=n, n_drugs=d[d.cls == c].pert_id.nunique()))
for poor, d in dp.groupby("poor"):
    for c, n in d.cls.value_counts().items():
        e = d[d.cls == c]
        S.append(dict(level="cohort_x_drug", poor_match=poor, cls=c, n=n, already_d102_absent=int(e.absent.sum()),
                      top10pct_not_absent=int(((e.pct >= 90) & ~e.absent).sum()), median_pct=round(e.pct.median(), 1)))
S = pd.DataFrame(S); S.to_csv(OUT / "size_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); print(S.to_string(index=False))
pm = r[r.poor_match & (r.cls == "line_only")].sort_values("spec_mech_pct", ascending=False)
print(pm[["group_id", "cmap_name", "target", "median_tpm", "line_tpm", "n_lines", "source", "spec_mech_pct", "target_absent"]].head(15).to_string(index=False))
tm = r[r.poor_match & (r.cls == "tumor_only")].sort_values("median_tpm", ascending=False)
print(tm[["group_id", "cmap_name", "target", "median_tpm", "line_tpm", "n_lines", "source", "spec_mech_pct"]].head(12).to_string(index=False))
