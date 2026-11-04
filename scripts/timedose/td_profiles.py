#!/usr/bin/env python
"""D88 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): time- and dose-resolved drug x line profiles.
Runs in `celligner`. Genes and gene order as scripts/79c (drug_line_profiles_genes.tsv; level-5 genes matched by symbol).
Arms (per drug x line pair of the 79c index):
  time pairs (6 h and 24 h hiq signatures):  A0 (79c profile), T6 (mean of all 6 h signatures), T24 (all 24 h)
  dose pairs (24 h, >= 3 doses, lowest dose with TAS >= 0.4 below the top dose):  A0, DL (signatures at the lowest
    active dose), DH (signatures at the highest dose)
Versions: raw, and specific = each arm's line vector (mean of that arm's profiles in the line, unit length) projected
out, as 79c (A0 specific = the 79c specific profile). Proliferation score = mean of HALLMARK E2F + G2M genes.
Outputs: data/processed/lincs/timedose/<comparison>_<arm>_<version>.f32 (profiles x genes, float32, C order),
         <comparison>_index.tsv (pert_id, ModelID, cmap_name, n_sig per arm, active_6h), <comparison>_prolif.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; OUT = P / "timedose"; OUT.mkdir(exist_ok=True)
TAS = 0.4

s = pd.read_csv(P / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
s = s[s.ModelID.notna() & (s.pert_type == "trt_cp")].copy()
s["t"] = s.pert_time.round().astype(int); s["dose"] = s.pert_dose.round(3)
X = np.load(P / "level5_bing_hiq.npy", mmap_mode="r")
lg = pd.read_csv(P / "level5_bing_genes.tsv", sep="\t")
pg = pd.read_csv(P / "drug_line_profiles_genes.tsv", sep="\t")
pos = pd.Series(np.arange(len(lg)), index=lg.gene_symbol)
assert pos.index.is_unique and pg.symbol.isin(pos.index).all()
cols = pos[pg.symbol].values
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); idx["row"] = np.arange(len(idx))
ng = len(pg)
A0 = {"raw": np.fromfile(P / "drug_line_profiles_bing.f32", dtype=np.float32).reshape(len(idx), ng),
      "specific": np.fromfile(P / "drug_line_profiles_bing_specific.f32", dtype=np.float32).reshape(len(idx), ng)}
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt")}
pm = pg.symbol.isin(set(gmt["HALLMARK_E2F_TARGETS"]) | set(gmt["HALLMARK_G2M_CHECKPOINT"])).values

def mean_sig(ix):
    return np.asarray(X[np.sort(ix)][:, cols], dtype=np.float64).mean(0)

def specific(M, lines):
    S = M.copy()
    for l in np.unique(lines):
        k = lines == l
        g = S[k].mean(0); g /= np.linalg.norm(g)
        S[k] = S[k] - np.outer(S[k] @ g, g)
    return S

by = s.groupby(["pert_id", "ModelID"])
pairs = {k: v for k, v in by.indices.items()}
time_rows, dose_rows = [], []
for r in idx.itertuples():
    ix = pairs.get((r.pert_id, r.ModelID))
    if ix is None:
        continue
    d = s.iloc[ix]
    if {6, 24} <= set(d.t):
        time_rows.append((r.row, r.pert_id, r.ModelID, d.cmap_name.iloc[0], d.index[d.t == 6].values, d.index[d.t == 24].values,
                          float(d[d.t == 6].tas.max())))
    e = d[d.t == 24]
    doses = np.sort(e.dose.unique())
    if len(doses) >= 3 and (e.tas >= TAS).any():
        lo = e[e.tas >= TAS].dose.min()
        if lo < doses[-1]:
            dose_rows.append((r.row, r.pert_id, r.ModelID, d.cmap_name.iloc[0], e.index[e.dose == lo].values, e.index[e.dose == doses[-1]].values, lo, doses[-1]))

def emit(name, rows, arms, extra):
    I = pd.DataFrame([x[:4] for x in rows], columns=["row", "pert_id", "ModelID", "cmap_name"])
    for k, v in extra.items():
        I[k] = v
    lines = I.ModelID.values
    prol = {}
    for arm, sel in arms.items():
        raw = A0["raw"][I.row.values].astype(np.float64) if sel is None else np.vstack([mean_sig(x[sel]) for x in rows])
        spec = A0["specific"][I.row.values].astype(np.float64) if sel is None else specific(raw, lines)
        for ver, M in [("raw", raw), ("specific", spec)]:
            M.astype(np.float32).tofile(OUT / f"{name}_{arm}_{ver}.f32")
            prol[f"{arm}_{ver}"] = M[:, pm].mean(1)
        if sel is not None:
            I[f"n_sig_{arm}"] = [len(x[sel]) for x in rows]
    I.to_csv(OUT / f"{name}_index.tsv", sep="\t", index=False)
    pd.DataFrame(prol).to_csv(OUT / f"{name}_prolif.tsv", sep="\t", index=False, float_format="%.5g")
    small = I.groupby("ModelID").size()
    print(f"{name}: {len(I)} pairs, {I.pert_id.nunique()} drugs, {I.ModelID.nunique()} lines "
          f"(pairs in lines with < 10 pairs: {int(small[small < 10].sum())})")

emit("time", time_rows, {"A0": None, "T6": 4, "T24": 5}, {"tas_6h": [x[6] for x in time_rows]})
emit("dose", dose_rows, {"A0": None, "DL": 4, "DH": 5}, {"dose_low": [x[6] for x in dose_rows], "dose_high": [x[7] for x in dose_rows]})
