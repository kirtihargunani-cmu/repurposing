#!/usr/bin/env python
"""D90 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): mechanism (MECH) drug x line profiles, aligned to the
79c index and genes so 79d / 79e / 79i run on them unchanged (runs in `celligner`).
Per 79c pair: 6 h profile at the lowest 6 h dose with TAS >= 0.4 if active at 6 h ("early_low"); else the 24 h profile
at the lowest active dose if >= 3 doses at 24 h and that dose is below the top ("late_low"); else the 79c profile
("timing_unknown"). Profile = mean of the hiq level-5 signatures at that time and dose. Line-specific version as 79c,
with each line's vector from the MECH raw profiles.
Outputs data/processed/lincs/d90/mech_{raw,specific}.f32 (profiles x genes, C order), results/lincs/d90/mech_source.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; OUT = P / "d90"; RES = ROOT / "results/lincs/d90"
OUT.mkdir(exist_ok=True); RES.mkdir(parents=True, exist_ok=True)
TAS = 0.4
s = pd.read_csv(P / "level5_bing_hiq_sigs.tsv.gz", sep="\t")
s = s[s.ModelID.notna() & (s.pert_type == "trt_cp")].copy()
s["t"] = s.pert_time.round().astype(int); s["dose"] = s.pert_dose.round(3)
X = np.load(P / "level5_bing_hiq.npy", mmap_mode="r")
lg = pd.read_csv(P / "level5_bing_genes.tsv", sep="\t"); pg = pd.read_csv(P / "drug_line_profiles_genes.tsv", sep="\t")
cols = pd.Series(np.arange(len(lg)), index=lg.gene_symbol)[pg.symbol].values
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t")
ng = len(pg)
M = np.fromfile(P / "drug_line_profiles_bing.f32", dtype=np.float32).reshape(len(idx), ng).astype(np.float64)
pairs = s.groupby(["pert_id", "ModelID"]).indices
src, info = [], []
for i, r in enumerate(idx.itertuples()):
    ix = pairs.get((r.pert_id, r.ModelID)); kind, pick, dose = "timing_unknown", None, np.nan
    if ix is not None:
        d = s.iloc[ix]
        e6 = d[(d.t == 6) & (d.tas >= TAS)]
        if len(e6):
            dose = e6.dose.min(); pick = d[(d.t == 6) & (d.dose == dose)].index.values; kind = "early_low"
        else:
            e = d[d.t == 24]; doses = np.sort(e.dose.unique())
            if len(doses) >= 3 and (e.tas >= TAS).any() and e[e.tas >= TAS].dose.min() < doses[-1]:
                dose = e[e.tas >= TAS].dose.min(); pick = e.index[e.dose == dose].values; kind = "late_low"
    if pick is not None:
        M[i] = np.asarray(X[np.sort(pick)][:, cols], dtype=np.float64).mean(0)
    src.append(kind); info.append(dose)
M.astype(np.float32).tofile(OUT / "mech_raw.f32")
S = M.copy()
for ix in idx.groupby("ModelID").indices.values():
    g = S[ix].mean(0); g /= np.linalg.norm(g); S[ix] = S[ix] - np.outer(S[ix] @ g, g)
S.astype(np.float32).tofile(OUT / "mech_specific.f32")
o = idx.assign(source=src, dose=info)
o.to_csv(RES / "mech_source.tsv", sep="\t", index=False)
print(o.source.value_counts().to_string())
