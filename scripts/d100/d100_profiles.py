#!/usr/bin/env python
"""D100 (rule fixed 2026-11-03, docs/DECISIONS.md): early / late drug x line profile sets with the 79c layout.
Timed pairs (6 h and 24 h hiq signatures; D88 time arm, data/processed/lincs/timedose/time_{T6,T24}_raw.f32): early = all
6 h signatures, late = all 24 h signatures; other pairs keep the 79c aggregated raw profile. Line-specific version: each
line's vector (mean of the arm's raw profiles in the line, unit length) projected out, as 79c.
Outputs data/processed/lincs/d100/{early,late}_{raw,specific}.f32, timed_pairs.tsv"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; OUT = P / "d100"; OUT.mkdir(exist_ok=True)
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); ng = len(pd.read_csv(P / "drug_line_profiles_genes.tsv", sep="\t"))
agg = np.fromfile(P / "drug_line_profiles_bing.f32", dtype=np.float32).reshape(len(idx), ng)
T = pd.read_csv(P / "timedose/time_index.tsv", sep="\t")
for arm, f in [("early", "T6"), ("late", "T24")]:
    X = np.fromfile(P / f"timedose/time_{f}_raw.f32", dtype=np.float32).reshape(len(T), ng)
    M = agg.astype(np.float64).copy(); M[T.row.values] = X
    M.astype(np.float32).tofile(OUT / f"{arm}_raw.f32")
    for ix in idx.groupby("ModelID").indices.values():
        g = M[ix].mean(0); g /= np.linalg.norm(g); M[ix] = M[ix] - np.outer(M[ix] @ g, g)
    M.astype(np.float32).tofile(OUT / f"{arm}_specific.f32")
    print(arm, "timed pairs replaced:", len(T))
T[["row", "pert_id", "ModelID"]].to_csv(OUT / "timed_pairs.tsv", sep="\t", index=False)
