#!/usr/bin/env python
"""Exploratory (2026-11-03): a lagged ARACNe across drugs, adapted from TimeDelay-ARACNE (Zoppoli et al. 2010).
TimeDelay-ARACNE needs one long time series; LINCS has many drugs with 2 time points (6 h, 24 h). Here the samples are
the drugs of one line, the regulators are measured at 6 h and the targets are each landmark gene's late change: its
24 h value minus what its own 6 h value predicts (OLS over the line's drugs), so persistence of the early response
is removed and a link means "early regulator state predicts what changes later".
Variants: mrna = the 81 landmark TFs' own 6 h z (as TimeDelay-ARACNE uses expression); act = 6 h regulator activity
from the tumor ARACNe network of the line's lineage (scripts/timedose/early_late.R, E signatures), which includes
regulators not measured on L1000 (ESR1, FOXA1, AR).
Inputs: data/processed/lincs/timedose/time_{T6,T24}_specific.f32 + time_index.tsv; early_late/E_<net>.f32.
Outputs: data/processed/lincs/pnet_lag_{mrna,act}/<LINE>/all.exp.tsv and regulators.txt (one per variant);
         run ARACNe with PNET_DIR=data/processed/lincs/pnet_lag_<v> scripts/pnet/pnet_run_aracne.sh
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
P = ROOT / "data/processed/lincs"; TD = P / "timedose"
LINES = [("MCF7", "ACH-000019", "BRCA"), ("PC3", "ACH-000090", "PRAD"), ("A375", "ACH-000219", "SKCM")]
pg = pd.read_csv(P / "drug_line_profiles_genes.tsv", sep="\t")
I = pd.read_csv(TD / "time_index.tsv", sep="\t"); n = len(I)
lm = pd.read_csv(P / "pnet/MCF7_specific/all.exp.tsv", sep="\t", usecols=["gene"]).gene
regs = (P / "pnet/regulators.txt").read_text().split()
gi = pd.Series(np.arange(len(pg)), index=pg.gene_id)
E = np.fromfile(TD / "time_T6_specific.f32", dtype=np.float32).reshape(n, len(pg))
L = np.fromfile(TD / "time_T24_specific.f32", dtype=np.float32).reshape(n, len(pg))
idm = pd.read_csv(ROOT / "results/step5/id_map.tsv", sep="\t"); idm = idm[idm.primary_for_hgnc == True]
s2e = dict(zip(idm.symbol, idm.gene_id))
for v in ["mrna", "act"]:
    (P / f"pnet_lag_{v}").mkdir(exist_ok=True)
act_regs = set()
for line, mid, net in LINES:
    k = np.where(I.ModelID.values == mid)[0]
    e, l = E[k][:, gi[lm].values].astype(np.float64), L[k][:, gi[lm].values].astype(np.float64)
    ec, lc = e - e.mean(0), l - l.mean(0)
    b = (ec * lc).sum(0) / np.maximum((ec * ec).sum(0), 1e-12)
    res = lc - ec * b                                   # late change not predicted by the gene's own 6 h value
    samples = [f"{p}" for p in I.pert_id.values[k]]
    tgt = pd.DataFrame(res.T, index=lm.values, columns=samples)
    # mrna variant
    r6 = pd.DataFrame(e[:, [list(lm).index(r) for r in regs]].T, index=["R6_" + r for r in regs], columns=samples)
    x = pd.concat([r6, tgt]); x.index.name = "gene"
    d = P / "pnet_lag_mrna" / line; d.mkdir(exist_ok=True); x.to_csv(d / "all.exp.tsv", sep="\t", float_format="%.5f")
    # act variant: 6 h activity from the lineage network
    rn = (TD / "early_late" / f"{net}.regulators.txt").read_text().split()
    A = np.fromfile(TD / "early_late" / f"E_{net}.f32", dtype=np.float32).reshape(n, len(rn))[k]
    ids = ["A6_" + (r if r.startswith("ENSG") else s2e.get(r, r)) for r in rn]
    ra = pd.DataFrame(A.T, index=ids, columns=samples); ra = ra[~ra.index.duplicated()]
    x = pd.concat([ra, tgt]); x.index.name = "gene"
    d = P / "pnet_lag_act" / line; d.mkdir(exist_ok=True); x.to_csv(d / "all.exp.tsv", sep="\t", float_format="%.5f")
    act_regs |= set(ra.index)
    print(f"{line}: {len(k)} drugs; late change kept {np.mean((res ** 2).sum(0) / np.maximum((lc ** 2).sum(0), 1e-12)):.2f} of the 24 h variance; "
          f"regulators mrna {len(regs)}, act {len(ra)}")
(P / "pnet_lag_mrna/regulators.txt").write_text("\n".join("R6_" + r for r in regs) + "\n")
(P / "pnet_lag_act/regulators.txt").write_text("\n".join(sorted(act_regs)) + "\n")
