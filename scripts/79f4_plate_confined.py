#!/usr/bin/env python
"""D109a plate-confinement flag (user 2026-11-06), drug level.

Compound plate = first token of the LINCS sig_id (e.g. REP.A007). plate_confined when (a) every active signature
(is_hiq, TAS >= 0.4) of the drug comes from one compound plate and (b) >= 1 other compound plate tested the drug at doses
overlapping the active range (0.5 x lowest .. 2 x highest active dose) with >= 3 signatures, none active, median TAS <
0.5 x the median TAS of the active signatures. Output results/lincs/d109/plate_confined.tsv; 79f3's
evidence_flags_v2.tsv gets the column plate_confined.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
s = pd.read_csv(ROOT / "data/raw/lincs2020/siginfo_beta.txt", sep="\t", low_memory=False,
                usecols=["sig_id", "pert_id", "cmap_name", "pert_dose", "pert_dose_unit", "tas", "is_hiq", "pert_type"])
s = s[(s.pert_type == "trt_cp") & s.pert_dose_unit.astype(str).str.lower().isin(["um", "nm"])].copy()
s["dose"] = pd.to_numeric(s.pert_dose, errors="coerce")
s.loc[s.pert_dose_unit.astype(str).str.lower().eq("nm"), "dose"] /= 1000
s["plate"] = s.sig_id.str.split("_").str[0]
s["active"] = (s.is_hiq == 1) & (s.tas >= 0.4)
act_drugs = set(s.pert_id[s.active])
rows = []
for p, d in s[s.pert_id.isin(act_drugs)].groupby("pert_id"):
    a = d[d.active]
    plates = a.plate.unique()
    res = dict(pert_id=p, cmap_name=d.cmap_name.iloc[0], n_active=len(a), n_active_plates=len(plates), plate_confined=False, active_plate="", other_plate="")
    if len(plates) == 1:
        lo, hi = 0.5 * a.dose.min(), 2 * a.dose.max()
        o = d[(d.plate != plates[0]) & d.dose.between(lo, hi)]
        for op, od in o.groupby("plate"):
            if len(od) >= 3 and not od.active.any() and od.tas.median() < 0.5 * a.tas.median():
                res.update(plate_confined=True, active_plate=plates[0], other_plate=op)
                break
        res["active_plate"] = plates[0]
    rows.append(res)
pc = pd.DataFrame(rows)
(RES / "d109").mkdir(exist_ok=True)
pc.to_csv(RES / "d109/plate_confined.tsv", sep="\t", index=False)
print(f"drugs with an active signature {len(pc)}; single active plate {int((pc.n_active_plates == 1).sum())}; plate_confined {int(pc.plate_confined.sum())}")
for n in ["metronidazole", "pentobarbital", "doxorubicin", "temsirolimus", "tofacitinib", "erlotinib", "TG-101348", "digoxin"]:
    x = pc[pc.cmap_name == n]
    print(f"  {n:14}", x[["n_active", "n_active_plates", "plate_confined", "active_plate", "other_plate"]].values.tolist())
ef = pd.read_csv(RES / "evidence_flags_v2.tsv", sep="\t", low_memory=False).drop(columns=["plate_confined"], errors="ignore")
ef["plate_confined"] = ef.pert_id.map(pc.set_index("pert_id").plate_confined).fillna(False).astype(bool)
ef.to_csv(RES / "evidence_flags_v2.tsv", sep="\t", index=False)
f = pd.read_csv(RES / "fda_ranking_v2_consensus.tsv", sep="\t")
sp = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False, usecols=["group_id", "pert_id", "cmap_name"])
t = f[f.fda_rank <= 10].merge(sp, on=["group_id", "cmap_name"]).merge(ef[["group_id", "pert_id", "plate_confined"]], on=["group_id", "pert_id"])
print(f"FDA top-10 rows {len(t)}; plate_confined {int(t.plate_confined.sum())}: {sorted(t[t.plate_confined].drug.unique())}")
