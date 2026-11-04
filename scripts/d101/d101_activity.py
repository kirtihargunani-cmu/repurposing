#!/usr/bin/env python
"""D101 (rule fixed 2026-11-03, docs/DECISIONS.md): drug-side ESR1 / AR activity from PROGENy Estrogen / Androgen.
In every production ARACNe network that scores ESR1 or AR, the regulator's column is replaced by the PROGENy pathway
score of the same profile (D99, data/processed/lincs/d99/progeny_specific.tsv.gz), z-scored over all profiles and
multiplied by the SD of the regulator's ARACNe activity in that network. Other networks / regulators are links to
production. CollecTRI unchanged (links). Output data/processed/lincs/drug_viper_specific_{v2_tf,collectri}_d101/"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; P = ROOT / "data/processed/lincs"
idx = pd.read_csv(P / "drug_line_profiles_index.tsv", sep="\t"); cond = idx.pert_id + "|" + idx.ModelID
pr = pd.read_csv(P / "d99/progeny_specific.tsv.gz", sep="\t").set_index("condition").reindex(cond)
assert pr.notna().all().all(), "PROGENy score missing for some profiles"
PW = {"ESR1": "Estrogen", "AR": "Androgen"}
z = {r: ((pr[p] - pr[p].mean()) / pr[p].std()).values for r, p in PW.items()}
src, dst = P / "drug_viper_specific_v2_tf", P / "drug_viper_specific_v2_tf_d101"; dst.mkdir(exist_ok=True)
rows = []
for f in sorted(src.glob("*.f32")):
    regs = (src / f"{f.stem}.regulators.txt").read_text().split()
    hit = [r for r in PW if r in regs]
    t = dst / f.name
    if t.exists() or t.is_symlink():
        t.unlink()
    if hit:
        A = np.fromfile(f, dtype=np.float32).reshape(len(idx), len(regs)).copy()
        for r in hit:
            j = regs.index(r); sd = A[:, j].std(); old = A[:, j].copy(); A[:, j] = z[r] * sd
            rows.append(dict(net=f.stem, regulator=r, sd=round(float(sd), 3), corr_old_new=round(float(np.corrcoef(old, A[:, j])[0, 1]), 3)))
        A.tofile(t)
    else:
        t.symlink_to(f.resolve())
    rt = dst / f"{f.stem}.regulators.txt"
    if not rt.exists():
        rt.symlink_to((src / f"{f.stem}.regulators.txt").resolve())
cd = P / "drug_viper_specific_collectri_d101"; cd.mkdir(exist_ok=True)
for f in (P / "drug_viper_specific_collectri").iterdir():
    if not (cd / f.name).exists():
        (cd / f.name).symlink_to(f.resolve())
out = pd.DataFrame(rows); (ROOT / "results/lincs/d101").mkdir(parents=True, exist_ok=True)
out.to_csv(ROOT / "results/lincs/d101/replaced_columns.tsv", sep="\t", index=False); print(out.to_string(index=False))
