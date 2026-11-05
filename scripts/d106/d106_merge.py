#!/usr/bin/env python
"""D106: merge the 12-group rerun of 61 into the quality_correction tables (other groups' rows kept from the backup)
and compare before vs after (patient concordance with the group MRs; PC 1 checks; per-regulator agreement of patient NES).
Outputs results/step6/aracne_project_v2_tf/quality_correction/{checks,patient_concordance,quality_index}.tsv (merged),
results/step6/aracne_project_v2_tf/d106_compare.tsv"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
Q = ROOT / "results/step6/aracne_project_v2_tf/quality_correction"; QB = Q.with_name("quality_correction_pre_d106")
B = ROOT / "data/processed/step6/aracne_project_v2_tf_patient"; BB = B.with_name(B.name + "_pre_d106")
G = list(pd.read_csv(ROOT / "config/group_networks.tsv", sep="\t").group_id)
for f in ["checks", "patient_concordance", "quality_index"]:
    new, old = pd.read_csv(Q / f"{f}.tsv", sep="\t"), pd.read_csv(QB / f"{f}.tsv", sep="\t")
    if set(new.group_id) <= set(G):   # first run after 61: merge; a rerun of this script finds the merged table
        pd.concat([old[~old.group_id.isin(G)], new], ignore_index=True).sort_values("group_id").to_csv(Q / f"{f}.tsv", sep="\t", index=False)
v2 = lambda k: k[k.version.str.startswith("v2")]   # each group also has a v1 reference row
co = v2(pd.read_csv(QB / "patient_concordance.tsv", sep="\t")).set_index("group_id")
cn = v2(pd.read_csv(Q / "patient_concordance.tsv", sep="\t")).set_index("group_id")
pc1 = lambda k: v2(k)[v2(k).pc == 1].set_index("group_id")
ko, kn = pc1(pd.read_csv(QB / "checks.tsv", sep="\t")), pc1(pd.read_csv(Q / "checks.tsv", sep="\t"))
rows = []
for g in G:
    o = pd.read_csv(BB / f"{g}.patient_nes.tsv.gz", sep="\t", index_col=0)
    n = pd.read_csv(B / f"{g}.patient_nes.tsv.gz", sep="\t", index_col=0)
    regs, pts = o.index.intersection(n.index), o.columns.intersection(n.columns)
    r = [np.corrcoef(o.loc[x, pts], n.loc[x, pts])[0, 1] for x in regs]
    rows.append(dict(group_id=g, n_patients=len(pts), regulators_old=len(o), regulators_new=len(n), regulators_shared=len(regs),
                     median_r_old_vs_new=round(float(np.nanmedian(r)), 3),
                     concordance_old=co.median_rho.get(g, np.nan), concordance_new=cn.median_rho.get(g, np.nan),
                     frac_hit_old=co.median_frac_hit.get(g, np.nan), frac_hit_new=cn.median_frac_hit.get(g, np.nan),
                     pc1_r2_quality_old=ko.r2_quality.get(g, np.nan), pc1_r2_quality_new=kn.r2_quality.get(g, np.nan),
                     pc1_r2_purity_old=ko.r2_purity.get(g, np.nan), pc1_r2_purity_new=kn.r2_purity.get(g, np.nan),
                     pc1_r2_plate_old=ko.r2_plate.get(g, np.nan), pc1_r2_plate_new=kn.r2_plate.get(g, np.nan),
                     eff_dim_old=ko.effective_dim.get(g, np.nan), eff_dim_new=kn.effective_dim.get(g, np.nan)))
C = pd.DataFrame(rows)
C.to_csv(Q.parent / "d106_compare.tsv", sep="\t", index=False, float_format="%.3f")
pd.set_option("display.width", 250)
print(C.round(3).to_string(index=False))
print("checks rows:", len(pd.read_csv(Q / "checks.tsv", sep="\t")), "| groups:", pd.read_csv(Q / "checks.tsv", sep="\t").group_id.nunique())
