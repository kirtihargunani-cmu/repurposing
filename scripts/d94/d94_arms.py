#!/usr/bin/env python
"""D94 (rule fixed 2026-11-03, docs/DECISIONS.md): inputs of the LUM and SPLIT arms.
SPLIT = the D93 map and Step 6 results (results/step6/aracne_project_v2_tf_d93). LUM = the same with BRCA-LumA / -LumB on
BRCA_LUM, their Step 6 rows from the D92a run (results/step6/aracne_project_v2_tf_d92b; same seed rule).
Outputs: results/step5/d94_{lum,split}_group_net.tsv, results/step6/aracne_project_v2_tf_d94{lum,split}/{group_mr,regulator_filter}.tsv
"""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
S6 = ROOT / "results/step6"
m = pd.read_csv(ROOT / "results/step5/d93_group_net.tsv", sep="\t")
lum = m.copy(); lum.loc[lum.group_id.isin(["BRCA-LumA", "BRCA-LumB"]), "net"] = "BRCA_LUM"
m.to_csv(ROOT / "results/step5/d94_split_group_net.tsv", sep="\t", index=False)
lum.to_csv(ROOT / "results/step5/d94_lum_group_net.tsv", sep="\t", index=False)
for arm in ["split", "lum"]:
    d = S6 / f"aracne_project_v2_tf_d94{arm}"; d.mkdir(exist_ok=True)
    for f in ["group_mr.tsv", "regulator_filter.tsv"]:
        a = pd.read_csv(S6 / "aracne_project_v2_tf_d93" / f, sep="\t")
        if arm == "lum":
            b = pd.read_csv(S6 / "aracne_project_v2_tf_d92b" / f, sep="\t")
            lu = ["BRCA-LumA", "BRCA-LumB"]
            a = pd.concat([a[~a.group_id.isin(lu)], b[b.group_id.isin(lu)]], ignore_index=True)
        a.to_csv(d / f, sep="\t", index=False)
    print(arm, "groups:", a.group_id.nunique())
