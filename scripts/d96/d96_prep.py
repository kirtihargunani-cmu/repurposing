#!/usr/bin/env python
"""D96 (rule fixed 2026-11-03, docs/DECISIONS.md): map and override folders for the 85-tumor arm.
Groups >= 85 tumors in the 8 split cohorts on their own networks; Step 6 rows from the D93 run (same networks and seeds)
except PCPG-PHEO / UCEC-EEC (own-tumor networks, run _d96new). Subgroups: LGG-noncodel rows from strata_d94, PHEO rows
from strata_d96. Outputs results/step5/d96_group_net.tsv, results/step6/aracne_project_v2_tf_d96[ / _strata_d96x]/."""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent; S6 = ROOT / "results/step6"
M = {"BRCA-LumA": "BRCA_LUMA", "BRCA-LumB": "BRCA_LUMB", "BRCA-Basal": "BRCA_BAS", "CESC-SCC": "CESC_SCC", "ESCA-AD": "ESCA_AD",
     "ESCA-SCC": "ESCA_SCC", "LGG-IDHmut-codel": "LGG_IDHMUT_CODEL", "LGG-IDHmut-noncodel": "LGG_IDHMUT_NONCODEL", "LGG-IDHwt": "LGG_IDHWT",
     "PCPG-PHEO": "PCPG_PHEO", "SARC-LMS": "SARC_LMS", "UCEC-EEC": "UCEC_EEC"}
pd.DataFrame({"group_id": list(M), "net": list(M.values())}).to_csv(ROOT / "results/step5/d96_group_net.tsv", sep="\t", index=False)
new = ["PCPG-PHEO", "UCEC-EEC"]
od = S6 / "aracne_project_v2_tf_d96"; od.mkdir(exist_ok=True)
for f in ["group_mr.tsv", "regulator_filter.tsv"]:
    a = pd.read_csv(S6 / "aracne_project_v2_tf_d93" / f, sep="\t"); b = pd.read_csv(S6 / "aracne_project_v2_tf_d96new" / f, sep="\t")
    pd.concat([a[a.group_id.isin(M) & ~a.group_id.isin(new)], b[b.group_id.isin(new)]]).to_csv(od / f, sep="\t", index=False)
sd = S6 / "aracne_project_v2_tf_strata_d96x"; sd.mkdir(exist_ok=True)
a = pd.read_csv(S6 / "aracne_project_v2_tf_strata_d94/group_mr.tsv", sep="\t"); b = pd.read_csv(S6 / "aracne_project_v2_tf_strata_d96/group_mr.tsv", sep="\t")
pd.concat([a[a.group_id.str.startswith("LGG")], b[b.group_id.str.startswith("PHEO")]]).to_csv(sd / "group_mr.tsv", sep="\t", index=False)
p1 = pd.concat([pd.read_csv(ROOT / "results/lincs/prolif_regulators_v2_tf_d93extra.tsv", sep="\t"),
                pd.read_csv(ROOT / "results/lincs/prolif_regulators_v2_tf_d96extra.tsv", sep="\t")])
p1.to_csv(ROOT / "results/lincs/prolif_regulators_v2_tf_d96all.tsv", sep="\t", index=False)
print("d96 groups:", len(M))
