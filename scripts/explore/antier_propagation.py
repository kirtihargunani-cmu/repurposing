"""Exploratory (2026-11-03): how the anti-estrogen misread propagates into the reversal score. Anti-estrogens in MCF7
(BRCA tumor ARACNe v2 TF network): per regulator contribution -w_i a_i / ||w|| to rev for BRCA-LumA / LumB / Basal
(w = tumor-vs-normal NES, all tested regulators, as 79e SCORE=weighted), split into proliferation-linked regulators
(79h P1 + 79i P2, BRCA) and the rest; for the early (6 h), 24 h, late-only and aggregated signatures (D89)."""
import numpy as np, pandas as pd
P = "data/processed/lincs/"
I = pd.read_csv(P + "timedose/time_index.tsv", sep="\t"); er = open(P + "timedose/early_late/BRCA.regulators.txt").read().split()
ae = ["tamoxifen", "fulvestrant", "raloxifene", "toremifene", "bazedoxifene", "clomifene"]
k = np.where((I.ModelID == "ACH-000019") & I.cmap_name.isin(ae))[0]
S = {s: pd.Series(np.fromfile(P + f"timedose/early_late/{s}_BRCA.f32", dtype=np.float32).reshape(len(I), len(er))[k].mean(0), index=er) for s in ["E", "L", "R", "A0"]}
mr = pd.read_csv("results/step6/aracne_project_v2_tf/group_mr.tsv", sep="\t").pivot(index="regulator", columns="group_id", values="nes")
pl = pd.read_csv("results/lincs/prolif_regulators_v2_tf.tsv", sep="\t"); p2 = pd.read_csv("results/lincs/prolif_regulators_p2_v2_tf_d85.tsv", sep="\t")
prol = set(pl[(pl.net == "BRCA") & pl.prolif_linked].regulator) | set(p2[(p2.net == "BRCA") & p2.p2].regulator)
rows = []
for s, lab in [("E", "early 6 h"), ("L", "24 h"), ("R", "late-only"), ("A0", "aggregated (current)")]:
    for g in ["BRCA-LumA", "BRCA-LumB", "BRCA-Basal"]:
        m = mr[g].dropna(); c = m.index.intersection(S[s].index); w = m[c]; x = S[s][c]
        cont = -(w * x) / np.sqrt((w ** 2).sum()); isp = c.isin(prol)
        rows.append((lab, g, round(cont.sum(), 2), round(cont[isp].sum(), 2), round(cont[~isp].sum(), 2), round(np.corrcoef(x, w)[0, 1], 3),
                     ", ".join(f"{r} {cont[r]:+.2f}" for r in cont[~isp].sort_values().index[:5])))
o = pd.DataFrame(rows, columns=["signature", "group", "rev", "rev_prolif_regs", "rev_other_regs", "r_drug_vs_tumor", "most_negative_other"])
o.to_csv("results/lincs/early_late/antier_propagation.tsv", sep="\t", index=False); print(o.to_string(index=False))
