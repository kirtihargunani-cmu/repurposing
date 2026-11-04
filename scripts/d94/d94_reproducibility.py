#!/usr/bin/env python
"""D94 network uncertainty (rule fixed 2026-11-03): Jaccard of the edge sets consolidated from bootstrap halves 1-50 and
51-100 (scripts/d94/d94_reproducibility.sh). Cohort gate: median reproducibility of the cohort's cell-type networks
(per arm) >= 80% of its pooled network's. Output results/lincs/d94/reproducibility.tsv, cohort_gate.tsv"""
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
H = ROOT / "data/processed/step5/aracne_v2_tf_halves"; OUT = ROOT / "results/lincs/d94"; OUT.mkdir(parents=True, exist_ok=True)
def edges(f):
    d = pd.read_csv(f, sep="\t", usecols=[0, 1]); return set(zip(d.iloc[:, 0], d.iloc[:, 1]))
rows = []
for d in sorted(H.iterdir()):
    a, b = edges(d / "h1/network.txt"), edges(d / "h2/network.txt")
    full = sum(1 for _ in open(ROOT / f"data/processed/step5/aracne_v2_tf/{d.name}/network.txt")) - 1
    rows.append(dict(net=d.name, edges_full=full, edges_h1=len(a), edges_h2=len(b), jaccard=round(len(a & b) / len(a | b), 4)))
R = pd.DataFrame(rows); R.to_csv(OUT / "reproducibility.tsv", sep="\t", index=False)
gate = []
for arm in ["lum", "split"]:
    m = pd.read_csv(ROOT / f"results/step5/d94_{arm}_group_net.tsv", sep="\t")
    m["cohort"] = m.group_id.str.split("-").str[0]
    for c, g in m.groupby("cohort"):
        nets = sorted(set(g.net)); cj = R.set_index("net").jaccard.reindex(nets)
        pj = R.set_index("net").jaccard.get(c)
        gate.append(dict(arm=arm, cohort=c, networks=";".join(nets), median_cell_type=round(cj.median(), 4), pooled=pj,
                         ratio=round(cj.median() / pj, 3), passes=bool(cj.median() >= 0.8 * pj)))
G = pd.DataFrame(gate); G.to_csv(OUT / "cohort_gate.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); print(R.to_string(index=False)); print(); print(G.to_string(index=False))
