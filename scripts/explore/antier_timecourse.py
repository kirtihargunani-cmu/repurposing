"""Exploratory (2026-11-03): anti-estrogens in MCF7 by time point: estrogen-response vs luminal-identity genes (level-5 z, hiq)."""
import numpy as np, pandas as pd
P = "data/processed/lincs/"
s = pd.read_csv(P + "level5_bing_hiq_sigs.tsv.gz", sep="\t"); g = pd.read_csv(P + "level5_bing_genes.tsv", sep="\t")
X = np.load(P + "level5_bing_hiq.npy", mmap_mode="r")
ae = ["fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene", "ZK-164015", "ICI-182780"]
m = s[(s.cell_iname == "MCF7") & s.cmap_name.isin(ae)]
canon = ["GREB1", "MYB", "XBP1", "PGR", "TFF1", "CA12", "SIAH2", "PDZK1", "CCND1", "IGFBP4", "NRIP1"]
lum = ["GATA3", "FOXA1", "ERBB4", "ANXA9", "HEBP1", "MINDY1", "SCCPDH", "SPDEF", "TBX3", "AR"]
gi = {x: i for i, x in enumerate(g.gene_symbol)}
for t, d in m.groupby("pert_time"):
    sub = np.asarray(X[d.index.values])
    for name, gs in [("estrogen-response", canon), ("luminal identity", lum)]:
        gs = [x for x in gs if x in gi]; v = sub[:, [gi[x] for x in gs]]
        print(f"{t:>4}h n={len(d):3d} {name:18s} mean z {v.mean():+.2f}  " + " ".join(f"{x}{v[:, k].mean():+.1f}" for k, x in enumerate(gs)))
