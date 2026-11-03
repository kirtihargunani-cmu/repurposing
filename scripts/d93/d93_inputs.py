#!/usr/bin/env python
"""D93 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): ARACNe inputs for cell-type networks in the 8 split
cohorts. Each of the 41 groups in BRCA, CESC, ESCA, LGG, PCPG, SARC, TGCT, UCEC with > 50 tumors gets its own network;
groups with <= 50 tumors join the network (among those) whose tumor centroid (log2 TPM + 1 on the genes shared by the
cohort inputs, centred over all their tumors) correlates best with theirs, and that network is built with them.
Samples = the group's tumors in the cohort input (data/processed/step5/aracne_inputs_v2/<COHORT>.exp.tsv); inputs as
scripts/42 (v2 gene set, TPM >= 1 in >= 10% of the network's tumors, nonzero variance, TF-only regulator list).
Networks that already exist with exactly the same tumors (BRCA_BAS, BRCA_HER2 from D92a) are reused.
Outputs: data/processed/step5/aracne_inputs_v2/<NET>.*, results/step5/d93_networks.tsv, results/step5/d93_group_net.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
IN = ROOT / "data/processed/step5/aracne_inputs_v2"
COH = ["BRCA", "CESC", "ESCA", "LGG", "PCPG", "SARC", "TGCT", "UCEC"]
MIN_N = 50
ok = set(pd.read_csv(ROOT / "results/step4_v2/signature_summary.tsv", sep="\t").query("status == 'ok'").group_id)
cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
sub_of = cov.drop_duplicates("sample").set_index("sample").group_id
X = {c: pd.read_csv(IN / f"{c}.exp.tsv", sep="\t", index_col=0) for c in COH}
samples = {}
for c, x in X.items():
    g = sub_of.reindex(x.columns)
    for gid in sorted(set(g.dropna()) & ok):
        samples[gid] = (c, list(x.columns[(g == gid).values]))
name = lambda gid: gid.replace("BRCA-LumA", "BRCA_LUMA").replace("BRCA-LumB", "BRCA_LUMB").replace("BRCA-Basal", "BRCA_BAS") \
    .replace("BRCA-Her2", "BRCA_HER2").replace("-", "_").upper()
big = {gid: v for gid, v in samples.items() if len(v[1]) > MIN_N}
small = {gid: v for gid, v in samples.items() if len(v[1]) <= MIN_N}
genes = sorted(set.intersection(*(set(x.index) for x in X.values())))
allx = pd.concat([X[c].loc[genes, s] for c, s in samples.values()], axis=1)
allx = allx.sub(allx.mean(axis=1), axis=0)
cen = {gid: allx[s].mean(axis=1) for gid, (c, s) in samples.items()}
members = {gid: [gid] for gid in big}
assign = []
for gid in small:
    r = {b: float(np.corrcoef(cen[gid], cen[b])[0, 1]) for b in big}
    best = max(r, key=r.get); members[best].append(gid)
    assign.append(dict(group_id=gid, n=len(samples[gid][1]), joins=best, r_best=round(r[best], 3),
                       runner_up=sorted(r, key=r.get)[-2], r_runner_up=round(sorted(r.values())[-2], 3)))
    print(f"{gid} ({len(samples[gid][1])}) joins {best} (r = {r[best]:.3f}; next {sorted(r, key=r.get)[-2]} {sorted(r.values())[-2]:.3f})")
pc = set(pd.read_csv(ROOT / "results/step3/genes_protein_coding_v23_v2.tsv", sep="\t").gene_id)
regs = pd.read_csv(ROOT / "results/step5/aracne_regulators_v2.tsv", sep="\t")
need = {}
for b, gs in members.items():
    net = name(b) + ("_X" if len(gs) > 1 else "")
    samp = [s for g in gs for s in samples[g][1]]
    f = IN / f"{net}.exp.tsv"
    if f.exists() and sorted(pd.read_csv(f, sep="\t", nrows=0).columns[1:]) == sorted(samp):
        need[net] = (gs, samp, "reused")
    else:
        need[net] = (gs, samp, "build")
todo = sorted({s for n, (gs, samp, st) in need.items() if st == "build" for s in samp})
tpm_f = ROOT / "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
header = pd.read_csv(tpm_f, sep="\t", nrows=0).columns
chunks = []
for ch in pd.read_csv(tpm_f, sep="\t", usecols=[header[0]] + todo, chunksize=5000, index_col=0):
    ch.index = ch.index.str.split(".").str[0]; chunks.append(ch[ch.index.isin(pc)])
x = np.log2(np.power(2.0, pd.concat(chunks)) - 0.001 + 1).clip(lower=0)
rows, gmap = [], []
for net, (gs, samp, st) in need.items():
    if st == "build":
        m = x[samp]; m = m[((m >= 1).mean(axis=1) >= 0.10) & (m.std(axis=1) > 0)]
        m.round(4).to_csv(IN / f"{net}.exp.tsv", sep="\t", index_label="gene")
        r = regs[regs.gene_id.isin(m.index)]
        r.gene_id.to_csv(IN / f"{net}.regulators.txt", index=False, header=False)
        r.loc[r.is_lambert_tf, "gene_id"].to_csv(IN / f"{net}.regulators_tf.txt", index=False, header=False)
        ng = m.shape[0]
    else:
        ng = sum(1 for _ in open(IN / f"{net}.exp.tsv")) - 1
    rows.append(dict(net=net, groups=";".join(gs), n_tumors=len(samp), n_genes=ng, status=st))
    gmap += [(g, net) for g in gs]
nets = pd.DataFrame(rows); nets.to_csv(ROOT / "results/step5/d93_networks.tsv", sep="\t", index=False)
pd.DataFrame(assign).to_csv(ROOT / "results/step5/d93_small_groups.tsv", sep="\t", index=False)
pd.DataFrame(gmap, columns=["group_id", "net"]).to_csv(ROOT / "results/step5/d93_group_net.tsv", sep="\t", index=False)
print(nets.to_string(index=False))
