#!/usr/bin/env python
"""TRACE-ARACNe stage 5: is the footprint stable ACROSS ARACNe networks?

D137 showed the footprint is a property of the network: CollecTRI and ARACNe nominate
nearly disjoint gene sets for the same target and give opposite answers. With CollecTRI
retired (user, 2026-11-17), the question that decides whether ARACNe-only is viable is
whether that instability was CollecTRI-vs-ARACNe, or is internal to ARACNe as well.

Method: rebuild the SAME footprint for the same target and the same drug x line response,
once per ARACNe cohort network, changing nothing but the network. Replay the identical L5
computation. If the answer swings across ARACNe networks as much as it swung between
CollecTRI and ARACNe, the footprint is not a property of the biology under any network and
ARACNe-only does not repair it.

Output: results/lincs/trace_aracne/network_stability.tsv
"""
import os, sys
from collections import defaultdict
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/lincs/trace_aracne"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "8"
HOP_CAP_NON_TF = 2
CASES = [
    dict(key="fulvestrant_ESR1", drug="fulvestrant", target="ESR1", ModelID="ACH-000019",
         line="MCF7", orient="lower", home="BRCA_LUMA"),
    dict(key="vemurafenib_BRAF", drug="vemurafenib", target="BRAF", ModelID="ACH-000219",
         line="A375", orient="lower", home="SKCM"),
]


def main():
    genes = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_genes.tsv", sep="\t")
    idx = pd.read_csv(ROOT / "data/processed/lincs/drug_line_profiles_index.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "results/lincs/drug_annotations.tsv", sep="\t")
    P = pd.read_parquet(ROOT / "results/lincs/d128/a1/target_tf_paths.parquet")
    P = P[P.network == "omnipath"]
    E = pd.read_csv(ROOT / "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz", sep="\t")
    E = E[E.in_lincs_gene_space]
    M = np.memmap(ROOT / "data/processed/lincs/drug_line_profiles_bing_specific.f32",
                  dtype="float32", mode="r", shape=(len(idx), len(genes)))
    gpos = {g: i for i, g in enumerate(genes.gene_id)}
    rows, footprints = [], {}

    for case in CASES:
        action = 1 if case["orient"] == "raise" else -1
        pid = set(annot[annot.cmap_name == case["drug"]].pert_id)
        pr = int(idx[(idx.ModelID == case["ModelID"]) & (idx.pert_id.isin(pid))]
                 .sort_values("n_sig", ascending=False).index[0])
        line_rows = idx.index[idx.ModelID == case["ModelID"]].to_numpy()
        pp_all = P[P.target_gene == case["target"]]
        for net, En in E.groupby("network"):
            regs = set(En.regulator)
            is_reg = case["target"] in regs
            cap = 0 if is_reg else HOP_CAP_NON_TF
            res = pp_all[(pp_all.hops <= cap) & (pp_all.path_status == "resolved")
                         & pp_all.tf.isin(regs)]
            if not len(res):
                rows.append(dict(case=case["key"], network=net, n_regulators=0, n_signed=0))
                continue
            acc = defaultdict(list)
            memb = {t: d for t, d in En[En.regulator.isin(set(res.tf))].groupby("regulator")}
            for t in res.itertuples():
                sub = memb.get(t.tf)
                if sub is None:
                    continue
                for g, m in zip(sub.gene_id, sub.tfmode):
                    acc[g].append(int(t.path_sign) * (1 if m >= 0 else -1))
            sign = {g: (v[0] * action if len(set(v)) == 1 else 0) for g, v in acc.items()}
            sign = {g: s for g, s in sign.items() if s != 0 and g in gpos}
            if not sign:
                rows.append(dict(case=case["key"], network=net, n_regulators=int(len(res)), n_signed=0))
                continue
            cols = np.array([gpos[g] for g in sign])
            e = np.array(list(sign.values()), dtype=float)
            z = np.asarray(M[pr][cols], dtype=float)
            mean = float((np.tanh(z / 2) * e).mean())
            Z = np.asarray(M[line_rows][:, cols], dtype=float)
            means = (np.tanh(Z / 2) * e).mean(axis=1)
            rows.append(dict(case=case["key"], network=net, home=net == case["home"],
                             n_regulators=int(len(res)), n_signed=len(sign),
                             concordance=round(float((np.sign(z) == e).mean()), 4),
                             mean_signed_fz=round(mean, 6),
                             generic_percentile=round(float((means < mean).mean() * 100), 1)))
            footprints[(case["key"], net)] = sign

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "network_stability.tsv", sep="\t", index=False)
    for case in CASES:
        d = df[(df.case == case["key"]) & (df.n_signed > 0)]
        print(f"\n=== {case['key']}  ({case['line']}), {len(d)} ARACNe networks with a footprint ===")
        print(f"  n_signed        median {d.n_signed.median():.0f}  range {d.n_signed.min()}-{d.n_signed.max()}")
        print(f"  concordance     median {100*d.concordance.median():.1f}%  "
              f"range {100*d.concordance.min():.1f}-{100*d.concordance.max():.1f}%")
        print(f"  mean signed f(z) median {d.mean_signed_fz.median():+.4f}  "
              f"range {d.mean_signed_fz.min():+.4f} to {d.mean_signed_fz.max():+.4f}  "
              f"-> POSITIVE in {int((d.mean_signed_fz > 0).sum())} of {len(d)} networks")
        print(f"  generic pct     median {d.generic_percentile.median():.1f}  "
              f"range {d.generic_percentile.min():.1f}-{d.generic_percentile.max():.1f}")
        h = d[d.home]
        if len(h):
            print(f"  HOME network {h.network.iloc[0]}: n={int(h.n_signed.iloc[0])}, "
                  f"conc {100*float(h.concordance.iloc[0]):.1f}%, mean {float(h.mean_signed_fz.iloc[0]):+.4f}, "
                  f"pct {float(h.generic_percentile.iloc[0]):.1f}")
        # pairwise footprint agreement between networks
        ks = [k for k in footprints if k[0] == case["key"]]
        jac, agree = [], []
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                a, b = footprints[ks[i]], footprints[ks[j]]
                sa, sb = set(a), set(b)
                inter = sa & sb
                jac.append(len(inter) / len(sa | sb))
                if inter:
                    agree.append(np.mean([a[g] == b[g] for g in inter]))
        print(f"  pairwise footprint Jaccard  median {np.median(jac):.3f}  "
              f"(10th-90th {np.percentile(jac,10):.3f}-{np.percentile(jac,90):.3f})")
        print(f"  sign agreement on shared genes median {100*np.median(agree):.1f}%")
    print("\nwrote", OUT / "network_stability.tsv")


if __name__ == "__main__":
    sys.exit(main())
