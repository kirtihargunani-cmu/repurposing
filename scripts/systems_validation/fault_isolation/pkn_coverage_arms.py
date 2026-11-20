#!/usr/bin/env python
"""Would a bigger prior-knowledge network actually help? (user, 2026-11-20, before pricing IPA)

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

D143 recorded that only ~33% of tested ARACNe regulators are nodes in our signed PKN and called the
gap a study bias. That was measured on the TF-only set. On the FULL set actually scored in D158
(4,233 regulators over 7 networks) coverage is 54.1%, and 84.6% of the missing regulators are named,
studied genes - so the gap looks like our strict signed-consensus filter and OmniPath's signalling
scope, not absence of literature. That is the case for a larger knowledge base such as IPA's.

Before paying for one, this measures what a larger PKN buys, using sources we already hold:

  curated        the production PKN: omnipath_signed_directed (consensus sign+direction,
                 >= 2 resources OR >= 2 PubMed ids)
  relaxed        omnipath_raw, consensus direction and exactly one consensus sign, NO evidence floor
  permissive     omnipath_raw, directed, sign from is_stimulation / is_inhibition (non-consensus ok)
  curated+signor union of the curated arm with the SIGNOR arm
  relaxed+signor union of the relaxed arm with SIGNOR

Reported per arm: edges, nodes, % activating, coverage of the 4,233 scored regulators, coverage of
the 387 annotated drug-target genes, and - the number that actually matters - how many
(drug target -> scored regulator) pairs are connected by a directed path of <= 3 hops. Coverage of
nodes is necessary but not sufficient; a target that cannot REACH a disease regulator is untraceable
however well curated it is.

The activation balance is reported alongside, because D145 found the ARACNe causal map returns an
EMPTY CARNIVAL solution at 96-97% activating edges. A larger PKN that is also more activation-biased
may buy coverage it cannot use.

Output: results/systems_validation/pkn/pkn_coverage_arms.{tsv,json}
"""
import json, os
from pathlib import Path
import numpy as np, pandas as pd
import scipy.sparse as sp

ROOT = Path("/home/kirti/projects/repurposing")
RAW = ROOT / "data/raw/omnipath"
OUT = ROOT / "results/systems_validation/pkn"
VD = ROOT / "data/processed/lincs/drug_viper_specific_v2_condfull"
MAX_HOPS = 3
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "4"


def load_raw(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d.rename(columns={"source": "source_uniprot", "target": "target_uniprot",
                          "source_genesymbol": "source", "target_genesymbol": "target"})
    d = d[d.source.notna() & d.target.notna()]
    d = d[d.source != d.target]
    d = d[~d.source.astype(str).str.contains("_") & ~d.target.astype(str).str.contains("_")]
    return d


def sign_consensus(d):
    m = d.consensus_direction.astype(bool) & (
        d.consensus_stimulation.astype(bool) ^ d.consensus_inhibition.astype(bool))
    x = d[m].copy()
    x["interaction"] = np.where(x.consensus_stimulation.astype(bool), 1, -1)
    return x[["source", "target", "interaction"]]


def sign_permissive(d):
    m = d.is_directed.astype(bool) & (d.is_stimulation.astype(bool) ^ d.is_inhibition.astype(bool))
    x = d[m].copy()
    x["interaction"] = np.where(x.is_stimulation.astype(bool), 1, -1)
    return x[["source", "target", "interaction"]]


def dedup(e):
    return e.drop_duplicates(["source", "target"]).reset_index(drop=True)


def reach(edges, srcs, dsts, max_hops):
    """how many (src, dst) pairs are joined by a directed path of <= max_hops."""
    nodes = pd.Index(pd.unique(pd.concat([edges.source, edges.target])))
    pos = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    A = sp.csr_matrix((np.ones(len(edges), np.int8),
                       (edges.source.map(pos).values, edges.target.map(pos).values)), shape=(n, n))
    si = [pos[s] for s in srcs if s in pos]
    di = np.array([pos[t] for t in dsts if t in pos])
    if not si or not len(di):
        return 0, 0, 0
    F = sp.csr_matrix((np.ones(len(si), np.int8), (np.arange(len(si)), si)), shape=(len(si), n))
    seen = F.copy()
    R = F.copy()
    for _ in range(max_hops):
        R = (R @ A)
        R.data[:] = 1
        seen = seen + R
        seen.data[:] = 1
    sub = seen[:, di]
    per_src = np.asarray((sub > 0).sum(1)).ravel()
    return int(per_src.sum()), int((per_src > 0).sum()), len(si)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    regs = set()
    for f in VD.glob("*.regulators.txt"):
        regs |= set(f.read_text().split())
    u = pd.read_csv(ROOT / "results/lincs/target_role/d121/units.tsv", sep="\t",
                    dtype={"tid": str, "action_flag": str})
    tgts = {g for x in u.genes.dropna() for g in str(x).split("|")}
    print(f"scored ARACNe regulators: {len(regs):,}   annotated drug-target genes: {len(tgts):,}")
    print(f"possible target x regulator pairs: {len(tgts) * len(regs):,}\n")

    op = load_raw(RAW / "omnipath_raw.tsv.gz")
    # what does each gate of the production filter cost? (user question, 2026-11-20)
    def nres(x):
        return x.fillna("").astype(str).map(lambda v: len({y for y in v.split(";") if y}))
    o2 = op.copy()
    o2["n_resources"] = nres(o2.sources); o2["n_references"] = nres(o2.references)
    g_dir = o2.consensus_direction.astype(bool)
    g_sgn = o2.consensus_stimulation.astype(bool) ^ o2.consensus_inhibition.astype(bool)
    g_ev = (o2.n_resources >= 2) | (o2.n_references >= 2)
    tot = len(o2)
    print("=== what each gate of the production filter costs (OmniPath raw) ===")
    print(f"  raw interactions (genes only, no self-loops/complexes) : {tot:>8,}")
    print(f"  gate 1  consensus_direction                            : {int(g_dir.sum()):>8,} kept "
          f"({100*g_dir.mean():>4.1f}%)")
    print(f"  gate 2  + exactly one consensus sign                   : {int((g_dir&g_sgn).sum()):>8,} kept "
          f"({100*(g_dir&g_sgn).mean():>4.1f}%)")
    print(f"  gate 3  + >=2 resources OR >=2 PubMed ids              : {int((g_dir&g_sgn&g_ev).sum()):>8,} kept "
          f"({100*(g_dir&g_sgn&g_ev).mean():>4.1f}%)")
    drop3 = (g_dir & g_sgn & ~g_ev)
    print(f"\n  gate 3 alone DELETES {int(drop3.sum()):,} signed, directed edges "
          f"({100*drop3.sum()/max((g_dir&g_sgn).sum(),1):.1f}% of what survives gates 1-2)")
    lost = set(o2[drop3].source) | set(o2[drop3].target)
    kept = set(o2[g_dir & g_sgn & g_ev].source) | set(o2[g_dir & g_sgn & g_ev].target)
    print(f"  nodes ONLY reachable via gate-3-deleted edges: {len(lost - kept):,}")
    print(f"  of those, ARACNe regulators we actually score: {len((lost - kept) & regs):,}")
    print(f"  of those, annotated drug-target genes        : {len((lost - kept) & tgts):,}\n")
    sg = load_raw(RAW / "signor_raw.tsv.gz")
    cur = pd.read_csv(OUT / "signed_pkn.tsv", sep="\t")[["source", "interaction", "target"]]
    sg_c = sign_consensus(sg)

    arms = {
        "curated (production)": dedup(cur),
        "relaxed": dedup(sign_consensus(op)),
        "permissive": dedup(sign_permissive(op)),
        "curated+signor": dedup(pd.concat([cur, sg_c], ignore_index=True)),
        "relaxed+signor": dedup(pd.concat([sign_consensus(op), sg_c], ignore_index=True)),
    }

    rows = []
    hdr = (f"{'arm':<22}{'edges':>9}{'nodes':>8}{'%act':>7}{'regulators':>12}{'targets':>10}"
           f"{'reachable pairs':>17}{'targets w/ >=1':>15}")
    print(hdr); print("-" * len(hdr))
    for name, e in arms.items():
        nodes = set(e.source) | set(e.target)
        act = 100 * float((e.interaction == 1).mean())
        rc, tc = len(regs & nodes), len(tgts & nodes)
        pairs, tgt_any, n_src = reach(e, sorted(tgts), sorted(regs), MAX_HOPS)
        rows.append(dict(arm=name, edges=len(e), nodes=len(nodes), pct_activating=round(act, 1),
                         regulators_covered=rc, regulators_pct=round(100 * rc / len(regs), 1),
                         targets_covered=tc, targets_pct=round(100 * tc / len(tgts), 1),
                         reachable_pairs=pairs, targets_with_any=tgt_any, max_hops=MAX_HOPS))
        print(f"{name:<22}{len(e):>9,}{len(nodes):>8,}{act:>6.1f}%"
              f"{rc:>8,} {100*rc/len(regs):>4.1f}%{tc:>6,} {100*tc/len(tgts):>3.0f}%"
              f"{pairs:>17,}{tgt_any:>12,}/{n_src}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "pkn_coverage_arms.tsv", sep="\t", index=False)
    json.dump(rows, open(OUT / "pkn_coverage_arms.json", "w"), indent=2)
    base = df.iloc[0]
    print(f"\nvs the production arm:")
    for _, r in df.iloc[1:].iterrows():
        print(f"  {r.arm:<22} edges x{r.edges/base.edges:>4.1f}  regulators +{r.regulators_covered-base.regulators_covered:>5,} "
              f"({r.regulators_pct - base.regulators_pct:+.1f} pts)  reachable pairs x{r.reachable_pairs/max(base.reachable_pairs,1):>4.1f}"
              f"  activating {r.pct_activating - base.pct_activating:+.1f} pts")
    print(f"\n-> {OUT/'pkn_coverage_arms.tsv'}")


if __name__ == "__main__":
    main()
