"""D128 shared graph code: signed routes from a target gene to regulators, and the gene footprint they imply.

Factored out of `a0_audit.py` (D128-A0 / A0.1), which proved it on the full input set; A1 and A2 import it rather
than restating it. Nothing here reads a drug response or any outcome: a footprint depends only on the annotated
target, the signed signaling prior and the regulon network.
"""
import pickle
from collections import defaultdict  # noqa: F401
from pathlib import Path

import numpy as np

import config as C

# Signs are bitmasks, not sets: 1 = {+1}, 2 = {-1}, 3 = both (ambiguous). The first version used Python sets inside
# the BFS and a pandas filter per reached TF, which made the loop single-threaded and compute-bound (see D128-A0
# result: 43 min without finishing). Masks, one cached BFS per gene and pre-filtered regulon gene lists replace it.
PLUS, MINUS, BOTH = 1, 2, 3


def mul(m1, m2):
    """Product of two sign masks."""
    if m1 == BOTH or m2 == BOTH:
        return BOTH
    return PLUS if m1 == m2 else MINUS


def build_graph(edges):
    g = {}
    for s, d, sign, disp in edges[["src", "dst", "sign", "sign_disputed"]].itertuples(index=False):
        g.setdefault(s, []).append((d, BOTH if bool(disp) else (PLUS if int(sign) == 1 else MINUS)))
    return g


def reach_regulators(graph, start, max_hops, regs):
    """BFS from start, keeping only regulators: {reg: (shortest hop, sign mask over all shortest paths)}.

    The mask is over *shortest* paths, so a walk at the largest hop arm also answers every smaller arm (filter on hop).
    """
    dist, mask = {start: 0}, {start: PLUS}
    out = {}
    if start in regs:
        out[start] = (0, PLUS)
    frontier = [start]
    for h in range(1, max_hops + 1):
        nxt = {}
        for u in frontier:
            mu = mask[u]
            for v, me in graph.get(u, ()):
                if v in dist and dist[v] < h:
                    continue
                s = mul(mu, me)
                nxt[v] = nxt.get(v, 0) | s
        for v, s in nxt.items():
            if v not in dist:
                dist[v], mask[v] = h, s
                if v in regs:
                    out[v] = (h, s)
            elif dist[v] == h:                      # reached again at the same (shortest) hop
                mask[v] |= s
                if v in regs:
                    out[v] = (h, mask[v])
        frontier = [v for v in nxt if dist.get(v) == h]
    return out


class Reach:
    """Per-gene reachable-regulator cache, persisted per (network, source checksum, max hops) as the plan asks."""

    def __init__(self, name, graph, regs, max_hops, checksum, cache_dir):
        self.name, self.graph, self.regs, self.max_hops = name, graph, regs, max_hops
        # Identity for the footprint memo: the NAME alone is not enough, because two different graphs may share one.
        self.key = f"{name}:{checksum[:12]}:h{max_hops}"
        self.path = Path(cache_dir) / f"reach_{name}_{checksum[:12]}_h{max_hops}.pkl"
        self.mem, self.hits, self.misses = {}, 0, 0
        if self.path.exists():
            self.mem = pickle.loads(self.path.read_bytes())

    def __call__(self, gene):
        r = self.mem.get(gene)
        if r is None:
            self.misses += 1
            r = self.mem[gene] = reach_regulators(self.graph, gene, self.max_hops, self.regs)
        else:
            self.hits += 1
        return r

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_bytes(pickle.dumps(self.mem, protocol=4))
        return self.path


def unit_footprint(genes, reach, reg_genes, max_hops, _cache={}):
    """Counts for one target unit. Expected gene sign up to the drug's action = path sign x tfmode.

    `reg_genes[tf] = (gene ids, tfmodes)`, already restricted to the LINCS gene space, so no pandas in this loop.
    Memoised on (network, unit gene set, hop arm): units share target genes across drugs.
    """
    key = (getattr(reach, "key", None) or reach.name, frozenset(genes), max_hops)
    if key in _cache:
        return _cache[key]
    tf_hop, tf_mask = {}, {}
    for g in genes:
        for t, (h, m) in reach(g).items():
            if h > max_hops:
                continue
            if t not in tf_hop or h < tf_hop[t]:
                tf_hop[t], tf_mask[t] = h, m
            elif h == tf_hop[t]:
                tf_mask[t] |= m
    gene_exp, gene_any, per_tf = {}, set(), {}
    for t, h in tf_hop.items():
        gl = reg_genes.get(t)
        if gl is None:
            continue
        ids, modes = gl
        gene_any.update(ids)
        if tf_mask[t] != BOTH:
            p = 1 if tf_mask[t] == PLUS else -1
            for gid, mode in zip(ids, modes):
                e = p * mode
                prev = gene_exp.get(gid)
                gene_exp[gid] = e if prev is None else (prev if prev == e else 0)   # 0 = conflicting expectations
            per_tf[t] = set(ids)
    resolved = {g for g, s in gene_exp.items() if s}
    conflicted = {g for g, s in gene_exp.items() if not s}
    hub = max((len(v & resolved) for v in per_tf.values()), default=0)
    hop0 = [t for t, h in tf_hop.items() if h == 0 and t in per_tf]
    hop0_resolved = set().union(*[per_tf[t] for t in hop0]) & resolved if hop0 else set()
    out = dict(
        gene_exp=gene_exp,
        tf_hop=tf_hop, n_tfs={h: sum(1 for x in tf_hop.values() if x <= h) for h in C.HOP_ARMS},
        n_anchored_tfs=len(tf_hop), n_resolved_tfs=sum(1 for m in tf_mask.values() if m != BOTH),
        n_ambiguous_tfs=sum(1 for m in tf_mask.values() if m == BOTH),
        genes_any=gene_any, genes_resolved=resolved, genes_conflicted=conflicted, hop0_resolved=hop0_resolved,
        hub_fraction=(hub / len(resolved)) if resolved else np.nan,
        hub_tf=max(per_tf, key=lambda t: len(per_tf[t] & resolved)) if per_tf and resolved else "")
    _cache[key] = out
    return out


def unit_subsets(genes, reach, reg_genes, regs, max_hops):
    """D128-A0.1: type and floor PER GENE, then take the union.

    A unit gene that is itself a CollecTRI regulator anchors its own regulon at zero hops (the `tf` subset); a gene
    that is not anchors through the signed signaling prior at 0-`max_hops` (the `signaling` subset). The earlier rule
    typed the whole unit `tf` as soon as any one gene was a regulator and then floored the unit on zero-hop genes
    alone, which floored palbociclib (CCND1|CDK4) at 0 zero-hop genes although CDK4 resolved 1,758 through signaling.
    """
    tf_genes = [g for g in genes if g in regs]
    sig_genes = [g for g in genes if g not in regs]
    tf_fp = unit_footprint(tf_genes, reach, reg_genes, 0) if tf_genes else None
    sig_fp = unit_footprint(sig_genes, reach, reg_genes, max_hops) if sig_genes else None
    kind = "mixed" if tf_genes and sig_genes else ("tf" if tf_genes else "signaling")
    return kind, tf_fp, sig_fp


def merge_subsets(tf_fp, sig_fp):
    """Union of the two subsets, each unique gene counted once (plan 3.4); disagreement between them is conflicting."""
    exp = {}
    for fp in (tf_fp, sig_fp):
        if fp is None:
            continue
        for g, e in fp["gene_exp"].items():
            exp[g] = e if g not in exp else (exp[g] if exp[g] == e else 0)
    resolved = {g for g, e in exp.items() if e}
    any_genes = set().union(*[fp["genes_any"] for fp in (tf_fp, sig_fp) if fp is not None]) if (tf_fp or sig_fp) else set()
    return any_genes, resolved, {g for g, e in exp.items() if not e}


def floors(kind, n_res, n_res_tfs, n_hop0, mult):
    """Floor for one subset: `tf` on its own regulon genes, `signaling` on resolved genes and anchored TFs."""
    fg, ft = C.FLOOR_GENES * mult, C.FLOOR_TFS * mult
    if kind == "tf":
        return n_hop0 >= fg
    return n_res >= fg and n_res_tfs >= ft


def unit_floor(tf_fp, sig_fp, mult, space=None):
    """D128-A0.1: the unit clears if EITHER subset clears its own floor. Returns (passes, which)."""
    def n(fp, key):
        if fp is None:
            return 0
        s = fp[key]
        return len(s & space) if space is not None else len(s)
    tf_ok = tf_fp is not None and floors("tf", 0, 0, n(tf_fp, "genes_resolved"), mult)
    sig_ok = sig_fp is not None and floors("signaling", n(sig_fp, "genes_resolved"), sig_fp["n_resolved_tfs"], 0, mult)
    which = "both" if (tf_ok and sig_ok) else ("tf" if tf_ok else ("signaling" if sig_ok else ""))
    return (tf_ok or sig_ok), which


