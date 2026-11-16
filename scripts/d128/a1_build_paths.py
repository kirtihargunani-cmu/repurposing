#!/usr/bin/env python
"""D128 A1: signed route index from drug target units to regulators to genes (plan Section 4 A1; D128-A1).

Drug-agnostic by construction: every table here is a function of the annotated target genes, the signed signaling
prior (D126 OmniPath cache, SIGNOR as a sensitivity arm) and the regulon network (CollecTRI, P3). No drug response,
no disease contrast and no outcome label is read, so nothing in A1 can depend on the direction a drug happens to
move. Expected gene signs are given UP TO the drug's action: the action sign is applied later, in A2.

Outputs results/lincs/d128/a1/:
  regulon_edges.parquet          network, regulator, gene_id, sign, confidence, source
  target_tf_paths.parquet        target_gene, network, tf, hops, path_sign, path_status (resolved / ambiguous)
  target_gene_footprints.parquet unit_id, gene_id, subset (tf / signaling / both), expected_sign_pre_action,
                                 status, feature_space
  path_ambiguity.tsv             per target gene and unit: ambiguity and hub concentration
  a1_run.json                    manifest: git SHA, input checksums, counts, timings
"""
import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import config as C  # noqa: E402
import d128io as IO  # noqa: E402
from graph import BOTH, MINUS, PLUS, Reach, merge_subsets, unit_subsets  # noqa: E402

A1 = C.OUT / "a1"
T0 = time.time()
LOG = []


def log(msg):
    line = f"[{time.time() - T0:7.1f}s] {msg}"
    LOG.append(line); print(line, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=C.THREADS)
    a = ap.parse_args()
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(a.threads)
    A1.mkdir(parents=True, exist_ok=True)

    genes = IO.profile_genes()
    gene_space = set(genes.gene_id)
    land = set(genes.gene_id[genes.feature_space == "landmark"])
    col = IO.collectri(C.COLLECTRI)
    regs = set(col.regulator)
    ck = col[col.target.isin(gene_space)]
    reg_genes = {r: (tuple(d.target), tuple(d.tfmode.astype(int))) for r, d in ck.groupby("regulator")}

    # --- regulon_edges.parquet: the signed, weighted TF -> gene relation (immutable between contexts) --------------
    E = col.rename(columns={"target": "gene_id", "tfmode": "sign", "likelihood": "confidence"}).copy()
    E["network"] = C.NETWORK
    E["source"] = "collectri_regulons_core"
    E["in_lincs_gene_space"] = E.gene_id.isin(gene_space)
    E["feature_space"] = np.where(E.gene_id.isin(land), "landmark",
                                  np.where(E.in_lincs_gene_space, "best inferred", "not measured"))
    E[["network", "regulator", "gene_id", "sign", "confidence", "source", "in_lincs_gene_space", "feature_space"]] \
        .to_parquet(A1 / "regulon_edges.parquet", index=False)
    log(f"regulon_edges: {len(E)} edges, {E.regulator.nunique()} regulators, "
        f"{int(E.in_lincs_gene_space.sum())} in the LINCS gene space")

    # --- target_tf_paths.parquet: every target gene of every unit, both signaling priors ---------------------------
    U = pd.read_csv(C.UNITS, sep="\t", dtype={"action_flag": str, "tid": str}).fillna({"action_flag": ""})
    U["genes_str"] = U.genes
    U["genes"] = U.genes.str.split("|")
    U["unit_id"] = U.cmap_name + "::" + U.tid.fillna("annot:" + U.genes_str)
    dup = U.unit_id.duplicated(keep=False)
    if dup.any():
        U.loc[dup, "unit_id"] = U.unit_id[dup] + "#" + U[dup].groupby("unit_id").cumcount().add(1).astype(str)
    targets = sorted({g for gs in U.genes for g in gs})

    REACH = {}
    for name, path in (("omnipath", C.OMNIPATH), ("signor", C.SIGNOR)):
        e = IO.signaling(path)
        from graph import build_graph
        REACH[name] = Reach(name, build_graph(e), regs, max(C.HOP_ARMS), IO.sha256(path), C.OUT / "a0/cache")

    rows = []
    for t in targets:
        for name, R in REACH.items():
            for tf, (h, m) in R(t).items():
                rows.append((t, name, tf, h, {PLUS: 1, MINUS: -1, BOTH: 0}[m],
                             "ambiguous" if m == BOTH else "resolved", tf in reg_genes,
                             len(reg_genes.get(tf, ((), ()))[0])))
    P = pd.DataFrame(rows, columns=["target_gene", "network", "tf", "hops", "path_sign", "path_status",
                                    "tf_has_regulon", "tf_regulon_genes"])
    P.to_parquet(A1 / "target_tf_paths.parquet", index=False)
    om = P[P.network == "omnipath"]
    log(f"target_tf_paths: {len(P)} rows for {len(targets)} target genes; omnipath {len(om)} "
        f"({int((om.path_status == 'ambiguous').sum())} ambiguous, "
        f"{int((om.hops == 0).sum())} zero-hop anchors)")

    # --- target_gene_footprints.parquet: expected sign per (unit, gene), before the drug's action ------------------
    fprows, amb = [], []
    for r in U.itertuples(index=False):
        kind, tf_fp, sig_fp = unit_subsets(r.genes, REACH["omnipath"], reg_genes, regs, C.MAX_HOPS)
        g_any, g_res, g_conf = merge_subsets(tf_fp, sig_fp)
        exp = {}
        for sub, fp in (("tf", tf_fp), ("signaling", sig_fp)):
            if fp is None:
                continue
            for g, e in fp["gene_exp"].items():
                if g in exp and exp[g][0] != e:
                    exp[g] = (0, exp[g][1] + "+" + sub)
                elif g in exp:
                    # D128-A2.1: a gene reached by BOTH subsets with the SAME expected sign is labelled `both`,
                    # not by whichever subset was walked first. The label is what the engine floors on (P10 as
                    # frozen by D128-A0.1 counts such a gene in each subset's own floor), so it has to be lossless.
                    exp[g] = (e, "both")
                else:
                    exp[g] = (e, sub)
        for g, (e, sub) in exp.items():
            fprows.append((r.unit_id, r.cmap_name, kind, sub, g, e,
                           "resolved" if e else "conflicted",
                           "landmark" if g in land else "best inferred"))
        for g in g_any - set(exp):      # reached only through an ambiguous TF: counted, never signed
            fprows.append((r.unit_id, r.cmap_name, kind, "ambiguous_tf", g, 0, "ambiguous_route",
                           "landmark" if g in land else "best inferred"))
        n_amb_tf = (tf_fp["n_ambiguous_tfs"] if tf_fp else 0) + (sig_fp["n_ambiguous_tfs"] if sig_fp else 0)
        n_anch = (tf_fp["n_anchored_tfs"] if tf_fp else 0) + (sig_fp["n_anchored_tfs"] if sig_fp else 0)
        hub = sig_fp if sig_fp is not None else tf_fp
        amb.append(dict(unit_id=r.unit_id, cmap_name=r.cmap_name, unit_kind=kind, genes=r.genes_str,
                        n_anchored_tfs=n_anch, n_ambiguous_tfs=n_amb_tf,
                        # D128-A2.1: the signalling subset's resolved-TF count, which the frozen P10 floor needs
                        # (>= FLOOR_TFS anchored TFs for that subset); A0's coverage.tsv carries the same number as
                        # `n_sig_subset_tfs`.
                        n_sig_subset_tfs=(sig_fp["n_resolved_tfs"] if sig_fp else 0),
                        n_tf_subset_resolved=(len(tf_fp["genes_resolved"]) if tf_fp else 0),
                        n_sig_subset_resolved=(len(sig_fp["genes_resolved"]) if sig_fp else 0),
                        tf_ambiguity_fraction=round(n_amb_tf / n_anch, 3) if n_anch else np.nan,
                        n_genes_any=len(g_any), n_resolved=len(g_res), n_conflicted=len(g_conf),
                        n_ambiguous_route=len(g_any - set(exp)),
                        gene_ambiguity_fraction=round((len(g_conf) + len(g_any - set(exp))) / len(g_any), 3) if g_any else np.nan,
                        hub_tf=hub["hub_tf"] if hub else "", hub_fraction=round(hub["hub_fraction"], 3)
                        if hub and hub["hub_fraction"] == hub["hub_fraction"] else np.nan))
    F = pd.DataFrame(fprows, columns=["unit_id", "cmap_name", "unit_kind", "subset", "gene_id",
                                      "expected_sign_pre_action", "status", "feature_space"])
    F.to_parquet(A1 / "target_gene_footprints.parquet", index=False)
    assert not F.duplicated(["unit_id", "gene_id"]).any(), "a gene appears twice in one unit's footprint"
    A = pd.DataFrame(amb)
    A.to_csv(A1 / "path_ambiguity.tsv", sep="\t", index=False)
    log(f"target_gene_footprints: {len(F)} (unit, gene) rows, {F.unit_id.nunique()} units; "
        f"{int((F.status == 'resolved').sum())} resolved, {int((F.status == 'conflicted').sum())} conflicted, "
        f"{int((F.status == 'ambiguous_route').sum())} via an ambiguous route")
    log(f"ambiguity: median TF-ambiguity fraction {A.tf_ambiguity_fraction.median():.3f}, "
        f"median gene-ambiguity fraction {A.gene_ambiguity_fraction.median():.3f}")

    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=C.ROOT, capture_output=True, text=True).stdout.strip()
    man = dict(stage="D128-A1", git_sha=sha, frozen=C.FROZEN, network=C.NETWORK, max_hops=C.MAX_HOPS,
               inputs={str(p.relative_to(C.ROOT)): IO.sha256(p) for p in (C.COLLECTRI, C.OMNIPATH, C.SIGNOR, C.UNITS)},
               counts=dict(regulon_edges=len(E), target_genes=len(targets), tf_paths=len(P), footprint_rows=len(F),
                           units=len(A), resolved=int((F.status == "resolved").sum()),
                           conflicted=int((F.status == "conflicted").sum()),
                           ambiguous_route=int((F.status == "ambiguous_route").sum())),
               log=LOG, seconds=round(time.time() - T0, 1))
    (A1 / "a1_run.json").write_text(json.dumps(man, indent=1, default=str))
    log("done")


if __name__ == "__main__":
    main()
