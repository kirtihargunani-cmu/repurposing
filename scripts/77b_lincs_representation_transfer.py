#!/usr/bin/env python3
"""Paired gene/activity/hierarchical LINCS cross-line transfer benchmark.

Consumes *the exact* pairs and specific profiles from 77_lincs_extrapolation_test.py.
No new held-out splits, perturbation filter, or dose/time aggregation are made here.

Gene profiles are the original specific Level-5 z-signatures. Signed, weighted
TF->target edges define a linear regulator *proxy*; this is NOT VIPER. An actual
precomputed VIPER matrix can be supplied, with a (drug,line) index and TF IDs.

Regulon-target metric: signed, confidence/overlap-weighted within-regulon Pearson r
between observed and transferred changes, aggregated across regulators. Regulon
membership is NOT dropped and shared genes are overlap-downweighted. This measures
transfer fidelity, NOT disease rescue or a model trained on the hierarchy.
"""

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_TAUS = (10, 25, 50)


def read_frame(path):
    path = Path(path)
    return pd.read_csv(path, sep="\t" if not str(path).endswith(".csv") else ",")


def read_ids(path, label):
    """An ordered one-column TSV/CSV/TXT, or numpy array; symbols as strings."""
    path = Path(path)
    if path.suffix == ".npy":
        a = np.load(path, allow_pickle=False)
        if a.ndim != 1:
            raise ValueError(f"{label}: expected 1D numpy array")
        return list(map(str, a))
    if path.suffix in (".txt", ".list"):
        return [s.strip() for s in path.read_text().splitlines() if s.strip()]
    f = read_frame(path)
    if f.shape[1] != 1:
        raise ValueError(f"{label}: expected one column (ordered to match .npy columns)")
    return f.iloc[:, 0].astype(str).tolist()


def read_gene_ids(path):
    """Read script 76's ordered gene metadata, preserving L1000 column order.

    Script 76 writes level5_landmark_genes.tsv in the same row order used to
    index the 978 columns. Prefer gene_symbol to match VIPER regulator targets.
    Also supports a single-column list via read_ids().
    """
    path = Path(path)
    if path.suffix in (".txt", ".list", ".npy"):
        return read_ids(path, "gene IDs")
    frame = read_frame(path)
    if "gene_symbol" in frame.columns:
        result = frame["gene_symbol"].astype(str).tolist()
    elif frame.shape[1] == 1:
        result = frame.iloc[:, 0].astype(str).tolist()
    else:
        raise ValueError("Gene metadata must include gene_symbol or have one column")
    if not result or any(not x or x.lower() == "nan" for x in result):
        raise ValueError("Missing gene symbols in ordered gene metadata")
    return result


def choose_col(cols, choices, required=True):
    lookup = {str(c).lower(): c for c in cols}
    for c in choices:
        if c.lower() in lookup:
            return lookup[c.lower()]
    if required:
        raise ValueError(f"Expected column among {choices}; found {list(cols)}")
    return None


def build_network(edges_file, genes, min_targets=5):
    """Generate signed linear TF matrix B and overlap-aware target-moment matrices.

    Required network TSV: regulator, target, and sign (or tfmode or mor).
    Optional: confidence/likelihood. If signed 'weight' is provided, abs(weight)
    acts as confidence. No sign is inferred from drug/disease expression.
    """
    f = read_frame(edges_file)
    tf_col = choose_col(f.columns, ["regulator", "source", "tf", "regulator_name"])
    tg_col = choose_col(f.columns, ["target", "gene", "target_gene", "target_name"])
    sg_col = choose_col(f.columns, ["sign", "mor", "tfmode", "mode", "weight"])
    conf_col = choose_col(f.columns, ["confidence", "likelihood", "edge_confidence", "support"], required=False)
    f = f[[tf_col, tg_col, sg_col] + ([conf_col] if conf_col else [])].copy()
    f.columns = ["regulator", "target", "sign"] + (["confidence"] if conf_col else [])
    f["regulator"] = f["regulator"].astype(str)
    f["target"] = f["target"].astype(str)
    f["sign"] = pd.to_numeric(f["sign"], errors="raise")
    if conf_col:
        f["confidence"] = pd.to_numeric(f["confidence"], errors="raise").abs()
    else:
        f["confidence"] = f["sign"].abs() if str(sg_col).lower() == "weight" else 1.0
    f["sign"] = np.sign(f["sign"])
    if len(set(genes)) != len(genes):
        raise ValueError("Duplicate gene labels; matrices require unique ordered gene IDs")
    known = set(genes)
    f = f[f.target.isin(known)].copy()
    f = f[(f.sign != 0) & (f.confidence > 0) & np.isfinite(f.confidence)].copy()
    # Contradictory duplicate edges cannot be safely collapsed; reject them.
    if f.duplicated(["regulator", "target"]).any():
        conflicting = f.groupby(["regulator", "target"])["sign"].nunique()
        if (conflicting > 1).any():
            raise ValueError("Network contains duplicate TF-target edges with conflicting signs")
        f = f.sort_values("confidence", ascending=False).drop_duplicates(["regulator", "target"])
    counts = f.groupby("regulator").size()
    keep = counts[counts >= min_targets].index
    f = f[f.regulator.isin(keep)].copy()
    if f.empty:
        raise ValueError("No edges remain after gene matching and min-target threshold. Check gene IDs.")
    regulators = sorted(f.regulator.unique())
    reg_idx = {r: i for i, r in enumerate(regulators)}
    gene_idx = {g: i for i, g in enumerate(genes)}
    rr = f.regulator.map(reg_idx).to_numpy()
    cc = f.target.map(gene_idx).to_numpy()
    signs = f.sign.to_numpy(dtype=float)
    conf = f.confidence.to_numpy(dtype=float)
    shape = (len(regulators), len(genes))
    # Smooth TF activity proxy; L2-normalize per regulon (not a VIPER NES).
    B = sparse.csr_matrix((signs * conf, (rr, cc)), shape=shape, dtype=float)
    norm = np.sqrt(np.asarray(B.multiply(B).sum(axis=1)).ravel())
    B = sparse.diags(1 / norm) @ B
    # Shared genes should not silently count multiple times across regulators.
    degree = np.bincount(cc, minlength=len(genes))
    q = conf / degree[cc]
    W = sparse.csr_matrix((q, (rr, cc)), shape=shape, dtype=float)
    Ws = sparse.csr_matrix((q * signs, (rr, cc)), shape=shape, dtype=float)
    total = np.asarray(W.sum(axis=1)).ravel()
    # Per-regulon moments are mean weighted by confidence, corrected for overlap.
    coverage = pd.DataFrame({
        "regulator": regulators,
        "n_landmark_targets": np.bincount(rr, minlength=len(regulators)),
        "weight_total_overlap_corrected": total,
    })
    return regulators, B.tocsr(), W.tocsr(), Ws.tocsr(), total, coverage


def pearson_rows(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim == 1:
        x, y = x[None, :], y[None, :]
    x = x - x.mean(axis=1, keepdims=True)
    y = y - y.mean(axis=1, keepdims=True)
    den = np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1)
    return np.divide((x * y).sum(axis=1), den, out=np.full(len(x), np.nan), where=den > 1e-12)


def intra_regulon_correlations(truth, pred, W, Ws, total):
    """Signed confidence-weighted Pearson across targets *within each TF*.

    Vectorized sparse moment calculations avoid a slow Python TF-by-TF loop.
    Signed target contributions (sign*x) are used to retain the mode of regulation.
    """
    truth, pred = np.asarray(truth, float), np.asarray(pred, float)
    u = np.asarray(Ws @ truth).ravel() / total
    v = np.asarray(Ws @ pred).ravel() / total
    uu = np.asarray(W @ (truth * truth)).ravel() / total - u * u
    vv = np.asarray(W @ (pred * pred)).ravel() / total - v * v
    uv = np.asarray(W @ (truth * pred)).ravel() / total - u * v
    den = np.sqrt(np.maximum(uu, 0) * np.maximum(vv, 0))
    return np.divide(uv, den, out=np.full_like(den, np.nan), where=den > 1e-12)


def stable_null_seed(pair_key, seed):
    h = hashlib.sha256(f"{seed}/{pair_key[0]}/{pair_key[1]}".encode()).digest()
    return int.from_bytes(h[:8], "little")


def make_predictions(i, drug, line, idx, lookup, drug_rows, line_rows, S, dist, taus, seed):
    others = [j for j in drug_rows[drug] if j != i]
    if not others:
        raise ValueError(f"No source cell line for pair {drug} x {line}")
    source_lines = idx.ModelID.to_numpy()[others]
    dd = dist.loc[line, source_lines].to_numpy(dtype=float)
    if not np.all(np.isfinite(dd)):
        raise ValueError(f"Nonfinite Celligner distance for {drug} x {line}")
    order = np.argsort(dd)
    nearest_idx = others[order[0]]
    p = {"nearest": (np.array([nearest_idx], dtype=int), np.array([1.0])),
         "mean": (np.array(others, dtype=int), np.ones(len(others)) / len(others))}
    for tau in taus:
        weights = np.exp(-(dd - dd.min()) / tau)
        weights /= weights.sum()
        p[f"weighted_tau{tau:g}"] = (np.array(others, dtype=int), weights)
    pool = [j for j in line_rows[source_lines[order[0]]] if idx.pert_id.iat[j] != drug]
    if pool:
        rng = np.random.default_rng(stable_null_seed((drug, line), seed))
        p["null"] = (np.array([rng.choice(pool)], dtype=int), np.array([1.0]))
    return p


def load_external_activities(matrix_path, index_path, regulator_ids_path, idx, regulators, layout="C"):
    """Align external (preferably VIPER-on-specific) activities by exact pair keys."""
    ids = read_ids(regulator_ids_path, "regulator IDs")
    ext_idx = read_frame(index_path)
    if str(matrix_path).endswith(".f32"):
        M = np.memmap(matrix_path, dtype="float32", mode="r",
                      shape=(len(ext_idx), len(ids)), order=layout)
    else:
        M = np.load(matrix_path, mmap_mode="r")
    if not {"pert_id", "ModelID"}.issubset(ext_idx.columns):
        raise ValueError("External activity index must contain pert_id and ModelID")
    if M.shape != (len(ext_idx), len(ids)):
        raise ValueError("External activity matrix shape disagrees with activity index / regulator IDs")
    if ext_idx.duplicated(["pert_id", "ModelID"]).any() or len(set(ids)) != len(ids):
        raise ValueError("Duplicate external activity row keys or regulator IDs")
    rows = pd.MultiIndex.from_frame(ext_idx[["pert_id", "ModelID"]].astype(str))
    wanted = pd.MultiIndex.from_frame(idx[["pert_id", "ModelID"]].astype(str))
    row_idx = rows.get_indexer(wanted)
    if (row_idx < 0).any():
        raise ValueError(f"External activities missing {(row_idx < 0).sum()} original drug-line rows")
    col_idx = pd.Index(ids).get_indexer(regulators)
    if (col_idx < 0).any():
        raise ValueError(f"External activities missing {(col_idx < 0).sum()} network regulators")
    return np.asarray(M[np.ix_(row_idx, col_idx)], dtype=float)


def measure(truth, pred, B, W, Ws, total, observed_activity, predicted_activity):
    gene_r = pearson_rows(truth, pred)[0]
    reg_r = pearson_rows(observed_activity, predicted_activity)[0]
    per_tf = intra_regulon_correlations(truth, pred, W, Ws, total)
    finite = np.isfinite(per_tf)
    direction = np.sign(observed_activity) == np.sign(predicted_activity)
    confident_direction = (np.abs(observed_activity) > 1e-10) & (np.abs(predicted_activity) > 1e-10)
    valid = finite & confident_direction
    return {
        "gene_r": float(gene_r),
        "activity_r": float(reg_r),
        "regulon_target_median_r": float(np.median(per_tf[finite])) if finite.any() else np.nan,
        "regulon_target_mean_r": float(np.mean(per_tf[finite])) if finite.any() else np.nan,
        "regulon_target_frac_positive": float(np.mean(per_tf[finite] > 0)) if finite.any() else np.nan,
        "activity_direction_agreement": float(np.mean(direction[confident_direction])) if confident_direction.any() else np.nan,
        "joint_regulon_fraction": float(np.mean((per_tf[valid] > 0.3) & direction[valid])) if valid.any() else np.nan,
        "n_regulons_valid": int(finite.sum()),
    }


def run(args):
    root = Path(args.root).resolve()
    proc = root / "data/processed/lincs"
    res = root / "results/lincs"
    pairs_file = Path(args.pairs) if args.pairs else res / "extrapolation_lolo_pairs.tsv"
    prof_file = Path(args.profiles) if args.profiles else proc / "specific_profiles.npy"
    idx_file = Path(args.index) if args.index else proc / "active_profiles_index.tsv.gz"
    dist_file = Path(args.dist) if args.dist else root / "data/processed/celligner/dist_line_line.tsv.gz"
    outdir = Path(args.output_dir) if args.output_dir else res
    outdir.mkdir(parents=True, exist_ok=True)

    S = np.load(prof_file, mmap_mode="r")
    idx = read_frame(idx_file)
    pairs = read_frame(pairs_file)
    dist = pd.read_csv(dist_file, sep="\t", index_col=0)
    gene_file = Path(args.genes) if args.genes else proc / "level5_landmark_genes.tsv"
    genes = read_gene_ids(gene_file)
    if S.shape != (len(idx), len(genes)):
        raise ValueError(f"Mismatch: S.shape={S.shape}, idx={len(idx)}, landmark genes={len(genes)}")
    for name in ("pert_id", "ModelID"):
        if name not in idx or name not in pairs:
            raise ValueError(f"Missing {name} from original index/pair table")
        idx[name] = idx[name].astype(str)
        pairs[name] = pairs[name].astype(str)
    dist.index = dist.index.astype(str)
    dist.columns = dist.columns.astype(str)
    if idx.duplicated(["pert_id", "ModelID"]).any() or pairs.duplicated(["pert_id", "ModelID"]).any():
        raise ValueError("Input contains duplicate drug-line keys")
    lookup = {(r.pert_id, r.ModelID): i for i, r in idx.iterrows()}
    missing = [(r.pert_id, r.ModelID) for r in pairs.itertuples(index=False)
               if (r.pert_id, r.ModelID) not in lookup]
    if missing:
        raise ValueError(f"Original LOLO pairs not found in source index: {missing[:5]}")
    regulators, B, W, Ws, total, coverage = build_network(args.edges, genes, args.min_targets)
    if args.activity_matrix:
        if not (args.activity_index and args.activity_ids):
            raise ValueError("External VIPER requires --activity-index and --activity-ids")
        A = load_external_activities(args.activity_matrix, args.activity_index, args.activity_ids, idx, regulators,
                                     args.activity_layout)
        activity_label = args.activity_label or "external_activity_matrix"
        warnings.warn("Verify external activity matrix was derived from the SAME specific-profile basis. "
                      "VIPER on raw L1000 profiles is not directly comparable to this benchmark.")
    else:
        # Linear target-weighted signed projection, NOT VIPER.
        A = np.asarray(S @ B.T, dtype=float)
        activity_label = "linear_signed_projection_not_VIPER"
    if not np.all(np.isfinite(A)):
        raise ValueError("Activity matrix has missing/nonfinite values")
    drug_rows = {k: list(v) for k, v in idx.groupby("pert_id").indices.items()}
    line_rows = {k: list(v) for k, v in idx.groupby("ModelID").indices.items()}
    taus = DEFAULT_TAUS if args.taus is None else tuple(args.taus)

    out = []
    for number, pair in enumerate(pairs.itertuples(index=False), 1):
        i = lookup[(pair.pert_id, pair.ModelID)]
        pred_defs = make_predictions(i, pair.pert_id, pair.ModelID, idx, lookup, drug_rows,
                                     line_rows, S, dist, taus, args.seed)
        target = np.asarray(S[i], dtype=float)
        target_a = A[i]
        for method, (source, weights) in pred_defs.items():
            pred = np.asarray(weights @ np.asarray(S[source], dtype=float), dtype=float)
            if args.activity_matrix:
                pred_a = weights @ A[source]
            else:
                pred_a = np.asarray(B @ pred).ravel()
            metrics = measure(target, pred, B, W, Ws, total, target_a, pred_a)
            out.append({"ModelID": pair.ModelID, "pert_id": pair.pert_id,
                        "method": method, "activity_representation": activity_label,
                        "n_source_lines": len(drug_rows[pair.pert_id]) - 1,
                        **metrics})
        if args.max_pairs and number >= args.max_pairs:
            break
        if args.progress_every and number % args.progress_every == 0:
            print(f"evaluated {number:,}/{len(pairs):,} held-out pairs", flush=True)
    detail = pd.DataFrame(out)
    if detail.empty:
        raise ValueError("No evaluation records")
    # Invariance check: original nearest and mean gene scores must match.
    check = detail.pivot(index=["ModelID", "pert_id"], columns="method", values="gene_r")
    compare = pairs.set_index(["ModelID", "pert_id"])
    for m, old_col in (("nearest", "r_nearest"), ("mean", "r_mean"),
                       *((f"weighted_tau{tau:g}", f"r_weighted_tau{tau:g}") for tau in taus)):
        if m in check and old_col in compare:
            common = check[m].dropna().index.intersection(compare.index)
            delta = np.abs(check.loc[common, m].to_numpy() - compare.loc[common, old_col].to_numpy())
            delta = delta[np.isfinite(delta)]
            if len(delta) and np.max(delta) > args.gene_tolerance:
                raise AssertionError(
                    f"Existing gene benchmark is NOT reproduced ({old_col}, max |delta|={np.max(delta):.5g}). "
                    "Check profile/index pairing before interpreting the new metrics.")
            print(f"Original {old_col} max abs difference: {np.max(delta) if len(delta) else np.nan:.6f}")
    detail.to_csv(outdir / "extrapolation_representation_pairs.tsv.gz", sep="\t", index=False,
                  compression="gzip", float_format="%.6g")
    aggs = {k: (k, "median") for k in ("gene_r", "activity_r", "regulon_target_median_r",
                                             "regulon_target_mean_r", "regulon_target_frac_positive",
                                             "activity_direction_agreement", "joint_regulon_fraction", "n_regulons_valid")}
    summary = detail.groupby("method", sort=False).agg(n_pairs=("pert_id", "size"), **aggs).reset_index()
    summary.to_csv(outdir / "extrapolation_representation_summary.tsv", sep="\t", index=False, float_format="%.6g")
    by_line = detail.groupby(["ModelID", "method"], sort=False).agg(n_pairs=("pert_id", "size"), **aggs).reset_index()
    by_line.to_csv(outdir / "extrapolation_representation_by_line.tsv", sep="\t", index=False, float_format="%.6g")
    coverage.to_csv(outdir / "extrapolation_regulon_coverage.tsv", sep="\t", index=False)
    manifest = dict(
        original_pair_count=len(pairs), evaluated_pair_count=detail[["pert_id", "ModelID"]].drop_duplicates().shape[0],
        methods=list(detail.method.unique()), activity_representation=activity_label,
        n_landmark_genes=len(genes), n_regulons=len(regulators), n_edges=int(W.nnz),
        min_targets=args.min_targets, seed=args.seed, source_files={
            "pairs": str(pairs_file), "profiles": str(prof_file), "index": str(idx_file),
            "dist": str(dist_file), "edges": str(args.edges), "genes": str(gene_file)},
        caveats=["Does not introduce new held-out splits", "Per-regulon target r is signed and overlap-adjusted",
                 "Linear projection is not VIPER", "Hierarchy metrics evaluate, not train, a transfer predictor",
                 "Random-null drug is deterministically reselected; not the original 77 draw",
                 "L1000 level-5 signatures are z-scores, NOT additive post-treatment expression"])
    (outdir / "extrapolation_representation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))
    return detail, summary


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=str(Path(__file__).resolve().parent.parent), help="Repository root")
    p.add_argument("--edges", required=True, help="Signed TF-target TSV (regulator,target,sign,confidence)")
    p.add_argument("--genes", help="Ordered landmark-gene metadata (default: script 76 level5_landmark_genes.tsv)")
    p.add_argument("--pairs", help="Default results/lincs/extrapolation_lolo_pairs.tsv")
    p.add_argument("--profiles", help="Default data/processed/lincs/specific_profiles.npy")
    p.add_argument("--index", help="Default data/processed/lincs/active_profiles_index.tsv.gz")
    p.add_argument("--dist", help="Default data/processed/celligner/dist_line_line.tsv.gz")
    p.add_argument("--output-dir", help="Default results/lincs")
    p.add_argument("--min-targets", type=int, default=5)
    p.add_argument("--taus", nargs="*", type=float, help="Softmax temperatures; default 10 25 50")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--progress-every", type=int, default=1000)
    p.add_argument("--gene-tolerance", type=float, default=0.002,
                   help="Max allowed discrepancy from original rounded r_nearest/r_mean")
    p.add_argument("--max-pairs", type=int, help="Debug only; preserve all pairs for results")
    p.add_argument("--activity-matrix", help="Optional precomputed drug-line VIPER/activity .npy (same baseline)")
    p.add_argument("--activity-index", help="TSV row keys for activity matrix (pert_id,ModelID)")
    p.add_argument("--activity-ids", help="Ordered regulator column IDs for activity matrix")
    p.add_argument("--activity-label", help="Report label for external activity matrix")
    p.add_argument("--activity-layout", choices=["C", "F"], default="C",
                   help="Memory layout of raw .f32 external activity matrix (R usually writes F)")
    return p.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
