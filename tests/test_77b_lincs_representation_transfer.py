"""Synthetic unit tests. No private/large LINCS data are needed."""
import importlib.util
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


MODULE = Path(__file__).resolve().parents[1] / "scripts/77b_lincs_representation_transfer.py"
spec = importlib.util.spec_from_file_location("lincs77b", MODULE)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def synthetic_repo(tmp_path):
    root = tmp_path / "repo"
    proc = root / "data/processed/lincs"
    res = root / "results/lincs"
    cell = root / "data/processed/celligner"
    for d in (proc, res, cell):
        d.mkdir(parents=True)
    idx = pd.DataFrame({"pert_id": ["D1", "D1", "D1", "D2", "D2"],
                        "ModelID": ["A", "B", "C", "A", "B"]})
    idx.to_csv(proc / "active_profiles_index.tsv.gz", sep="\t", index=False)
    S = np.array([[1.0, -2, 3.0, 0.5, -0.6, 0.8],
                  [1.1, -1.9, 3.2, 0.4, -0.8, 0.6],
                  [-1.0, 2.0, -3.0, -0.5, 0.6, -0.8],
                  [-1.2, 1.7, 0.4, -1.0, 2.0, -2.0],
                  [-1.1, 1.8, 0.3, -0.8, 1.8, -2.1]], np.float32)
    np.save(proc / "specific_profiles.npy", S)
    dist = pd.DataFrame([[0, 10, 50], [10, 0, 45], [50, 45, 0]],
                        index=["A", "B", "C"], columns=["A", "B", "C"])
    dist.to_csv(cell / "dist_line_line.tsv.gz", sep="\t")
    pairs = pd.DataFrame([
        {"pert_id": "D1", "ModelID": "A", "r_nearest": mod.pearson_rows(S[0], S[1])[0],
         "r_mean": mod.pearson_rows(S[0], S[[1, 2]].mean(0))[0]},
        {"pert_id": "D1", "ModelID": "B", "r_nearest": mod.pearson_rows(S[1], S[0])[0],
         "r_mean": mod.pearson_rows(S[1], S[[0, 2]].mean(0))[0]},
        {"pert_id": "D1", "ModelID": "C", "r_nearest": mod.pearson_rows(S[2], S[1])[0],
         "r_mean": mod.pearson_rows(S[2], S[[0, 1]].mean(0))[0]},
        {"pert_id": "D2", "ModelID": "A", "r_nearest": mod.pearson_rows(S[3], S[4])[0],
         "r_mean": mod.pearson_rows(S[3], S[4])[0]},
        {"pert_id": "D2", "ModelID": "B", "r_nearest": mod.pearson_rows(S[4], S[3])[0],
         "r_mean": mod.pearson_rows(S[4], S[3])[0]},
    ])
    pairs.to_csv(res / "extrapolation_lolo_pairs.tsv", sep="\t", index=False, float_format="%.4f")
    edge = pd.DataFrame({
        "regulator": ["TF1"] * 3 + ["TF2"] * 3 + ["TF3"] * 3,
        "target": ["G1", "G2", "G3", "G2", "G4", "G5", "G3", "G5", "G6"],
        "sign": [1, 1, -1, -1, 1, 1, 1, -1, 1],
        "confidence": [0.9, 0.8, 0.7, 0.3, 0.9, 0.5, 0.8, 0.9, 0.4],
    })
    edge.to_csv(root / "edges.tsv", sep="\t", index=False)
    (root / "genes.txt").write_text("\n".join([f"G{i}" for i in range(1, 7)]) + "\n")
    return root


def args(root, **kw):
    o = dict(root=str(root), edges=str(root / "edges.tsv"), genes=str(root / "genes.txt"),
             pairs=None, profiles=None, index=None, dist=None, output_dir=None, min_targets=3,
             taus=None, seed=0, progress_every=0, gene_tolerance=0.002, max_pairs=None,
             activity_matrix=None, activity_index=None, activity_ids=None, activity_label=None,
             activity_layout="C")
    o.update(kw)
    return Namespace(**o)


def test_reuses_exact_pairs_and_gene_baseline(tmp_path):
    root = synthetic_repo(tmp_path)
    df, summ = mod.run(args(root))
    assert df[["ModelID", "pert_id"]].drop_duplicates().shape[0] == 5
    assert set(df.method) == {"nearest", "mean", "weighted_tau10", "weighted_tau25", "weighted_tau50", "null"}
    original = pd.read_csv(root / "results/lincs/extrapolation_lolo_pairs.tsv", sep="\t")
    got = df[df.method == "nearest"].set_index(["pert_id", "ModelID"]).gene_r
    for p in original.itertuples(index=False):
        assert np.isclose(got.loc[p.pert_id, p.ModelID], p.r_nearest, atol=0.002)
    assert summ.n_pairs.min() == 5
    assert (root / "results/lincs/extrapolation_representation_manifest.json").exists()
    assert df.activity_representation.eq("linear_signed_projection_not_VIPER").all()
    assert df.regulon_target_median_r.notna().any()


def test_signs_count_and_overlaps_are_correct(tmp_path):
    root = synthetic_repo(tmp_path)
    regs, B, W, Ws, total, coverage = mod.build_network(root / "edges.tsv", [f"G{i}" for i in range(1, 7)], 3)
    assert len(regs) == 3
    assert coverage.n_landmark_targets.tolist() == [3, 3, 3]
    # G2 belongs to two TFs; its fractional weight is confidence/2 for TF1.
    assert np.isclose(W[0, 1], 0.8 / 2)
    assert np.isclose(Ws[1, 1], -0.3 / 2)
    truth = np.array([1, -2, 3, 0.5, -0.6, 0.8], float)
    per_tf = mod.intra_regulon_correlations(truth, truth, W, Ws, total)
    np.testing.assert_allclose(per_tf, 1.0)
    pred = -truth
    per_tf_opposite = mod.intra_regulon_correlations(truth, pred, W, Ws, total)
    np.testing.assert_allclose(per_tf_opposite, -1.0)


def test_detects_incorrect_old_benchmark_pairing(tmp_path):
    root = synthetic_repo(tmp_path)
    p = root / "results/lincs/extrapolation_lolo_pairs.tsv"
    df = pd.read_csv(p, sep="\t")
    df.loc[0, "r_nearest"] = -0.8
    df.to_csv(p, sep="\t", index=False)
    with pytest.raises(AssertionError, match="NOT reproduced"):
        mod.run(args(root))


def test_external_viper_alignment_is_explicit(tmp_path):
    root = synthetic_repo(tmp_path)
    idx = pd.read_csv(root / "data/processed/lincs/active_profiles_index.tsv.gz", sep="\t")
    idx = idx.iloc[::-1].reset_index(drop=True)
    idx.to_csv(root / "activity_idx.tsv", sep="\t", index=False)
    (root / "reg_ids.txt").write_text("TF3\nTF2\nTF1\n")
    v = np.random.default_rng(2).normal(size=(len(idx), 3))
    np.save(root / "activity.npy", v)
    df, _ = mod.run(args(root, activity_matrix=str(root / "activity.npy"),
                         activity_index=str(root / "activity_idx.tsv"),
                         activity_ids=str(root / "reg_ids.txt"), activity_label="VIPER_specific"))
    assert df.activity_representation.eq("VIPER_specific").all()
    assert df.activity_r.notna().any()


def test_rejects_conflicting_edge_signs(tmp_path):
    root = synthetic_repo(tmp_path)
    p = root / "edges.tsv"
    x = pd.read_csv(p, sep="\t")
    x = pd.concat([x, x.iloc[[0]].assign(sign=-1)], ignore_index=True)
    x.to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError, match="conflicting"):
        mod.build_network(p, [f"G{i}" for i in range(1, 7)], min_targets=3)


def test_explicit_temperatures_preserve_old_weighted_column_names(tmp_path):
    root = synthetic_repo(tmp_path)
    df, _ = mod.run(args(root, taus=[10.0, 25.0]))
    assert set(df.method) == {"nearest", "mean", "weighted_tau10", "weighted_tau25", "null"}


def test_reads_script76_ordered_gene_metadata_by_default(tmp_path):
    root = synthetic_repo(tmp_path)
    genes = [f"G{i}" for i in range(1, 7)]
    pd.DataFrame({"gene_id": [str(i) for i in range(1, 7)],
                  "gene_symbol": genes, "ensembl_id": [f"ENSG{i}" for i in range(1, 7)]}).to_csv(
        root / "data/processed/lincs/level5_landmark_genes.tsv", sep="\t", index=False)
    df, summary = mod.run(args(root, genes=None))
    assert len(df) > 0
    assert (summary.n_pairs == 5).all()
