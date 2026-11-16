"""D128 io checks (plan Section 6.1 test_io.py): byte size / orientation, gene order, duplicates, signs.

Synthetic fixtures only; no LINCS, R, credentials or GPU needed.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/d128"))
import d128io as IO  # noqa: E402


def _write_matrix(tmp_path, rows, cols):
    m = np.arange(rows * cols, dtype=np.float32).reshape(rows, cols)
    p = tmp_path / "m.f32"
    m.tofile(p)
    return p, m


def test_profiles_profile_major(tmp_path):
    p, m = _write_matrix(tmp_path, 3, 5)
    x = IO.profiles(p, 3, 5)
    assert x.shape == (3, 5)
    np.testing.assert_array_equal(x[1], m[1])          # row = profile, column = gene


def test_profiles_wrong_shape_fails(tmp_path):
    p, _ = _write_matrix(tmp_path, 3, 5)
    with pytest.raises(SystemExit, match="bytes"):
        IO.profiles(p, 4, 5)


def test_index_duplicates_fail(tmp_path):
    p = tmp_path / "idx.tsv"
    pd.DataFrame(dict(pert_id=["a", "a"], ModelID=["L1", "L1"], n_sig=[1, 2])).to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit, match="duplicated"):
        IO.profile_index(p)


def test_genes_order_kept_and_space_checked(tmp_path):
    p = tmp_path / "g.tsv"
    g = pd.DataFrame(dict(gene_id=["ENSG3", "ENSG1", "ENSG2"], symbol=["C", "A", "B"],
                          feature_space=["landmark", "best inferred", "landmark"]))
    g.to_csv(p, sep="\t", index=False)
    assert IO.profile_genes(p).gene_id.tolist() == ["ENSG3", "ENSG1", "ENSG2"]   # never sorted
    g.loc[0, "feature_space"] = "imputed"
    g.to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit, match="feature_space"):
        IO.profile_genes(p)


def test_collectri_sign_and_confidence(tmp_path):
    p = tmp_path / "c.tsv"
    pd.DataFrame(dict(regulator=["ESR1"], target=["ENSG1"], tfmode=[0], likelihood=[1.0])).to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit, match="tfmode"):
        IO.collectri(p)
    pd.DataFrame(dict(regulator=["ESR1"], target=["ENSG1"], tfmode=[1], likelihood=[0.0])).to_csv(p, sep="\t", index=False)
    with pytest.raises(SystemExit, match="likelihood"):
        IO.collectri(p)
