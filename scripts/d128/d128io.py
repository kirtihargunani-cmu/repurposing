"""D128: read-only loaders with strict shape / orientation checks (plan Section 2.1-2.2, E14).

Nothing here writes. Every loader fails loudly (SystemExit with the dataset named) instead of guessing.
"""
import gzip
import hashlib
import os
from pathlib import Path

import numpy as np
import pandas as pd

from config import PROFILE_GENES, PROFILE_INDEX


class InputError(SystemExit):
    pass


def sha256(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def profile_index(path=PROFILE_INDEX):
    idx = pd.read_csv(path, sep="\t")
    need = {"pert_id", "ModelID", "n_sig"}
    if not need <= set(idx.columns):
        raise InputError(f"{path}: missing columns {need - set(idx.columns)}")
    dup = idx.duplicated(["pert_id", "ModelID"])
    if dup.any():
        raise InputError(f"{path}: {int(dup.sum())} duplicated (pert_id, ModelID) rows; aggregated profiles must be unique")
    return idx.reset_index(drop=True)


def profile_genes(path=PROFILE_GENES):
    g = pd.read_csv(path, sep="\t")
    need = {"gene_id", "symbol", "feature_space"}
    if not need <= set(g.columns):
        raise InputError(f"{path}: missing columns {need - set(g.columns)}")
    if g.gene_id.duplicated().any():
        raise InputError(f"{path}: duplicated gene_id {g.gene_id[g.gene_id.duplicated()].tolist()[:5]}")
    bad = set(g.feature_space) - {"landmark", "best inferred"}
    if bad:
        raise InputError(f"{path}: unknown feature_space values {bad}")
    return g.reset_index(drop=True)   # order is the matrix column order: never sort


def profiles(path, n_rows, n_cols):
    """Profile-major float32 matrix (rows = profile index, cols = gene list), memory-mapped read-only."""
    path = Path(path)
    size = os.path.getsize(path)
    if size != n_rows * n_cols * 4:
        raise InputError(f"{path}: {size} bytes != {n_rows} rows x {n_cols} cols x 4 bytes (wrong layout or index)")
    return np.memmap(path, dtype=np.float32, mode="r", shape=(n_rows, n_cols))


def collectri(path):
    r = pd.read_csv(path, sep="\t")
    need = {"regulator", "target", "tfmode", "likelihood"}
    if not need <= set(r.columns):
        raise InputError(f"{path}: missing columns {need - set(r.columns)}")
    if not r.tfmode.isin([-1, 1]).all():
        raise InputError(f"{path}: tfmode outside +-1 ({sorted(r.tfmode.unique())[:5]})")
    if not (np.isfinite(r.likelihood) & (r.likelihood > 0)).all():
        raise InputError(f"{path}: non-finite or non-positive likelihood")
    return r


def signaling(path):
    e = pd.read_csv(path, sep="\t", usecols=["source_genesymbol", "target_genesymbol", "sign", "sign_disputed"])
    if not e["sign"].isin([-1, 1]).all():
        raise InputError(f"{path}: sign outside +-1")
    return e.rename(columns={"source_genesymbol": "src", "target_genesymbol": "dst"})


def depmap_matrix(path, genes):
    """DepMap gene matrix (rows ModelID, columns 'SYMBOL (entrez)'): the requested symbols only."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split(",")
    # DepMap leaves the ModelID column unnamed in some files (pandas calls it "Unnamed: 0"), and a usecols *callable*
    # with index_col=0 indexes the filtered frame, which silently takes a data column as the index. Select by
    # position instead: column 0 is always ModelID.
    keep = {i: c.split(" (")[0] for i, c in enumerate(header) if i and c.split(" (")[0] in genes}
    missing = set(genes) - set(keep.values())
    m = pd.read_csv(path, index_col=0, usecols=[0, *keep], low_memory=False)
    m.columns = [keep[i] for i in sorted(keep)]
    m.index.name = "ModelID"
    m.attrs["missing_genes"] = sorted(missing)
    return m


def de_table(proc_dir, unit, contrast, cols=("gene_id", "symbol", "logFC", "z", "adj.P.Val")):
    p = Path(proc_dir) / f"{unit}__{contrast}.tsv.gz"
    if not p.exists():
        return None
    d = pd.read_csv(p, sep="\t", usecols=list(cols))
    if d.gene_id.duplicated().any():
        raise InputError(f"{p}: duplicated gene_id")
    return d
