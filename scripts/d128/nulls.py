"""D128 A2 nulls (plan 3.5): seeded controls that preserve the structures which would otherwise explain a score.

Every null is scored through the SAME functions as the real target, so a difference cannot come from a different
code path. No therapeutic label enters any of these operations.

  degree_matched_targets   another target unit of similar footprint size and TF count: tests whether the score only
                           reflects how much of the transcriptome a target can reach.
  shuffled_edge_signs      the footprint's expected signs permuted, keeping which genes are in it: tests whether the
                           signed structure matters or only the gene set.
  gene_label_permutation   responses permuted across genes, matched by feature space: tests whether the specific
                           gene-to-response pairing matters.
"""
import numpy as np

import config as C


def rng(seed_offset=0):
    return np.random.default_rng(C.SEED + seed_offset)


def shuffled_edge_signs(expected, seed_offset=0):
    """Permute the expected signs over the same genes (preserves the sign balance and the gene set)."""
    e = np.asarray(expected, dtype=float).copy()
    rng(seed_offset).shuffle(e)
    return e


def gene_label_permutation(z, landmark_mask=None, seed_offset=0):
    """Permute responses across genes, within feature space when given (landmark and inferred differ in spread)."""
    z = np.asarray(z, dtype=float).copy()
    r = rng(seed_offset)
    if landmark_mask is None:
        r.shuffle(z)
        return z
    lm = np.asarray(landmark_mask, bool)
    for mask in (lm, ~lm):
        if mask.sum() > 1:
            idx = np.flatnonzero(mask)
            z[idx] = z[r.permutation(idx)]
    return z


def degree_matched_units(unit_id, candidates, n_genes, n_tfs, k=20, tol=0.25, seed_offset=0, exclude_genes=None):
    """Pick up to k other units whose footprint size and TF count are within `tol` of this one.

    `candidates` is a DataFrame with unit_id, n_sign_resolved, n_resolved_tfs and genes. Matching on both size and TF
    count avoids a null that is simply smaller than the real footprint. `exclude_genes` drops any candidate sharing a
    target gene with the real unit: without it the pool fills with other drugs against the SAME target (every ESR1
    drug has an identical footprint), whose scores are identical by construction, and the null becomes the observation.
    """
    c = candidates[candidates.unit_id != unit_id]
    if exclude_genes is not None and "genes" in c.columns:
        ex = set(exclude_genes)
        c = c[~c.genes.fillna("").apply(lambda g: bool(ex & set(str(g).split("|"))))]
    lo, hi = n_genes * (1 - tol), n_genes * (1 + tol)
    tlo, thi = max(0, n_tfs * (1 - tol)), n_tfs * (1 + tol) + 1
    m = c[(c.n_sign_resolved >= lo) & (c.n_sign_resolved <= hi)
          & (c.n_resolved_tfs >= tlo) & (c.n_resolved_tfs <= thi)]
    if len(m) > k:
        m = m.iloc[rng(seed_offset).choice(len(m), k, replace=False)]
    return m


def summarise(values, observed):
    """Null median, spread and the observed value's empirical percentile among them."""
    v = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
    if not len(v):
        return dict(n_null=0, null_median=np.nan, null_sd=np.nan, percentile=np.nan, z=np.nan)
    sd = float(v.std(ddof=1)) if len(v) > 1 else np.nan
    return dict(n_null=len(v), null_median=float(np.median(v)), null_sd=sd,
                percentile=float((v < observed).mean() * 100) if np.isfinite(observed) else np.nan,
                z=float((observed - v.mean()) / sd) if np.isfinite(observed) and sd and np.isfinite(sd) else np.nan)


# --- fast paths (D128-A3) -------------------------------------------------------------------------------------------
# A3 scores the whole comparison population, so a null that calls score_mechanism once per draw is not affordable
# (16 M calls). These helpers compute the SAME number by the SAME permutations and are only used where the shortcut is
# exact: every expected sign non-zero and every response finite, so `usable` is the whole gene set and the weighted
# mean of score.py reduces to a plain mean. tanh is element-wise, so permuting f(z) equals f(permuted z).
# tests/d128/test_scores.py asserts element-wise equality with score_mechanism. Otherwise the caller must use the
# per-draw path.

def exact_shortcut_ok(expected, z):
    e = np.asarray(expected, dtype=float)
    return bool(len(e)) and bool(np.all(e != 0)) and bool(np.all(np.isfinite(np.asarray(z, dtype=float))))


def m_null_sign_shuffle(expected, fz, n):
    e = np.asarray(expected, dtype=float)
    fz = np.asarray(fz, dtype=float)
    return np.array([float((shuffled_edge_signs(e, i) * fz).mean()) for i in range(n)])


def m_null_gene_label(expected, fz, landmark_mask, n):
    e = np.asarray(expected, dtype=float)
    fz = np.asarray(fz, dtype=float)
    return np.array([float((e * gene_label_permutation(fz, landmark_mask, i)).mean()) for i in range(n)])
