"""D128 A2: the M and R scores as pure functions (plan 3.2, 3.3, 3.4; Section 0 frozen values in config.py).

Both take typed arrays and return a result object carrying the score, its coverage and every per-gene contribution,
so a score can always be explained gene by gene. Neither reads a file, a label or a global.

Conventions, all from the frozen Section 0:
  f(z)  = tanh(z / F_DIVISOR)                       bounded drug response, P9
  M     = sum_g w_g e_g f(z_g) / sum_g w_g          e_g = action x path x tfmode, P7 (w = 1)
  R     = - sum_g w_g s_g f(z_g) / sum_g w_g        s_g = sign(logFC), P8 (w = 1 over adj.P < 0.05 genes)
A gene counts ONCE per target unit (plan 3.4). A conflicting or ambiguous gene is excluded from the signed numerator
and counted in coverage. An empty or under-covered gene set returns NA with a reason: **a score of 0 means the
evidence cancelled, NA means there was none**, and the two are never merged.
"""
from dataclasses import dataclass, field

import numpy as np

import config as C


@dataclass
class Score:
    """One score with everything needed to explain or audit it."""
    value: float                       # NaN when not scoreable
    status: str                        # "scored" | "abstained"
    abstain_reason: str = ""
    n_genes: int = 0                   # genes in the signed numerator
    n_landmark: int = 0
    n_inferred: int = 0
    n_excluded_unresolved: int = 0     # conflicted or ambiguous-route genes in the footprint
    n_missing_response: int = 0        # footprint genes absent from the profile
    weight_sum: float = 0.0
    hub_fraction: float = np.nan       # largest single TF's share of the signed weight
    contributions: dict = field(default_factory=dict)   # gene_id -> signed contribution

    @property
    def scored(self):
        return self.status == "scored"


def bounded(z, divisor=None):
    """f(z) = tanh(z / divisor), bounded to [-1, 1] (P9). NaN propagates; it never becomes 0."""
    return np.tanh(np.asarray(z, dtype=float) / (C.F_DIVISOR if divisor is None else divisor))


def _aggregate(gene_ids, expected, z, weights, landmark_mask, min_genes, kind_label,
               n_unresolved=0, tf_of_gene=None):
    """Shared core: signed weighted mean of expected x f(z) over unique, usable genes."""
    gene_ids = np.asarray(gene_ids)
    expected = np.asarray(expected, dtype=float)
    z = np.asarray(z, dtype=float)
    w = np.ones(len(gene_ids)) if weights is None else np.asarray(weights, dtype=float)
    lm = np.zeros(len(gene_ids), bool) if landmark_mask is None else np.asarray(landmark_mask, bool)

    if len(gene_ids) != len(set(gene_ids.tolist())):
        raise ValueError("a gene appears twice: unique-gene aggregation is violated (plan 3.4)")

    usable = (expected != 0) & np.isfinite(z) & np.isfinite(w) & (w > 0)
    n_missing = int((~np.isfinite(z)).sum())
    if not usable.any():
        return Score(np.nan, "abstained", "no_usable_genes", n_excluded_unresolved=n_unresolved,
                     n_missing_response=n_missing)
    if usable.sum() < min_genes:
        return Score(np.nan, "abstained", f"below_floor_{kind_label}:{int(usable.sum())}<{min_genes}",
                     n_genes=int(usable.sum()), n_excluded_unresolved=n_unresolved, n_missing_response=n_missing)

    fz = bounded(z[usable])
    e, ww, g = expected[usable], w[usable], gene_ids[usable]
    contrib = ww * e * fz
    total = float(ww.sum())
    value = float(contrib.sum() / total)
    if not -1.0 - 1e-9 <= value <= 1.0 + 1e-9:
        raise AssertionError(f"score {value} outside [-1, 1]")

    hub = np.nan
    if tf_of_gene is not None:
        by_tf = {}
        for gid, wt in zip(g, ww):
            by_tf[tf_of_gene.get(gid, "")] = by_tf.get(tf_of_gene.get(gid, ""), 0.0) + wt
        hub = max(by_tf.values()) / total if by_tf else np.nan

    return Score(float(np.clip(value, -1, 1)), "scored", "", n_genes=int(usable.sum()),
                 n_landmark=int(lm[usable].sum()), n_inferred=int((~lm[usable]).sum()),
                 n_excluded_unresolved=n_unresolved, n_missing_response=n_missing, weight_sum=total,
                 hub_fraction=hub, contributions=dict(zip(g.tolist(), contrib.tolist())))


def score_mechanism(gene_ids, expected_pre_action, action_sign, z, weights=None, landmark_mask=None,
                    min_genes=None, n_unresolved=0, tf_of_gene=None, action_flag=""):
    """M: does the observed response match the route from the annotated target? (plan 3.2)

    `expected_pre_action` is A1's sign (path x tfmode); `action_sign` is +1 raise / -1 lower. An ambiguous or assumed
    action gives no expectation at all under P11, so M abstains rather than guessing a direction.
    """
    if action_flag in C.M_AMBIGUOUS_FLAGS:
        return Score(np.nan, "abstained", f"ambiguous_action:{action_flag}")
    if action_sign not in (1, -1):
        return Score(np.nan, "abstained", "unknown_action_sign")
    e = np.asarray(expected_pre_action, dtype=float) * action_sign
    return _aggregate(gene_ids, e, z, weights, landmark_mask,
                      C.FLOOR_GENES if min_genes is None else min_genes, "M",
                      n_unresolved=n_unresolved, tf_of_gene=tf_of_gene)


def score_rescue(gene_ids, disease_sign, z, weights=None, landmark_mask=None, min_genes=None,
                 n_unresolved=0, tf_of_gene=None):
    """R: do the target-linked genes move against the disease direction? (plan 3.3)

    Positive R means opposing the disease direction. R carries no action sign, so a unit whose action is ambiguous
    can still be scored for R; the two coverages are tracked separately.
    """
    s = -np.asarray(disease_sign, dtype=float)      # the minus of the definition, applied once
    return _aggregate(gene_ids, s, z, weights, landmark_mask,
                      C.FLOOR_GENES if min_genes is None else min_genes, "R",
                      n_unresolved=n_unresolved, tf_of_gene=tf_of_gene)


def aggregate_units(scores, weights=None):
    """Mean over scoreable target units (P13, unweighted). Abstained units are excluded, never counted as 0."""
    ok = [s for s in scores if s.scored]
    if not ok:
        return Score(np.nan, "abstained", "no_scoreable_unit")
    w = np.ones(len(ok)) if weights is None else np.asarray(weights, dtype=float)[[i for i, s in enumerate(scores) if s.scored]]
    v = float(np.average([s.value for s in ok], weights=w))
    return Score(v, "scored", "", n_genes=sum(s.n_genes for s in ok), n_landmark=sum(s.n_landmark for s in ok),
                 n_inferred=sum(s.n_inferred for s in ok), weight_sum=float(w.sum()))
