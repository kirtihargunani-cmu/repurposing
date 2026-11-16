"""D128 A2 pass gate (plan 6.1 test_scores.py / test_abstain.py / test_nulls.py): hand-calculated M and R.

Each expected value below is computed by hand in its comment. The gate also requires: scores bounded to [-1, 1], a
score of zero distinguishable from NA, and no confident score for a unit with no genes.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/d128"))
import nulls as N  # noqa: E402
import score as S  # noqa: E402

G10 = [f"g{i}" for i in range(10)]


def test_bounded_transform_matches_tanh_and_is_bounded():
    assert S.bounded(0) == 0.0
    assert S.bounded(2) == pytest.approx(np.tanh(1.0))        # tanh(2/2) = tanh(1) = 0.7616
    assert -1 < S.bounded(1e6) <= 1 and -1 <= S.bounded(-1e6) < 0
    assert np.isnan(S.bounded(np.nan))                         # NaN never becomes 0


def test_m_perfect_agreement_is_one():
    """All 10 genes expected up (pre-action +1, action -1 -> expected -1) and all respond at z = -1e6.

    M = sum(1 * -1 * f(-1e6)) / 10 = sum(1 * -1 * -1) / 10 = 1.
    """
    m = S.score_mechanism(G10, [1] * 10, -1, [-1e6] * 10)
    assert m.scored and m.value == pytest.approx(1.0) and m.n_genes == 10


def test_m_perfect_disagreement_is_minus_one():
    m = S.score_mechanism(G10, [1] * 10, -1, [1e6] * 10)
    assert m.value == pytest.approx(-1.0)


def test_m_hand_calculated_mixed_case():
    """6 genes expected +1, 4 expected -1, action +1; z = +2 for all.

    f(2) = tanh(1) = 0.761594. M = (6*(+1) + 4*(-1)) * 0.761594 / 10 = 0.2 * 0.761594 = 0.1523188.
    """
    m = S.score_mechanism(G10, [1] * 6 + [-1] * 4, 1, [2.0] * 10)
    assert m.value == pytest.approx(0.2 * np.tanh(1.0), rel=1e-9)
    assert m.contributions["g0"] == pytest.approx(np.tanh(1.0))


def test_m_action_sign_flips_the_score():
    a = S.score_mechanism(G10, [1] * 10, 1, [2.0] * 10)
    b = S.score_mechanism(G10, [1] * 10, -1, [2.0] * 10)
    assert a.value == pytest.approx(-b.value)


def test_zero_score_is_not_na():
    """5 genes agree, 5 disagree: the evidence cancels. The score is 0.0 and SCORED, not NA."""
    m = S.score_mechanism(G10, [1] * 5 + [-1] * 5, 1, [2.0] * 10)
    assert m.scored and m.value == pytest.approx(0.0)
    assert not np.isnan(m.value)                               # distinguishable from an abstention
    nn = S.score_mechanism([], [], 1, [])
    assert nn.status == "abstained" and np.isnan(nn.value)


def test_no_genes_never_gets_a_confident_score():
    m = S.score_mechanism([], [], -1, [])
    assert m.status == "abstained" and m.abstain_reason == "no_usable_genes" and m.n_genes == 0


def test_below_floor_abstains_with_a_reason():
    m = S.score_mechanism(G10[:4], [1] * 4, -1, [-2.0] * 4)
    assert m.status == "abstained" and m.abstain_reason.startswith("below_floor_M:4<10")


def test_ambiguous_action_abstains_rather_than_assuming():
    m = S.score_mechanism(G10, [1] * 10, -1, [-2.0] * 10, action_flag="direction_assumed")
    assert m.status == "abstained" and "ambiguous_action" in m.abstain_reason


def test_conflicting_genes_excluded_from_numerator_but_counted():
    """Expected 0 = conflicting: excluded from the signed sum, still reported in coverage."""
    m = S.score_mechanism(G10, [1] * 10 + [], 1, [2.0] * 10, n_unresolved=7)
    assert m.n_genes == 10 and m.n_excluded_unresolved == 7
    m2 = S.score_mechanism(G10, [1] * 6 + [0] * 4, 1, [2.0] * 10)
    assert m2.n_genes == 6                                     # the four zeros do not enter


def test_missing_response_is_not_zero():
    z = [2.0] * 7 + [np.nan] * 3
    m = S.score_mechanism(G10, [1] * 10, 1, z, min_genes=5)
    assert m.n_genes == 7 and m.n_missing_response == 3
    assert m.value == pytest.approx(np.tanh(1.0))              # the 3 missing do not dilute toward 0


def test_duplicate_gene_raises():
    with pytest.raises(ValueError, match="twice"):
        S.score_mechanism(["g1", "g1"] + G10[:8], [1] * 10, 1, [2.0] * 10)


def test_weight_scaling_invariance():
    a = S.score_mechanism(G10, [1] * 10, 1, [2.0] * 10, weights=[1] * 10)
    b = S.score_mechanism(G10, [1] * 10, 1, [2.0] * 10, weights=[7] * 10)
    assert a.value == pytest.approx(b.value)


def test_r_opposes_the_disease_direction():
    """Disease up (+1) and the drug pushes down (z = -1e6): R = -(+1 * -1) = +1."""
    r = S.score_rescue(G10, [1] * 10, [-1e6] * 10)
    assert r.value == pytest.approx(1.0)
    r2 = S.score_rescue(G10, [1] * 10, [1e6] * 10)
    assert r2.value == pytest.approx(-1.0)


def test_r_hand_calculated():
    """7 genes up, 3 down; z = -2 for all. f(-2) = -tanh(1).

    R = -[7*(+1)*(-0.7616) + 3*(-1)*(-0.7616)]/10 = -[-5.331 + 2.285]/10 = 0.3046.
    """
    r = S.score_rescue(G10, [1] * 7 + [-1] * 3, [-2.0] * 10)
    assert r.value == pytest.approx(0.4 * np.tanh(1.0), rel=1e-9)


def test_r_scoreable_when_action_is_ambiguous():
    """R carries no action sign, so an ambiguous action blocks M but not R (separate coverages)."""
    m = S.score_mechanism(G10, [1] * 10, -1, [-2.0] * 10, action_flag="direction_assumed")
    r = S.score_rescue(G10, [1] * 10, [-2.0] * 10)
    assert m.status == "abstained" and r.scored


def test_hub_fraction_reported():
    tf = {g: ("HUB" if i < 8 else "OTHER") for i, g in enumerate(G10)}
    m = S.score_mechanism(G10, [1] * 10, 1, [2.0] * 10, tf_of_gene=tf)
    assert m.hub_fraction == pytest.approx(0.8)


def test_unit_aggregation_skips_abstentions():
    a = S.score_mechanism(G10, [1] * 10, -1, [-1e6] * 10)      # +1
    b = S.score_mechanism([], [], -1, [])                      # abstained
    agg = S.aggregate_units([a, b])
    assert agg.scored and agg.value == pytest.approx(1.0)      # not averaged with a 0
    assert S.aggregate_units([b, b]).status == "abstained"


# --- nulls ---------------------------------------------------------------------------------------------------
def test_null_permutations_are_seeded_and_reproducible():
    e = [1, -1] * 5
    assert list(N.shuffled_edge_signs(e)) == list(N.shuffled_edge_signs(e))
    z = np.arange(10.0)
    assert list(N.gene_label_permutation(z)) == list(N.gene_label_permutation(z))


def test_sign_shuffle_preserves_the_sign_balance():
    e = [1] * 6 + [-1] * 4
    s = N.shuffled_edge_signs(e)
    assert sorted(s) == sorted(np.array(e, float))


def test_gene_label_permutation_stays_within_feature_space():
    z = np.array([10.0, 11.0, 12.0, 1.0, 2.0, 3.0])
    lm = np.array([True] * 3 + [False] * 3)
    p = N.gene_label_permutation(z, lm)
    assert set(p[:3]) == {10.0, 11.0, 12.0} and set(p[3:]) == {1.0, 2.0, 3.0}


def test_degree_matched_units_match_on_size_and_tfs():
    c = pd.DataFrame(dict(unit_id=["u1", "u2", "u3", "u4"], n_sign_resolved=[100, 105, 1000, 95],
                          n_resolved_tfs=[10, 11, 90, 9]))
    m = N.degree_matched_units("u1", c, 100, 10)
    assert set(m.unit_id) == {"u2", "u4"}                       # u3 is an order of magnitude larger


def test_null_scoring_does_not_alter_the_real_target():
    e = np.array([1.0] * 10)
    before = e.copy()
    N.shuffled_edge_signs(e)
    assert list(e) == list(before)


def test_null_summary_percentile():
    s = N.summarise([0.0, 0.1, 0.2, 0.3], 0.25)
    assert s["n_null"] == 4 and s["percentile"] == pytest.approx(75.0)


# --- D128-A3: the vectorised null fast path must be the SAME number, not an approximation -------------------------
def test_fast_nulls_equal_score_mechanism():
    """nulls.m_null_* must reproduce the per-draw score_mechanism path exactly where the shortcut applies."""
    rng = np.random.default_rng(7)
    n = 60
    e = rng.choice([-1.0, 1.0], n)
    z = rng.normal(0, 2, n)
    lm = rng.random(n) < 0.3
    assert N.exact_shortcut_ok(e, z)
    fz = S.bounded(z)
    slow = [S.score_mechanism(list(range(n)), N.shuffled_edge_signs(e, i), 1, z, landmark_mask=lm).value
            for i in range(25)]
    assert np.allclose(N.m_null_sign_shuffle(e, fz, 25), slow, atol=1e-12)
    slow = [S.score_mechanism(list(range(n)), e, 1, N.gene_label_permutation(z, lm, i), landmark_mask=lm).value
            for i in range(25)]
    assert np.allclose(N.m_null_gene_label(e, fz, lm, 25), slow, atol=1e-12)


def test_fast_null_shortcut_refused_when_not_exact():
    """A zero expected sign or a missing response must send the caller back to the per-draw path."""
    assert not N.exact_shortcut_ok([1.0, 0.0, -1.0], [1.0, 1.0, 1.0])
    assert not N.exact_shortcut_ok([1.0, -1.0], [1.0, np.nan])
    assert not N.exact_shortcut_ok([], [])
