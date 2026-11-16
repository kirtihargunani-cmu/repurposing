"""D128-A3.2: the sign convention of M_disease, pinned against score.py's M and R.

M_disease is M's own formula fed the DISEASE contrast instead of the drug profile, so its sign is only meaningful
if it is read the same way round as M and R. The cases below are hand-calculated and fix the reading:

    M_disease > 0  the tumor already moved the way an inhibitor would move it -> the reference is target-ACTIVE
                   relative to the tumor -> a drug that works (M > 0) is scored NEGATIVELY by R  (reference_inverted)
    M_disease < 0  the tumor moved against the inhibited direction -> the drug reverses it toward normal (aligned)
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/d128"))
import score as S  # noqa: E402

G = [f"d{i}" for i in range(12)]
SUB = ["tf"] * 12


def _m_disease(expected_pre_action, action, z_disease):
    """Exactly what a3_reference_activity.py computes: M with the disease contrast in place of the drug profile."""
    return S.score_mechanism(G, expected_pre_action, action, z_disease, subset=SUB, n_signaling_tfs=0)


def test_m_disease_is_positive_when_the_tumor_already_looks_target_inhibited():
    """Target-driven genes (expected UP if the target is inhibited) are already UP in the tumor vs the reference.

    e = +1 on every gene after the action sign; disease z = +12 on every gene; f(12/2) = tanh(6) ~ 1.
    M_disease = mean(+1 * +1) = +1 -> the reference runs the pathway; R cannot be read here.
    """
    m = _m_disease([-1] * 12, -1, [12.0] * 12)          # pre-action -1 x action -1 = +1
    assert m.scored and m.value == pytest.approx(np.tanh(6.0))


def test_m_disease_is_negative_when_the_tumor_looks_target_active():
    """The same footprint, with the tumor moved the other way: the drug would reverse it toward normal."""
    m = _m_disease([-1] * 12, -1, [-12.0] * 12)
    assert m.scored and m.value == pytest.approx(-np.tanh(6.0))


def test_the_inverted_case_makes_a_working_drug_score_negative_on_R():
    """The whole point of the flag, in one assertion.

    Disease: every footprint gene UP in the tumor (sign(logFC) = +1), and an inhibitor is expected to push them UP
    too (e = +1) - that is M_disease > 0. A drug that does exactly what its target annotation says (z = +2 on every
    gene) therefore gets M > 0 and R < 0 on the SAME genes: R penalises the mechanistically correct drug.
    """
    e_pre, action = [-1] * 12, -1                       # e = +1
    md = _m_disease(e_pre, action, [12.0] * 12)
    m = S.score_mechanism(G, e_pre, action, [2.0] * 12, subset=SUB, n_signaling_tfs=0)
    r = S.score_rescue(G, [1] * 12, [2.0] * 12)         # disease sign +1 on every gene
    assert md.value > 0 and m.value > 0 and r.value < 0
    assert m.value == pytest.approx(-r.value)           # with e == disease sign, R is exactly -M

    # and the mirror image: when the tumor moved against the inhibited direction, M and R agree in sign
    md2 = _m_disease(e_pre, action, [-12.0] * 12)
    r2 = S.score_rescue(G, [-1] * 12, [2.0] * 12)
    assert md2.value < 0 and m.value > 0 and r2.value > 0


def test_m_disease_obeys_the_same_frozen_floor_as_m():
    """A unit the frozen P10 floor will not admit cannot carry a flag either."""
    g = [f"e{i}" for i in range(20)]
    m = S.score_mechanism(g, [1] * 20, -1, [10.0] * 20, subset=["signaling"] * 20, n_signaling_tfs=1)
    assert m.status == "abstained" and m.abstain_reason.startswith("below_floor_M_unit:")
