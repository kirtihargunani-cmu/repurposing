"""D128 A1 pass gate (plan 6.1 test_graph.py): every directed / sign / path case against a hand calculation.

Each expectation below is worked out by hand in its comment. Nothing here reads a drug response: A1 must not depend
on the direction a drug moves.
"""
import itertools
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/d128"))
import graph as G  # noqa: E402

REGS = {"TFA", "TFB", "TFC", "DIRECT"}
REG_GENES = {"TFA": (("ga1", "ga2"), (1, -1)), "TFB": (("gb1",), (1,)), "TFC": (("gc1",), (-1,)),
             "DIRECT": (("gd1", "gd2"), (1, -1))}


_N = itertools.count()


def R(edges, max_hops=2, regs=REGS):
    """One fixture graph. Each gets a unique `key`, so the footprint memo cannot confuse two fixtures."""
    g = G.build_graph(pd.DataFrame(edges, columns=["src", "dst", "sign", "sign_disputed"]))

    def reach(gene):
        return G.reach_regulators(g, gene, max_hops, regs)
    reach.name, reach.max_hops, reach.key = "fixture", max_hops, f"fixture-{next(_N)}"
    return reach


def test_direct_tf_anchor_zero_hops():
    """A target that IS a TF anchors at hop 0 with sign +1: its regulon signs pass through unchanged."""
    fp = G.unit_footprint(["DIRECT"], R([("X", "Y", 1, False)]), REG_GENES, 0)
    assert fp["tf_hop"] == {"DIRECT": 0}
    assert fp["gene_exp"] == {"gd1": 1, "gd2": -1}          # +1 * +1 = +1 ; +1 * -1 = -1


def test_inhibitor_to_activator_one_hop():
    """KIN -(-1)-> TFA, TFA activates ga1 (+1) and represses ga2 (-1).

    Expected before the drug's action: ga1 = -1 * +1 = -1, ga2 = -1 * -1 = +1.
    """
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFA", -1, False)]), REG_GENES, 2)
    assert fp["tf_hop"] == {"TFA": 1}
    assert fp["gene_exp"] == {"ga1": -1, "ga2": 1}


def test_two_hop_sign_multiplies():
    """KIN -(-1)-> MID -(-1)-> TFB(+1 on gb1): (-1)(-1)(+1) = +1."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "MID", -1, False), ("MID", "TFB", -1, False)]), REG_GENES, 2)
    assert fp["tf_hop"] == {"TFB": 2} and fp["gene_exp"] == {"gb1": 1}


def test_hop_cap_excludes_longer_route():
    """The same route at a 1-hop cap reaches nothing: the cap is a cap, not a shortcut."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "MID", -1, False), ("MID", "TFB", -1, False)], max_hops=1), REG_GENES, 1)
    assert fp["tf_hop"] == {} and fp["gene_exp"] == {}


def test_equal_length_opposite_paths_are_ambiguous():
    """Two 1-hop routes of opposite sign to TFB: ambiguous, so gb1 carries no signed expectation."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFB", 1, False), ("KIN", "TFB", -1, False)]), REG_GENES, 2)
    assert fp["n_ambiguous_tfs"] == 1 and fp["gene_exp"] == {}
    assert "gb1" in fp["genes_any"]                          # still counted in coverage


def test_shortest_path_wins_over_longer_opposite():
    """KIN -(+1)-> TFB at 1 hop, and a 2-hop route of opposite sign. Only the shortest counts, so TFB stays +1."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFB", 1, False), ("KIN", "MID", -1, False), ("MID", "TFB", 1, False)]),
                          REG_GENES, 2)
    assert fp["tf_hop"]["TFB"] == 1 and fp["gene_exp"] == {"gb1": 1}


def test_disputed_edge_is_ambiguous_not_a_choice():
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFB", 1, True)]), REG_GENES, 2)
    assert fp["n_ambiguous_tfs"] == 1 and fp["gene_exp"] == {}


def test_cycle_terminates_and_keeps_shortest():
    """A cycle must not loop: TFA at 1 hop keeps hop 1 despite the cycle back through MID."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFA", 1, False), ("TFA", "MID", 1, False), ("MID", "KIN", 1, False),
                                      ("MID", "TFA", 1, False)]), REG_GENES, 2)
    assert fp["tf_hop"]["TFA"] == 1 and fp["gene_exp"] == {"ga1": 1, "ga2": -1}


def test_shared_gene_counted_once_per_unit():
    """TFA and TFB both regulate gs with the SAME expected sign: one entry, not two."""
    rg = {"TFA": (("gs",), (1,)), "TFB": (("gs",), (1,))}
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFA", 1, False), ("KIN", "TFB", 1, False)]), rg, 2)
    assert fp["gene_exp"] == {"gs": 1}
    assert len(fp["genes_resolved"]) == 1


def test_shared_gene_with_opposite_expectations_is_conflicted():
    rg = {"TFA": (("gs",), (1,)), "TFB": (("gs",), (-1,))}
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFA", 1, False), ("KIN", "TFB", 1, False)]), rg, 2)
    assert fp["gene_exp"] == {"gs": 0}                       # 0 = conflicting, never a silent pick
    assert fp["genes_resolved"] == set() and fp["genes_conflicted"] == {"gs"}


def test_multi_gene_complex_is_one_unit():
    """A two-component complex contributes the union of its routes, each gene still counted once."""
    fp = G.unit_footprint(["KIN", "DIRECT"], R([("KIN", "TFB", 1, False)]), REG_GENES, 2)
    assert fp["gene_exp"] == {"gb1": 1, "gd1": 1, "gd2": -1}


def test_isolated_target_yields_nothing_not_zero():
    fp = G.unit_footprint(["ORPHAN"], R([("KIN", "TFB", 1, False)]), REG_GENES, 2)
    assert fp["tf_hop"] == {} and fp["genes_any"] == set()
    assert fp["n_anchored_tfs"] == 0                          # A2 must abstain, not score 0


@pytest.mark.parametrize("flip", [1, -1])
def test_signs_flip_with_the_prior_not_with_any_drug_response(flip):
    """A1 depends only on the prior: flipping every edge sign flips every expectation. No outcome enters."""
    fp = G.unit_footprint(["KIN"], R([("KIN", "TFA", flip, False)]), REG_GENES, 2)
    assert fp["gene_exp"] == {"ga1": flip, "ga2": -flip}
