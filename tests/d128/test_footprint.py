"""D128 A0 footprint counting (E01, E03, E04, E05): zero-hop TF anchor, sign composition, ambiguity, unique genes."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/d128"))
import a0_audit as A  # noqa: E402

EDGES = pd.DataFrame(dict(
    src=["KIN", "KIN", "MID", "KIN", "MID2", "KIN"],
    dst=["TF1", "MID", "TF2", "MID2", "TF3", "TF3"],
    sign=[1, -1, 1, 1, 1, -1],
    sign_disputed=[False, False, False, False, False, False]))
REGS = {"TF1", "TF2", "TF3", "ESR1"}
# regulon genes, already restricted to the measured gene space (as a0_audit precomputes them)
REG_GENES = {"TF1": (("g1", "g2"), (1, -1)), "TF2": (("g1",), (1,)), "TF3": (("g3",), (1,)),
             "ESR1": (("g4", "g5"), (1, -1))}


class FakeReach:
    """Stands in for the persisted Reach cache: same call signature, no disk."""

    def __init__(self, name, edges, max_hops=2):
        self.name, self.graph, self.max_hops = name, A.build_graph(edges), max_hops

    def __call__(self, gene):
        return A.reach_regulators(self.graph, gene, self.max_hops, REGS)


def test_zero_hop_tf_anchor():                      # E01: a TF target is its own anchor
    fp = A.unit_footprint(["ESR1"], FakeReach("t_zero", EDGES), REG_GENES, 2)
    assert fp["tf_hop"] == {"ESR1": 0}
    assert fp["genes_resolved"] == {"g4", "g5"} and fp["hop0_resolved"] == {"g4", "g5"}


def test_sign_composition_and_conflict():         # E03 composition; E05 g1 counted once; conflicting signs
    fp = A.unit_footprint(["KIN"], FakeReach("t_sign", EDGES), REG_GENES, 2)
    # TF1 at 1 hop (+), TF2 at 2 hops via MID (- x + = -) -> g1 expected + from TF1 and - from TF2 -> conflicted
    assert fp["tf_hop"]["TF1"] == 1 and fp["tf_hop"]["TF2"] == 2
    assert "g1" in fp["genes_conflicted"] and "g1" not in fp["genes_resolved"]
    assert "g2" in fp["genes_resolved"]


def test_opposite_paths_same_hop_are_ambiguous():  # E04: KIN->TF3 at 1 hop with both signs
    e = pd.concat([EDGES, pd.DataFrame(dict(src=["KIN"], dst=["TF3"], sign=[1], sign_disputed=[False]))])
    fp = A.unit_footprint(["KIN"], FakeReach("t_amb", e), REG_GENES, 2)
    assert fp["n_ambiguous_tfs"] >= 1
    assert "g3" in fp["genes_any"] and "g3" not in fp["genes_resolved"]   # kept in coverage, out of the signed set


def test_disputed_edge_is_ambiguous():
    e = EDGES.copy(); e.loc[0, "sign_disputed"] = True
    fp = A.unit_footprint(["KIN"], FakeReach("t_disp", e), REG_GENES, 1)
    assert fp["n_ambiguous_tfs"] >= 1


def test_hop_arm_filters_without_rewalking():
    r = FakeReach("t_hop", EDGES, max_hops=3)
    assert A.unit_footprint(["KIN"], r, REG_GENES, 1)["n_anchored_tfs"] == 2      # TF1, TF3 at 1 hop
    assert A.unit_footprint(["KIN"], r, REG_GENES, 2)["n_anchored_tfs"] == 3      # + TF2 at 2 hops


def test_sign_mask_product():
    assert A.mul(A.PLUS, A.MINUS) == A.MINUS and A.mul(A.MINUS, A.MINUS) == A.PLUS
    assert A.mul(A.BOTH, A.PLUS) == A.BOTH


# --- D128-A0.1: per-gene typing, floors per subset, union for the unit -------------------------------------------
def test_mixed_unit_clears_on_its_signaling_gene():
    """A unit of one regulator with a tiny regulon and one kinase must not be floored on the regulator alone."""
    regs = {"TFSMALL", "TF1", "TF2", "TF3"}
    reg_genes = dict(REG_GENES, TFSMALL=(("x1",), (1,)))      # 1 regulon gene: below the floor of 10
    big = {f"b{i}" for i in range(12)}
    reg_genes["TF1"] = (tuple(sorted(big)), tuple([1] * 12))   # reached from KIN at 1 hop
    r = FakeReach("t_mixed", EDGES, max_hops=2)
    r_regs = regs

    def reach(gene):
        return A.reach_regulators(r.graph, gene, 2, r_regs)
    reach.name, reach.max_hops = "t_mixed_fn", 2
    kind, tf_fp, sig_fp = A.unit_subsets(["TFSMALL", "KIN"], reach, reg_genes, regs, 2)
    assert kind == "mixed"
    assert len(tf_fp["genes_resolved"]) == 1                   # the tf subset alone is below the floor
    ok, which = A.unit_floor(tf_fp, sig_fp, 1.0)
    assert ok and which == "signaling"                         # the signaling subset carries it
    any_g, res, conf = A.merge_subsets(tf_fp, sig_fp)
    assert "x1" in res and big <= res                           # union, each gene once
    assert len(res) == len(set(res))


def test_union_never_loses_a_pure_subset_unit():
    regs = {"ESR1", "TF1", "TF2", "TF3"}
    r = FakeReach("t_pure", EDGES, max_hops=2)

    def reach(gene):
        return A.reach_regulators(r.graph, gene, 2, regs)
    reach.name, reach.max_hops = "t_pure_fn", 2
    kind, tf_fp, sig_fp = A.unit_subsets(["ESR1"], reach, REG_GENES, regs, 2)
    assert kind == "tf" and sig_fp is None
    assert A.unit_floor(tf_fp, sig_fp, 0.1)[1] == "tf"         # clears on its own regulon at a low floor
