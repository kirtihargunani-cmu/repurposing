#!/usr/bin/env python
"""Gate O0 — verify EVERY hand-computable fixture the spec asserts in O-MATH, against the
installed packages. A fixture whose stated answer is wrong is a spec defect, not a code bug.

Spec: docs/SPEC_systems_theoretic_mechanism_validation.md, O-MATH "Required hand-computable
regression tests". Env: sysval. Output: results/systems_validation/synthetic/gate_o0_fixtures.json
"""
import json, sys
from pathlib import Path
import numpy as np, networkx as nx, scipy

ROOT = Path("/home/kirti/projects/repurposing")
OUT = ROOT / "results/systems_validation/synthetic"
TOL = 1e-10
res, fails = [], 0


def ctrb(A, B):
    K = A.shape[0]
    return np.hstack([np.linalg.matrix_power(A, i) @ B for i in range(K)])


def obsv(A, C):
    K = A.shape[0]
    return np.vstack([C @ np.linalg.matrix_power(A, i) for i in range(K)])


def check(name, got, expected, note=""):
    global fails
    ok = got == expected
    fails += (not ok)
    res.append(dict(fixture=name, expected=expected, got=got, pass_=ok, note=note))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:52s} expected {expected!r}, got {got!r}"
          + (f"  — {note}" if note else ""))


print(f"numpy {np.__version__} | scipy {scipy.__version__} | networkx {nx.__version__}\n")
print("O-MATH-2: Kalman rank fixtures")
A = np.array([[0., 1.], [0., 0.]])
check("two-state chain, rank(Kc)", int(np.linalg.matrix_rank(ctrb(A, np.array([[0.], [1.]])))), 2)
check("two-state chain, rank(Ko)", int(np.linalg.matrix_rank(obsv(A, np.array([[1., 0.]])))), 2)
check("wrong actuator B=[1,0]^T, rank(Kc)", int(np.linalg.matrix_rank(ctrb(A, np.array([[1.], [0.]])))), 1)
check("wrong sensor C=[0,1], rank(Ko)", int(np.linalg.matrix_rank(obsv(A, np.array([[0., 1.]])))), 1)

print("\nO-MATH-1: structural matching fixtures")


def unmatched_right(edges, nodes):
    """Spec's own construction: bipartite j_left -> i_right for every permitted A[i,j], i.e. edge j->i."""
    B = nx.Graph()
    B.add_nodes_from([f"{n}L" for n in nodes], bipartite=0)
    B.add_nodes_from([f"{n}R" for n in nodes], bipartite=1)
    B.add_edges_from([(f"{j}L", f"{i}R") for (j, i) in edges])
    m = nx.bipartite.maximum_matching(B, top_nodes=[f"{n}L" for n in nodes])
    matched_right = {v for v in m if v.endswith("R")}
    return sorted(n for n in nodes if f"{n}R" not in matched_right), len(m) // 2


u, msize = unmatched_right([(1, 2), (2, 3)], [1, 2, 3])
check("directed path 1->2->3, unmatched right node(s)", u, [1], f"matching size {msize}")
u2, msize2 = unmatched_right([(1, 2), (2, 3), (3, 1)], [1, 2, 3])
check("directed cycle 1->2->3->1, unmatched right node(s)", u2, [],
      f"perfect matching size {msize2}; >=1 input still required (spec says so)")

print("\nO-MATH-3: reachability / conditioning fixtures")
# chain 1->2->3 dynamics, actuator only at node 3 (non-upstream): states 1,2 unreachable
A3 = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])   # A[i,j] = j->i
check("actuator at non-upstream node only, rank(Kc) < 3",
      bool(np.linalg.matrix_rank(ctrb(A3, np.array([[0.], [0.], [1.]]))) < 3), True,
      f"rank = {int(np.linalg.matrix_rank(ctrb(A3, np.array([[0.],[0.],[1.]]))))}")
check("actuator at upstream node, rank(Kc) == 3",
      int(np.linalg.matrix_rank(ctrb(A3, np.array([[1.], [0.], [0.]])))), 3)
# exact rank but poor conditioning
eps = 1e-8
Aw = np.array([[0., eps], [0., 0.]])
Kw = ctrb(Aw, np.array([[0.], [1.]]))
sv = np.linalg.svd(Kw, compute_uv=False)
check("weak coupling: full numeric rank yet tiny singular value",
      bool(np.linalg.matrix_rank(Kw) == 2 and sv.min() < 1e-6), True,
      f"singular values {sv.round(12).tolist()}")
# Gramian on an unreachable direction must NOT be inverted
from scipy.linalg import expm
B1 = np.array([[0.], [0.], [1.]])
T, n = 1.0, 200
ts = np.linspace(0, T, n)
Wc = sum(expm(A3 * t) @ B1 @ B1.T @ expm(A3.T * t) * (T / n) for t in ts)
check("Wc singular when directions unreachable (do not invert)",
      bool(np.linalg.matrix_rank(Wc, tol=1e-9) < 3), True,
      f"rank(Wc) = {int(np.linalg.matrix_rank(Wc, tol=1e-9))}")

print(f"\nGATE O0: {len(res) - fails}/{len(res)} fixtures behave as the spec states")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "gate_o0_fixtures.json").write_text(json.dumps(dict(
    gate="O0 (hand-computable O-MATH fixtures)", versions=dict(
        numpy=np.__version__, scipy=scipy.__version__, networkx=nx.__version__),
    n_fixtures=len(res), n_pass=len(res) - fails, results=res), indent=1, default=str))
print("wrote", OUT / "gate_o0_fixtures.json")
sys.exit(0)
