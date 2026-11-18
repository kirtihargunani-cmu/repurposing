#!/usr/bin/env python
"""Gate S0 — reproduce the published Angulo/Moog/Liu invariant exactly, symbolically and
numerically, and test the accessibility rank claim the spec makes.

Spec: SPEC_systems_theoretic_mechanism_validation.md S-MATH-1 / S-MATH-2.
Env: sysval. Output: results/systems_validation/synthetic/gate_s0_angulo.json
"""
import json, sys
from pathlib import Path
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp

OUT = Path("/home/kirti/projects/repurposing/results/systems_validation/synthetic")
res, fails = [], 0


def check(name, got, expected, note=""):
    global fails
    ok = (got == expected)
    fails += (not ok)
    res.append(dict(fixture=name, expected=expected, got=got, pass_=ok, note=note))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:55s} {got!r}" + (f"  — {note}" if note else ""))


print(f"sympy {sp.__version__} | numpy {np.__version__}\n")
x1, x2, x3, u1, u3 = sp.symbols("x1 x2 x3 u1 u3", positive=True)

# ---- S-MATH-1: the exact invariant, symbolically -------------------------------------------
print("S-MATH-1: symbolic invariant under x3-only control")
f = sp.Matrix([x1 * (-1 + x3), x2 * (1 - x3), x3 * (sp.Rational(-1, 2) + sp.Rational(3, 2) * x3)])
xi = x1 * x2
g3 = sp.Matrix([0, 0, 1])
xidot = (sp.Matrix([xi]).jacobian([x1, x2, x3]) * (f + g3 * u3))[0]
check("d(x1*x2)/dt with x3-only control simplifies to 0", sp.simplify(xidot) == 0, True,
      f"symbolic result: {sp.simplify(xidot)}")

g1 = sp.Matrix([1, 0, 0])
xidot_u1 = sp.simplify((sp.Matrix([xi]).jacobian([x1, x2, x3]) * (f + g3 * u3 + g1 * u1))[0])
check("with x1 actuation, d(x1*x2)/dt == x2*u1", sp.simplify(xidot_u1 - x2 * u1) == 0, True,
      f"symbolic result: {xidot_u1}")

# ---- S-MATH-1: numerically conserved under several controls --------------------------------
print("\nS-MATH-1: numerical conservation, x(0)=(2,3,1), product must stay 6")
def rhs(t, x, uf):
    a, b, c = x
    return [a * (-1 + c), b * (1 - c), c * (-0.5 + 1.5 * c) + uf(t)]

controls = {"u3=0": lambda t: 0.0, "u3=-0.5 const": lambda t: -0.5,
            "u3=pulse": lambda t: -2.0 if 0.5 < t < 0.8 else 0.0,
            "u3=sin": lambda t: 0.3 * np.sin(5 * t)}
worst = 0.0
for lab, uf in controls.items():
    s = solve_ivp(rhs, [0, 2.0], [2.0, 3.0, 1.0], args=(uf,), rtol=1e-11, atol=1e-13, dense_output=True)
    prod = s.y[0] * s.y[1]
    dev = float(np.max(np.abs(prod - 6.0)))
    worst = max(worst, dev)
    print(f"    {lab:16s} max |x1*x2 - 6| = {dev:.3e}")
check("invariant conserved under all four controls (<1e-6)", bool(worst < 1e-6), True,
      f"worst deviation {worst:.3e}")

# target (1,1,1) has product 1 != 6 -> must be declared unreachable under x3-only control
check("target (1,1,1) product differs from initial -> unreachable by x3 alone",
      bool(abs(1.0 * 1.0 - 6.0) > 1e-9), True, "product 1 vs conserved 6")

# with u1 the product DOES move
def rhs2(t, x):
    a, b, c = x
    return [a * (-1 + c) + 0.5, b * (1 - c), c * (-0.5 + 1.5 * c)]
s2 = solve_ivp(rhs2, [0, 1.0], [2.0, 3.0, 1.0], rtol=1e-11, atol=1e-13)
p_end = float(s2.y[0, -1] * s2.y[1, -1])
check("adding x1 actuation breaks the invariant", bool(abs(p_end - 6.0) > 1e-3), True,
      f"product moved 6 -> {p_end:.4f}")

# ---- S-MATH-2: Lie-bracket accessibility rank ----------------------------------------------
print("\nS-MATH-2: accessibility rank via Lie brackets (spec convention [f,g]=Dg·f - Df·g)")
X = sp.Matrix([x1, x2, x3])
def lie(a, b):
    return b.jacobian(X) * a - a.jacobian(X) * b
fg3 = sp.simplify(lie(f, g3))
M1 = sp.Matrix.hstack(g3, fg3)
r1 = M1.rank()
check("one actuator (x3 only): accessibility rank < 3", bool(r1 < 3), True,
      f"rank(span{{g3,[f,g3]}}) = {r1}")
# the invariant's differential must annihilate the accessible directions
dxi = sp.Matrix([[x2, x1, 0]])
ann = sp.simplify(dxi * M1)
check("d(xi) annihilates the one-actuator accessible directions", ann == sp.zeros(1, 2), True,
      f"dxi·[g3, [f,g3]] = {ann.T.tolist()}")
M2 = sp.Matrix.hstack(g1, g3, fg3)
r2 = sp.simplify(M2.subs({x1: 2, x2: 3, x3: 1})).rank()
check("two actuators (x1 and x3): locally full rank 3 when x2>0", r2, 3,
      "evaluated at x=(2,3,1)")

print(f"\nGATE S0: {len(res)-fails}/{len(res)} fixtures behave as the spec states")
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "gate_s0_angulo.json").write_text(json.dumps(dict(
    gate="S0 (Angulo invariant + accessibility)", versions=dict(sympy=sp.__version__, numpy=np.__version__),
    n_fixtures=len(res), n_pass=len(res)-fails, results=res), indent=1, default=str))
print("wrote", OUT / "gate_s0_angulo.json")
sys.exit(0)
