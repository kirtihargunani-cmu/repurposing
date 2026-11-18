#!/usr/bin/env python
"""Enumerate ALTERNATIVE OPTIMAL CARNIVAL solutions on open solvers, by no-good cuts.

INTERPRETER: /home/kirti/miniforge3/envs/sysval/bin/python   (needs only the stdlib)

WHY THIS EXISTS. The spec (F-MATH, SPEC section 4.8a) requires reporting causal ambiguity: how many
and how different are the explanations that fit the observations equally well. CARNIVAL exposes
solution pooling only through CPLEX, which is commercial and unavailable here; lpSolve has no pooling
options and CBC only a pool gap. Audit 2026-11-18 confirmed the consequence empirically: on a network
where two sources are provably interchangeable, CARNIVAL returned a single confident driver.

METHOD (standard integer-programming technique, no new mathematics). Let x be the vector of node
binaries in CARNIVAL's own LP. Solve; record the optimum z* and the incumbent assignment. Then add
the no-good cut

    sum_{i : x_i = 1} x_i  -  sum_{j : x_j = 0} x_j  <=  |{i : x_i = 1}| - 1

which forbids exactly that assignment and nothing else, and re-solve. Repeat. Every solution whose
objective equals z* (within --tol) is an ALTERNATIVE OPTIMUM. The first objective strictly worse than
z* proves the optimal set is exhausted, and enumeration stops.

Cuts are applied to NODE binaries (nU*, nD*) only, not edge binaries: two solutions that activate the
same nodes through different redundant edges are the same mechanistic explanation, and cutting on
edges would inflate the count with distinctions that carry no biological meaning.

LIMITS, stated plainly:
  * This enumerates optima of the LP CARNIVAL wrote. It cannot reveal explanations that CARNIVAL's
    own formulation excludes.
  * --max-solutions truncates; a run that hits the cap reports `exhausted: false` and its count is a
    LOWER BOUND, never "the number of explanations".
  * Each cut makes the next solve harder. Time grows with the number of solutions found.

Usage:
  enumerate_optima.py --lp FILE.lp --cbc /path/to/cbc [--max-solutions N] [--tol 1e-6] [--out J.json]
"""
import argparse, json, re, shutil, subprocess, sys, tempfile
from pathlib import Path

SECTION = re.compile(r"^\s*(Bounds|Binaries|Binary|Generals|Integers|End)\s*$", re.I)
BIN_START = re.compile(r"^\s*(Binaries|Binary)\s*$", re.I)
BIN_END = re.compile(r"^\s*(Bounds|Generals|Integers|Semi|SOS|End)\s*$", re.I)


def declared_binaries(lp_text: str) -> set:
    """The variables the LP actually declares binary.

    Needed because CARNIVAL also emits CONTINUOUS distance variables named nDs1, nDs2, ... which a
    naive startswith(("nU","nD")) filter captures. Cutting on a non-binary variable is not a no-good
    cut - it silently over-constrains the model (observed: the optimum jumped 0.8 -> 2.0 and
    enumeration stopped after one solution, reporting an ambiguous case as unique).
    """
    out, inside = set(), False
    for ln in lp_text.splitlines():
        if BIN_START.match(ln):
            inside = True
            continue
        if inside:
            if BIN_END.match(ln):
                break
            out.update(ln.split())
    return out


def solve(lp: Path, cbc: str, seconds: int):
    """Run CBC on an LP file; return (objective, {var: value}) or (None, {}) if not optimal."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "sol.txt"
        p = subprocess.run([cbc, str(lp), "seconds", str(seconds), "solve", "solu", str(out)],
                           capture_output=True, text=True)
        if not out.exists():
            return None, {}, p.stdout[-400:]
        lines = out.read_text().splitlines()
        if not lines:
            return None, {}, "empty solution file"
        head = lines[0]
        if "objective value" not in head.lower():
            return None, {}, head
        obj = float(head.rsplit("objective value", 1)[1].strip())
        vals = {}
        for ln in lines[1:]:
            f = ln.split()
            if len(f) >= 3:
                try:
                    vals[f[1]] = float(f[2])
                except ValueError:
                    pass
        return obj, vals, head


def add_cut(lp_text: str, ones, zeros, tag: str) -> str:
    """Insert a no-good cut forbidding exactly the assignment (ones=1, zeros=0)."""
    terms = " + ".join(ones) if ones else ""
    neg = " ".join(f"- {z}" for z in zeros)
    body = " ".join(x for x in (terms, neg) if x)
    cut = f"{tag}:\t {body} <= {len(ones) - 1}\n"
    lines = lp_text.splitlines(keepends=True)
    for i, ln in enumerate(lines):
        if SECTION.match(ln):
            return "".join(lines[:i]) + cut + "".join(lines[i:])
    raise RuntimeError("no Bounds/Binaries/End section found in LP")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lp", required=True)
    ap.add_argument("--cbc", default="/home/kirti/miniforge3/envs/sysval/bin/cbc")
    ap.add_argument("--max-solutions", type=int, default=25)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--seconds", type=int, default=300)
    ap.add_argument("--prefix", default="nU,nD", help="binary prefixes to cut on")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    prefixes = tuple(a.prefix.split(","))

    work = Path(tempfile.mkdtemp()) / "model.lp"
    shutil.copy(a.lp, work)
    text = work.read_text()
    binaries = declared_binaries(text)
    if not binaries:
        print("ERROR: no Binaries section found in the LP", file=sys.stderr)
        return 2
    cuttable = {b for b in binaries if b.startswith(prefixes)}
    if not cuttable:
        print(f"ERROR: no declared binary matches prefixes {prefixes}", file=sys.stderr)
        return 2
    print(f"{len(binaries)} declared binaries; cutting on {len(cuttable)} matching {prefixes}")

    z_star, sols, status = None, [], "exhausted"
    while len(sols) < a.max_solutions:
        work.write_text(text)
        obj, vals, msg = solve(work, a.cbc, a.seconds)
        if obj is None:
            status = "infeasible_or_failed" if not sols else "exhausted"
            break
        if z_star is None:
            z_star = obj
        elif obj > z_star + a.tol:
            status = "exhausted"          # proven: no further solution attains z*
            break
        nodes = {k: v for k, v in vals.items() if k in cuttable}
        ones = sorted(k for k, v in nodes.items() if v > 0.5)
        zeros = sorted(k for k, v in nodes.items() if v <= 0.5)
        sols.append(dict(index=len(sols) + 1, objective=obj, n_active_nodes=len(ones),
                         active=ones))
        if not ones:
            status = "exhausted"
            break
        text = add_cut(text, ones, zeros, f"nogood{len(sols)}")
    else:
        status = "truncated_at_max_solutions"

    # how much do the explanations actually differ?
    sets = [set(s["active"]) for s in sols]
    jac = []
    for i in range(len(sets)):
        for j in range(i + 1, len(sets)):
            u = sets[i] | sets[j]
            jac.append(len(sets[i] & sets[j]) / len(u) if u else 1.0)
    always = set.intersection(*sets) if sets else set()
    ever = set.union(*sets) if sets else set()
    res = dict(lp=str(a.lp), optimum=z_star, n_declared_binaries=len(binaries),
               n_cuttable_binaries=len(cuttable), n_alternative_optima=len(sols), status=status,
               exhausted=(status == "exhausted"),
               core_nodes_in_every_solution=sorted(always),
               nodes_in_some_solution=sorted(ever),
               n_core=len(always), n_union=len(ever),
               median_pairwise_jaccard=(sorted(jac)[len(jac) // 2] if jac else None),
               solutions=sols,
               method="no-good cuts on node binaries; CPLEX solution pooling unavailable",
               cut_prefixes=list(prefixes))
    print(f"optimum {z_star}  alternative optima: {len(sols)}  status: {status}")
    print(f"  nodes active in EVERY optimum : {len(always)}")
    print(f"  nodes active in SOME optimum  : {len(ever)}")
    if jac:
        print(f"  median pairwise Jaccard       : {res['median_pairwise_jaccard']:.3f}")
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        Path(a.out).write_text(json.dumps(res, indent=1))
        print("wrote", a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
