#!/usr/bin/env python
"""Per-stage runtime history and ETA for the formal pipeline (D187, 2026-11-22).

INTERPRETER: any; stdlib only, so it runs in either env.

User, 2026-11-22: *"let's make sure that we log the runtime of each major stage as well, it will be
useful to have estimates as we are most likely to package this pipeline as one that is hopefully
consistent, reproducible over other diseases."*

WHAT THIS IS FOR. Absolute seconds do not transfer to another disease set. A stage that takes 600 s
over 41 cohorts and 3,742 drugs tells you nothing about the same stage over 12 cohorts until it is
expressed PER UNIT OF WORK. So every stage is reported three ways:

    total       wall-clock seconds, as measured
    per_group   seconds / scored cohort      (stages whose work scales with cohorts)
    per_drug    ms / scored drug             (stages whose work scales with the compound set)

Which denominator is right differs per stage and is NOT guessed here - both are printed, and the
one that is stable across runs of different size is the one that generalises. With a single run
behind it, neither is yet established; this file is the instrument, not the answer.

SOURCE. `logs/BUILD_LOG.tsv`, which `scripts/runlog.sh` appends to with two rows per run (start,
end) and a `dur_s` on the end row. Append-only, so a start row with no end row IS a crash record
and is reported as such rather than being silently dropped from the medians.

KEYING. A stage is keyed by (script, env signature), not by script alone, because several scripts
run more than once per pipeline with different environments - `79e` for groups and strata, `79e3`
for the default and the ALL=1 pass. Where runlog does not record the distinguishing variable the
runs collapse into one key; those keys are flagged `AMBIGUOUS` rather than averaged silently.

    scripts/stage_timings.py            history table
    scripts/stage_timings.py --eta      estimated total for a full run_pipeline.sh
    scripts/stage_timings.py --write    also write results/lincs/pipeline_timings.tsv
"""
import csv
import re
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEDGER = ROOT / "logs/BUILD_LOG.tsv"
RUNNER = ROOT / "scripts/run_pipeline.sh"
OUT = ROOT / "results/lincs/pipeline_timings.tsv"

# env vars that distinguish two runs of the SAME script within one pipeline
KEYVARS = ("STRATA", "SFX", "SPLIT_SFX", "SFXS", "ALL", "NO_TARGETS", "EDGES")


def env_sig(env: str) -> str:
    """The subset of the recorded env that distinguishes stages of the same script."""
    got = dict(p.split("=", 1) for p in env.split() if "=" in p)
    return " ".join(f"{k}={got[k]}" for k in KEYVARS if k in got)


def read_ledger():
    """-> (done, crashed): done is {(script, sig): [dur_s]}, crashed is [(run_id, script)]."""
    if not LEDGER.exists():
        print(f"no ledger at {LEDGER}", file=sys.stderr)
        return {}, []
    starts, done, ends = {}, defaultdict(list), set()
    rows = list(csv.DictReader(LEDGER.open(), delimiter="\t"))
    for r in rows:
        key = (r["script"], env_sig(r.get("env", "") or ""))
        if r["event"] == "start":
            starts[r["run_id"]] = key
        elif r["event"] == "end":
            ends.add(r["run_id"])
            if r["status"] == "ok" and r.get("dur_s"):
                try:
                    done[key].append(int(r["dur_s"]))
                except ValueError:
                    pass
    crashed = [(rid, k[0]) for rid, k in starts.items() if rid not in ends]
    return done, crashed


def scale():
    """Current problem size, for the per-unit rates. Missing inputs -> None, never a guess."""
    s = {}
    f = ROOT / "results/lincs/drug_reversal_v2_tf.tsv"
    if f.exists():
        groups, drugs = set(), set()
        with f.open() as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for r in rd:
                groups.add(r.get("group_id"))
                drugs.add(r.get("pert_id"))
        s["n_groups"], s["n_drugs"] = len(groups), len(drugs)
    f = ROOT / "results/lincs/stratum_reversal_v2_tf.tsv"
    if f.exists():
        with f.open() as fh:
            s["n_strata"] = len({r.get("stratum") for r in csv.DictReader(fh, delimiter="\t")})
    return s


def stages():
    """The ordered stage list, parsed from run_pipeline.sh so there is ONE source of truth."""
    if not RUNNER.exists():
        return []
    text = re.sub(r"\\\n\s*", " ", RUNNER.read_text())   # join line continuations
    out = []
    for line in text.splitlines():
        m = re.match(r'\s*run\s+"([^"]+)"\s+(.*?)\s+--\s+(.*)$', line)
        if not m:
            continue
        label, cmd = m.group(1), m.group(3)
        script = next((Path(a).name for a in reversed(cmd.split())
                       if a.endswith((".py", ".R"))), None)
        sig = env_sig(" ".join(p for p in cmd.split() if "=" in p and not p.startswith("/")))
        if script:
            out.append((label, script, sig))
    return out


def main() -> int:
    done, crashed = read_ledger()
    sc = scale()
    ng, nd = sc.get("n_groups"), sc.get("n_drugs")

    print(f"scale: {sc if sc else 'unknown (core tables not built yet)'}\n")
    print(f"{'stage (label from run_pipeline.sh)':42s} {'n':>2s} {'median':>8s} {'max':>8s} "
          f"{'s/group':>8s} {'ms/drug':>8s}")
    print("-" * 84)

    # Stages of the same script that the ledger cannot tell apart. runlog.sh records the env it
    # sees in ITS OWN process, but run_pipeline.sh passes per-stage variables as `env K=V -- cmd`,
    # which reach only the child - so e.g. the groups and strata runs of 50b collapse to one key.
    # They are reported with the shared history and flagged, never silently averaged.
    per_script = defaultdict(list)
    for label, script, _ in stages():
        per_script[script].append(label)

    total, unknown = 0, []
    for label, script, sig in stages():
        v, ambiguous = done.get((script, sig)), False
        if not v:
            # The ledger has no run under this stage's signature. Fall back to the script's runs
            # regardless of signature - correct when the script appears once in the pipeline,
            # and an unresolved collapse when it appears more than once.
            cand = [d for (s2, _), d in done.items() if s2 == script]
            v = [x for d in cand for x in d] or None
            ambiguous = bool(v) and len(per_script[script]) > 1
        if not v:
            print(f"{label:42s} {'-':>2s} {'NO HISTORY':>8s}")
            unknown.append(label)
            continue
        med = st.median(v)
        total += med
        pg = f"{med/ng:8.1f}" if ng else f"{'-':>8s}"
        pd_ = f"{med/nd*1000:8.1f}" if nd else f"{'-':>8s}"
        tag = "  AMBIGUOUS (shared with the other stages of this script)" if ambiguous else ""
        print(f"{label:42s} {len(v):2d} {med:7.0f}s {max(v):7d}s {pg} {pd_}{tag}")

    print("-" * 84)
    h, m = divmod(int(total) // 60, 60)
    print(f"{'TOTAL of stages with history':42s}    {total:7.0f}s  ({h}h{m:02d}m)")
    if unknown:
        print(f"\n{len(unknown)} stage(s) with NO recorded runtime - the total above EXCLUDES them:")
        for u in unknown:
            print(f"    {u}")
    if crashed:
        print(f"\n{len(crashed)} run(s) with a start row and no end row (crashed, or still running):")
        for rid, s2 in crashed[-5:]:
            print(f"    {rid}  {s2}")

    if "--write" in sys.argv:
        OUT.parent.mkdir(parents=True, exist_ok=True)
        with OUT.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["label", "script", "env_sig", "n_runs", "median_s", "max_s",
                        "s_per_group", "ms_per_drug", "n_groups", "n_drugs"])
            for label, script, sig in stages():
                v = done.get((script, sig)) or done.get((script, ""))
                if not v:
                    w.writerow([label, script, sig, 0, "", "", "", "", ng or "", nd or ""])
                    continue
                med = st.median(v)
                w.writerow([label, script, sig, len(v), f"{med:.0f}", max(v),
                            f"{med/ng:.2f}" if ng else "", f"{med/nd*1000:.2f}" if nd else "",
                            ng or "", nd or ""])
        print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
