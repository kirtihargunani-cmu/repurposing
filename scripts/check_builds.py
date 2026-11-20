#!/usr/bin/env python
"""Audit build provenance: crashed runs, unrecorded outputs, and stale/mixed-age arms.

INTERPRETER: any python3 (stdlib only)

    python scripts/check_builds.py            # audit everything registered
    python scripts/check_builds.py --runs     # just the run ledger (crashes, died)
    python scripts/check_builds.py --strict   # exit 1 if any ERROR-level finding

Checks, in order of how badly they have bitten us:

  MIXED_AGE   files in one output dir span more than MAX_SPAN_H hours. This is the
              D162 signature: drug_viper_raw_v2_tf held 28 files from 11-02 while its
              sibling arm was rebuilt on 11-20, and 79e reads BOTH arms per network.
  NOT_REWRITTEN  the manifest records files the latest run did not touch -- a partial
              rebuild left older files in place under current-looking names.
  COUNT       file count != the run's --expect.
  UNRECORDED  files changed after the manifest was written, i.e. something wrote here
              outside a logged run. Provenance is no longer trustworthy.
  NO_MANIFEST registered output has no manifest at all (pre-dates this machinery).
  DIED        ledger has a start row with no end row -- the process was killed.
  CRASH       ledger end row with a non-zero exit code.

Registry: config/build_outputs.tsv (path, owner_script, expect, note). Anything with a
manifest is audited whether registered or not.
"""
import json, os, sys, csv
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/home/kirti/projects/repurposing")
LEDGER = ROOT / "logs/BUILD_LOG.tsv"
REGISTRY = ROOT / "config/build_outputs.tsv"
MAX_SPAN_H = 24.0          # files in one build should land within a day of each other
MANIFEST = ".build_manifest.json"

ERROR, WARN, INFO = "ERROR", "WARN", "INFO"
RANK = {ERROR: 0, WARN: 1, INFO: 2}


def ago(ts):
    d = datetime.now(timezone.utc) - datetime.fromisoformat(ts).astimezone(timezone.utc)
    return f"{d.days}d" if d.days else f"{d.seconds // 3600}h"


def registry():
    out = []
    if REGISTRY.exists():
        with open(REGISTRY) as f:
            for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
                out.append(r)
    return out


def check_runs():
    f = []
    if not LEDGER.exists():
        return [(INFO, "ledger", "no logs/BUILD_LOG.tsv yet -- nothing has been run under scripts/runlog.sh")]
    with open(LEDGER) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    starts = {r["run_id"]: r for r in rows if r["event"] == "start"}
    ends = {r["run_id"]: r for r in rows if r["event"] == "end"}
    for rid, s in starts.items():
        if rid not in ends:
            f.append((ERROR, f"run {rid}", f"DIED: {s['script']} started {ago(s['ts'])} ago, no end row "
                                           f"(killed/OOM/reboot). env='{s['env']}' log={s['log']}"))
    for rid, e in ends.items():
        if e["status"] != "ok":
            f.append((ERROR, f"run {rid}", f"CRASH: {e['script']} exit={e['exit_code']} after {e['dur_s']}s "
                                           f"env='{e['env']}' log={e['log']}"))
    return f


def check_output(path, expect=None, owner=None):
    p = ROOT / path if not os.path.isabs(path) else Path(path)
    f = []
    if not p.exists():
        return [(WARN, path, f"registered output does not exist (owner: {owner})")]
    mf = p / MANIFEST if p.is_dir() else Path(str(p) + MANIFEST)
    if p.is_dir():
        files = [q for q in p.iterdir() if q.is_file() and not q.name.startswith(".build_manifest")]
    else:
        files = [p]
    if not files:
        return [(WARN, path, "output is empty")]
    mt = {q.name: q.stat().st_mtime for q in files}
    span = (max(mt.values()) - min(mt.values())) / 3600
    if span > MAX_SPAN_H:
        old = sorted(mt, key=mt.get)[:4]
        new = sorted(mt, key=mt.get)[-1]
        f.append((ERROR, path, f"MIXED_AGE: {len(files)} files span {span:.0f}h "
                               f"({datetime.fromtimestamp(min(mt.values())):%m-%d} .. "
                               f"{datetime.fromtimestamp(max(mt.values())):%m-%d}). "
                               f"oldest: {', '.join(old)}{' ...' if len(files) > 4 else ''}; newest: {new}. "
                               f"A consumer reading all of these mixes builds."))
    if expect and str(expect).strip() and len(files) != int(expect):
        f.append((ERROR, path, f"COUNT: {len(files)} files, expected {expect}"))
    if not mf.exists():
        f.append((WARN, path, f"NO_MANIFEST: {len(files)} files with no provenance "
                              f"(newest {datetime.fromtimestamp(max(mt.values())):%Y-%m-%d %H:%M}). "
                              f"Re-run under scripts/runlog.sh, or register it as known-unprovenanced."))
        return f
    try:
        man = json.load(open(mf))
    except Exception as e:
        return f + [(ERROR, path, f"manifest unreadable: {e}")]
    nr = man.get("files_not_rewritten_by_this_run") or []
    if nr:
        f.append((ERROR, path, f"NOT_REWRITTEN: run {man.get('run_id')} left {len(nr)} pre-existing file(s) "
                               f"untouched: {', '.join(nr[:6])}{' ...' if len(nr) > 6 else ''}. Partial rebuild."))
    fin = man.get("finished_at")
    if fin:
        ft = datetime.fromisoformat(fin).timestamp()
        newer = [n for n, t in mt.items() if t > ft + 2]
        if newer:
            f.append((ERROR, path, f"UNRECORDED: {len(newer)} file(s) written after the manifest "
                                   f"({', '.join(sorted(newer)[:5])}{' ...' if len(newer) > 5 else ''}) -- "
                                   f"something wrote here outside a logged run."))
    if man.get("n_files") != len(files):
        f.append((WARN, path, f"manifest says {man.get('n_files')} files, found {len(files)}"))
    if not f:
        f.append((INFO, path, f"ok: {len(files)} files, run {man.get('run_id')}, "
                              f"{ago(fin) if fin else '?'} old, env='{man.get('env') or '-'}'"))
    return f


def main():
    args = sys.argv[1:]
    findings = check_runs()
    if "--runs" not in args:
        seen = set()
        for r in registry():
            findings += check_output(r["path"], r.get("expect"), r.get("owner_script"))
            seen.add(r["path"])
        for mf in ROOT.glob("**/" + MANIFEST):
            rel = str(mf.parent.relative_to(ROOT))
            if rel not in seen:
                findings += check_output(rel)
                seen.add(rel)

    findings.sort(key=lambda x: (RANK[x[0]], x[1]))
    n_err = sum(1 for lv, _, _ in findings if lv == ERROR)
    n_warn = sum(1 for lv, _, _ in findings if lv == WARN)
    for lv, where, msg in findings:
        if lv == INFO and "-v" not in args:
            continue
        print(f"[{lv:5}] {where}\n         {msg}")
    print(f"\n{n_err} error, {n_warn} warn, {len(findings)} checks"
          f"{'' if n_err else '  -- no stale or crashed builds'}")
    if "--strict" in args and n_err:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
