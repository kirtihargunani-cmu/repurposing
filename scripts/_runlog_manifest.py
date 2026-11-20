#!/usr/bin/env python
"""Write .build_manifest.json for one output of a runlog.sh run. Not called directly.

INTERPRETER: any python3 (stdlib only, run with -I)

For a directory: <dir>/.build_manifest.json describing every file in it.
For a file:      <name>.build_manifest.json beside it.

The manifest answers the question D162 could not: *which run wrote this, when, with
what settings, and did it replace something older?* `superseded` lists files that were
present before this run and that this run did NOT rewrite -- the stale-arm signature.
"""
import json, os, sys, hashlib
from datetime import datetime, timezone

ROOT = "/home/kirti/projects/repurposing"


def stamp(ts):
    return datetime.fromtimestamp(ts, timezone.utc).astimezone().isoformat(timespec="seconds")


def main():
    target = os.path.abspath(sys.argv[1])
    env = os.environ
    start = float(env.get("START_EPOCH", 0)) or None

    # what was there before this run
    pre = {}
    pf = env.get("PRE", "")
    if pf and os.path.exists(pf):
        for line in open(pf):
            p, sz, mt = line.rstrip("\n").split("\t")
            pre[os.path.abspath(p)] = (int(sz), float(mt))

    if os.path.isdir(target):
        names = sorted(n for n in os.listdir(target)
                       if os.path.isfile(os.path.join(target, n)) and not n.startswith(".build_manifest"))
        files = {n: dict(bytes=os.path.getsize(os.path.join(target, n)),
                         mtime=stamp(os.path.getmtime(os.path.join(target, n))))
                 for n in names}
        mf = os.path.join(target, ".build_manifest.json")
    else:
        n = os.path.basename(target)
        files = {n: dict(bytes=os.path.getsize(target), mtime=stamp(os.path.getmtime(target)))}
        mf = target + ".build_manifest.json"
        names = [n]

    # files this run did not rewrite: present before, mtime unchanged
    now_paths = {os.path.join(target, n) if os.path.isdir(target) else target: n for n in names}
    not_rewritten = []
    for p, n in now_paths.items():
        ap = os.path.abspath(p)
        if ap in pre and abs(os.path.getmtime(p) - pre[ap][1]) < 1e-6:
            not_rewritten.append(n)

    mtimes = [os.path.getmtime(os.path.join(target, n)) if os.path.isdir(target) else os.path.getmtime(target)
              for n in names]

    prev = None
    if os.path.exists(mf):
        try:
            old = json.load(open(mf))
            prev = dict(run_id=old.get("run_id"), finished_at=old.get("finished_at"),
                        n_files=old.get("n_files"), env=old.get("env"))
        except Exception:
            prev = dict(run_id="unreadable")

    man = dict(
        schema="repurposing/build_manifest/1",
        run_id=env.get("RUN_ID"),
        finished_at=stamp(datetime.now().timestamp()),
        duration_s=int(env.get("DUR", 0) or 0),
        script=env.get("SCRIPT"),
        cmd=env.get("CMD"),
        env=env.get("ENVSTR"),
        git_sha=env.get("GIT_SHA"),
        git_state=env.get("GIT_DIRTY"),
        note=env.get("NOTE") or None,
        expect=int(env["EXPECT"]) if env.get("EXPECT") else None,
        n_files=len(names),
        mtime_span_hours=round((max(mtimes) - min(mtimes)) / 3600, 2) if mtimes else 0.0,
        files_not_rewritten_by_this_run=sorted(not_rewritten),
        supersedes=prev,
        files=files,
    )
    with open(mf, "w") as f:
        json.dump(man, f, indent=1)
    warn = ""
    if not_rewritten:
        warn = f"  !! {len(not_rewritten)} pre-existing file(s) NOT rewritten by this run"
    print(f"runlog: manifest {os.path.relpath(mf, ROOT)}  n_files={len(names)}"
          f"  span={man['mtime_span_hours']}h{warn}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
