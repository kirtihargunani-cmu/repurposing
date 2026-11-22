#!/usr/bin/env python
"""Gate: do the builds on disk use the parameters we decided on? (D188, 2026-11-22)

INTERPRETER: any; stdlib + yaml.

`check_builds.py` answers "is this build STALE". It has never answered "did this build use the
thresholds we adopted", and that gap is what let the per-timepoint TAS floor be measured (D156),
analysed at length (D169), implemented (`scripts/tas_floor.py`) and then never switched on. The
value `"tas_anchor": null` was sitting in `drug_condition_profiles_run.json` the whole time.
Nothing compared it to anything.

This reads `config/adopted_parameters.yaml` and each build's own run.json, and reports:

    MISMATCH    the build recorded a value that disagrees with the adopted one
    IN_FORCE    a parameter marked design_only is actually active in a build  (the reverse error)
    UNRECORDED  a parameter is adopted but the build does not record it at all - an unverifiable
                threshold, which is the state that produced this script
    MISSING     the run.json named in the registry does not exist

Exit status is 1 if anything failed, so it can gate a chain the way check_builds.py does.

    scripts/check_params.py            report
    scripts/check_params.py --strict   exit 1 on UNRECORDED too (default: warn)
"""
import json
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    print("needs pyyaml", file=sys.stderr)
    sys.exit(2)

ROOT = Path(__file__).resolve().parent.parent
REG = ROOT / "config/adopted_parameters.yaml"


def main() -> int:
    strict = "--strict" in sys.argv
    reg = yaml.safe_load(REG.read_text())["parameters"]

    errors, warns, checks = [], [], 0
    for name, p in reg.items():
        status = p.get("status")
        records = p.get("records") or []
        if not records:
            if status == "adopted":
                warns.append((name, "UNRECORDED", "adopted but no build records it; cannot be verified"))
            continue
        for rec in records:
            f = ROOT / rec["file"]
            key = rec["key"]
            checks += 1
            if not f.exists():
                # `pending: true` means the recording code is in place but the build that writes
                # it has not been rerun yet (e.g. 79d, which is expensive and upstream of the
                # chain). That is a gap to close, not a reason to block every run - but it stays
                # visible as a warning so it cannot be forgotten, which is the whole point.
                if rec.get("pending"):
                    warns.append((name, "PENDING",
                                  f"{rec['file']} not written yet; appears on the next build of "
                                  f"{rec.get('writer', 'that stage')}"))
                else:
                    errors.append((name, "MISSING", f"{rec['file']} does not exist"))
                continue
            try:
                got = json.loads(f.read_text()).get(key, "__absent__")
            except Exception as e:                      # noqa: BLE001 - report, never crash the gate
                errors.append((name, "MISSING", f"{rec['file']} unreadable: {e}"))
                continue

            if got == "__absent__":
                (errors if status == "adopted" else warns).append(
                    (name, "UNRECORDED", f"{rec['file']} has no key '{key}'"))
                continue

            if status == "adopted":
                want = p.get("value")
                if got != want:
                    errors.append((name, "MISMATCH",
                                   f"{rec['file']}:{key} = {got!r}, adopted value is {want!r}"))
            elif status == "design_only":
                want = p.get("design_only_value")
                if got != want:
                    errors.append((name, "IN_FORCE",
                                   f"{rec['file']}:{key} = {got!r}, but {name} is design_only "
                                   f"(expected {want!r}) - a parameter is in force that was "
                                   f"never adopted"))

    for n, kind, msg in errors:
        print(f"[ERROR] {n:22s} {kind:10s} {msg}")
    for n, kind, msg in warns:
        print(f"[WARN ] {n:22s} {kind:10s} {msg}")
    print(f"\n{len(errors)} error, {len(warns)} warn, {checks} checks over {len(reg)} parameters")

    if errors:
        return 1
    if warns and strict:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
