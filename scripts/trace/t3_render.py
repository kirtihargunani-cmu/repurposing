#!/usr/bin/env python
"""Mechanism trace, stage 3: render the flat record and the prose.

Section 0 of docs/PROCEDURE_mechanism_trace.md: the trace is a DATA STRUCTURE and the narrative is a
rendering of it. This script writes nothing that is not already a field in links_L1_L4.json,
links_L5_L7.json or verdicts.json. It reads no data.
"""
import json
import textwrap
from pathlib import Path

import pandas as pd

OUT = Path("/home/kirti/projects/repurposing/results/lincs/trace")
FIELDS = ["expected", "observed", "verdict", "provenance", "independence", "concentration", "alternatives"]


def main():
    links = json.loads((OUT / "links_L1_L4.json").read_text()) + \
        json.loads((OUT / "links_L5_L7.json").read_text())
    verdicts = {v["pair"]: v for v in json.loads((OUT / "verdicts.json").read_text())}
    rows = []
    for l in sorted(links, key=lambda x: (x["pair"], int(x["link"][1:]))):
        r = dict(pair=l["pair"], link=l["link"], question=l["question"])
        for f in FIELDS:
            v = l[f]
            r[f] = " | ".join(v) if isinstance(v, list) else v
        r["subpredictions"] = " | ".join(
            f"{s['id']}: {'PASS' if s['passed'] else 'FAIL'} ({s['value']})" for s in l.get("subpredictions", []))
        rows.append(r)
    pd.DataFrame(rows).to_csv(OUT / "trace_records.tsv", sep="\t", index=False)

    w = textwrap.TextWrapper(width=98, subsequent_indent="  ")
    out = ["MECHANISM TRACES - first application of docs/PROCEDURE_mechanism_trace.md Section 10.",
           "Rendered from the records; every sentence below is a field. No new assertion is added here.", ""]
    for pair, v in verdicts.items():
        out += [f"=== {pair} ===",
                f"VERDICT: {v['verdict'].upper()} - break point {v['break_point']} "
                f"(links at the worst level: {', '.join(v['all_links_at_worst'])}).",
                f"Chain: " + "  ".join(f"{c['link']}:{c['verdict']}" for c in v["chain"]),
                f"Effective independent links (Section 3): "
                f"{v['independence']['effective_independent_links']} of {v['independence']['total_links']}; "
                f"links adding no new source: {', '.join(v['independence']['links_adding_no_new_source'])}.",
                ""]
        for l in [x for x in links if x["pair"] == pair]:
            out.append(f"-- {l['link']} [{l['verdict']}] {l['question']}")
            for f in FIELDS[:-1]:
                out += w.wrap(f"   {f}: {l[f]}")
            for a in l["alternatives"]:
                out += w.wrap(f"   alternative: {a}")
            for s in l.get("subpredictions", []):
                out += w.wrap(f"   sub-prediction {s['id']} {'PASS' if s['passed'] else 'FAIL'}: "
                              f"{s['statement']} -> {s['value']}")
            out.append("")
        out += ["ALTERNATIVES CARRIED INTO THE VERDICT (Section 6):"]
        for a in v["alternatives_carried"]:
            out += w.wrap(f"   [{a['link']}] {a['alternative']}")
        out += ["", "NOT DONE (Section 7): the five negative controls have not been run for this application, so",
                "the procedure's false-positive rate is unmeasured and neither trace counts as a positive.", ""]
    (OUT / "trace_narrative.txt").write_text("\n".join(out))
    print("wrote", OUT / "trace_records.tsv", "and", OUT / "trace_narrative.txt")


if __name__ == "__main__":
    main()
