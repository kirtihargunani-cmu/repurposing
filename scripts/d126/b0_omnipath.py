#!/usr/bin/env python
"""D126 B0: pull OmniPath, filter to signed and directed interactions, cache, checksum, register.

Rule: docs/SPEC_graph_metric_scoring.md B2, with the choices and the two corrections of D126-B0.

The prior is the *cached file*, not a version string: OmniPath is a service over moving upstream resources and has
no data release label (omnipathdb.org/about returns the server software version). So this script writes the raw
response and the filtered edge table under data/raw/omnipath/, checksums both, and registers the query URL, the
retrieval date, the row counts and the signed-and-directed subset size in config/datasets.yaml with a `review:` key
saying a re-pull is a new decision rather than a refresh.

Filter (B2), every threshold named here and nowhere else:
  consensus_direction, exactly one of consensus_stimulation / consensus_inhibition (B2's "consensus_sign"; the
  service has no such field -- D126-B0 (2)), and >= MIN_RESOURCES distinct resources or >= MIN_REFERENCES distinct
  PubMed ids. Edges where some source says stimulation and another inhibition are kept and marked sign_disputed.

Arms: the primary prior (datasets=omnipath) and the SIGNOR sensitivity arm (the same query, resources=SIGNOR).

Outputs: data/raw/omnipath/{<arm>_raw.tsv.gz, <arm>_signed_directed.tsv.gz}, results/lincs/d126/b0_summary.tsv,
         config/datasets.yaml entry `omnipath`.
"""
import gzip
import hashlib
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# --- config block: the pre-registered filter (B2) and the query (D126-B0) ---------------------------------------
MIN_RESOURCES = 2
MIN_REFERENCES = 2
BASE = "https://omnipathdb.org/interactions"
QUERY = dict(datasets="omnipath", organisms="9606", genesymbols="yes",
             fields="sources,references,curation_effort", format="tsv")
ARMS = {"omnipath": {}, "signor": {"resources": "SIGNOR"}}      # the SIGNOR sensitivity arm is the same query
USER_AGENT = "OmniPath-client"      # generic: no user or project identity is sent
ABOUT = "https://omnipathdb.org/about?format=json"

RAW = ROOT / "data/raw/omnipath"
OUT = ROOT / "results/lincs/d126"
REG = ROOT / "config/datasets.yaml"


def url_for(arm):
    q = dict(QUERY, **ARMS[arm])
    return BASE + "?" + "&".join(f"{k}={v}" for k, v in sorted(q.items()))


def fetch(url, dest):
    r = subprocess.run(["curl", "-sS", "--fail", "--compressed", "-A", USER_AGENT, url, "-o", str(dest)],
                       capture_output=True, text=True, timeout=1800)
    if r.returncode:
        sys.exit(f"fetch failed: {url}\n{r.stderr}")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def n_tokens(s, pmid=False):
    if not isinstance(s, str) or not s:
        return 0
    return len({p.split(":")[-1] if pmid else p for p in s.split(";") if p})


def filter_edges(x):
    """B2's signed-and-directed subset, with D126-B0 (2)'s reading of 'consensus_sign'."""
    for c in ["is_directed", "is_stimulation", "is_inhibition", "consensus_direction",
              "consensus_stimulation", "consensus_inhibition"]:
        x[c] = x[c].astype(bool)
    x["n_resources"] = [n_tokens(s) for s in x.sources]
    x["n_references"] = [n_tokens(s, pmid=True) for s in x.references]
    x["consensus_signed"] = x.consensus_stimulation ^ x.consensus_inhibition      # exactly one
    x["sign"] = pd.NA
    x.loc[x.consensus_signed & x.consensus_stimulation, "sign"] = 1
    x.loc[x.consensus_signed & x.consensus_inhibition, "sign"] = -1
    x["sign_disputed"] = x.is_stimulation & x.is_inhibition
    keep = (x.consensus_direction & x.consensus_signed
            & ((x.n_resources >= MIN_RESOURCES) | (x.n_references >= MIN_REFERENCES)))
    return x, keep


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    today = date.today().isoformat()
    about = subprocess.run(["curl", "-sS", "-A", USER_AGENT, ABOUT], capture_output=True, text=True,
                           timeout=120).stdout.strip()
    print("service:", about[:160])

    rows, reg = [], {}
    for arm in ARMS:
        url = url_for(arm)
        raw = RAW / f"{arm}_raw.tsv.gz"
        tmp = RAW / f"{arm}_raw.tsv"
        print(f"\n[{arm}] {url}")
        fetch(url, tmp)
        x = pd.read_csv(tmp, sep="\t")
        with gzip.open(raw, "wt") as f:                 # cache the response as retrieved
            x.to_csv(f, sep="\t", index=False)
        tmp.unlink()
        x, keep = filter_edges(x)
        e = x[keep].copy()
        cols = ["source", "target", "source_genesymbol", "target_genesymbol", "sign", "sign_disputed",
                "n_resources", "n_references", "curation_effort", "sources", "references"]
        out = RAW / f"{arm}_signed_directed.tsv.gz"
        e[cols].to_csv(out, sep="\t", index=False)
        nodes = len(set(e.source_genesymbol) | set(e.target_genesymbol))
        rows.append(dict(arm=arm, retrieved=today, url=url, raw_rows=len(x), kept_edges=int(keep.sum()),
                         nodes=nodes, activating=int((e.sign == 1).sum()), inhibiting=int((e.sign == -1).sum()),
                         sign_disputed=int(e.sign_disputed.sum()),
                         dropped_no_consensus_direction=int((~x.consensus_direction).sum()),
                         dropped_unsigned=int((~x.consensus_signed).sum()),
                         dropped_thin_support=int((x.consensus_direction & x.consensus_signed & ~keep).sum()),
                         raw_sha256=sha256(raw), edges_sha256=sha256(out)))
        reg[arm] = rows[-1]
        print(f"  raw {len(x)} -> signed and directed {int(keep.sum())} edges over {nodes} nodes "
              f"({rows[-1]['activating']} activating, {rows[-1]['inhibiting']} inhibiting, "
              f"{rows[-1]['sign_disputed']} sign_disputed)")

    s = pd.DataFrame(rows)
    s.to_csv(OUT / "b0_summary.tsv", sep="\t", index=False)
    (OUT / "b0_run.json").write_text(json.dumps(dict(service=about, retrieved=today, arms=reg,
                                                     min_resources=MIN_RESOURCES, min_references=MIN_REFERENCES),
                                                indent=1))

    # --- config/datasets.yaml entry (D126-B0 (1)): the cached file is the dependency -----------------------------
    entry = [
        "",
        "omnipath:  # D126 signed signaling prior (scripts/d126/b0_omnipath.py). The CACHED FILES are the dependency.",
        f"  release: no-upstream-release  # {about.split('.')[0].replace('Hello, this is the ', '')} is the server build, not a data version",
        f"  retrieved: {today}",
        "  access: open",
        "  review: |",
        "    OmniPath is a service over moving upstream resources and publishes no data release label, so a checksum",
        "    pins reproducibility but gives no way to notice drift. The analysis depends on the cached tables below,",
        "    not on the service. Re-pulling yields a different prior and is a new decision, not a refresh: any D126",
        "    number names the retrieval date above. Subset sizes are the signed-and-directed counts, not OmniPath's",
        "    headline interaction count.",
        f"  filter: consensus_direction and exactly one of consensus_stimulation / consensus_inhibition, and"
        f" (>= {MIN_RESOURCES} distinct resources or >= {MIN_REFERENCES} distinct PubMed ids)",
        "  files:",
    ]
    for r in rows:
        entry += [
            f"    {r['arm']}_raw:             {{url: \"{r['url']}\", rows: {r['raw_rows']}, sha256: {r['raw_sha256']}}}",
            f"    {r['arm']}_signed_directed: {{path: data/raw/omnipath/{r['arm']}_signed_directed.tsv.gz,"
            f" edges: {r['kept_edges']}, nodes: {r['nodes']}, sha256: {r['edges_sha256']}}}",
        ]
    txt = REG.read_text()
    if "\nomnipath:" in txt:
        sys.exit("config/datasets.yaml already has an `omnipath` entry; refusing to overwrite it (D126-B0 (1))")
    REG.write_text(txt.rstrip("\n") + "\n" + "\n".join(entry) + "\n")

    pd.set_option("display.width", 200)
    print("\n" + s.drop(columns=["url", "raw_sha256", "edges_sha256"]).to_string(index=False))
    print(f"\nregistered in {REG.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
