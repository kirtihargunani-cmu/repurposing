#!/usr/bin/env python
"""F-MATH 3: export a frozen, signed prior-knowledge network (PKN) for CARNIVAL.

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

Spec: "Export separate, frozen signed_pkn.tsv representations from the repo's curated signaling
sources. Preserve source and target gene symbols/IDs, direction, sign (including unknown/conflicted),
interaction provenance, degree, and tissue restrictions. ARACNe coexpression links must not silently
become directed biochemical causal interactions."

PRIMARY SOURCE = data/raw/omnipath/omnipath_signed_directed.tsv.gz, because that is what the rest of
the project already uses as its signalling prior (scripts/d126/b1_config.py PRIOR; scripts/d128/
config.py OMNIPATH, "primary signaling prior (P6)"). Building from omnipath_raw.tsv.gz instead would
give a DIFFERENT, larger and far more activation-biased network (70,389 edges, 90.5% activating,
against 15,860 edges and 70.3% activating here) and would silently diverge from D126/D128.

SIGNOR is exported as a SEPARATE arm (the spec's "separate, frozen representations"), matching
d128/config.py where SIGNOR is the P6 sensitivity prior. It is never merged into the OmniPath file.

Sign handling - nothing is guessed:
  activating / inhibiting   sign +1 / -1 from the curated consensus  -> the CARNIVAL file
  sign_disputed             carried per edge as a column, NOT dropped (1,184 of 15,860)
  unknown / conflicted      taken from omnipath_raw.tsv.gz (directed, neither or both stimulation
                            and inhibition) and written to pkn_excluded_edges.tsv, never signed
Complex nodes (OmniPath writes them as underscore-joined symbols, e.g. ABI1_EPS8_SOS1) are NOT genes
and cannot be matched to TF activity measurements. They are excluded from the CARNIVAL file and
recorded with reason `complex_node`, rather than being split or silently kept.

ARACNe is not read by this script at all.

Outputs (results/systems_validation/pkn/):
  signed_pkn.tsv / signor_signed_pkn.tsv   CARNIVAL-ready: source, interaction, target + provenance
  pkn_excluded_edges.tsv                   every excluded edge with its reason
  pkn_manifest.json                        counts, sha256 of inputs and outputs, git SHA
"""
import hashlib, json, os, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
RAW = ROOT / "data/raw/omnipath"
OUT = ROOT / "results/systems_validation/pkn"
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[v] = "4"
KEEP = ["source", "interaction", "target", "source_uniprot", "target_uniprot", "sign_disputed",
        "n_resources", "n_references", "curation_effort", "sources"]


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def build(src, label):
    """One curated signed-directed export -> (carnival_ready, excluded)."""
    d = pd.read_csv(src, sep="\t", low_memory=False)
    d = d.rename(columns={"source_genesymbol": "source", "target_genesymbol": "target",
                          "source": "source_uniprot", "target": "target_uniprot",
                          "sign": "interaction"})
    for c in ("n_resources", "n_references", "curation_effort", "sources", "sign_disputed"):
        if c not in d:
            d[c] = pd.NA
    d["sign_disputed"] = d.sign_disputed.fillna(False).astype(bool)
    d = d[d.source.notna() & d.target.notna()]
    exc = []
    self_loop = d.source == d.target
    exc.append(d[self_loop].assign(excluded_reason="self_loop"))
    d = d[~self_loop]
    cx = d.source.astype(str).str.contains("_") | d.target.astype(str).str.contains("_")
    exc.append(d[cx].assign(excluded_reason="complex_node"))
    d = d[~cx]
    d = (d.sort_values("curation_effort", ascending=False)
           .drop_duplicates(["source", "target", "interaction"]))
    opp = d.groupby(["source", "target"]).interaction.nunique() > 1
    bad = set(opp[opp].index)
    is_bad = [(a, b) in bad for a, b in zip(d.source, d.target)]
    exc.append(d[is_bad].assign(excluded_reason="opposite_signs_for_same_pair"))
    d = d[[not x for x in is_bad]]
    nodes = pd.Index(sorted(set(d.source) | set(d.target)))
    d["source_out_degree"] = d.source.map(d.groupby("source").size().reindex(nodes).fillna(0).astype(int))
    d["target_in_degree"] = d.target.map(d.groupby("target").size().reindex(nodes).fillna(0).astype(int))
    d["tissue_restriction"] = "none_recorded"
    d = d[KEEP + ["source_out_degree", "target_in_degree", "tissue_restriction"]]
    d = d.sort_values(["source", "target"]).reset_index(drop=True)
    e = pd.concat(exc, ignore_index=True) if exc else pd.DataFrame()
    if len(e):
        e = e[[c for c in KEEP if c in e] + ["excluded_reason"]].assign(prior=label)
    return d, e, nodes


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pkn, exc, nodes = build(RAW / "omnipath_signed_directed.tsv.gz", "omnipath")

    # unknown / conflicted DIRECTED edges, from the raw export: preserved, never signed
    raw = pd.read_csv(RAW / "omnipath_raw.tsv.gz", sep="\t", low_memory=False)
    raw = raw[raw.is_directed == True].copy()
    stim, inh = raw.is_stimulation == True, raw.is_inhibition == True
    raw["excluded_reason"] = None
    raw.loc[stim & inh, "excluded_reason"] = "sign_conflicted"
    raw.loc[~stim & ~inh, "excluded_reason"] = "sign_unknown"
    unk = raw[raw.excluded_reason.notna()].rename(
        columns={"source_genesymbol": "source", "target_genesymbol": "target",
                 "source": "source_uniprot", "target": "target_uniprot"})
    unk["interaction"] = 0
    unk["sign_disputed"] = False
    for c in ("n_resources", "n_references"):
        if c not in unk:
            unk[c] = pd.NA
    unk = unk[[c for c in KEEP if c in unk] + ["excluded_reason"]].assign(prior="omnipath")
    unk = unk[unk.source.notna() & unk.target.notna() & (unk.source != unk.target)]
    excluded = pd.concat([exc, unk], ignore_index=True).drop_duplicates(
        ["source", "target", "interaction", "excluded_reason"])

    f_pkn, f_exc = OUT / "signed_pkn.tsv", OUT / "pkn_excluded_edges.tsv"
    pkn.to_csv(f_pkn, sep="\t", index=False)
    excluded.to_csv(f_exc, sep="\t", index=False)
    outs = {f_pkn: pkn, f_exc: excluded}

    sig_src = RAW / "signor_signed_directed.tsv.gz"
    sig_n = None
    if sig_src.exists():
        sig, _, sig_nodes = build(sig_src, "signor")
        f_sig = OUT / "signor_signed_pkn.tsv"
        sig.to_csv(f_sig, sep="\t", index=False)
        outs[f_sig] = sig
        sig_n = dict(edges=int(len(sig)), nodes=int(len(sig_nodes)))

    hub = pkn.source.value_counts()
    man = dict(
        artifact="F-MATH 3 frozen signed prior-knowledge network",
        spec="docs/SPEC_systems_theoretic_mechanism_validation.md",
        primary_source="data/raw/omnipath/omnipath_signed_directed.tsv.gz",
        why_this_source="matches scripts/d126/b1_config.py PRIOR and scripts/d128/config.py OMNIPATH",
        git_sha=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                               text=True).stdout.strip(),
        aracne_used=False,
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in
                [RAW / "omnipath_raw.tsv.gz", RAW / "omnipath_signed_directed.tsv.gz"]
                + ([sig_src] if sig_src.exists() else [])},
        outputs={str(k.relative_to(ROOT)): sha256(k) for k in outs},
        omnipath=dict(signed_edges=int(len(pkn)), nodes=int(len(nodes)),
                      activating=int((pkn.interaction == 1).sum()),
                      inhibiting=int((pkn.interaction == -1).sum()),
                      sign_disputed=int(pkn.sign_disputed.sum()),
                      max_out_degree=int(hub.max()), top_hubs=hub.head(10).to_dict()),
        signor=sig_n,
        excluded=dict(excluded.excluded_reason.value_counts().items()),
        excluded_total=int(len(excluded)),
        tissue_restriction="none_recorded (these interactions are context-free)")
    (OUT / "pkn_manifest.json").write_text(json.dumps(man, indent=1, default=str))

    o = man["omnipath"]
    print(f"signed_pkn.tsv : {o['signed_edges']:,} edges, {o['nodes']:,} nodes "
          f"({o['activating']:,} activating, {o['inhibiting']:,} inhibiting, "
          f"{o['sign_disputed']:,} sign-disputed kept as a flag)")
    if sig_n:
        print(f"signor arm     : {sig_n['edges']:,} edges, {sig_n['nodes']:,} nodes (separate file)")
    print(f"excluded       : {man['excluded_total']:,} -> {man['excluded']}")
    print(f"max out-degree : {o['max_out_degree']} ({list(o['top_hubs'])[:5]})")
    print("wrote", OUT)


if __name__ == "__main__":
    sys.exit(main())
