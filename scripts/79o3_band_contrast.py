#!/usr/bin/env python
"""Trial-outcome contrast across FDA-rank bands, with a blind inter-pass agreement check (D176).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

D175 found `trial_outcomes_v2.tsv` covers fda_rank 1-10 only, so there was no low-rank control and
"do top-ranked drugs do better in trials" could not be asked. D176 fetched mid and bottom bands.

THE CONFOUND THIS SCRIPT EXISTS TO CHECK. The top band was classified in an earlier session; mid and
bottom are classified now. Any drift in how "positive" is judged between passes maps exactly onto the
contrast being measured. So 60 already-classified top-band trials were mixed, shuffled and unlabelled,
into the new reading queue. Agreement between their old and new labels is measured BEFORE the bands
are compared, and the comparison is only reported as sound if agreement holds.

Stage 1 - quote verification: every reading must carry a contiguous excerpt that occurs
    character-for-character in that trial's own abstracts (whitespace/case normalised). Failures are
    DISCARDED, never repaired.
Stage 2 - agreement: exact outcome match, binary (positive vs not) match, and Cohen's kappa over the
    60 blind re-reads.
Stage 3 - contrast: per (unit, drug) pair, positive rate conditional on >= 1 classified completed
    trial, by band, with Kruskal-Wallis and pairwise Mann-Whitney.

Reported only; nothing is gated or promoted.
"""
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, fisher_exact

ROOT = Path("/home/kirti/projects/repurposing")
D = ROOT / "results/lincs/d108"
norm = lambda s: re.sub(r"\s+", " ", str(s or "")).strip().lower()


def load_readings(paths):
    out = {}
    for f in paths:
        for line in open(f):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get("nct"):
                out[r["nct"]] = r
    return out


def main():
    q = {json.loads(l)["nct"]: json.loads(l) for l in open(D / "reading_queue.jsonl")}
    key = json.load(open(D / "blind_key.json"))
    blind_files = sorted(glob.glob(str(D / "readings/blind_chunk*.jsonl")))
    if not blind_files:
        print("no blind readings yet", file=sys.stderr)
        return 1
    old_files = [f for f in glob.glob(str(D / "readings/*.jsonl")) if "blind_chunk" not in f]
    new, old = load_readings(blind_files), load_readings(old_files)
    print(f"blind readings {len(new):,} over {len(blind_files)} files; prior readings {len(old):,}")

    # ---- stage 1: verbatim quote check
    def ok_quote(r):
        if r.get("outcome") in ("not_reported", "not_this_trial"):
            return True
        qt = norm(r.get("quote"))
        if len(qt.split()) < 5:
            return False
        hay = " ".join(norm(p.get("title", "") + " " + p.get("abstract", ""))
                       for p in q.get(r["nct"], {}).get("papers", []))
        return qt in hay
    bad = [n for n, r in new.items() if not ok_quote(r)]
    print(f"quote check: {len(new) - len(bad):,} pass, {len(bad):,} DISCARDED"
          + (f" (e.g. {bad[:3]})" if bad else ""))
    new = {n: r for n, r in new.items() if n not in bad}

    # ---- stage 2: agreement on the blind re-reads
    rr = [n for n in key["reread"] if n in new and n in old]
    print(f"\n=== inter-pass agreement, {len(rr)} blind re-reads of already-classified trials ===")
    if rr:
        a = pd.Series([old[n]["outcome"] for n in rr])
        b = pd.Series([new[n]["outcome"] for n in rr])
        exact = (a == b).mean()
        pos = lambda s: s.eq("positive")
        binacc = (pos(a) == pos(b)).mean()
        po, pe = binacc, (pos(a).mean() * pos(b).mean() + (1 - pos(a).mean()) * (1 - pos(b).mean()))
        kappa = (po - pe) / (1 - pe) if pe < 1 else float("nan")
        print(f"  exact outcome agreement : {exact:.3f} ({int((a==b).sum())}/{len(rr)})")
        print(f"  binary positive-vs-not  : {binacc:.3f}   Cohen's kappa {kappa:.3f}")
        print(f"  prior pass positive rate {pos(a).mean():.3f}  |  new pass {pos(b).mean():.3f}")
        print("\n  confusion (rows = prior pass, cols = new pass):")
        print(pd.crosstab(a, b).to_string())
        if abs(pos(a).mean() - pos(b).mean()) > 0.10 or kappa < 0.6:
            print("\n  !! AGREEMENT IS POOR - the band contrast below pools two passes that do not")
            print("     measure the same thing. Treat it as uninterpretable and re-read the top band.")

    # ---- stage 3: the contrast. prior labels for top, new labels for mid/bottom.
    t = pd.read_csv(D / "trials.tsv", sep="\t", low_memory=False)
    lab = {**old, **new}                      # new wins where both exist (it is the later read)
    t["outcome"] = t.nct.map(lambda n: lab.get(n, {}).get("outcome"))
    t["role"] = t.nct.map(lambda n: lab.get(n, {}).get("drug_role"))
    cl = t[t.outcome.notna() & ~t.outcome.isin(["not_this_trial"])
           & ~t.role.isin(["background", "comparator"])]
    g = (cl.assign(pos=cl.outcome.eq("positive"))
         .groupby(["unit", "drug", "band"]).agg(n_classified=("pos", "size"), n_pos=("pos", "sum"))
         .reset_index())
    g["pos_rate"] = g.n_pos / g.n_classified
    print(f"\n=== band contrast: positive rate given >= 1 classified completed trial ===")
    print(f"    (drug_role background/comparator excluded - the trial does not test the drug)")
    print(g.groupby("band").agg(pairs=("pos_rate", "size"), trials=("n_classified", "sum"),
                                mean_pos_rate=("pos_rate", "mean"),
                                any_positive=("n_pos", lambda v: (v > 0).mean())).round(3).to_string())
    bands = {b: g[g.band == b].pos_rate.values for b in ("top", "mid", "bottom") if (g.band == b).any()}
    if len(bands) >= 2:
        if len(bands) == 3:
            print(f"\n  Kruskal-Wallis across bands        : p = {kruskal(*bands.values()).pvalue:.3g}")
        for x, y in (("top", "bottom"), ("top", "mid"), ("mid", "bottom")):
            if x in bands and y in bands:
                p = mannwhitneyu(bands[x], bands[y], alternative="greater").pvalue
                print(f"  Mann-Whitney {x} > {y:<7}      : p = {p:.3g}")
        if "top" in bands and "bottom" in bands:
            tt, bb = g[g.band == "top"], g[g.band == "bottom"]
            tab = [[int(tt.n_pos.sum()), int(tt.n_classified.sum() - tt.n_pos.sum())],
                   [int(bb.n_pos.sum()), int(bb.n_classified.sum() - bb.n_pos.sum())]]
            odds, pf = fisher_exact(tab, alternative="greater")
            print(f"  Fisher on trial counts top vs bottom: OR = {odds:.2f}, p = {pf:.3g}  {tab}")
    g.to_csv(ROOT / "results/lincs/trial_band_contrast.tsv", sep="\t", index=False)
    print("\nwrote results/lincs/trial_band_contrast.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
