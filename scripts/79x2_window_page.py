#!/usr/bin/env python
"""Browsable page of the window-dimensioned drug ranking (D172; runs in `celligner`).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

One self-contained HTML file (data embedded), built from `drug_ranking_by_window<SFX>.tsv` (79k).
Every drug is shown with the TIMEPOINT WINDOW its score came from, because `79c` otherwise averages
6 h and 24 h into one number and the window disappears: clomifene in BRCA-LumA reads 98.0 early and
12.9 late, and pools to 46.7 (D170/D171).

Per cohort, one row per drug: the lead score (`layered_pct`, D122), the mechanism-specific
percentile, and the per-window reversal percentiles with the gap between them. Rows carry the
window badge early / late / both / unknown, and the gate flags that must be read WITH the scores
(D168): `target_absent` means every curated mechanism target is unexpressed, so a high reversal
percentile beside it is a suppressed false positive, not a lead.

Honest limits shown on the page itself, not only here:
  * only 29 of 41 cohorts have ANY window data; 12 have zero condition-level profiles in their
    matched cell lines, so `unknown` there is permanent, not pending
  * `n_conditions` is usually 1, so a large early-late gap is a flag for inspection, not an
    effect size
  * early and late percentiles are computed over DIFFERENT populations (the drugs having that
    window) and compare as rankings, not magnitudes

Where several LINCS compound ids share a `cmap_name`, the one with the most conditions is kept
(most evidence), and `n_ids` records how many were collapsed.

Inputs:  results/lincs/drug_ranking_by_window<SFX>.tsv (79k)
         scripts/templates/window_rankings.html (placeholder __DATA__)
Output:  results/lincs/drug_rankings_by_window<SFX>.html
         published copy (private): https://claude.ai/artifact/Loe6nCcucsB9m7TksExQ8j
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
RES = ROOT / "results/lincs"
SFX = os.environ.get("SPLIT_SFX", "_v2_tf")
TPL = ROOT / "scripts/templates/window_rankings.html"


def num(v, nd=1):
    """JSON-safe rounded float, None for missing - the template renders None as a dash."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if not np.isfinite(f) else round(f, nd)


def main():
    src = RES / f"drug_ranking_by_window{SFX}.tsv"
    if not src.exists():
        print(f"missing {src} - run 79k_window_long.py first", file=sys.stderr)
        return 1
    d = pd.read_csv(src, sep="\t", low_memory=False)
    tb = lambda v: v.astype(str).isin(["True", "TRUE", "1"])
    d = d[tb(d.passes_gate)].copy()

    # one row per (cohort, drug name, window): keep the id with the most conditions
    d["_n"] = pd.to_numeric(d.n_conditions, errors="coerce").fillna(0)
    d = d.sort_values("_n", ascending=False)
    ids = d.groupby(["group_id", "cmap_name"]).pert_id.nunique().rename("n_ids")
    d = d.drop_duplicates(["group_id", "cmap_name", "window"])

    # wide: one record per (cohort, drug) carrying both windows
    piv = d.pivot_table(index=["group_id", "cmap_name"], columns="window",
                        values=["window_score_pct", "n_conditions"], aggfunc="first")
    flat = d.drop_duplicates(["group_id", "cmap_name"]).set_index(["group_id", "cmap_name"])
    rows = {}
    for (g, nm), r in flat.iterrows():
        def pv(k, w):
            try:
                return piv.loc[(g, nm), (k, w)]
            except KeyError:
                return np.nan
        e, l_ = pv("window_score_pct", "early"), pv("window_score_pct", "late")
        p = pv("window_score_pct", "unknown")
        ce, cl = pv("n_conditions", "early"), pv("n_conditions", "late")
        win = ("both" if pd.notna(e) and pd.notna(l_) else
               "early" if pd.notna(e) else "late" if pd.notna(l_) else "unknown")
        gap = num(e - l_) if pd.notna(e) and pd.notna(l_) else None
        rows.setdefault(g, []).append([
            nm,                                   # 0 drug
            num(r.get("layered_pct")),            # 1 lead score
            num(r.get("spec_mech_pct")),          # 2 mechanism-specific pct
            num(e), num(l_), num(p),              # 3,4,5 early / late / unknown-timepoint (pooled score) pct
            gap,                                  # 6 early - late
            win,                                  # 7 badge
            num(ce, 0), num(cl, 0),               # 8,9 conditions per window
            int(ids.get((g, nm), 1)),             # 10 ids collapsed
            bool(tb(pd.Series([r.get("target_absent")])).iloc[0]),   # 11
            bool(tb(pd.Series([r.get("prolif_driven")])).iloc[0]),   # 12
            (str(r.get("evidence_targets"))[:90] if pd.notna(r.get("evidence_targets")) else ""),  # 13
        ])
    for g in rows:
        rows[g].sort(key=lambda x: (x[1] is None, -(x[1] or 0)))

    cov = {g: sum(1 for r in v if r[7] != "unknown") for g, v in rows.items()}

    # D192: per-cohort approval RECOVERY, per window, matched against the pooled score on the SAME
    # rows. Optional - the page still builds without it. It is NOT accuracy: label_phase == 0 is
    # not a negative, so a low value is consistent with discovery. See 79ag's docstring.
    acc = {}
    af = RES / "window_accuracy.tsv"
    if af.exists():
        A = pd.read_csv(af, sep="\t")
        for r in A.itertuples():
            acc.setdefault(r.cohort, {})[r.window] = [
                round(float(r.auc_window), 3), round(float(r.se), 3),
                round(float(r.auc_pooled_same_rows), 3), round(float(r.delta), 3),
                int(r.n_pos), int(r.n)]
    meta = {
        "accuracy": acc,
        "cohorts": sorted(rows, key=lambda g: (-cov[g] / max(1, len(rows[g])), g)),
        "n_rows": int(sum(len(v) for v in rows.values())),
        "n_cohorts": len(rows),
        "n_with_window": int(sum(cov.values())),
        "n_both": int(sum(1 for v in rows.values() for r in v if r[7] == "both")),
        "coverage": {g: [cov[g], len(rows[g])] for g in rows},
    }
    data = json.dumps({"meta": meta, "rows": rows}, separators=(",", ":"))
    page = TPL.read_text().replace("__DATA__", data, 1)
    out = RES / f"drug_rankings_by_window{SFX}.html"
    out.write_text(page)
    print(f"79x2: {meta['n_rows']:,} (cohort, drug) rows over {meta['n_cohorts']} cohorts; "
          f"{meta['n_with_window']:,} carry a real window, {meta['n_both']:,} have both")
    print(f"      {len(page) / 1e6:.2f} MB -> results/lincs/{out.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
