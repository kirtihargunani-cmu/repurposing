#!/usr/bin/env python
"""One row per (disease, drug, TIMEPOINT WINDOW) — no unlabelled averaged number (D171).

INTERPRETER: /home/kirti/miniforge3/envs/celligner/bin/python

User, 2026-11-20: *"don't even return a pooled score, if a drug works against a target for a disease
at an early timepoint, return the drug with a column for time point saying early."*

WHAT THIS CHANGES. `79c` collapses every drug x line to one profile by averaging over doses AND
times, so a drug gets one score per disease and the window it came from is invisible. Where the two
windows disagree that average means neither: clomifene in BRCA-LumA reads **97.9 early vs 12.9
late** and pools to 46.7 (D170). This script makes the window an explicit dimension of the output.

    window_coverage == both        -> TWO rows, window = early and window = late
    window_coverage == early_only  -> ONE row,  window = early
    window_coverage == late_only   -> ONE row,  window = late
    window_coverage == unknown     -> ONE row,  window = unknown (D180; was "pooled"), with window_unknown_reason

**`unknown` (formerly `pooled`) is a label, not a window.** It means "averaged over whatever timepoints existed, window
unknown" and is the honest answer where condition-level data does not exist. It is kept rather than
dropped because dropping it would delete most of the atlas: only **29 of 41 cohorts** have any
window data, and **12 have zero condition-level profiles in their matched cell lines** - no window
can ever be assigned there from this dataset. Nine cohorts (BRCA x4, PRAD, SKCM, LUSC, LUAD,
COAD/READ) carry nearly all of it.

READ `window_score_pct` WITHIN ITS WINDOW. early and late percentiles are computed over different
populations (the drugs having that window), so they compare as rankings, not as magnitudes. And the
window score comes from the CONDITION-level pipeline while `rev`/`rev_pct` come from the CELL-level
one: same formula (79e:285), differently aggregated inputs, so they are not expected to match even
for a single-window drug.

`n_conditions` is usually 1 (D155 run.json median), so a large early-late gap is a flag for
inspection, not a measured effect size.

Env: SPLIT_SFX (default _v2_tf). Output: results/lincs/drug_ranking_by_window<SFX>.tsv
"""
import os
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/home/kirti/projects/repurposing")
RES = ROOT / "results/lincs"
WIN = ROOT / "results/systems_validation/integrated/window_reversal.tsv.gz"
SFX = os.environ.get("SPLIT_SFX", "_v2_tf")
KEEP = ["group_id", "pert_id", "cmap_name", "passes_gate", "gate_reason", "target_absent",
        "evidence_targets", "layered_pct", "spec_mech_pct", "rev", "rev_pct", "prolif_driven"]


def main():
    sp = RES / f"drug_reversal_split{SFX}.tsv"
    if not sp.exists():
        print(f"missing {sp}", file=sys.stderr)
        return 1
    d = pd.read_csv(sp, sep="\t", low_memory=False)
    w = pd.read_csv(WIN, sep="\t", low_memory=False).rename(columns={"cohort": "group_id"})
    w = w.drop_duplicates(["group_id", "pert_id"])
    nc = w.set_index(["group_id", "pert_id"])[["n_conditions_early", "n_conditions_late"]]

    base = d[[c for c in KEEP if c in d.columns]].copy()
    out = []
    for side in ("early", "late"):
        m = d.window_coverage.isin(["both", f"{side}_only"]) & d[f"rev_{side}"].notna()
        if not m.any():
            continue
        r = base[m].copy()
        r["window"] = side
        r["window_score"] = d.loc[m, f"rev_{side}"].values
        r["window_score_pct"] = d.loc[m, f"rev_{side}_pct"].values
        r["window_agreement"] = d.loc[m, "window_agreement"].values
        r["window_gap_pct"] = d.loc[m, "window_gap_pct"].values
        key = pd.MultiIndex.from_arrays([r.group_id, r.pert_id])
        r["n_conditions"] = nc[f"n_conditions_{side}"].reindex(key).values
        out.append(r)
    m = d.window_coverage.isna() | d.window_coverage.eq("unknown")
    r = base[m].copy()
    r["window"] = "unknown"
    r["window_unknown_reason"] = d.loc[m, "window_unknown_reason"].values if "window_unknown_reason" in d else pd.NA
    r["window_score"] = d.loc[m, "rev"].values
    r["window_score_pct"] = d.loc[m, "rev_pct"].values if "rev_pct" in d else pd.NA
    for c in ("window_agreement", "window_gap_pct", "n_conditions"):
        r[c] = pd.NA
    out.append(r)

    o = pd.concat(out, ignore_index=True)
    o = o.sort_values(["group_id", "cmap_name", "window"]).reset_index(drop=True)
    f = RES / f"drug_ranking_by_window{SFX}.tsv"
    o.to_csv(f, sep="\t", index=False)

    print(f"79k: {len(d):,} split rows -> {len(o):,} (disease, drug, window) rows")
    print(o.window.value_counts().to_string())
    tb = lambda v: v.astype(str).isin(["True", "TRUE", "1"])
    g = o[tb(o.passes_gate)]
    print(f"\ngate-passing rows {len(g):,}; drugs with BOTH windows: "
          f"{g[g.window.isin(['early','late'])].groupby(['group_id','pert_id']).size().eq(2).sum():,}")
    print(f"cohorts with any real window: {o[o.window.ne('unknown')].group_id.nunique()} of {o.group_id.nunique()}")
    print(f"wrote {f.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
