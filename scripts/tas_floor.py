"""Per-timepoint activity floor (D156). Shared by every script that filters on TAS.

SPEC: docs/SPEC_per_timepoint_tas_floor.md

TAS (transcriptional activity score) measures response MAGNITUDE, and early responses are genuinely
smaller - median TAS 0.264 at 6 h against 0.289 at 24 h over is_hiq trt_cp signatures. Reproducibility
is handled separately by is_hiq. A single absolute floor therefore demands that a 6 h response be as
large as a 24 h one to count as active, and discards 76% of 6 h signatures (median TAS 0.219,
clustered just under 0.4).

The rule: the floor sits at the same QUANTILE of each timepoint's own TAS distribution. "Active" then
means the same thing relative to what is achievable at that timepoint. Anchoring on the status quo -
0.4 is the 0.645 quantile at 24 h - raises the 48 h floor as well as lowering the early ones.

Usage:
    from tas_floor import active_mask
    act = sig[active_mask(sig)]          # honours TAS_ANCHOR; default = uniform 0.4

Environment:
    TAS_ANCHOR   unset -> uniform TAS_MIN (production, unchanged)
                 a float in (0,1) -> per-timepoint floor at that quantile
    TAS_MIN      the uniform floor, default 0.4
"""
import os
import numpy as np
import pandas as pd

TAS_MIN = float(os.environ.get("TAS_MIN", "0.4"))
TAS_ANCHOR = os.environ.get("TAS_ANCHOR", "")
MIN_SIGS_PER_TIME = 200          # below this a timepoint keeps the uniform floor


def floors(sig, anchor=None):
    """Per-timepoint floor from the signatures' own TAS distribution. Returns {pert_time: floor}."""
    a = float(anchor if anchor is not None else TAS_ANCHOR)
    out = {}
    for t, g in sig.groupby("pert_time"):
        out[t] = float(g.tas.quantile(a)) if len(g) >= MIN_SIGS_PER_TIME else TAS_MIN
    return out


def active_mask(sig, anchor=None, verbose=True):
    """Boolean mask of ACTIVE signatures. is_hiq is NOT applied here - callers already do that."""
    a = anchor if anchor is not None else (TAS_ANCHOR or None)
    if not a:
        return sig.tas >= TAS_MIN
    f = floors(sig, a)
    thr = sig.pert_time.map(f).fillna(TAS_MIN)
    m = sig.tas >= thr
    if verbose:
        base = sig.tas >= TAS_MIN
        print(f"TAS_ANCHOR={a}: per-timepoint floors "
              + ", ".join(f"{t:g}h {v:.3f}" for t, v in sorted(f.items()) if len(str(t)) < 8)
              + f" | active {int(m.sum()):,} vs {int(base.sum()):,} at a uniform {TAS_MIN}", flush=True)
    return m
