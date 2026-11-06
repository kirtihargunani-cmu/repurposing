"""Join the D109 evidence flags (79f3, results/lincs/evidence_flags_v2.tsv) onto a cohort x drug table.

WEAK = flags that count as weak evidence (hidden on the pages, excluded from 79p / 79o): dose_inconsistent,
inactive_matched_line, low_exposure. generic_mechanism is shown only.
"""
from pathlib import Path

import pandas as pd

RES = Path(__file__).resolve().parent.parent / "results/lincs"
FLAGS = ["dose_inconsistent", "plate_confined", "inactive_matched_line", "low_exposure", "generic_mechanism"]
# D109 result: rule (1) dose_inconsistent not adopted (flags 45% of judged profiles, incl. doxorubicin / temsirolimus;
# does not separate plate artifacts from high-dose cytotoxicity); computed and kept in the table, not used.
WEAK = ["plate_confined", "inactive_matched_line", "low_exposure"]   # plate_confined: D109a (replaces the dose rule)


def apply_evidence(df, unit="group_id"):
    f = RES / "evidence_flags_v2.tsv"
    df = df.drop(columns=[c for c in FLAGS + ["exposure_ratio", "weak_d109"] if c in df])
    if not f.exists():
        for c in FLAGS:
            df[c] = False
        df["weak_d109"] = False
        return df
    e = pd.read_csv(f, sep="\t", low_memory=False)
    if "exposure_ratio" not in e:
        e["exposure_ratio"] = float("nan")
    for c in FLAGS:
        if c not in e:
            e[c] = False
    e = e.drop_duplicates(["group_id", "pert_id"]).set_index(["group_id", "pert_id"])[FLAGS + ["exposure_ratio"]]
    df = df.join(e, on=[unit, "pert_id"])
    for c in FLAGS:
        df[c] = df[c].astype(str).str.lower().eq("true")
    df["dose_inconsistent"] = False   # D109 rule (1) not adopted
    df["weak_d109"] = df[WEAK].any(axis=1)
    return df
