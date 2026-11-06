"""Shared selection of FDA-approved drug rows per disease and subgroup (used by 79o MODE=fda_all and 79x).

Rows: drugs that pass the target gate (79b) and are matched to a Drugs@FDA ingredient (79p, fda_matches.tsv; marketed
or discontinued). Several LINCS ids of one ingredient -> the best on the cohort's lead score (diseases: spec_mech, then
spec, rev; subgroups: spec_sib_mech, then spec_mech, rev), as 79p does. Unlike 79p, flagged (incl. D107 local_route) and non-reversing drugs are
kept (the page shows them, weak evidence hidden by default). `drug` = FDA ingredient with salt words removed (79p rule).
"""
from pathlib import Path

import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from route_flag import apply_route  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SALTS = {"HYDROCHLORIDE", "DIHYDROCHLORIDE", "HCL", "SODIUM", "DISODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC",
         "MESYLATE", "DIMESYLATE", "MALEATE", "TARTRATE", "BITARTRATE", "CITRATE", "SULFATE", "BISULFATE", "PHOSPHATE",
         "ACETATE", "DIACETATE", "SUCCINATE", "FUMARATE", "HEMIFUMARATE", "HYDROBROMIDE", "BROMIDE", "CHLORIDE", "TOSYLATE",
         "BESYLATE", "MALATE", "LACTATE", "GLUCONATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "HYDRATE", "ANHYDROUS",
         "TROMETHAMINE", "MEGLUMINE", "ESYLATE", "NAPSYLATE", "PAMOATE", "STEARATE", "NITRATE", "OXALATE", "SALICYLATE",
         "HYCLATE", "DECANOATE", "PROPIONATE", "VALERATE", "DIPROPIONATE", "BENZOATE", "CYPIONATE", "ENANTHATE"}


def truthy(v):
    return str(v) in ("True", "true", "1", "1.0")


def drug_name(ingredient):
    """79p display name (SALTS identical to 79p; a shorter list split names, e.g. "Nintedanib esylate", D107): salt words removed, capitalized"""
    return " ".join(w for w in str(ingredient).split() if w.upper() not in SALTS).capitalize()


def fda_rows():
    """(diseases, subgroups) DataFrames from the split tables, FDA drugs only, one row per cohort x ingredient;
    subgroups have group_id = stratum and parent_group."""
    fm = pd.read_csv(RES / "fda_matches.tsv", sep="\t").drop_duplicates("pert_id").set_index("pert_id").fda_ingredient
    g = pd.read_csv(RES / "drug_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False)
    s = pd.read_csv(RES / "stratum_reversal_split_v2_consensus.tsv", sep="\t", low_memory=False).rename(columns={"stratum": "group_id"})
    out = []
    for df, order in ((g, ["spec_mech", "spec", "rev"]), (s, ["spec_sib_mech", "spec_mech", "rev"])):
        df = df[df.passes_gate.map(truthy) & df.pert_id.isin(fm.index)].copy()
        df["fda_ingredient"] = df.pert_id.map(fm)
        df["drug"] = df.fda_ingredient.map(drug_name)
        df = df.sort_values(order, ascending=False).drop_duplicates(["group_id", "fda_ingredient"])
        out.append(apply_route(df))   # D107 local_route flag (kept, shown on the pages)
    return out[0], out[1]
