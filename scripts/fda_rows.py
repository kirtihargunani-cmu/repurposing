"""Shared selection of FDA-approved drug rows per disease and subgroup (used by 79o MODE=fda_all and 79x).

Rows: drugs that pass the target gate (79b) and are matched to a Drugs@FDA ingredient (79p, fda_matches.tsv; marketed
or discontinued). Several LINCS ids of one ingredient -> the best on the cohort's lead score (diseases: spec_mech, then
spec, rev; subgroups: spec_sib_mech, then spec_mech, rev), among the compounds 79p would rank (eligible(); D115), else
among all. Unlike 79p, flagged (incl. D107 local_route) and non-reversing drugs are
kept (the page shows them, weak evidence hidden by default). `drug` = FDA ingredient with salt words removed (79p rule).
"""
from pathlib import Path

import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from route_flag import apply_route  # noqa: E402
from evidence_flags import apply_evidence  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
SALTS = {"HYDROCHLORIDE", "DIHYDROCHLORIDE", "HCL", "SODIUM", "DISODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC",
         "MESYLATE", "DIMESYLATE", "MALEATE", "TARTRATE", "BITARTRATE", "CITRATE", "SULFATE", "BISULFATE", "PHOSPHATE",
         "ACETATE", "DIACETATE", "SUCCINATE", "FUMARATE", "HEMIFUMARATE", "HYDROBROMIDE", "BROMIDE", "CHLORIDE", "TOSYLATE",
         "BESYLATE", "MALATE", "LACTATE", "GLUCONATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "HYDRATE", "ANHYDROUS",
         "TROMETHAMINE", "MEGLUMINE", "ESYLATE", "NAPSYLATE", "PAMOATE", "STEARATE", "NITRATE", "OXALATE", "SALICYLATE",
         "HYCLATE", "DECANOATE", "PROPIONATE", "VALERATE", "DIPROPIONATE", "BENZOATE", "CYPIONATE", "ENANTHATE",
         # D112: salts / esters found missing, and griseofulvin particle-size labels
         "DITOSYLATE", "CAMSYLATE", "ISETHIONATE", "XINAFOATE", "PIVALATE",
         "MICROSIZE", "ULTRAMICROSIZE", "MICROCRYSTALLINE", "ULTRAMICROCRYSTALLINE"}
# D112: LINCS cmap_name -> Drugs@FDA ingredient (INN vs USAN spellings; reviewed one by one)
FDA_ALIAS = {"tioguanine": "THIOGUANINE", "cyclosporin-a": "CYCLOSPORINE", "noretynodrel": "NORETHYNODREL", "clomifene": "CLOMIPHENE"}
# D112: LINCS cmap_names never matched to Drugs@FDA (wrong ChEMBL / Hub names): MRE-269 is selexipag's metabolite, not treprostinil
FDA_EXCLUDE = {"MRE-269"}


def truthy(v):
    return str(v) in ("True", "true", "1", "1.0")


def drug_name(ingredient):
    """79p display name (SALTS identical to 79p; a shorter list split names, e.g. "Nintedanib esylate", D107): salt words removed, capitalized"""
    return " ".join(w.strip(",") for w in str(ingredient).split() if w.strip(",").upper() not in SALTS).strip(" ,").capitalize()   # D112: "GRISEOFULVIN, ULTRAMICROSIZE"


def eligible(df, unit="group_id"):
    """D115: the 79p eligibility test, one definition for 79p and fda_rows. Needs apply_route / apply_evidence columns."""
    t = lambda c: df[c].map(truthy) if c in df else pd.Series(False, index=df.index)
    rs = df["receptor_status"] if "receptor_status" in df else pd.Series("", index=df.index)
    return (df.passes_gate.map(truthy) & ~t("unmatched_lines") & ~t("low_transfer") & (df.rev > 0) & ~t("target_absent")
            & ~rs.eq("receptor_mismatch") & ~t("local_route") & ~t("weak_d109"))   # D102c / D104c / D107 / D109


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
        df = apply_evidence(apply_route(df))   # D107 / D109 flags (kept, shown on the pages)
        df["eligible"] = eligible(df)
        # D115: best eligible compound of the ingredient (= the one 79p ranks); else the best compound overall, with its flags
        df = df.sort_values(["eligible"] + order, ascending=False).drop_duplicates(["group_id", "fda_ingredient"])
        out.append(df)
    return out[0], out[1]
