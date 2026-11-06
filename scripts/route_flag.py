"""Route flag (D107, user 2026-11-06): drugs with no systemic formulation cannot reach a tumor at a useful dose.

local_route = True when a drug has no approved way into the bloodstream:
  FDA-matched drugs (79p, fda_matches.tsv): every Drugs@FDA product of the ingredient (any marketing status, alone or in
    combination) has only local routes. Routes are the part of Products.txt `Form` after ";" (comma-separated).
  Other drugs: ChEMBL topical and not oral / parenteral (79a route_topical_only).
Systemic routes (SYSTEMIC below): oral and other routes that deliver to the blood by design (sublingual, buccal,
transmucosal, transdermal, rectal, enteral) and every injection route that reaches the circulation. Everything else is
local: topical, ophthalmic, otic, nasal, inhalation, vaginal, dental, intravitreal, intrathecal, intravesical,
intraperitoneal, intra-articular, irrigation, implantation, etc. Nasal and inhaled products of systemic drugs exist,
but a drug approved ONLY nasally / inhaled has no established systemic dose, which is what the flag is about.
Rectal counts as systemic (diazepam, acetaminophen); a rectal-only local drug (mesalamine enema) is missed, accepted.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results/lincs"
FDA = ROOT / "data/raw/drugsatfda/extracted"
SYSTEMIC = {"ORAL", "ORAL-28", "ORAL-21", "ORAL-20", "ORAL SUSPENSION", "ORALLY DISINTEGRATING", "ORAL AND TABLET",
            "SUBLINGUAL", "BUCCAL", "TRANSMUCOSAL", "TRANSDERMAL", "RECTAL", "ENTERAL",
            "INJECTION", "INJECTABLE", "INTRAVENOUS", "IV (INFUSION)", "IM-IV", "INTRAMUSCULAR", "SUBCUTANEOUS",
            "SUBCUTANEOUS LYOPHILIZED POWER", "INTRA-ARTERIAL", "INTRAOSSEOUS", "INTRADERMAL"}


def routes(form):
    if not isinstance(form, str) or ";" not in form:
        return set()
    return {r.strip().upper() for r in form.split(";", 1)[1].split(",") if r.strip()}


def fda_routes(norm):
    """normalized ingredient key -> set of routes over all Drugs@FDA products (norm = 79p's name normalizer)"""
    p = pd.read_csv(FDA / "Products.txt", sep="\t", dtype=str, encoding="latin-1", on_bad_lines="skip")
    p["ingredient"] = p.ActiveIngredient.str.split(";")
    p = p.explode("ingredient")
    p["key"] = p.ingredient.map(norm)
    p["r"] = p.Form.map(routes)
    p = p[p.key != ""]
    return p.groupby("key").r.agg(lambda s: set().union(*s))


def is_local(rs):
    """True if routes are known and none is systemic; None if no route is known"""
    return None if not rs else not (rs & SYSTEMIC)


def local_route_by_pert():
    """pert_id -> (local_route bool, route source, routes) from fda_matches.tsv (79p) and drug_annotations.tsv (79a)"""
    ann = pd.read_csv(RES / "drug_annotations.tsv", sep="\t", low_memory=False).drop_duplicates("pert_id").set_index("pert_id")
    out = pd.DataFrame(index=ann.index)
    out["local_route"] = ann.route_topical_only.astype(str).str.lower().eq("true")
    out["route_source"] = "chembl"
    out["routes"] = ""
    f = RES / "fda_matches.tsv"
    if f.exists():
        m = pd.read_csv(f, sep="\t").drop_duplicates("pert_id").set_index("pert_id")
        if "local_route" in m:
            k = m.index[m.local_route.notna()]
            out.loc[out.index.intersection(k), "local_route"] = m.loc[out.index.intersection(k), "local_route"].astype(str).str.lower().eq("true")
            out.loc[out.index.intersection(k), "route_source"] = "drugsatfda"
            out.loc[out.index.intersection(k), "routes"] = m.loc[out.index.intersection(k), "fda_routes"].fillna("")
    return out


BLADDER = ("TCGA-BLCA", "BLCA-")   # D107 exception: intravesical reaches bladder tumors


def apply_route(df, unit="group_id"):
    """adds local_route (bool), route_source, routes to a cohort x drug table with pert_id; D107 bladder exception"""
    r = local_route_by_pert()
    df = df.drop(columns=[c for c in ("local_route", "route_source", "routes") if c in df]).join(r, on="pert_id")
    df["local_route"] = df.local_route.fillna(False).astype(bool)
    bl = df[unit].astype(str).str.startswith(BLADDER) & df.routes.fillna("").str.contains("INTRAVESICAL")
    df.loc[bl, "local_route"] = False
    return df
