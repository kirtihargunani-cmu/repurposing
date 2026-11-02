#!/usr/bin/env python3
"""Export the real numbers the animated deck needs into presentation/scene_data.json.

Every number shown on screen should come from this file, never typed by hand. Run from anywhere:
    python3 presentation/export_scene_data.py
Reads only git-tracked results on the celligner branch (standard library only). On main, where the drug-side v2
results are not present, use the committed scene_data.json or run this from the celligner worktree.
"""
import csv
import json
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LINCS = REPO / "results" / "lincs"
CELL = REPO / "results" / "celligner"


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x, nd=3):
    try:
        return round(float(x), nd)
    except (TypeError, ValueError):
        return None


NEEDED = [CELL / "disease_line_weights_long.tsv"] + [LINCS / n for n in (
    "worked_example_regulators.tsv", "split_controls_v2_consensus.tsv", "drug_reversal_summary_v2_consensus.tsv",
    "data_gaps_v2.tsv", "fda_ranking_v2_consensus.tsv", "fda_ranking_strata_v2_consensus.tsv")]
missing = [str(p.relative_to(REPO)) for p in NEEDED if not p.exists()]
if missing:
    raise SystemExit("Missing drug-side results (they live on the celligner branch):\n  " + "\n  ".join(missing)
                     + "\nKeep the committed scene_data.json, or run this script from ~/projects/repurposing-celligner.")

out = {}

# Step 1: Celligner line weights for the example disease (top 8 lines).
w = [r for r in rows(CELL / "disease_line_weights_long.tsv") if r["group_id"] == "TCGA-LUAD"]
w.sort(key=lambda r: -float(r["weight"]))
out["celligner_luad_top_lines"] = [
    {"line": r["CellLineName"], "weight": f(r["weight"], 4), "lineage": r["OncotreeLineage"],
     "in_tissue_key": r["in_key"] == "True"}
    for r in w[:8]
]

# Step 6 worked example: per-regulator contributions (top 12 each way + totals).
we = rows(LINCS / "worked_example_regulators.tsv")
examples = {}
for key in sorted({(r["group_id"], r["cmap_name"]) for r in we}):
    sub = [r for r in we if (r["group_id"], r["cmap_name"]) == key]
    sub.sort(key=lambda r: -float(r["contribution"]))
    examples[f"{key[0]}|{key[1]}"] = {
        "rev": f(sum(float(r["contribution"]) for r in sub)),
        "n_regulators": len(sub),
        "n_positive": sum(float(r["contribution"]) > 0 for r in sub),
        "top_positive": [
            {"regulator": r["regulator"], "disease_nes": f(r["weight_used"]), "drug_activity": f(r["drug_activity"]),
             "contribution": f(r["contribution"]), "prolif_linked": r["prolif_linked"] == "True"}
            for r in sub[:12]
        ],
        "top_negative": [
            {"regulator": r["regulator"], "disease_nes": f(r["weight_used"]), "drug_activity": f(r["drug_activity"]),
             "contribution": f(r["contribution"])}
            for r in sub[-6:]
        ],
    }
out["worked_examples"] = examples

# Known-drug controls (percentiles).
out["controls"] = [
    {k: (r[k] if k in ("unit", "expected", "drug") else f(r[k], 1))
     for k in ("unit", "expected", "drug", "rev", "spec", "spec_mech", "prolif_share")}
    for r in rows(LINCS / "split_controls_v2_consensus.tsv")
]

# Evidence status per disease (D78).
summ = rows(LINCS / "drug_reversal_summary_v2_consensus.tsv")
out["evidence_status_counts"] = dict(Counter(r["evidence_status"] for r in summ))
out["evidence_status"] = [
    {"group": r["group_id"], "status": r["evidence_status"], "best_line": r["best_line"],
     "best_line_match": f(r["best_line_match"], 4)}
    for r in summ
]

# Data gaps (C24).
gaps = rows(LINCS / "data_gaps_v2.tsv")
out["data_gap_counts"] = dict(Counter(r["data_gap_step"] or "pass" for r in gaps))
out["data_gaps"] = [{"group": r["group_id"], "step": r["data_gap_step"] or "pass"} for r in gaps]

# Top 3 FDA drugs per disease and subgroup by spec_mech / spec_sib_mech.
def top3(path, key):
    d = {}
    for r in rows(path):
        if int(float(r["fda_rank"])) <= 3:
            d.setdefault(r[key], []).append(r["fda_ingredient"].title())
    return d

out["fda_top3"] = top3(LINCS / "fda_ranking_v2_consensus.tsv", "group_id")
out["fda_top3_subgroups"] = top3(LINCS / "fda_ranking_strata_v2_consensus.tsv", "stratum")

# Fixed numbers from the docs (disease side on main, drug side spec); kept here so scenes read one file.
out["pipeline_facts"] = {
    "disease_groups": 49, "patients": 10211, "genes_v2": 17784, "signature_groups": 41,
    "aracne_networks": 29, "subgroups": 14, "lincs_drugs": 3058, "fda_matched_drugs": 527,
    "generic_drugs_by_rev": 32, "low_transfer_share_pct": 11.5, "knockdown_tfs": 336,
}

dest = Path(__file__).resolve().parent / "scene_data.json"
dest.write_text(json.dumps(out, indent=1))
print(f"wrote {dest} ({dest.stat().st_size // 1024} KB)")
