#!/usr/bin/env python
"""Where good-quality data run out, per disease group and subgroup (user 2026-11-01, C24; LLM notes R3).

Rule (user): if there is no good-quality, lineage-matched data, report "no good-quality data at <step>" instead of a
drug list. Lineage-matched = the cell line's DepMap OncotreeCode is in the group's lineage codes
(config/group_lineages.tsv), whatever its organ.
Checks, in pipeline order; the first that fails is the group's data gap:
  D1 cell-line matching (DepMap)   at least one lineage-matched DepMap line exists
  D2 drug data (LINCS)             at least one lineage-matched line has good-quality compound profiles
                                   (is_hiq and TAS >= 0.4, the drug side's active signatures; drug_line_profiles_index.tsv)
  D1b line resemblance (Celligner) one of those lines resembles the disease: line match >= 0.1 (D76 cutoff). Line match
                                   here = the line's mean Celligner weight over the group's patients, weights computed
                                   over all lineage-matched lines (softmax, same tau as 73), divided by the top lineage
                                   line's weight
  D2b enough drugs (LINCS)         the matched lines (D1b) together carry >= MIN_COMPOUNDS compounds (10; user 2026-11-01)
Reported for every group: n lineage lines in DepMap, those with LINCS data (name, compounds, line match), the best one.
Outputs: results/lincs/data_gaps_v2.tsv (groups) and data_gaps_strata_v2.tsv (Step 10 subgroups; parent's lineage,
         the subgroup's own patients). Columns data_gap_step ("" = none), good_quality_data (bool), data_gap_note.
This report does not change the rankings; the layer must show data_gap_note instead of a drug list when
good_quality_data is False.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data/processed"
RES = ROOT / "results/lincs"
LINE_MATCH_MIN = 0.1  # D76
MIN_COMPOUNDS = int(os.environ.get("MIN_COMPOUNDS", "10"))  # compounds on matched lineage lines (user 2026-11-01)
TOPK = 10             # tau rule of 71 / 73
STEP = {"D1": "drug-side Step 1, cell-line matching (DepMap)", "D2": "drug-side Step 2, drug data (LINCS)",
        "D1b": "drug-side Step 1, cell-line matching (Celligner resemblance)",
        "D2b": "drug-side Step 2, drug data (LINCS: too few compounds in matched lines)"}


def main():
    lin = pd.read_csv(ROOT / "config/group_lineages.tsv", sep="\t", comment="#").set_index("group_id")
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    dist = pd.read_csv(PROC / "celligner/dist_tumor_model.tsv.gz", sep="\t", index_col=0)
    models = models[models.index.isin(dist.columns) & (models.OncotreePrimaryDisease != "Non-Cancerous")]
    dist = dist[models.index]
    ds = np.sort(dist.values, axis=1)
    tau = float(np.median(ds[:, TOPK - 1] - ds[:, 0]) / np.log(10))
    qc = pd.read_csv(ROOT / "results/celligner/patient_qc.tsv", sep="\t", index_col=0)
    prof = pd.read_csv(PROC / "lincs/drug_line_profiles_index.tsv", sep="\t")
    ncpd = prof.groupby("ModelID").pert_id.nunique()

    def assess(name, lineage_row, samples, parent):
        codes = [c for c in str(lineage_row.codes).split(";") if c and c != "nan"]
        L = list(models.index[models.OncotreeCode.isin(codes)])
        r = dict(group_id=name, parent=parent, lineage=lineage_row.lineage, n_lineage_lines=len(L))
        withdata = [m for m in L if m in ncpd.index]
        r["n_lineage_lines_lincs"] = len(withdata)
        r["n_compounds_lineage"] = int(sum(ncpd[m] for m in withdata))
        s = [x for x in samples if x in dist.index]
        if L and s:
            d = dist.loc[s, L].values
            w = np.exp(-(d - d.min(axis=1, keepdims=True)) / tau)
            wm = pd.Series((w / w.sum(axis=1, keepdims=True)).mean(axis=0), index=L)
            lm = wm / wm.max()
            top = wm.idxmax()
            r["top_lineage_line"] = f"{models.CellLineName[top]} ({models.OncotreeCode[top]})"
        else:
            lm = pd.Series(dtype=float)
            r["top_lineage_line"] = ""
        r["lincs_lineage_lines"] = "; ".join(
            f"{models.CellLineName[m]} ({models.OncotreeCode[m]}, {ncpd[m]} compounds, match {lm.get(m, 0):.3f})"
            for m in sorted(withdata, key=lambda m: -lm.get(m, 0)))
        matched = [m for m in withdata if lm.get(m, 0) >= LINE_MATCH_MIN]
        r["n_matched_lines"] = len(matched)
        r["n_compounds_matched"] = int(prof[prof.ModelID.isin(matched)].pert_id.nunique())
        # best line: the matched line with most compounds; without a matched line, the closest line with data
        best = (max(matched, key=lambda m: ncpd[m]) if matched else
                max(withdata, key=lambda m: lm.get(m, 0)) if withdata else None)
        r["best_line"] = models.CellLineName[best] if best else ""
        r["best_line_match"] = round(float(lm.get(best, 0)), 4) if best else 0.0
        r["best_line_compounds"] = int(ncpd[best]) if best else 0
        if not L:
            gap, why = "D1", f"DepMap has no {lineage_row.lineage} cell line"
        elif not withdata:
            gap, why = "D2", (f"none of the {len(L)} {lineage_row.lineage} lines in DepMap has good-quality LINCS compound "
                              f"profiles")
        elif r["best_line_match"] < LINE_MATCH_MIN:
            gap, why = "D1b", (f"the {lineage_row.lineage} lines with LINCS data do not resemble this disease (best "
                               f"{r['best_line']}, line match {r['best_line_match']:.3f} < {LINE_MATCH_MIN})")
        elif r["n_compounds_matched"] < MIN_COMPOUNDS:
            gap, why = "D2b", (f"the matched {lineage_row.lineage} lines carry only {r['n_compounds_matched']} compounds "
                               f"(minimum {MIN_COMPOUNDS})")
        else:
            gap, why = "", ""
        r["data_gap_step"] = gap
        r["good_quality_data"] = gap == ""
        r["data_gap_note"] = (f"No good-quality data at {STEP[gap]}: {why}. A lineage-matched drug list cannot be reported."
                              if gap else f"Lineage-matched data: {r['n_compounds_matched']} compounds in {r['n_matched_lines']} "
                                          f"matched line(s); best {r['best_line']} ({r['best_line_compounds']} compounds, "
                                          f"line match {r['best_line_match']:.2f}).")
        return r

    rows = [assess(g, lin.loc[g], qc.index[qc.group_id == g], g) for g in lin.index]
    out = pd.DataFrame(rows)
    out.to_csv(RES / "data_gaps_v2.tsv", sep="\t", index=False)

    sp = pd.read_csv(ROOT / "results/step10_v2/stratum_patients.tsv", sep="\t")
    srows = [assess(st, lin.loc[d.group_id.iloc[0]], d["sample"].tolist(), d.group_id.iloc[0])
             for st, d in sp.groupby("stratum") if d.group_id.iloc[0] in lin.index]
    sout = pd.DataFrame(srows)
    sout.to_csv(RES / "data_gaps_strata_v2.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    for o in (out, sout):
        print(o[["group_id", "n_lineage_lines", "n_lineage_lines_lincs", "best_line", "best_line_compounds",
                 "best_line_match", "n_matched_lines", "n_compounds_matched", "data_gap_step"]].to_string(index=False))
    print(f"\ngroups with good-quality lineage-matched data: {int(out.good_quality_data.sum())} of {len(out)}; "
          f"subgroups: {int(sout.good_quality_data.sum())} of {len(sout)}")


if __name__ == "__main__":
    main()
