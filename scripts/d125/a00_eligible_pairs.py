#!/usr/bin/env python
"""D125 A0.0: eligible-pair gate for A0b (docs/SPEC_graph_metric_scoring.md section A0.0; D125-A0 (1)).

Counts the drug x line pairs (79c profiles) with >= A00_MIN_REPLICATES active level-5 signatures (is_hiq, TAS >= 0.4: the
signatures 79c averages) and the lines they fall in. Fewer than A00_MIN_PAIRS pairs, or fewer than A00_MIN_LINES lines ->
A0b is confined to the qualifying lines and reported as a line-level bound; else A0b runs at cohort level. Level-4
replicate counts (siginfo distil_ids; profiles not on disk) and the A0a strata ((dose, time) cells) are reported per pair.
Outputs: results/lincs/d125/a0/eligible_pairs.tsv (per pair), eligible_pairs_summary.tsv (gate and per-line counts),
         eligible_pairs_gate.json (the verdict A0b reads).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
# --- config block (values in a0_config.py only) ---
from a0_config import A00_MIN_REPLICATES, A00_MIN_PAIRS, A00_MIN_LINES, A0A_STRATUM, A0A_MIN_STRATA, TAS_MIN  # noqa: E402
from a0_common import OUT, ROOT, active_signatures, production_index, init_work, write_json, log  # noqa: E402


def main():
    init_work()
    act = active_signatures()
    idx = production_index()
    g = act.groupby("col")
    per = idx.assign(
        n_active_l5=g.size().reindex(idx.col).values,
        n_l4_replicates=g.n_l4.sum().reindex(idx.col).values,
        n_dose_time_cells=g.apply(lambda d: len(d[list(A0A_STRATUM)].drop_duplicates())).reindex(idx.col).values)
    names = act.drop_duplicates("col").set_index("col")[["cmap_name", "cell_iname"]]
    per = per.join(names, on="col")
    per["a0a_eligible"] = per.n_dose_time_cells >= A0A_MIN_STRATA
    per["a0b_eligible"] = per.n_active_l5 >= A00_MIN_REPLICATES
    per["l4_ge_min"] = per.n_l4_replicates >= A00_MIN_REPLICATES     # reported only (level-4 reading of A0.0)
    per = per[["col", "pert_id", "cmap_name", "ModelID", "cell_iname", "n_active_l5", "n_l4_replicates",
               "n_dose_time_cells", "a0a_eligible", "a0b_eligible", "l4_ge_min"]]
    per.to_csv(OUT / "eligible_pairs.tsv", sep="\t", index=False)

    e = per[per.a0b_eligible]
    n_pairs, n_lines = len(e), e.ModelID.nunique()
    cohort_level = n_pairs >= A00_MIN_PAIRS and n_lines >= A00_MIN_LINES
    l4 = per[per.l4_ge_min]
    lines = (per.groupby(["ModelID", "cell_iname"])
             .agg(profiles=("col", "size"), a0b_pairs=("a0b_eligible", "sum"), a0a_pairs=("a0a_eligible", "sum"),
                  l4_ge_min_pairs=("l4_ge_min", "sum"))
             .reset_index().sort_values("a0b_pairs", ascending=False))
    head = pd.DataFrame([
        dict(item="definition", value=f"active level-5 signatures (is_hiq, TAS >= {TAS_MIN}) per 79c pair >= {A00_MIN_REPLICATES}"),
        dict(item="profiles (79c pairs)", value=len(per)),
        dict(item="active signatures", value=int(per.n_active_l5.sum())),
        dict(item="eligible pairs", value=n_pairs),
        dict(item="eligible lines", value=n_lines),
        dict(item="gate", value=f">= {A00_MIN_PAIRS} pairs in >= {A00_MIN_LINES} lines"),
        dict(item="A0b scope", value="cohort-level" if cohort_level else "line-level bound (qualifying lines only)"),
        dict(item="reported: pairs with >= 6 level-4 replicates", value=len(l4)),
        dict(item="reported: their lines", value=l4.ModelID.nunique()),
        dict(item="A0a pairs (>= 2 dose x time cells)", value=int(per.a0a_eligible.sum())),
    ])
    with open(OUT / "eligible_pairs_summary.tsv", "w") as fh:
        head.to_csv(fh, sep="\t", index=False)
        fh.write("\n")
        lines.to_csv(fh, sep="\t", index=False)
    write_json(OUT / "eligible_pairs_gate.json", dict(
        n_pairs=n_pairs, n_lines=n_lines, min_pairs=A00_MIN_PAIRS, min_lines=A00_MIN_LINES,
        min_replicates=A00_MIN_REPLICATES, cohort_level=cohort_level,
        qualifying_lines=sorted(e.ModelID.unique().tolist())))
    log(f"A0.0: {n_pairs} pairs with >= {A00_MIN_REPLICATES} active signatures in {n_lines} lines -> "
        f"A0b {'cohort-level' if cohort_level else 'line-level bound'} (level-4 reading: {len(l4)} pairs, "
        f"{l4.ModelID.nunique()} lines); A0a pairs {int(per.a0a_eligible.sum())}")
    print(lines.head(12).to_string(index=False))


if __name__ == "__main__":
    main()
