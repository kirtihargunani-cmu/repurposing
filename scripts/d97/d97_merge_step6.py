#!/usr/bin/env python
"""D97: merge the Step 6 runs of the switched groups (RES_TAG=_d97run) into the production tables: rows of the run's
groups replace production rows in every results table; per-group files in the processed folder are copied over."""
import shutil, sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent.parent
for base in ["aracne_project_v2_tf", "aracne_project_v2_tf_strata"]:
    run_r, prod_r = ROOT / f"results/step6/{base}_d97run", ROOT / f"results/step6/{base}"
    run_p, prod_p = ROOT / f"data/processed/step6/{base}_d97run", ROOT / f"data/processed/step6/{base}"
    if not run_r.exists():
        continue
    for f in sorted(run_r.glob("*.tsv")):
        new = pd.read_csv(f, sep="\t") if f.stat().st_size > 1 else pd.DataFrame()
        if new.empty or "group_id" not in new or not (prod_r / f.name).exists():
            print(f"{base}/{f.name}: skipped (empty or no group_id)"); continue
        old = pd.read_csv(prod_r / f.name, sep="\t")
        out = pd.concat([old[~old.group_id.isin(new.group_id)], new], ignore_index=True)
        out.to_csv(prod_r / f.name, sep="\t", index=False)
        print(f"{base}/{f.name}: replaced {old.group_id.isin(new.group_id).sum()} rows of {new.group_id.nunique()} groups")
    if run_p.exists():
        for f in run_p.iterdir():
            shutil.copy2(f, prod_p / f.name); print(f"copied {base}/{f.name}")
