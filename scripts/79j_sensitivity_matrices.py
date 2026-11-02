#!/usr/bin/env python
"""Drug x cell-line sensitivity matrices from PRISM and GDSC, mapped to LINCS drugs and DepMap lines (runs in `celligner`).

Input: data/raw/prism_20q2, data/raw/gdsc_8.5 (scripts/09b). Four datasets, kept separate:
  prism_sec   PRISM Repurposing 20Q2 secondary screen, 8-dose curves: kill = 1 - AUC (AUC of fitted viability)
  prism_pri   PRISM Repurposing 20Q2 primary screen, one dose (2.5 uM): kill = - replicate-collapsed log2 fold change
  gdsc1/2     GDSC release 8.5 fitted dose response: kill = 1 - AUC (AUC of fitted viability over the tested range)
Higher kill = more sensitive. Several entries per drug x line (screens, doses, GDSC drug IDs of one compound) are
averaged. Lines: PRISM depmap_id = ModelID; GDSC SANGER_MODEL_ID / COSMIC_ID -> DepMap Model.csv. Lines failing
PRISM STR profiling and non-cancerous models are dropped.
Drugs -> LINCS pert_id: PRISM by Broad ID (first 13 characters = LINCS pert_id), else by name; GDSC by name (DRUG_NAME
and SYNONYMS, normalized: lower case, letters and digits only) against LINCS cmap_name.
Outputs: data/processed/sensitivity/<dataset>.tsv.gz (pert_id x ModelID kill, NaN = not screened),
         results/lincs/sensitivity_validation/overlap.tsv (per dataset: drugs, lines, mapped to LINCS, scored drugs)
"""
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
OUT = ROOT / "data/processed/sensitivity"
VAL = ROOT / "results/lincs/sensitivity_validation"
PRISM = RAW / "prism_20q2/prism-repurposing-20q2-"
norm = lambda s: re.sub(r"[^a-z0-9]", "", str(s).lower())


def lincs_names():
    c = pd.read_csv(RAW / "lincs2020/compoundinfo_beta.txt", sep="\t", usecols=["pert_id", "cmap_name"])
    scored = set(pd.read_csv(ROOT / "results/lincs/drug_reversal.tsv", sep="\t", usecols=["pert_id"]).pert_id)
    c = c.drop_duplicates()
    c["key"] = c.cmap_name.map(norm)
    return c, scored


def to_pert(df, names, idcol=None):
    """Long table with a 'name' (and optionally a Broad-ID) column -> rows per LINCS pert_id."""
    by_id = pd.Series(dtype=object)
    if idcol:
        df["pert_id"] = df[idcol].astype(str).str[:13].where(df[idcol].astype(str).str[:13].isin(names.pert_id))
    else:
        df["pert_id"] = np.nan
    miss = df[df.pert_id.isna()].drop(columns="pert_id")
    hit = df[df.pert_id.notna()]
    keyed = miss.assign(key=miss.name.map(norm)).merge(names[["key", "pert_id"]], on="key")
    return pd.concat([hit, keyed.drop(columns="key")], ignore_index=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    VAL.mkdir(parents=True, exist_ok=True)
    names, scored = lincs_names()
    models = pd.read_csv(RAW / "depmap_24q4/Model.csv")
    cancer = set(models.ModelID[models.OncotreePrimaryDisease != "Non-Cancerous"])
    mats, ov = {}, []

    # PRISM secondary: 1 - AUC
    s = pd.read_csv(f"{PRISM}secondary-screen-dose-response-curve-parameters.csv", low_memory=False)
    if "passed_str_profiling" in s:
        s = s[s.passed_str_profiling.astype(str).str.upper() == "TRUE"]
    s = s.rename(columns={"depmap_id": "ModelID"}).dropna(subset=["auc"])
    s["kill"] = 1 - s.auc
    mats["prism_sec"] = (s, "broad_id")

    # PRISM primary: - log fold change at 2.5 uM
    lfc = pd.read_csv(f"{PRISM}primary-screen-replicate-collapsed-logfold-change.csv", index_col=0)
    ti = pd.read_csv(f"{PRISM}primary-screen-replicate-collapsed-treatment-info.csv")
    cl = pd.read_csv(f"{PRISM}primary-screen-cell-line-info.csv")
    ok = set(cl.depmap_id[cl.passed_str_profiling.astype(str).str.upper() == "TRUE"]) if "passed_str_profiling" in cl \
        else set(cl.depmap_id)
    lfc.index = lfc.index.map(cl.set_index("row_name").depmap_id)
    lfc = lfc[lfc.index.isin(ok)]
    p = lfc.rename_axis("ModelID").reset_index().melt("ModelID", var_name="column_name", value_name="lfc").dropna()
    p = p.merge(ti[["column_name", "broad_id", "name"]], on="column_name")
    p["kill"] = -p.lfc
    mats["prism_pri"] = (p, "broad_id")

    # GDSC1 / GDSC2: 1 - AUC (xlsx read through R readxl; the celligner env has no Excel reader)
    comp = pd.read_csv(RAW / "gdsc_8.5/screened_compounds_rel_8.5.csv")
    syn = comp.assign(name=comp.SYNONYMS.fillna("").str.split(",")).explode("name")[["DRUG_ID", "name"]]
    syn = pd.concat([comp[["DRUG_ID", "DRUG_NAME"]].rename(columns={"DRUG_NAME": "name"}), syn])
    syn["name"] = syn.name.str.strip()
    syn = syn[syn.name != ""].drop_duplicates()
    sanger = models.dropna(subset=["SangerModelID"]).drop_duplicates("SangerModelID").set_index("SangerModelID").ModelID
    cosmic = models.dropna(subset=["COSMICID"]).drop_duplicates("COSMICID").set_index("COSMICID").ModelID
    for ds in ["GDSC1", "GDSC2"]:
        tsv = OUT / f"{ds}_raw.tsv"
        if not tsv.exists():
            subprocess.run([str(Path.home() / "miniforge3/envs/repurposing/bin/Rscript"), "-e",
                            f'd <- readxl::read_excel("{RAW}/gdsc_8.5/{ds}_fitted_dose_response_27Oct23.xlsx"); '
                            f'data.table::fwrite(d, "{tsv}", sep = "\\t")'], check=True)
        g = pd.read_csv(tsv, sep="\t")
        g["ModelID"] = g.SANGER_MODEL_ID.map(sanger).fillna(g.COSMIC_ID.map(cosmic))
        g = g.dropna(subset=["ModelID", "AUC"]).merge(syn, on="DRUG_ID")
        g["kill"] = 1 - g.AUC
        mats[ds.lower()] = (g, None)

    for ds, (d, idcol) in mats.items():
        n_drugs_all = d[idcol].astype(str).str[:13].nunique() if idcol else d.DRUG_ID.nunique()
        n_lines_all = d.ModelID.nunique()
        m = to_pert(d[[c for c in ["ModelID", "name", "kill", idcol] if c]].copy(), names, idcol)
        m = m[m.ModelID.isin(cancer)]
        # several names of one GDSC drug can hit the same pert_id: average per drug x line once
        X = m.groupby(["pert_id", "ModelID"]).kill.mean().unstack()
        X.to_csv(OUT / f"{ds}.tsv.gz", sep="\t", float_format="%.4f")
        ov.append(dict(dataset=ds, drugs=n_drugs_all, lines=n_lines_all, lincs_drugs_mapped=len(X),
                       lines_cancer_mapped=X.shape[1], scored_drugs_mapped=len(set(X.index) & scored),
                       entries=int(X.notna().sum().sum())))
        print(ov[-1], flush=True)
    pd.DataFrame(ov).to_csv(VAL / "overlap.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
