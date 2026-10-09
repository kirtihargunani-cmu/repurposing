#!/usr/bin/env python
"""Drugs available per tissue and per disease: measured directly in the tissue's LINCS lines, and borrowable as the
cross-line consensus of drugs active elsewhere (scripts/77_lincs_extrapolation_test.py).

Tissue = DepMap OncotreeLineage of LINCS lines mapped to DepMap (results/lincs/cell_coverage.tsv); cancer lines by
default, non-cancerous LINCS lines (HA1E kidney, MCF 10A breast) counted in separate columns. Counts from siginfo
metadata: tested = any compound signature; high_quality = is_hiq; active = is_hiq and TAS >= 0.4.
Borrowable = not active in the tissue, but active in >= k cancer lines elsewhere (k = 1, 3, 5). Leave-one-line-out
accuracy of the cross-line mean by number of lines: 1 r 0.30, 2 0.38, 3-4 0.46, 5-8 0.53, 9+ 0.57 (ceiling 0.55).
Diseases (41 Step 4 signature groups) use their strict tissue key (config/celligner_tissue_keys.tsv).
Outputs: results/lincs/tissue_drug_counts.tsv, disease_drug_counts.tsv
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/lincs2020"
RES = ROOT / "results/lincs"
TAS_MIN, KS = 0.4, [1, 3, 5]


def counts(s, act_lines_per_drug, own_lines):
    """tested / high-quality / active drug sets in lines `own_lines`, plus borrowable counts"""
    x = s[s.ModelID.isin(own_lines)]
    a = set(x.loc[x.active, "pert_id"])
    out = dict(lincs_lines=x.ModelID.nunique(), drugs_tested=x.pert_id.nunique(),
               drugs_high_quality=x.loc[x.is_hiq == 1, "pert_id"].nunique(), drugs_active=len(a))
    # active lines outside the tissue
    elsewhere = act_lines_per_drug[~act_lines_per_drug.ModelID.isin(own_lines)].groupby("pert_id").ModelID.nunique()
    for k in KS:
        b = set(elsewhere[elsewhere >= k].index) - a
        out[f"borrowable_ge{k}lines"] = len(b)
        out[f"total_ge{k}lines"] = len(a) + len(b)
    return out


def main():
    models = pd.read_csv(ROOT / "data/raw/depmap_24q4/Model.csv", index_col=0)
    cov = pd.read_csv(RES / "cell_coverage.tsv", sep="\t").dropna(subset=["ModelID"])
    s = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                    usecols=["pert_id", "pert_type", "cell_iname", "is_hiq", "tas"])
    s = s[s.pert_type == "trt_cp"].merge(cov[["cell_iname", "ModelID"]], on="cell_iname")
    s = s.join(models[["OncotreeLineage", "OncotreePrimaryDisease"]], on="ModelID")
    s["active"] = (s.is_hiq == 1) & (s.tas >= TAS_MIN)
    s["cancer"] = s.OncotreePrimaryDisease != "Non-Cancerous"
    sc = s[s.cancer]
    act_lines = sc.loc[sc.active, ["pert_id", "ModelID"]].drop_duplicates()
    print(f"active drugs in cancer lines: {act_lines.pert_id.nunique()}; in >= 3 lines "
          f"{(act_lines.groupby('pert_id').size() >= 3).sum()}, >= 5 {(act_lines.groupby('pert_id').size() >= 5).sum()}")

    cancer_models = models[models.OncotreePrimaryDisease != "Non-Cancerous"]
    rows = []
    for lin, m in cancer_models.groupby("OncotreeLineage"):
        r = dict(tissue=lin, depmap_lines=len(m), **counts(sc, act_lines, set(m.index)))
        nc = s[~s.cancer & (s.OncotreeLineage == lin)]
        r["noncancer_lincs_lines"] = ", ".join(sorted(set(nc.cell_iname)))
        r["noncancer_drugs_active"] = nc.loc[nc.active, "pert_id"].nunique()
        rows.append(r)
    t = pd.DataFrame(rows).sort_values("drugs_active", ascending=False)
    t.to_csv(RES / "tissue_drug_counts.tsv", sep="\t", index=False)

    sig = pd.read_csv(ROOT / "results/step4/signature_summary.tsv", sep="\t")
    groups = sorted(set(sig.loc[sig.status == "ok", "group_id"]))
    keys = pd.read_csv(ROOT / "config/celligner_tissue_keys.tsv", sep="\t", comment="#").set_index("group_id")
    drows = []
    for g in groups:
        k = keys.loc[g]
        m = cancer_models[cancer_models.OncotreeLineage.isin(str(k.lineages).split(";"))]
        if isinstance(k.diseases, str) and k.diseases:
            m = m[m.OncotreePrimaryDisease.isin(k.diseases.split(";"))]
        drows.append(dict(group_id=g, key=k.key, **counts(sc, act_lines, set(m.index))))
    d = pd.DataFrame(drows)
    d.to_csv(RES / "disease_drug_counts.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(t.drop(columns="noncancer_lincs_lines").to_string(index=False))
    print(t.loc[t.noncancer_lincs_lines != "", ["tissue", "noncancer_lincs_lines", "noncancer_drugs_active"]].to_string(index=False))
    print(d.to_string(index=False))


if __name__ == "__main__":
    main()
