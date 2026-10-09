#!/usr/bin/env python
"""LINCS 2020 coverage: which cell lines were profiled with which compounds, linked to the Celligner matrices.

Sources: data/raw/lincs2020 (scripts/07_download_lincs.sh). Compound signatures only (pert_type trt_cp; shRNA /
CRISPR / ORF / ligand signatures are counted but not used). A signature is LINCS's consensus of replicates for one
compound x dose x time x line; is_hiq marks the CMap high-quality set (replicate agreement).
Cell line IDs: LINCS cell_iname -> DepMap ModelID by CCLE name, then Cellosaurus RRID, then stripped name.
Coverage tiers per line: compounds with >= 1 high-quality signature (n_cpd_hiq).
Disease coverage (Celligner hybrid weights, results/celligner/disease_line_weights_long.tsv): share of each
disease's line weight on LINCS-profiled lines, and per compound the weight of the disease's lines it was tested in.
Unprofiled weighted lines: distance to the nearest profiled line (data/processed/celligner/dist_line_line.tsv.gz).
Outputs: results/lincs/{cell_coverage, compounds, disease_coverage, disease_line_coverage}.tsv
         data/processed/lincs/{compound_line.tsv.gz, disease_compound_coverage.tsv.gz}
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/lincs2020"
DEPMAP = ROOT / "data/raw/depmap_24q4"
RES = ROOT / "results/lincs"
PROC = ROOT / "data/processed/lincs"
TIERS = [1, 100, 1000]  # minimum compounds (high-quality) for a line to count as profiled


def strip(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def map_lines(cells, models):
    ccle = models.reset_index().dropna(subset=["CCLEName"]).drop_duplicates("CCLEName").set_index("CCLEName").ModelID
    rrid = models.reset_index().dropna(subset=["RRID"]).drop_duplicates("RRID").set_index("RRID").ModelID
    st = models.reset_index().assign(s=models.StrippedCellLineName.map(strip).values)
    st = st[~st.s.duplicated(keep=False)].set_index("s").ModelID
    out = []
    for _, c in cells.iterrows():
        mid, how = None, None
        if isinstance(c.ccle_name, str) and c.ccle_name in ccle.index:
            mid, how = ccle[c.ccle_name], "ccle_name"
        elif isinstance(c.cellosaurus_id, str) and c.cellosaurus_id in rrid.index:
            mid, how = rrid[c.cellosaurus_id], "rrid"
        elif strip(c.cell_iname) in st.index:
            mid, how = st[strip(c.cell_iname)], "name"
        out.append((c.cell_iname, mid, how))
    return pd.DataFrame(out, columns=["cell_iname", "ModelID", "map_method"])


def main():
    RES.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    models = pd.read_csv(DEPMAP / "Model.csv", index_col=0)
    cells = pd.read_csv(RAW / "cellinfo_beta.txt", sep="\t", na_values=['""', ""], keep_default_na=False)
    sig = pd.read_csv(RAW / "siginfo_beta.txt", sep="\t", low_memory=False,
                      usecols=["sig_id", "pert_id", "cmap_name", "pert_type", "cell_iname", "pert_dose", "pert_time",
                               "is_hiq", "is_exemplar_sig", "tas"])
    cpd = pd.read_csv(RAW / "compoundinfo_beta.txt", sep="\t", na_values=['""'], keep_default_na=False)

    m = map_lines(cells, models)
    cp = sig[sig.pert_type == "trt_cp"]
    per = cp.groupby("cell_iname").agg(n_sig=("sig_id", "size"), n_cpd=("pert_id", "nunique"),
                                       n_sig_hiq=("is_hiq", "sum"))
    per["n_cpd_hiq"] = cp[cp.is_hiq == 1].groupby("cell_iname").pert_id.nunique()
    per["n_times"] = cp.groupby("cell_iname").pert_time.nunique()
    other = sig[sig.pert_type.str.startswith("trt_") & (sig.pert_type != "trt_cp")].groupby("cell_iname").size()
    cov = cells[["cell_iname", "cell_type", "cell_lineage", "primary_disease", "subtype", "ccle_name"]] \
        .merge(m, on="cell_iname", how="left").merge(per, left_on="cell_iname", right_index=True, how="left")
    cov["n_sig_genetic_other"] = cov.cell_iname.map(other)
    cov = cov.fillna({c: 0 for c in ["n_sig", "n_cpd", "n_sig_hiq", "n_cpd_hiq", "n_times", "n_sig_genetic_other"]})
    cov = cov.join(models[["CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]], on="ModelID")
    cov = cov.sort_values("n_cpd_hiq", ascending=False)
    cov.to_csv(RES / "cell_coverage.tsv", sep="\t", index=False)

    # compound x line (high-quality signatures), and per-compound summary
    cl = cp.groupby(["pert_id", "cell_iname"]).agg(n_sig=("sig_id", "size"), n_sig_hiq=("is_hiq", "sum"),
                                                   max_tas=("tas", "max")).reset_index()
    cl = cl.merge(m[["cell_iname", "ModelID"]], how="left")
    cl.to_csv(PROC / "compound_line.tsv.gz", sep="\t", index=False)
    hq = cl[cl.n_sig_hiq > 0]
    names = cp.drop_duplicates("pert_id").set_index("pert_id").cmap_name
    comp = pd.DataFrame({"cmap_name": names})
    comp["n_lines"] = cl.groupby("pert_id").cell_iname.nunique()
    comp["n_lines_hiq"] = hq.groupby("pert_id").cell_iname.nunique()
    comp["n_lines_hiq_depmap"] = hq.dropna(subset=["ModelID"]).groupby("pert_id").ModelID.nunique()
    comp["target"] = cpd.groupby("pert_id").target.apply(lambda x: ";".join(sorted(set(x.dropna()) - {""})))
    comp["moa"] = cpd.groupby("pert_id").moa.apply(lambda x: ";".join(sorted(set(x.dropna()) - {""})))
    comp = comp.fillna({"n_lines": 0, "n_lines_hiq": 0, "n_lines_hiq_depmap": 0}).sort_values("n_lines_hiq", ascending=False)
    comp.to_csv(RES / "compounds.tsv", sep="\t", index_label="pert_id")

    # disease coverage via Celligner weights
    w = pd.read_csv(ROOT / "results/celligner/disease_line_weights_long.tsv", sep="\t")
    prof = cov.dropna(subset=["ModelID"]).groupby("ModelID").n_cpd_hiq.max()
    w["lincs_n_cpd_hiq"] = w.ModelID.map(prof).fillna(0).astype(int)
    dist = pd.read_csv(ROOT / "data/processed/celligner/dist_line_line.tsv.gz", sep="\t", index_col=0)
    covered = [x for x in prof[prof >= 100].index if x in dist.columns]
    near = dist[covered].copy()
    for x in covered:
        near.loc[x, x] = np.inf  # a profiled line's own entry does not count as its neighbor
    w["nearest_profiled"] = w.ModelID.map(near.idxmin(axis=1)).map(models.CellLineName)
    w["dist_nearest_profiled"] = w.ModelID.map(near.min(axis=1)).round(2)
    w.to_csv(RES / "disease_line_coverage.tsv", sep="\t", index=False)

    rows = []
    cpd_rows = []
    for gid, g in w.groupby("group_id"):
        r = dict(group_id=gid, key=g.key.iloc[0])
        for t in TIERS:
            r[f"weight_on_lines_ge{t}cpd"] = round(g.loc[g.lincs_n_cpd_hiq >= t, "weight"].sum(), 3)
            r[f"n_lines_ge{t}cpd"] = int((g.lincs_n_cpd_hiq >= t).sum())
        top = g[g.lincs_n_cpd_hiq >= 100].sort_values("weight", ascending=False).head(3)
        r["top_profiled_lines"] = ", ".join(f"{a} {b:.3f} ({c} cpd)" for a, b, c in
                                            zip(top.CellLineName, top.weight, top.lincs_n_cpd_hiq))
        unp = g[g.lincs_n_cpd_hiq < 100]
        r["unprofiled_weighted_mean_dist"] = round(np.average(unp.dist_nearest_profiled, weights=unp.weight), 1) \
            if unp.weight.sum() > 0 else np.nan
        # per compound: weight of this disease's lines where it has a high-quality signature
        gw = g.set_index("ModelID").weight
        cw = hq.dropna(subset=["ModelID"]).assign(w=lambda d: d.ModelID.map(gw)).dropna(subset=["w"])
        cw = cw.groupby("pert_id").w.sum()
        r["n_cpd_any_profiled_line"] = int((cw > 0).sum())
        r["n_cpd_weight_ge_0.10"] = int((cw >= 0.10).sum())
        r["n_cpd_weight_ge_0.25"] = int((cw >= 0.25).sum())
        r["max_cpd_weight"] = round(cw.max(), 3) if len(cw) else 0.0
        rows.append(r)
        cpd_rows.append(cw.rename("disease_weight_covered").reset_index().assign(group_id=gid))
    dc = pd.DataFrame(rows).sort_values("weight_on_lines_ge100cpd")
    dc.to_csv(RES / "disease_coverage.tsv", sep="\t", index=False)
    pd.concat(cpd_rows).merge(comp[["cmap_name", "moa"]], left_on="pert_id", right_index=True, how="left") \
        .round(4).to_csv(PROC / "disease_compound_coverage.tsv.gz", sep="\t", index=False)

    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
    print(f"LINCS cell lines: {len(cells)}; with compound signatures: {(cov.n_sig > 0).sum()}; "
          f">= 100 hiq compounds: {(cov.n_cpd_hiq >= 100).sum()}; >= 1000: {(cov.n_cpd_hiq >= 1000).sum()}")
    print("mapping:", m.map_method.value_counts(dropna=False).to_dict())
    print(f"compounds: {cp.pert_id.nunique()} with signatures; {hq.pert_id.nunique()} with a high-quality one")
    print(cov[cov.n_cpd_hiq >= 100][["cell_iname", "CellLineName", "OncotreeLineage", "n_sig", "n_cpd", "n_cpd_hiq",
                                     "map_method"]].to_string(index=False))
    print(dc.drop(columns="top_profiled_lines").to_string(index=False))


if __name__ == "__main__":
    main()
