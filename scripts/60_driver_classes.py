#!/usr/bin/env python
"""Step 7 support: how often each group's patients carry kinase-pathway vs TF-driver alterations.

Sources: MC3 PanCanAtlas MAF (FILTER == PASS) and GISTIC2 thresholded gene-level calls (Xena TCGA hub).
Patients: TCGA primary_per_patient tumors of each included group; denominators are patients present in
the MAF (mutation calls) or in GISTIC (copy number). TARGET groups are not covered by MC3/GISTIC.

An alteration counts only if it is likely functional:
  oncogene        hotspot missense / in-frame (listed below) or high-level amplification (GISTIC = 2)
  tumor suppressor truncating mutation (nonsense, frameshift, splice) or deep deletion (GISTIC = -2)
Classes
  kinase_pathway  RTK / RAS-MAPK / PI3K-AKT signaling activated by mutation or amplification
  tf_driver       transcription factor altered by mutation / amplification / loss
  epigenetic_tsg  non-TF tumor suppressors and chromatin regulators
Output: results/step7/driver_classes.tsv (group x class / gene percentages)
"""
import gzip
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw/tcga_pancan"
OUT = ROOT / "results/step7"

# gene -> (class, role, hotspot protein positions or None = any listed variant type)
DRIVERS = {
    # kinase / RTK / RAS-MAPK / PI3K
    "BRAF": ("kinase_pathway", "onc", {600, 594, 597, 601, 469, 464, 466}),
    "KRAS": ("kinase_pathway", "onc", {12, 13, 61, 146, 117}),
    "NRAS": ("kinase_pathway", "onc", {12, 13, 61}),
    "HRAS": ("kinase_pathway", "onc", {12, 13, 61}),
    "EGFR": ("kinase_pathway", "onc", {858, 861, 719, 790, 289, 598} | set(range(745, 753))),
    "ERBB2": ("kinase_pathway", "onc", {310, 755, 777, 842}),
    "PIK3CA": ("kinase_pathway", "onc", {542, 545, 546, 1047, 88, 345, 420, 1043}),
    "AKT1": ("kinase_pathway", "onc", {17}),
    "MET": ("kinase_pathway", "onc", {1010, 1250}),
    "FGFR1": ("kinase_pathway", "onc", set()), "FGFR2": ("kinase_pathway", "onc", {252, 253, 549}),
    "FGFR3": ("kinase_pathway", "onc", {248, 249, 372, 375, 650}),
    "KIT": ("kinase_pathway", "onc", set(range(550, 592)) | {816, 822}),
    "PDGFRA": ("kinase_pathway", "onc", {842}), "MAP2K1": ("kinase_pathway", "onc", {57, 99, 124, 203}),
    "CDK4": ("kinase_pathway", "onc", {24}), "CCND1": ("kinase_pathway", "onc", set()),
    "NF1": ("kinase_pathway", "tsg", None), "PTEN": ("kinase_pathway", "tsg", None),
    # transcription factors
    "TP53": ("tf_driver", "tsg_tp53", None), "MYC": ("tf_driver", "onc", set()),
    "MYCN": ("tf_driver", "onc", set()), "SOX2": ("tf_driver", "onc", set()), "TP63": ("tf_driver", "onc", set()),
    "NKX2-1": ("tf_driver", "onc", set()), "AR": ("tf_driver", "onc", {878, 875, 702}),
    "FOXA1": ("tf_driver", "onc", set(range(247, 272))), "GATA3": ("tf_driver", "tsg", None),
    "CTNNB1": ("tf_driver", "onc", set(range(32, 46))), "NFE2L2": ("tf_driver", "onc", set(range(24, 35)) | set(range(75, 83))),
    "SOX9": ("tf_driver", "tsg", None), "ESR1": ("tf_driver", "onc", {536, 537, 538}),
    # epigenetic / non-TF tumor suppressors
    "VHL": ("epigenetic_tsg", "tsg", None), "PBRM1": ("epigenetic_tsg", "tsg", None),
    "SETD2": ("epigenetic_tsg", "tsg", None), "BAP1": ("epigenetic_tsg", "tsg", None),
    "KEAP1": ("epigenetic_tsg", "tsg", None), "CDKN2A": ("epigenetic_tsg", "tsg", None),
    "ARID1A": ("epigenetic_tsg", "tsg", None), "KMT2D": ("epigenetic_tsg", "tsg", None),
    "APC": ("epigenetic_tsg", "tsg", None), "SMAD4": ("epigenetic_tsg", "tsg", None),
    "RB1": ("epigenetic_tsg", "tsg", None), "STK11": ("epigenetic_tsg", "tsg", None),
}
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation"}
HOT = {"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(ROOT / "results/step3/covariates_toil.tsv", sep="\t", low_memory=False)
    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t")
    inc = set(groups.loc[groups.status == "include", "group_id"])
    pts = cov[(cov.role == "tumor") & cov.primary_per_patient & cov.group_id.isin(inc) & (cov.source == "TCGA")]
    pts = pts[["group_id", "patient", "sample"]].drop_duplicates()

    maf = pd.read_csv(RAW / "mc3.v0.2.8.PUBLIC.maf.gz", sep="\t", usecols=[0, 8, 15, 36, 53, 108],
                      dtype=str, low_memory=False)
    maf.columns = ["gene", "vclass", "barcode", "hgvsp", "pos", "filter"]
    profiled_mut = set(maf.barcode.str[:12])
    maf = maf[(maf["filter"] == "PASS") & maf.gene.isin(DRIVERS)]
    maf["patient"] = maf.barcode.str[:12]
    maf["p1"] = pd.to_numeric(maf.pos.str.split("/").str[0].str.split("-").str[0], errors="coerce")

    def functional(r):
        cls, role, hot = DRIVERS[r.gene]
        if role.startswith("tsg"):
            return r.vclass in TRUNC or (role == "tsg_tp53" and r.vclass in HOT)  # TP53 missense are drivers
        return r.vclass in HOT and hot is not None and r.p1 in hot
    mut = maf[maf.apply(functional, axis=1)][["patient", "gene"]].drop_duplicates()

    with gzip.open(RAW / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz", "rt") as fh:
        cn = pd.read_csv(fh, sep="\t", index_col=0)
    cn = cn.loc[cn.index.intersection(list(DRIVERS))]
    cn.columns = cn.columns.str[:12]
    cn = cn.loc[:, ~cn.columns.duplicated()]
    profiled_cn = set(cn.columns)
    amp = [(p, g) for g in cn.index for p in cn.columns[cn.loc[g] == 2] if DRIVERS[g][1] == "onc"]
    dele = [(p, g) for g in cn.index for p in cn.columns[cn.loc[g] == -2] if DRIVERS[g][1].startswith("tsg")]
    cna = pd.DataFrame(amp + dele, columns=["patient", "gene"])

    alt = pd.concat([mut, cna]).drop_duplicates()
    alt["cls"] = alt.gene.map(lambda g: DRIVERS[g][0])
    rows = []
    for gid, d in pts.groupby("group_id"):
        p = set(d.patient) & profiled_mut & profiled_cn
        if len(p) < 20:
            continue
        a = alt[alt.patient.isin(p)]
        row = {"group_id": gid, "n_profiled": len(p)}
        for c in ["kinase_pathway", "tf_driver", "epigenetic_tsg"]:
            row[f"pct_{c}"] = round(100 * a.loc[a.cls == c, "patient"].nunique() / len(p), 1)
        top = (a.groupby("gene").patient.nunique() / len(p) * 100).sort_values(ascending=False)
        row["top_drivers"] = ", ".join(f"{g} {v:.0f}%" for g, v in top.head(6).items())
        rows.append(row)
    out = pd.DataFrame(rows)
    mr = pd.read_csv(ROOT / "results/step6/collectri/fdr10_group_summary.tsv", sep="\t")[["group_id", "pass"]]
    out = out.merge(mr, on="group_id", how="left").rename(columns={"pass": "collectri_fdr10_regulators"})
    out.sort_values("collectri_fdr10_regulators").to_csv(OUT / "driver_classes.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 70)
    print(out.sort_values("collectri_fdr10_regulators").to_string(index=False))
    zero = out[out.collectri_fdr10_regulators == 0]; some = out[out.collectri_fdr10_regulators > 0]
    for c in ["pct_kinase_pathway", "pct_tf_driver", "pct_epigenetic_tsg"]:
        print(f"{c}: median {zero[c].median():.0f}% in groups with 0 regulators vs {some[c].median():.0f}% in groups with >= 1")


if __name__ == "__main__":
    main()
