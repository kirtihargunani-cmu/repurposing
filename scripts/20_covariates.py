#!/usr/bin/env python
"""Step 3: covariate tables.

Outputs (results/step3/):
  covariates_toil.tsv       one row per (sample, group) from results/step1/sample_groups_toil.tsv
  covariates_beataml.tsv    one row per BeatAML RNA-seq sample in the expression matrix
  genes_protein_coding_v23.tsv   protein-coding genes from the GENCODE v23 GTF (matches Toil)

primary_per_patient: within each (group, role, patient) exactly one sample is flagged, chosen by
SAMPLE_PRIORITY (diagnosis before metastasis before recurrence/post-treatment). GTEx donors that
contribute several tissues to one group keep the tissue listed first in that group's gtex_tissues.
Run scripts/21_estimate.R first; ESTIMATE scores are merged if results/step3/estimate_toil.tsv exists.
"""
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
OUT = ROOT / "results/step3"

# TCGA/TARGET sample-type code -> (class, priority); lower priority wins
SAMPLE_TYPES = {
    "01": ("primary", 1), "09": ("primary", 1), "03": ("primary", 2), "05": ("new_primary", 3),
    "06": ("metastatic", 4), "07": ("metastatic", 4),
    "02": ("recurrent", 5), "04": ("recurrent", 5), "40": ("recurrent", 5),
    "41": ("post_treatment", 6), "42": ("post_treatment", 6),
    "11": ("normal_adjacent", 1), "14": ("normal_adjacent", 2),
}


def code(sample):
    return sample.split("-")[3][:2]


def num(s):
    return pd.to_numeric(s, errors="coerce")


def target_clinical():
    frames = []
    for f in sorted((RAW / "target_clinical").glob("*ClinicalData*.xlsx")):
        x = pd.read_excel(f)
        if "Age at Diagnosis in Days" not in x:
            continue
        x = x.rename(columns={"TARGET USI": "patient"})
        x["age"] = num(x["Age at Diagnosis in Days"]) / 365.25
        x["sex_clin"] = x["Gender"].str.lower()
        x["blast_pct"] = num(x.get("Bone marrow leukemic blast percentage (%)"))
        x["tumor_pct"] = num(x.get("Percent Tumor"))
        x["protocol"] = x["Protocol"].astype(str)
        x["cell_of_origin"] = x.get("Cell of Origin")
        frames.append(x[["patient", "age", "sex_clin", "blast_pct", "tumor_pct", "protocol", "cell_of_origin"]])
    return pd.concat(frames).drop_duplicates("patient").set_index("patient")


def all_blasts():
    """TARGET-ALL blast %: mixed fractions, percents, '>70' and 'yes many' -> numeric percent or NaN."""
    b = pd.read_excel(RAW / "target_clinical/TARGET_ALL_Blast_Data_20230727.xlsx")
    raw = b["Blast %"].astype(str).str.strip().str.lower()
    val = num(raw.str.lstrip(">"))
    val = val.where(val > 1, val * 100)
    return pd.DataFrame({"patient": b["TARGET USI"].values, "blast_pct": val.values,
                         "blast_pct_raw": b["Blast %"].astype(str).values}).drop_duplicates("patient").set_index("patient")


def tcga_plates():
    f = sorted((RAW / "tcga_gdc").glob("rnaseq_aliquots_*.tsv"))[-1]
    a = pd.read_csv(f, sep="\t", dtype=str).iloc[:, 0].dropna()
    a = a[a.str.startswith("TCGA")].sort_values()
    parts = a.str.split("-")
    df = pd.DataFrame({"sample": a.str[:15].values, "plate": parts.str[5].values, "center": parts.str[6].values})
    n = df.groupby("sample").plate.nunique().rename("n_plates")
    return df.drop_duplicates("sample").set_index("sample").join(n)


def protein_coding():
    rows = []
    with gzip.open(RAW / "gencode/gencode.v23.annotation.gtf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene" or 'gene_type "protein_coding"' not in f[8]:
                continue
            gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            rows.append(dict(gene_id=gid.split(".")[0], gene_id_versioned=gid,
                             symbol=re.search(r'gene_name "([^"]+)"', f[8]).group(1), chrom=f[0]))
    return pd.DataFrame(rows)


def toil_covariates():
    sg = pd.read_csv(ROOT / "results/step1/sample_groups_toil.tsv", sep="\t", dtype=str)
    groups = pd.read_csv(ROOT / "results/step1/disease_groups_toil.tsv", sep="\t", dtype=str).fillna("")
    with gzip.open(RAW / "solid_toil/TcgaTargetGTEX_phenotype.txt.gz", "rt", encoding="latin-1") as fh:
        ph = pd.read_csv(fh, sep="\t", dtype=str)
    ph.columns = ["sample", "category", "disease", "site", "sample_type", "sex", "study"]
    ph = ph.set_index("sample")

    d = sg.copy()
    d["source"] = d["sample"].str.split("-").str[0].replace({"GTEX": "GTEx"})
    d["sex"] = d["sample"].map(ph.sex).str.lower().replace({"": np.nan})
    is_gtex = d.source == "GTEx"
    d["sample_type_code"] = np.where(is_gtex, "GTEx", d["sample"].map(lambda s: code(s) if not s.startswith("GTEX") else ""))
    d["sample_type_class"] = d.sample_type_code.map(lambda c: SAMPLE_TYPES.get(c, ("other", 9))[0])
    d.loc[is_gtex, "sample_type_class"] = "normal_gtex"
    d["_prio"] = d.sample_type_code.map(lambda c: SAMPLE_TYPES.get(c, ("other", 9))[1]).astype(int)

    # GTEx: donor = first two fields; priority = tissue order within the group's gtex_tissues
    d.loc[is_gtex, "patient"] = d.loc[is_gtex, "sample"].str.split("-").str[:2].str.join("-")
    tissue_order = {g.group_id: {t: i for i, t in enumerate(g.gtex_tissues.split(";"))} for g in groups.itertuples()}
    gt = d.loc[is_gtex]
    d.loc[is_gtex, "_prio"] = [tissue_order[g].get(t.removeprefix("GTEx "), 99) for g, t in zip(gt.group_id, gt.sample_type)]

    d = d.sort_values(["group_id", "role", "patient", "_prio", "sample"])
    d["primary_per_patient"] = ~d.duplicated(["group_id", "role", "patient"])
    chosen = d[d.primary_per_patient & (d.role == "tumor")]
    d["patient_has_untreated_tumor"] = pd.Series(pd.NA, index=d.index, dtype="boolean")
    untreated = chosen.set_index(["group_id", "patient"]).sample_type_class.isin(["primary", "new_primary", "metastatic"])
    tum = d.role == "tumor"
    d.loc[tum, "patient_has_untreated_tumor"] = [untreated.get((g, p)) for g, p in zip(d.loc[tum, "group_id"], d.loc[tum, "patient"])]

    # age
    clin = pd.read_csv(RAW / "tcga_pancan/Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t", dtype=str)
    tcga_age = num(clin.drop_duplicates("_PATIENT").set_index("_PATIENT").age_at_initial_pathologic_diagnosis)
    tc = target_clinical()
    gs = pd.read_csv(RAW / "gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
    gtex_age = gs.AGE.map(lambda a: np.mean([int(v) for v in a.split("-")]) + 0.5)
    d["age"] = np.select([d.source == "TCGA", d.source == "TARGET", is_gtex],
                         [d.patient.map(tcga_age), d.patient.map(tc.age), d.patient.map(gtex_age)], np.nan)
    d["age_is_binned"] = is_gtex
    d["sex"] = d.sex.fillna(d.patient.map(tc.sex_clin).where(d.source == "TARGET"))
    d["age"] = d.age.round(1)

    # site / batch / RNA quality
    plates = tcga_plates()
    sa = pd.read_csv(RAW / "gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t", dtype=str,
                     usecols=["SAMPID", "SMCENTER", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMNABTCH"])
    sa = sa.set_index("SAMPID")
    d["site"] = np.select([d.source == "TCGA", is_gtex, d.source == "TARGET"],
                          ["TSS-" + d["sample"].str.split("-").str[1], d["sample"].map(sa.SMCENTER),
                           d.patient.map(tc.protocol)], None)
    d["batch"] = np.select([d.source == "TCGA", is_gtex],
                           ["plate-" + d["sample"].map(plates.plate).fillna("NA"), d["sample"].map(sa.SMGEBTCH)], None)
    d["tcga_n_plates"] = d["sample"].map(plates.n_plates)
    d["gtex_rin"] = num(d["sample"].map(sa.SMRIN))
    d["gtex_ischemic_min"] = num(d["sample"].map(sa.SMTSISCH))

    # purity-related
    ab = pd.read_csv(RAW / "tcga_pancan/TCGA_mastercalls.abs_tables_JSedit.fixed.txt", sep="\t", dtype=str)
    ab = ab[ab["call status"] == "called"].drop_duplicates("array").set_index("array")
    d["absolute_purity"] = num(d["sample"].map(ab.purity))
    laml = pd.read_csv(RAW / "tcga_pancan/LAML_clinicalMatrix", sep="\t", dtype=str).set_index("sampleID")
    laml_blast = num(laml.lab_procedure_bone_marrow_blast_cell_outcome_percent_value)
    ab_all = all_blasts()
    tcont = pd.read_excel(RAW / "target_clinical/TARGET_AML_ClinicalData_Discovery_and_Validation_Tumor_Content_and_RIN_Supplement_20230720.xlsx")
    tcont["sample"] = tcont["TARGET Sample ID "].str[:19]
    tcont = tcont.drop_duplicates("sample").set_index("sample")
    diag = d.sample_type_class == "primary"  # clinical blasts/tumor % refer to the diagnostic sample
    d["blast_pct"] = np.select(
        [d.cohort == "TCGA-LAML", d.cohort.str.startswith("TARGET-AML") & diag, d.cohort == "TARGET-ALL"],
        [d["sample"].map(laml_blast), d.patient.map(tc.blast_pct), d.patient.map(ab_all.blast_pct).where(diag)], np.nan)
    d["blast_pct_raw"] = np.where((d.cohort == "TARGET-ALL") & diag, d.patient.map(ab_all.blast_pct_raw), None)
    d["tumor_content_pct"] = np.select(
        [d.cohort.str.startswith("TARGET-AML"), (d.cohort == "TARGET-NBL") & diag],
        [num(d["sample"].map(tcont["Tumor content (%)"])), d.patient.map(tc.tumor_pct)], np.nan)
    d["target_rin"] = num(d["sample"].map(tcont["RIN of  RNA"]))

    est = OUT / "estimate_toil.tsv"
    if est.exists():
        d = d.merge(pd.read_csv(est, sep="\t"), on="sample", how="left")
    else:
        print("WARNING: ESTIMATE scores not found; run scripts/21_estimate.R")

    # cross-check: TARGET-ALL lineage call vs clinical cell of origin
    xc = d[(d.cohort == "TARGET-ALL") & (d.role == "tumor")]
    print("TARGET-ALL group vs clinical cell of origin:\n",
          pd.crosstab(xc.group_id, xc.patient.map(tc.cell_of_origin).fillna("NA")).to_string())

    cols = ["sample", "patient", "group_id", "cohort", "role", "source", "sample_type", "sample_type_code",
            "sample_type_class", "primary_per_patient", "patient_has_untreated_tumor", "sex", "age",
            "age_is_binned", "site", "batch", "tcga_n_plates", "gtex_rin", "gtex_ischemic_min", "target_rin",
            "absolute_purity", "estimate_stromal", "estimate_immune", "estimate_score", "blast_pct",
            "blast_pct_raw", "tumor_content_pct"]
    return d[[c for c in cols if c in d]].sort_values(["group_id", "role", "sample"])


def beataml_covariates():
    c = pd.read_excel(RAW / "aml_beataml2/beataml_wv1to4_clinical.xlsx")
    expr_cols = pd.read_csv(RAW / "aml_beataml2/beataml_waves1to4_norm_exp_dbgap.txt", sep="\t", nrows=0).columns
    c = c[c.dbgap_rnaseq_sample.isin(expr_cols)].copy()
    drug = pd.read_csv(RAW / "aml_beataml2/beataml_probit_curve_fits_v4_dbgap.txt", sep="\t",
                       usecols=["dbgap_rnaseq_sample"]).dropna()
    stage_prio = {"Initial Diagnosis": 1, "Residual": 2, "Unknown": 3, "Relapse": 4, "Remission": 5}
    spec_prio = {"Bone Marrow Aspirate": 1, "Peripheral Blood": 2, "Leukapheresis": 3}
    c["_p"] = (c.diseaseStageAtSpecimenCollection.map(stage_prio).fillna(3) * 10
               + c.specimenType.map(spec_prio).fillna(3))
    c = c.sort_values(["dbgap_subject_id", "_p", "dbgap_rnaseq_sample"])
    c["primary_per_patient"] = ~c.duplicated("dbgap_subject_id")
    out = pd.DataFrame({
        "sample": c.dbgap_rnaseq_sample, "patient": c.dbgap_subject_id.astype(str),
        "group_id": "BeatAML", "role": "tumor", "source": "BeatAML",
        "specimen_type": c.specimenType, "disease_stage": c.diseaseStageAtSpecimenCollection,
        "primary_per_patient": c.primary_per_patient,
        "is_denovo": c.isDenovo, "is_relapse": c.isRelapse, "dx_at_acquisition": c.dxAtSpecimenAcquisition,
        "sex": c.consensus_sex.str.lower(), "age": num(c.ageAtSpecimenAcquisition),
        "site": "center-" + c.centerID.astype(str), "batch": c.cohort,
        "blast_pct_bm": num(c["%.Blasts.in.BM"]), "blast_pct_pb": num(c["%.Blasts.in.PB"]),
        "has_drug_response": c.dbgap_rnaseq_sample.isin(drug.dbgap_rnaseq_sample),
    })
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pc = protein_coding()
    pc.to_csv(OUT / "genes_protein_coding_v23.tsv", sep="\t", index=False)
    toil = toil_covariates()
    toil.to_csv(OUT / "covariates_toil.tsv", sep="\t", index=False)
    beat = beataml_covariates()
    beat.to_csv(OUT / "covariates_beataml.tsv", sep="\t", index=False)
    print(f"protein-coding genes: {len(pc)}")
    print(f"Toil rows: {len(toil)}; flagged primary_per_patient: {toil.primary_per_patient.sum()}")
    print(f"BeatAML RNA samples: {len(beat)}; subjects: {beat.patient.nunique()}")


if __name__ == "__main__":
    main()
