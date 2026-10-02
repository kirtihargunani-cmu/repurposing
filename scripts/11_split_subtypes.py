#!/usr/bin/env python
"""Step 1b: split mixed-histology Toil cohorts into disease groups and assign samples.

Inputs:  config/toil_disease_curation.tsv, config/subtype_splits.tsv,
         TCGA PanCanAtlas clinical + subtypes, Toil phenotype and TPM (TARGET-ALL lineage only).
Outputs: results/step1/disease_groups_toil.tsv   one row per disease group (split or whole cohort)
         results/step1/sample_groups_toil.tsv    sample -> group, role (tumor/adjacent_normal/gtex_normal)
         results/step1/target_all_lineage.tsv    B/T marker scores used to split TARGET-ALL
A group is used downstream only if it has >= MIN_PATIENTS tumor patients.
"""
import gzip
import importlib.util
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data/raw"
OUT = ROOT / "results/step1"
MIN_PATIENTS = 30
TUMOR_TYPES = re.compile(r"Tumor|Cancer|Metastatic|New Primary")

# B vs T lineage markers (Ensembl, unversioned); values are log2(TPM+0.001)
B_MARKERS = {"ENSG00000177455": "CD19", "ENSG00000196092": "PAX5", "ENSG00000105369": "CD79A"}
T_MARKERS = {"ENSG00000198851": "CD3E", "ENSG00000167286": "CD3D", "ENSG00000160654": "CD3G"}
B_MIN = 2.0  # max(CD19, PAX5) at or above this -> B lineage
T_MIN = 5.0  # otherwise mean CD3 at or above this -> T lineage

spec = importlib.util.spec_from_file_location("ground", ROOT / "scripts/10_ground_toil.py")
ground = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ground)


def patient(sample):
    return "-".join(sample.split("-")[:3])


def all_lineage(samples):
    """Score TARGET-ALL samples for B vs T lineage from the Toil TPM matrix."""
    tpm = RAW / "solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
    header = pd.read_csv(tpm, sep="\t", nrows=0).columns
    cols = [header[0]] + [s for s in samples if s in header]
    markers = {**B_MARKERS, **T_MARKERS}
    hits = []
    for chunk in pd.read_csv(tpm, sep="\t", usecols=cols, chunksize=5000):
        gene = chunk[header[0]].str.split(".").str[0]
        hits.append(chunk[gene.isin(markers)].assign(gene=gene[gene.isin(markers)].map(markers)))
    x = pd.concat(hits).drop(columns=header[0]).set_index("gene").T.copy()
    x["B_max"] = x[["CD19", "PAX5"]].max(axis=1)
    x["T_mean"] = x[["CD3E", "CD3D", "CD3G"]].mean(axis=1)
    x["lineage"] = "unassigned"
    x.loc[x.T_mean >= T_MIN, "lineage"] = "T"
    x.loc[x.B_max >= B_MIN, "lineage"] = "B"
    x.index.name = "sample"
    return x


def main():
    cur = pd.read_csv(ROOT / "config/toil_disease_curation.tsv", sep="\t", dtype=str).fillna("")
    splits = pd.read_csv(ROOT / "config/subtype_splits.tsv", sep="\t", dtype=str).fillna("")
    grounding = pd.read_csv(OUT / "disease_grounding_toil.tsv", sep="\t", dtype=str).fillna("")
    mondo = ground.parse_obo(RAW / "mondo/mondo.obo", "MONDO:")
    oncotree = {x["code"]: x for x in json.loads((RAW / "oncotree/tumorTypes.json").read_text())}
    with gzip.open(RAW / "solid_toil/TcgaTargetGTEX_phenotype.txt.gz", "rt", encoding="latin-1") as fh:
        pheno = pd.read_csv(fh, sep="\t", dtype=str).fillna("")
    pheno.columns = ["sample", "category", "disease", "site", "sample_type", "sex", "study"]
    gtex = pheno[(pheno.study == "GTEX") & (pheno.sample_type == "Normal Tissue")]

    # patient-level labels for splitting
    clin = pd.read_csv(RAW / "tcga_pancan/Survival_SupplementalTable_S1_20171025_xena_sp", sep="\t", dtype=str)
    histology = clin.drop_duplicates("_PATIENT").set_index("_PATIENT").histological_type
    sub = pd.read_csv(RAW / "tcga_pancan/TCGASubtype.20170308.tsv.gz", sep="\t", dtype=str)
    sub["patient"] = sub.sampleID.map(patient)
    # take the tumor sample's call (TCGA sample codes 01-09); the table also labels adjacent normals (11)
    sub = sub[sub.sampleID.str.split("-").str[3].str[:2].astype(int) < 10]
    subtype = (sub.dropna(subset=["Subtype_Selected"]).sort_values("sampleID")
               .drop_duplicates("patient").set_index("patient").Subtype_Selected)

    # samples per cohort
    samples = []
    for r in cur.itertuples(index=False):
        coh = pheno[(pheno.study == r.toil_study) & (pheno.category == r.toil_category)]
        for s in coh.itertuples(index=False):
            role = ("tumor" if TUMOR_TYPES.search(s.sample_type)
                    else "adjacent_normal" if s.sample_type == "Solid Tissue Normal" else None)
            if role:
                samples.append(dict(sample=s.sample, patient=patient(s.sample), cohort=r.cohort,
                                    role=role, sample_type=s.sample_type))
    samples = pd.DataFrame(samples)

    lineage = all_lineage(samples.loc[samples.cohort == "TARGET-ALL", "sample"].tolist())
    labels = {"histological_type": histology, "Subtype_Selected": subtype}

    def assign(row, rules):
        if rules.split_on.iloc[0] == "lineage":
            value = lineage.lineage.get(row["sample"], "")
        else:
            value = labels[rules.split_on.iloc[0]].get(row["patient"], "")
        value = "" if pd.isna(value) else str(value)
        for rule in rules.itertuples(index=False):
            if value and re.fullmatch(rule.match_regex, value):
                return rule.group_id
        return f"{row['cohort']}-unassigned"

    split_cohorts = set(splits.cohort)
    samples["group_id"] = [
        assign(row, splits[splits.cohort == row["cohort"]]) if row["cohort"] in split_cohorts else row["cohort"]
        for row in samples.to_dict("records")
    ]

    # group definitions: split groups from config, whole cohorts otherwise
    base = grounding.set_index("cohort")
    groups = []
    for r in splits.itertuples(index=False):
        groups.append(dict(group_id=r.group_id, cohort=r.cohort, group_label=r.group_label,
                           split_on=r.split_on, mondo_id=r.mondo_id, oncotree=r.oncotree,
                           gtex_tissues=r.gtex_tissues, gtex_match=r.gtex_match, gtex_sex=r.gtex_sex,
                           min_normals=r.min_normals, arm=base.loc[r.cohort, "arm"], notes=r.notes))
    for c in base.index.difference(split_cohorts):
        b = base.loc[c]
        groups.append(dict(group_id=c, cohort=c, group_label=b.mondo_label, split_on="",
                           mondo_id=b.mondo_id, oncotree=b.oncotree, gtex_tissues=b.gtex_tissues,
                           gtex_match=b.gtex_match, gtex_sex="", min_normals="", arm=b.arm, notes=b.notes))
    for c in sorted(split_cohorts):
        groups.append(dict(group_id=f"{c}-unassigned", cohort=c, group_label="no subtype label",
                           split_on="", mondo_id="", oncotree="", gtex_tissues="", gtex_match="",
                           gtex_sex="", min_normals="", arm=base.loc[c, "arm"], notes=""))
    groups = pd.DataFrame(groups).fillna("")

    # validate identifiers
    errors = []
    gtex_tissues = set(gtex.category)
    for g in groups.itertuples(index=False):
        if g.mondo_id and (g.mondo_id not in mondo or mondo[g.mondo_id]["obsolete"]):
            errors.append(f"{g.group_id}: MONDO {g.mondo_id} missing or obsolete")
        for c in filter(None, g.oncotree.split(";")):
            if c not in oncotree:
                errors.append(f"{g.group_id}: OncoTree {c} not in pinned release")
        for t in filter(None, g.gtex_tissues.split(";")):
            if t not in gtex_tissues:
                errors.append(f"{g.group_id}: GTEx tissue '{t}' not in Toil phenotype")
    if errors:
        print("VALIDATION FAILED:\n  " + "\n  ".join(errors), file=sys.stderr)
        sys.exit(1)

    # GTEx normals per group
    gtex_rows = []
    for g in groups.itertuples(index=False):
        for t in filter(None, g.gtex_tissues.split(";")):
            pool = gtex[gtex.category == t]
            if g.gtex_sex:  # e.g. female-only breast normals
                pool = pool[pool.sex.str.lower() == g.gtex_sex]
            for s in pool["sample"]:
                gtex_rows.append(dict(sample=s, patient=patient(s), cohort=g.cohort, role="gtex_normal",
                                      sample_type=f"GTEx {t}", group_id=g.group_id))
    sample_groups = pd.concat([samples, pd.DataFrame(gtex_rows)], ignore_index=True)

    # counts and inclusion
    tum = samples[samples.role == "tumor"]
    groups["mondo_label"] = groups.mondo_id.map(lambda m: mondo[m]["name"] if m else "")
    groups["n_tumor"] = groups.group_id.map(tum.group_id.value_counts()).fillna(0).astype(int)
    groups["n_tumor_patients"] = groups.group_id.map(tum.groupby("group_id").patient.nunique()).fillna(0).astype(int)
    adj = samples[samples.role == "adjacent_normal"]
    groups["n_adjacent_normal"] = groups.group_id.map(adj.group_id.value_counts()).fillna(0).astype(int)
    groups["n_gtex_normal"] = groups.group_id.map(
        sample_groups[sample_groups.role == "gtex_normal"].group_id.value_counts()).fillna(0).astype(int)
    reason = pd.Series("", index=groups.index)
    reason[groups.notes.str.startswith("EXCLUDE")] = "excluded by curation"
    reason[groups.group_id.str.endswith("-unassigned")] = "no subtype label"
    reason[(reason == "") & (groups.n_tumor_patients < MIN_PATIENTS)] = f"< {MIN_PATIENTS} tumor patients"
    groups["status"] = reason.map(lambda x: "exclude" if x else "include")
    groups["exclude_reason"] = reason

    cols = ["group_id", "cohort", "group_label", "split_on", "mondo_id", "mondo_label", "oncotree",
            "gtex_tissues", "gtex_match", "gtex_sex", "min_normals", "arm", "n_tumor", "n_tumor_patients", "n_adjacent_normal",
            "n_gtex_normal", "status", "exclude_reason", "notes"]
    groups = groups[cols].sort_values(["status", "cohort", "group_id"], ascending=[False, True, True])
    groups.to_csv(OUT / "disease_groups_toil.tsv", sep="\t", index=False)
    sample_groups.sort_values(["group_id", "role", "sample"]).to_csv(
        OUT / "sample_groups_toil.tsv", sep="\t", index=False)
    lineage.round(3).to_csv(OUT / "target_all_lineage.tsv", sep="\t")

    inc = groups[groups.status == "include"]
    print(f"{len(groups)} groups; {len(inc)} included "
          f"({inc.n_tumor.sum()} tumors, {inc.n_tumor_patients.sum()} patients)")


if __name__ == "__main__":
    main()
