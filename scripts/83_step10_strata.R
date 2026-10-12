#!/usr/bin/env Rscript
# Step 10 input: the patient strata that get their own disease signature and drug matching (user 2026-10-12).
#   TCGA-COAD  two groups from Step 8 (modules10, K = 2), named by CMS (Step 9b): cluster 2 = "CMS1/3-like"
#              (MSI-like; CMS1 + CMS3), cluster 1 = "CMS2/4-like" (CMS2 + CMS4). Replicated in GSE39582 and CPTAC-COAD.
#   TCGA-READ  three Step 8 clusters (modules10, K = 3) as separate subgroups, each matched to drugs independently;
#              tentative (plate-associated, not replicable in the colon-only cohorts).
#   BRCA-*     no clusters; continuous axes from 82_breast_axes.R (kept PCs; patient scores in breast_axes_scores.tsv).
# Outputs: results/step8/step10_strata.tsv (one row per stratum or axis), step10_patient_strata.tsv (COAD / READ patients)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); res_dir <- "results/step8"
asg <- fread(file.path(res_dir, "assignments.tsv"))[feature_set == "modules10"]
lab <- data.table(group_id = c("TCGA-COAD", "TCGA-COAD", "TCGA-READ", "TCGA-READ", "TCGA-READ"), K = c(2, 2, 3, 3, 3),
  cluster = c(2L, 1L, 1L, 2L, 3L),
  stratum = c("COAD-CMS1/3-like", "COAD-CMS2/4-like", "READ-C1", "READ-C2", "READ-C3"),
  description = c("MSI-like; CMS1 (immune/MSI) + CMS3 (metabolic); holds 43 of 48 MSI tumors",
                  "CMS2 (canonical WNT/MYC) + CMS4 (mesenchymal)",
                  "high ILF3/PRMT1/E2F2 proliferation-type modules; CMS mixed", "mostly CMS4 / CMS2", "CMS3-enriched"),
  status = c("replicated (GSE39582, CPTAC-COAD)", "replicated (GSE39582, CPTAC-COAD)", rep("tentative (plate-associated; no rectal replication cohort)", 3)))
pat <- asg[lab, on = .(group_id, K, cluster), nomatch = 0][, .(group_id, stratum, sample)]
stopifnot(!anyDuplicated(pat$sample))
lab[, n_patients := pat[, .N, by = stratum][match(lab$stratum, stratum), N]]
ax <- fread(file.path(res_dir, "breast_axes_pca.tsv"))[keep == TRUE]
an <- fread(file.path(res_dir, "breast_axes_annotation.tsv"))
ax <- an[ax, on = .(group_id, axis)]
br <- ax[, .(group_id, kind = "continuous axis", stratum = paste0(group_id, "-", axis), n_patients,
             description = sprintf("%.0f%% of module variance; %s", 100 * var_explained, sub(";.*", "", top_hallmarks)),
             status = "TCGA only (not yet replicated)")]
strata <- rbind(lab[, .(group_id, kind = "subgroup", stratum, n_patients, description, status)], br)
fwrite(strata, file.path(res_dir, "step10_strata.tsv"), sep = "\t")
fwrite(pat, file.path(res_dir, "step10_patient_strata.tsv"), sep = "\t")
print(strata[, .(stratum, n_patients, description)])
