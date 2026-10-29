#!/usr/bin/env Rscript
# Step 10 v2: the patient subgroups that get their own signature, group / patient VIPER and drug matching (D72).
# Kept (D71 rule, user 2026-10-28): the five groups passing 84_subgroup_consistency.R on the TF arm, K = 2, modules10.
#   TCGA-BLCA, TCGA-HNSC   clusters from the run without the recurring module (results/step8_v2_tf_no_c8; same splits
#                          as with it, ARI 0.87 / 0.91), so the subgroups do not rest on the C8 module (D26 open)
#   TCGA-COAD, TCGA-ACC, LGG-IDHmut-noncodel   results/step8_v2_tf
# Names from the PanCanAtlas subtypes of the cluster members (Subtype_Selected / Subtype_mRNA; checked 2026-10-28).
# v1 strata (COAD CMS-like, READ C1-C3, breast axes) are not carried over: READ and BRCA-LumA clusters did not reproduce.
# Outputs: results/step10_v2/strata.tsv (one row per subgroup), results/step10_v2/stratum_patients.tsv (stratum, group_id, sample)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
out <- "results/step10_v2"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
lab <- data.table(
  group_id = rep(c("TCGA-BLCA", "TCGA-HNSC", "TCGA-COAD", "TCGA-ACC", "LGG-IDHmut-noncodel"), each = 2),
  source   = rep(c("results/step8_v2_tf_no_c8", "results/step8_v2_tf_no_c8", "results/step8_v2_tf", "results/step8_v2_tf",
                   "results/step8_v2_tf"), each = 2),
  cluster  = rep(1:2, 5),
  stratum  = c("BLCA-basal-squamous-like", "BLCA-luminal-like", "HNSC-basal-mesenchymal-like", "HNSC-atypical-classical-like",
               "COAD-CIN-like", "COAD-MSI-like", "ACC-proliferative", "ACC-low-proliferation",
               "LGG-noncodel-LGr3", "LGG-noncodel-LGr2"),
  description = c("all 31 TCGA mRNA cluster III (basal-squamous) tumors", "TCGA mRNA clusters I / II (luminal-papillary)",
                  "86 of 87 Basal, 67 of 75 Mesenchymal", "46 of 67 Atypical (HPV-enriched), Classical mixed",
                  "148 of 171 GI.CIN", "45 of 48 GI.HM-indel (MSI); as v1 COAD-CMS1/3-like",
                  "all 25 steroid-phenotype + proliferation (+ 4 of 4 low + proliferation)", "steroid-phenotype high / low without proliferation",
                  "LGr3-dominant (140 of 208 LGr3)", "21 of 22 LGr2"))
pat <- rbindlist(lapply(split(lab, by = "source"), function(l) {
  a <- fread(file.path(l$source[1], "assignments.tsv"))[feature_set == "modules10" & K == 2]
  a[l, on = .(group_id, cluster), nomatch = 0][, .(stratum, group_id, sample)]
}))
stopifnot(!anyDuplicated(pat$sample), uniqueN(pat$stratum) == nrow(lab))
lab[, n_patients := pat[, .N, by = stratum][match(lab$stratum, stratum), N]]
fwrite(lab, file.path(out, "strata.tsv"), sep = "\t")
fwrite(pat, file.path(out, "stratum_patients.tsv"), sep = "\t")
print(lab[, .(stratum, group_id, n_patients, source = basename(source))])
