#!/usr/bin/env Rscript
# Step 10 v2: the patient subgroups that get their own signature, group / patient VIPER and drug matching (D72).
# Kept (D79, user 2026-10-28, option B): every group passing the D71 rule on the clustering WITHOUT the recurring module
#   (results/step8_v2_tf_no_c8; the module treated as a nuisance axis, D26): TCGA-BLCA, TCGA-HNSC, LGG-IDHmut-noncodel,
#   PCPG-PHEO, TCGA-THCA, TGCT-SEM, all K = 2, modules10; plus TCGA-COAD as an explicit exception (prediction strength
#   0.748 < 0.75, D71 correction) from results/step8_v2_tf: MSI-like vs CIN-like, the TCGA subtype split (ARI 0.42),
#   replicated in GSE39582 and CPTAC-COAD in v1. TCGA-ACC dropped (not significant without the module; was in D72).
# Names from the PanCanAtlas subtypes of the cluster members (Subtype_Selected / mRNA / DNA methylation; 2026-10-28).
# v1 strata (COAD CMS-like, READ C1-C3, breast axes) are not carried over: READ and BRCA-LumA clusters did not reproduce.
# Outputs: results/step10_v2/strata.tsv (one row per subgroup), results/step10_v2/stratum_patients.tsv (stratum, group_id, sample)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
out <- "results/step10_v2"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
NC8 <- "results/step8_v2_tf_no_c8"
lab <- data.table(
  group_id = rep(c("TCGA-BLCA", "TCGA-HNSC", "LGG-IDHmut-noncodel", "PCPG-PHEO", "TCGA-THCA", "TGCT-SEM", "TCGA-COAD"), each = 2),
  source   = c(rep(NC8, 12), rep("results/step8_v2_tf", 2)),
  cluster  = rep(1:2, 7),
  stratum  = c("BLCA-basal-squamous-like", "BLCA-luminal-like", "HNSC-basal-mesenchymal-like", "HNSC-atypical-classical-like",
               "LGG-noncodel-LGr2", "LGG-noncodel-LGr3", "PHEO-pseudohypoxia-Wnt", "PHEO-kinase-signaling",
               "THCA-BRAF-like-classical", "THCA-RAS-like-follicular", "TGCT-SEM-C1", "TGCT-SEM-C2",
               "COAD-CIN-like", "COAD-MSI-like"),
  description = c("all 31 TCGA mRNA cluster III (basal-squamous) tumors", "TCGA mRNA clusters I / II (luminal-papillary)",
                  "86 of 87 Basal, 67 of 75 Mesenchymal", "46 of 67 Atypical (HPV-enriched), Classical mixed",
                  "19 of 22 LGr2; = v2 LGG-noncodel-LGr2 of D72 (88 of 89)", "LGr3-dominant (141 of 208 LGr3)",
                  "all 36 Pseudohypoxia, 16 of 22 Wnt-altered", "54 of 64 Kinase signaling, Cortical admixture",
                  "THCA mRNA 3-5, classical methylation, 30 of 34 tall cell (BRAF-like)", "129 of 133 THCA mRNA 1, 110 of 114 follicular methylation (RAS-like)",
                  "no TCGA subtype labels", "no TCGA subtype labels",
                  "148 of 171 GI.CIN (exception, D79)", "45 of 48 GI.HM-indel (MSI); exception, D79"))
pat <- rbindlist(lapply(split(lab, by = "source"), function(l) {
  a <- fread(file.path(l$source[1], "assignments.tsv"))[feature_set == "modules10" & K == 2]
  a[l, on = .(group_id, cluster), nomatch = 0][, .(stratum, group_id, sample)]
}))
stopifnot(!anyDuplicated(pat$sample), uniqueN(pat$stratum) == nrow(lab))
lab[, n_patients := pat[, .N, by = stratum][match(lab$stratum, stratum), N]]
fwrite(lab, file.path(out, "strata.tsv"), sep = "\t")
fwrite(pat, file.path(out, "stratum_patients.tsv"), sep = "\t")
print(lab[, .(stratum, group_id, n_patients, source = basename(source))])
