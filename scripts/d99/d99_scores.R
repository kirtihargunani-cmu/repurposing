#!/usr/bin/env Rscript
# D99 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): PROGENy pathway scores (decoupleR run_mlm, human top 100,
# cached data/raw/progeny/progeny_human_top100.tsv) for
#   drug profiles: production 79c line-specific (and raw), D88 arms (time_T6 / time_T24 / dose_DL, line-specific)
#   disease signatures: Step 4 v2 tumor-vs-normal z (41 groups)
# and DoRothEA A-B TF activity (decoupleR run_ulm, results/step5/regulons_dorothea_ABC.tsv) on the production profiles.
# Outputs data/processed/lincs/d99/{progeny_<set>.tsv.gz, dorothea_ab_specific.tsv.gz}, results/lincs/d99/progeny_disease.tsv
suppressPackageStartupMessages({library(data.table); library(decoupleR)})
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
od <- "data/processed/lincs/d99"; rd <- "results/lincs/d99"; dir.create(od, FALSE, TRUE); dir.create(rd, FALSE, TRUE)
net <- fread("data/raw/progeny/progeny_human_top100.tsv")
pg <- fread("data/processed/lincs/drug_line_profiles_genes.tsv"); idx <- fread("data/processed/lincs/drug_line_profiles_index.tsv")
keep <- !duplicated(pg$symbol)
score <- function(M, nw, fun = run_mlm, mor = "weight") {
  r <- fun(mat = M, network = nw, .source = "source", .target = "target", .mor = mor, minsize = 5)
  dcast(as.data.table(r)[statistic %in% c("mlm", "ulm")], condition ~ source, value.var = "score")
}
rdmat <- function(f, n, cols) { M <- matrix(readBin(f, "numeric", size = 4, n = nrow(pg) * n), nrow = nrow(pg))[keep, ]
  rownames(M) <- pg$symbol[keep]; colnames(M) <- cols; M }
cols <- paste(idx$pert_id, idx$ModelID, sep = "|")
for (s in c("specific", "raw")) {
  f <- file.path("data/processed/lincs", if (s == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32")
  x <- score(rdmat(f, nrow(idx), cols), net); fwrite(x, file.path(od, sprintf("progeny_%s.tsv.gz", s)), sep = "\t")
  message(s, ": ", nrow(x), " profiles")
}
for (a in list(c("time", "T6"), c("time", "T24"), c("dose", "DL"))) {
  I <- fread(sprintf("data/processed/lincs/timedose/%s_index.tsv", a[1]))
  x <- score(rdmat(sprintf("data/processed/lincs/timedose/%s_%s_specific.f32", a[1], a[2]), nrow(I), paste(I$pert_id, I$ModelID, sep = "|")), net)
  fwrite(x, file.path(od, sprintf("progeny_%s.tsv.gz", a[2])), sep = "\t"); message(a[2], ": ", nrow(x))
}
z <- fread("results/step4_v2/z_tumor_vs_gtex.tsv"); za <- fread("results/step4_v2/z_tumor_vs_adjacent.tsv")   # each group's own contrast (as 50_viper)
ct <- fread("results/step4_v2/signature_summary.tsv")[status == "ok", .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
za <- za[match(z$symbol, za$symbol)]
Z <- sapply(ct$group_id, function(g) if (ct[group_id == g, contrast] == "gtex") z[[g]] else za[[g]]); rownames(Z) <- z$symbol
Z <- Z[!duplicated(rownames(Z)) & !is.na(rownames(Z)), ]; Z[is.na(Z)] <- 0
fwrite(score(Z, net), file.path(rd, "progeny_disease.tsv"), sep = "\t")
do <- fread("results/step5/regulons_dorothea_ABC.tsv")[confidence %in% c("A", "B"), .(source, target, mor)]
x <- score(rdmat("data/processed/lincs/drug_line_profiles_bing_specific.f32", nrow(idx), cols), do, run_ulm, "mor")
fwrite(x, file.path(od, "dorothea_ab_specific.tsv.gz"), sep = "\t")
message("D99 SCORES DONE")
