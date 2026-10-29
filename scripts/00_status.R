#!/usr/bin/env Rscript
# What is current for the active build? Run before any analysis.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); source("R/paths.R")
b <- build()
cat(sprintf("active build: %s\n  %s\n\n", b$name, b$description))
keys <- c("signatures", "regulons", "step6", "patient_nes", "step8", "step8_features")
for (k in keys) {
  p <- b[[k]]
  ex <- file.exists(p)
  extra <- ""
  if (ex && dir.exists(p)) extra <- sprintf("%d files", length(list.files(p)))
  if (ex && !dir.exists(p)) extra <- sprintf("%.1f MB", file.size(p) / 1e6)
  cat(sprintf("  %-16s %-5s %-50s %s\n", k, if (ex) "ok" else "MISSING", p, extra))
}
for (f in c(p_group_mr(b), p_regulator_filter(b))) {
  cat(sprintf("  %-16s %-5s %-50s %s\n", basename(f), if (file.exists(f)) "ok" else "MISSING", f,
              if (file.exists(f)) sprintf("%d rows", nrow(fread(f, select = 1L))) else ""))
}
g <- tryCatch(p_groups(b), error = function(e) character())
cat(sprintf("\n  groups with ok signature: %d\n", length(g)))
