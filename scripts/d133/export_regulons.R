#!/usr/bin/env Rscript
# D133: export the regulons VIPER actually scored, restricted to the LINCS gene universe exactly as
# scripts/79d_drug_viper.R does it (names(tfmode) %in% rownames(P), P rows = drug_line_profiles_genes.tsv$gene_id).
# AUDIT ONLY: reads production inputs, writes nothing production. Output: one TSV per network + collectri.
# Usage: Rscript scripts/d133/export_regulons.R <out_dir> <net1,net2,...>
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; nets <- strsplit(args[2], ",")[[1]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
genes <- fread("data/processed/lincs/drug_line_profiles_genes.tsv")
univ <- genes$gene_id
cat(sprintf("LINCS gene universe: %d\n", length(univ)))

dump <- function(reglist, tag) {
  rows <- rbindlist(lapply(names(reglist), function(r) {
    x <- reglist[[r]]; k <- names(x$tfmode) %in% univ
    if (!any(k)) return(NULL)
    data.table(regulator = r, target = names(x$tfmode)[k], tfmode = as.numeric(x$tfmode)[k],
               likelihood = as.numeric(x$likelihood)[k])
  }))
  f <- file.path(out, paste0(tag, ".tsv.gz"))
  fwrite(rows, f, sep = "\t")
  cat(sprintf("%-12s %6d regulators, %8d regulator-target rows -> %s\n", tag, uniqueN(rows$regulator), nrow(rows), f))
}

arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
cat(sprintf("ARACNe rds: %d networks\n", length(arac)))
for (nt in nets) { stopifnot(nt %in% names(arac)); dump(arac[[nt]], paste0("aracne_", nt)) }
rm(arac); gc()
dump(readRDS("results/step5/regulons_core.rds")[["collectri"]], "collectri")
