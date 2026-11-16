#!/usr/bin/env Rscript
# D128-A3.2 (reference-activity flag): per-GTEx-tissue mean expression, for the reference-tissue specificity
# diagnostic tau. Read-only with respect to every production input; writes one small table under the D128 sandbox.
#
#   tau(gene, cohort) = mean log2-CPM in the cohort's GTEx reference tissue - median over all GTEx tissues
#
# Samples: exactly the GTEx samples Step 4a kept (results/step4_v2/gtex_reference.tsv, use == TRUE), so the tissue
# means describe the same reference the disease contrast was computed against. Normalisation: TMM over those samples,
# log2-CPM with prior count 2, as elsewhere in the project (scripts/32_gtex_reference.R).
#
# Output: results/lincs/d128/a3/gtex_tissue_logcpm.tsv.gz  (gene_id, <one column per tissue>)
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../.."))
setwd(root)
OUT <- "results/lincs/d128/a3"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

g <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
s <- intersect(g$sample, colnames(counts))
cat(sprintf("GTEx samples kept by step 4a: %d; present in the count matrix: %d; tissues: %d\n",
            nrow(g), length(s), uniqueN(g[sample %in% s]$tissue)))
g <- g[sample %in% s]
y <- calcNormFactors(DGEList(counts[, g$sample]))
e <- cpm(y, log = TRUE, prior.count = 2)
rm(counts, y); invisible(gc())

tis <- sort(unique(g$tissue))
M <- sapply(tis, function(t) rowMeans(e[, g[tissue == t]$sample, drop = FALSE]))
D <- data.table(gene_id = sub("\\..*", "", rownames(e)))
D <- cbind(D, as.data.table(round(M, 4)))
fwrite(D, file.path(OUT, "gtex_tissue_logcpm.tsv.gz"), sep = "\t")
cat(sprintf("wrote %s: %d genes x %d tissues\n", file.path(OUT, "gtex_tissue_logcpm.tsv.gz"), nrow(D), length(tis)))
