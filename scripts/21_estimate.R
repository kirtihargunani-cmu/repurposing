#!/usr/bin/env Rscript
# Step 3: ESTIMATE stromal / immune scores for every Toil sample used in Step 1 groups.
# ssGSEA is rank-based within a sample, so log2(TPM+0.001) input is fine.
# Platform "illumina" -> scores only; ESTIMATE's purity formula is calibrated on Affymetrix,
# so tumor purity is not reported here (see ABSOLUTE cross-check in 20_covariates.py).
# Output: results/step3/estimate_toil.tsv
suppressPackageStartupMessages({library(data.table); library(estimate); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
n_cores <- 48
test <- nzchar(Sys.getenv("ESTIMATE_TEST"))  # quick check: 40 samples, 4 cores, no output file

samples <- unique(fread("results/step1/sample_groups_toil.tsv")$sample)
if (test) { samples <- head(samples, 40); n_cores <- 4 }
tpm_f <- "data/raw/solid_toil/TcgaTargetGtex_rsem_gene_tpm.gz"
hdr <- names(fread(cmd = paste("zcat", tpm_f, "| head -1")))
keep <- c(hdr[1], intersect(samples, hdr))
message(length(keep) - 1, " of ", length(samples), " samples found in TPM matrix")

x <- fread(cmd = paste("zcat", tpm_f), select = keep)
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
setnames(x, 1, "id")
x[, symbol := pm$gene[match(id, pm$id)]]
data(common_genes)
x <- x[symbol %in% common_genes$GeneSymbol]
x[, id := NULL]
x <- x[, lapply(.SD, mean), by = symbol]  # collapse duplicate symbols
message(nrow(x), " of ", nrow(common_genes), " ESTIMATE genes present")

tmp <- file.path(tempdir(), "estimate"); dir.create(tmp, showWarnings = FALSE)
cols <- setdiff(names(x), "symbol")
chunks <- split(cols, cut(seq_along(cols), n_cores, labels = FALSE))
res <- mclapply(seq_along(chunks), function(i) {
  gct <- file.path(tmp, sprintf("in_%02d.gct", i)); out <- file.path(tmp, sprintf("out_%02d.gct", i))
  m <- x[, c("symbol", chunks[[i]]), with = FALSE]
  writeLines(c("#1.2", paste(nrow(m), length(chunks[[i]]), sep = "\t")), gct)
  fwrite(cbind(m[, .(NAME = symbol, Description = symbol)], m[, -1]), gct, sep = "\t", append = TRUE,
         col.names = TRUE)
  capture.output(estimateScore(gct, out, platform = "illumina"))
  s <- fread(out, skip = 2, header = TRUE)
  setnames(s, c("NAME", "Description", chunks[[i]]))  # estimate mangles '-' in sample names
  s <- melt(s[, -2], id.vars = "NAME", variable.name = "sample")
  dcast(s, sample ~ NAME, value.var = "value")
}, mc.cores = n_cores)
res <- rbindlist(res)
setnames(res, c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
         c("estimate_stromal", "estimate_immune", "estimate_score"))
if (test) { print(head(res)); print(summary(res)); quit(save = "no") }
dir.create("results/step3", showWarnings = FALSE, recursive = TRUE)
fwrite(res, "results/step3/estimate_toil.tsv", sep = "\t")
message("ESTIMATE DONE: ", nrow(res), " samples")
