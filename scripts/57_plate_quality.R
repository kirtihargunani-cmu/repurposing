#!/usr/bin/env Rscript
# Diagnostic: what is the TCGA "plate" signal? No pipeline outputs change.
# In TCGA, a plate is a processing batch at the Biospecimen Core Resource (BCR): samples are shipped from the tissue
# source site (TSS = collecting hospital), RNA is extracted and plated, then sequenced at one center (UNC). So plate
# mixes collection site, collection/processing time and RNA quality. For each group (purity-corrected ARACNe
# patient NES, PCs 1-3):
#   r2_plate, r2_tss              R^2 of the PC on plate / on TSS alone (chance level = (levels - 1) / (n - 1))
#   r2_plate_within_tss           extra R^2 of plate after TSS is in the model (plate effect inside one hospital)
#   r2_tss_within_plate           extra R^2 of TSS after plate
#   cramers_v_plate_tss           how strongly plates and hospitals overlap (1 = each plate from one hospital)
#   r2_length_index               R^2 of the PC on an RNA-quality proxy: per-sample mean z of long genes minus
#                                 short genes (top vs bottom genomic-span quintile, protein-coding, expressed).
#                                 Degraded RNA loses long-transcript signal, so this index tracks RNA integrity.
#   r2_length_index_plate         how much of the length index itself is explained by plate
# Output: results/step6/aracne_project/diagnostics/plate_quality.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
NES_DIR <- "data/processed/step6/aracne_project_purity"; N_PC <- 3
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))

r2f <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  X <- model.matrix(~ f); 1 - sum(qr.resid(qr(X), y)^2) / sum((y - mean(y))^2) }
r2_extra <- function(y, base, add) {  # R^2 gained by adding `add` to a model that already has `base`
  b <- droplevels(factor(base)); a <- droplevels(factor(add)); if (nlevels(a) < 2) return(NA_real_)
  X0 <- if (nlevels(b) > 1) model.matrix(~ b) else matrix(1, length(y)); X1 <- cbind(X0, model.matrix(~ a)[, -1, drop = FALSE])
  tss <- sum((y - mean(y))^2); (sum(qr.resid(qr(X0), y)^2) - sum(qr.resid(qr(X1), y)^2)) / tss
}
cramer <- function(a, b) { t <- table(a, b); chi <- suppressWarnings(chisq.test(t, correct = FALSE)$statistic)
  sqrt(chi / (sum(t) * (min(dim(t)) - 1))) }

run <- function(f) {
  gid <- sub("\\.patient_nes\\.tsv\\.gz$", "", basename(f))
  x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator
  meta <- cov[match(rownames(A), sample)]
  lc <- cpm(calcNormFactors(DGEList(counts[, rownames(A)])), log = TRUE, prior.count = 1)
  lc <- lc[apply(lc, 1, median) >= 1 & rownames(lc) %in% names(span), ]
  z <- (lc - rowMeans(lc)) / apply(lc, 1, sd)
  q <- quantile(span[rownames(z)], c(0.2, 0.8))
  lidx <- colMeans(z[span[rownames(z)] >= q[2], ]) - colMeans(z[span[rownames(z)] <= q[1], ])
  pc <- prcomp(A, center = TRUE); ve <- pc$sdev^2 / sum(pc$sdev^2)
  plate <- meta$batch; tss <- meta$site; n <- nrow(A)
  rbindlist(lapply(1:N_PC, function(k) { y <- pc$x[, k]
    data.table(group_id = gid, pc = k, var_explained = round(ve[k], 3), n = n,
      n_plates = uniqueN(plate), n_tss = uniqueN(tss),
      r2_plate = round(r2f(y, plate), 3), chance_plate = round((uniqueN(plate) - 1) / (n - 1), 3),
      r2_tss = round(r2f(y, tss), 3), chance_tss = round((uniqueN(tss) - 1) / (n - 1), 3),
      r2_plate_within_tss = round(r2_extra(y, tss, plate), 3), r2_tss_within_plate = round(r2_extra(y, plate, tss), 3),
      cramers_v_plate_tss = round(if (uniqueN(plate) > 1 && uniqueN(tss) > 1) cramer(plate, tss) else NA_real_, 2),
      r2_length_index = round(cor(y, lidx)^2, 3), r2_length_index_plate = round(r2f(lidx, plate), 3))
  }))
}
files <- Sys.glob(file.path(NES_DIR, "*.patient_nes.tsv.gz"))
res <- rbindlist(mclapply(files, function(f) tryCatch(run(f), error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL }),
                          mc.cores = 12, mc.preschedule = FALSE))
fwrite(res[order(group_id, pc)], "results/step6/aracne_project/diagnostics/plate_quality.tsv", sep = "\t")
message("PLATE QUALITY DONE: ", uniqueN(res$group_id), " groups")
