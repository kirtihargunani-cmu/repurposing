#!/usr/bin/env Rscript
# Step 6c: purity correction of patient-level regulator activity (input to Step 8), all groups.
#   For each group, each regulator's patient NES is regressed on a purity covariate (linear, per group); the
#   residual plus the regulator's mean is kept, so only the part of the score that tracks purity is removed
#   (method reviewed in scripts/54_patient_nes_diagnostics.R; adopted 2026-10-07). Group-level msviper results
#   are unchanged (the Step 4 signatures are already ESTIMATE-adjusted).
#   Covariate: ESTIMATE score, except ABSOLUTE_GROUPS, which use ABSOLUTE purity (DNA-based). In LGG, ESTIMATE
#   agrees poorly with ABSOLUTE (rho -0.28 to -0.33) and ESTIMATE correction left ABSOLUTE-linked variation in
#   (PI decision 2026-10-08). Patients without the covariate are dropped from the corrected matrix (counted).
# REGULON=aracne_project (default) | collectri | dorothea_ABC
# Outputs: data/processed/step6/<REGULON>_purity/<group>.patient_nes.tsv.gz (corrected, same layout as input)
#          results/step6/<REGULON>/purity_correction/summary.tsv     per group: variance removed, purity PC
#            before/after (R^2 on ESTIMATE and on ABSOLUTE purity, which is DNA-based and independent of the
#            RNA), effective dimensionality, PC1 plate R^2, ESTIMATE-ABSOLUTE agreement
#          results/step6/<REGULON>/purity_correction/regulators.tsv  per group x regulator: R^2 on ESTIMATE, sign
suppressPackageStartupMessages({library(data.table); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
REGULON <- Sys.getenv("REGULON", "aracne_project"); N_PC <- 5
ABSOLUTE_GROUPS <- c("LGG-IDHmut-codel", "LGG-IDHmut-noncodel", "LGG-IDHwt")
in_dir <- file.path("data/processed/step6", REGULON); out_dir <- file.path("data/processed/step6", paste0(REGULON, "_purity"))
res_dir <- file.path("results/step6", REGULON, "purity_correction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]

r2f <- function(y, f) { X <- model.matrix(~ f); 1 - sum(qr.resid(qr(X), y)^2) / sum((y - mean(y))^2) }
r2c <- function(y, v) { ok <- !is.na(v); if (sum(ok) < 20) NA_real_ else cor(y[ok], v[ok])^2 }
pcs <- function(M) { p <- prcomp(M, center = TRUE); list(x = p$x[, 1:N_PC], ve = p$sdev^2 / sum(p$sdev^2)) }

run <- function(f) {
  gid <- sub("\\.patient_nes\\.tsv\\.gz$", "", basename(f))
  x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator
  meta <- cov[match(rownames(A), sample)]
  covariate <- if (gid %in% ABSOLUTE_GROUPS) "absolute_purity" else "estimate_score"
  ok <- !is.na(meta[[covariate]]); A <- A[ok, , drop = FALSE]; meta <- meta[ok]
  est <- meta[[covariate]]; mu <- colMeans(A)
  q <- qr(cbind(1, est - mean(est)))
  R <- qr.resid(q, A) + matrix(mu, nrow(A), ncol(A), byrow = TRUE)
  fwrite(data.table(regulator = colnames(R), round(t(R), 3)), file.path(out_dir, basename(f)), sep = "\t")

  ss0 <- colSums(sweep(A, 2, mu)^2); r2reg <- 1 - colSums(sweep(R, 2, mu)^2) / ss0
  regs <- data.table(group_id = gid, regulator = colnames(A), covariate = covariate, r2_purity = round(r2reg, 3),
                     direction = sign(qr.coef(q, A)[2, ]))
  b <- pcs(A); a <- pcs(R); plate <- factor(meta$batch)
  purity_pc <- which.max(apply(b$x, 2, r2c, v = est))
  summ <- data.table(group_id = gid, covariate = covariate, n_patients = nrow(A), n_dropped_no_covariate = sum(!ok),
    n_absolute = sum(!is.na(meta$absolute_purity)),
    rho_estimate_absolute = if (sum(!is.na(meta$absolute_purity)) >= 20)
      round(cor(meta$estimate_score, meta$absolute_purity, use = "complete", method = "spearman"), 2) else NA_real_,
    variance_removed = round(1 - sum(apply(R, 2, var)) / sum(apply(A, 2, var)), 3),
    frac_regulators_r2_gt_0.5 = round(mean(r2reg > 0.5), 3),
    top_purity_regulators = paste(head(colnames(A)[order(-r2reg)], 8), collapse = ","),
    purity_pc_before = purity_pc, purity_pc_var_before = round(b$ve[purity_pc], 3),
    purity_pc_r2_covariate_before = round(r2c(b$x[, purity_pc], est), 3),
    max_pc_r2_covariate_after = round(max(apply(a$x, 2, r2c, v = est)), 3),
    max_pc_r2_estimate_before = round(max(apply(b$x, 2, r2c, v = meta$estimate_score)), 3),
    max_pc_r2_estimate_after = round(max(apply(a$x, 2, r2c, v = meta$estimate_score)), 3),
    max_pc_r2_absolute_before = round(max(apply(b$x, 2, r2c, v = meta$absolute_purity)), 3),
    max_pc_r2_absolute_after = round(max(apply(a$x, 2, r2c, v = meta$absolute_purity)), 3),
    effective_dim_before = round(sum(b$ve)^2 / sum(b$ve^2), 2), effective_dim_after = round(sum(a$ve)^2 / sum(a$ve^2), 2),
    pc1_var_after = round(a$ve[1], 3),
    pc1_r2_plate_after = if (nlevels(plate) >= 2) round(r2f(a$x[, 1], plate), 3) else NA_real_,
    plate_r2_chance = if (nlevels(plate) >= 2) round((nlevels(plate) - 1) / (nrow(A) - 1), 3) else NA_real_)
  list(summ = summ, regs = regs)
}
files <- Sys.glob(file.path(in_dir, "*.patient_nes.tsv.gz"))
res <- mclapply(files, function(f) tryCatch(run(f), error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL }),
                mc.cores = 12, mc.preschedule = FALSE)
summ <- rbindlist(lapply(res, `[[`, "summ"))[order(-variance_removed)]
fwrite(summ, file.path(res_dir, "summary.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "regs"))[order(group_id, -r2_purity)], file.path(res_dir, "regulators.tsv"), sep = "\t")
message(sprintf("PURITY CORRECTION DONE (%s): %d of %d groups", REGULON, nrow(summ), length(files)))
