#!/usr/bin/env Rscript
# Diagnostics on patient-level regulator activity (input to Step 8); no pipeline outputs are changed.
#  1. Plate (TCGA batch): is the plate pattern already in raw expression, or only in the regulator scores?
#     For raw tumor logCPM (TMM, genes with median logCPM >= 1), ARACNe patient NES and CollecTRI patient NES:
#       pc_plate     R^2 of PCs 1-5 on plate vs a shuffled-plate null (1,000 shuffles, plate sizes kept;
#                    plates with many small levels inflate R^2, the null shows the chance level)
#       feat_plate   fraction of features (genes / regulators) with plate F-test p < 0.01 (chance = 0.01)
#  2. ESTIMATE: each regulator's ARACNe patient NES regressed on ESTIMATE score (residual + mean kept);
#     variance removed, most purity-driven regulators, PCA before/after (plate, ESTIMATE, ABSOLUTE purity).
#  3. Redundancy of the full ARACNe regulator set: modules of regulators whose patient NES correlate
#     |r| >= 0.7 (average-linkage on 1 - |r|), and regulators tracking a proliferation score
#     (mean z of E2F1, FOXM1, MYBL2 NES).
# Outputs: results/step6/aracne_project/diagnostics/{pc_plate,feat_plate,estimate_pca,estimate_regulators,
#          redundancy}.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(54)
GROUPS <- c("BRCA-Basal", "BRCA-Her2", "BRCA-LumA", "BRCA-LumB", "TCGA-COAD", "TCGA-READ")
N_SHUF <- 1000; N_PC <- 5; MOD_R <- 0.7; PROLIF <- c("E2F1", "FOXM1", "MYBL2")
out <- "results/step6/aracne_project/diagnostics"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
regs_meta <- fread("results/step5/aracne_regulators.tsv")
read_nes <- function(regulon, gid) {
  x <- fread(sprintf("data/processed/step6/%s/%s.patient_nes.tsv.gz", regulon, gid))
  m <- as.matrix(x[, -1]); rownames(m) <- x$regulator; t(m)  # patients x regulators
}
r2 <- function(y, f) { X <- model.matrix(~ f); 1 - sum(qr.resid(qr(X), y)^2) / sum((y - mean(y))^2) }
r2_cont <- function(y, v) { ok <- !is.na(v); if (sum(ok) < 10) NA_real_ else cor(y[ok], v[ok])^2 }
feat_p <- function(M, f) {  # F-test p per column of M (samples x features) on factor f
  X <- model.matrix(~ f); q <- qr(X); df1 <- q$rank - 1; df2 <- nrow(M) - q$rank
  rss1 <- colSums(qr.resid(q, M)^2); rss0 <- colSums(sweep(M, 2, colMeans(M))^2)
  pf(((rss0 - rss1) / df1) / (rss1 / df2), df1, df2, lower.tail = FALSE)
}
pca <- function(M) { p <- prcomp(M, center = TRUE, scale. = FALSE); list(x = p$x[, 1:N_PC], ve = p$sdev^2 / sum(p$sdev^2)) }
eff_dim <- function(ve) sum(ve)^2 / sum(ve^2)

run <- function(gid) {
  meta <- cov[group_id == gid]
  A <- read_nes("aracne_project", gid); C <- read_nes("collectri", gid)
  s <- intersect(intersect(rownames(A), meta$sample), colnames(counts))
  meta <- meta[match(s, sample)]; A <- A[s, ]; C <- C[s, ]
  plate <- factor(meta$batch)
  lc <- cpm(calcNormFactors(DGEList(counts[, s])), log = TRUE, prior.count = 1)
  lc <- t(lc[apply(lc, 1, median) >= 1, ])

  # 1. plate
  shuf <- replicate(N_SHUF, sample(plate), simplify = FALSE)
  mats <- list(raw_expression = lc, aracne_nes = A, collectri_nes = C)
  pc_plate <- rbindlist(lapply(names(mats), function(nm) {
    p <- pca(mats[[nm]])
    rbindlist(lapply(1:N_PC, function(k) {
      obs <- r2(p$x[, k], plate); nul <- vapply(shuf, function(f) r2(p$x[, k], f), numeric(1))
      data.table(group_id = gid, data = nm, pc = k, var_explained = round(p$ve[k], 3), r2_plate = round(obs, 3),
                 r2_chance_median = round(median(nul), 3), r2_chance_95 = round(quantile(nul, 0.95), 3),
                 p_shuffle = (sum(nul >= obs) + 1) / (N_SHUF + 1),
                 r2_estimate = round(r2_cont(p$x[, k], meta$estimate_score), 3),
                 r2_absolute = round(r2_cont(p$x[, k], meta$absolute_purity), 3))
    }))
  }))
  feat_plate <- rbindlist(lapply(names(mats), function(nm) {
    p <- feat_p(mats[[nm]], plate)
    data.table(group_id = gid, data = nm, n_features = length(p), frac_p01 = round(mean(p < 0.01), 3),
               frac_fdr05 = round(mean(p.adjust(p, "BH") < 0.05), 3))
  }))
  feat_plate[, `:=`(n_patients = length(s), n_plates = nlevels(plate))]

  # 2. ESTIMATE regression on ARACNe NES
  est <- meta$estimate_score; X <- cbind(1, est - mean(est))
  fit <- qr(X); Ares <- qr.resid(fit, A) + matrix(colMeans(A), nrow(A), ncol(A), byrow = TRUE)
  reg_r2 <- 1 - colSums((Ares - matrix(colMeans(A), nrow(A), ncol(A), byrow = TRUE))^2) / colSums(sweep(A, 2, colMeans(A))^2)
  est_reg <- data.table(group_id = gid, regulator = colnames(A), r2_estimate = round(reg_r2, 3),
                        slope_sign = sign(qr.coef(fit, A)[2, ]))
  est_pca <- rbindlist(lapply(c("before", "after"), function(w) {
    M <- if (w == "before") A else Ares; p <- prcomp(M, center = TRUE, scale. = FALSE); ve <- p$sdev^2 / sum(p$sdev^2)
    rbindlist(lapply(1:N_PC, function(k) data.table(group_id = gid, when = w, pc = k, var_explained = round(ve[k], 3),
      effective_dim = round(eff_dim(ve), 2), r2_plate = round(r2(p$x[, k], plate), 3),
      r2_estimate = round(r2_cont(p$x[, k], est), 3), r2_absolute = round(r2_cont(p$x[, k], meta$absolute_purity), 3),
      top_regulators = paste(names(sort(abs(p$rotation[, k]), decreasing = TRUE))[1:6], collapse = ","))))
  }))
  est_pca[, total_var_removed := round(1 - sum(apply(Ares, 2, var)) / sum(apply(A, 2, var)), 3)]

  # 3. redundancy
  R <- cor(A); hc <- hclust(as.dist(1 - abs(R)), method = "average"); mod <- cutree(hc, h = 1 - MOD_R)
  ms <- sort(table(mod), decreasing = TRUE)
  pro <- intersect(PROLIF, colnames(A))
  pscore <- rowMeans(scale(A[, pro, drop = FALSE]))
  rp <- cor(A, pscore)[, 1]
  pmod <- as.integer(names(sort(table(mod[pro]), decreasing = TRUE))[1])
  red <- data.table(group_id = gid, n_regulators = ncol(A), n_modules = length(ms), n_singletons = sum(ms == 1),
                    largest_module = as.integer(ms[1]), prolif_module_size = sum(mod == pmod),
                    n_track_prolif_07 = sum(abs(rp) >= MOD_R), frac_track_prolif_07 = round(mean(abs(rp) >= MOD_R), 3),
                    tf_in_prolif_module = sum(mod == pmod & colnames(A) %in% regs_meta[is_lambert_tf == TRUE, symbol]),
                    kinase_in_prolif_module = sum(mod == pmod & colnames(A) %in% regs_meta[is_kinase == TRUE, symbol]),
                    modules_ge10 = sum(ms >= 10), regs_in_modules_ge10 = sum(ms[ms >= 10]))
  list(pc_plate = pc_plate, feat_plate = feat_plate, est_reg = est_reg, est_pca = est_pca, red = red)
}
res <- mclapply(GROUPS, function(g) tryCatch(run(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                mc.cores = 6, mc.preschedule = FALSE)
pick <- function(k) rbindlist(lapply(res, `[[`, k))
fwrite(pick("pc_plate"), file.path(out, "pc_plate.tsv"), sep = "\t")
fwrite(pick("feat_plate"), file.path(out, "feat_plate.tsv"), sep = "\t")
fwrite(pick("est_pca"), file.path(out, "estimate_pca.tsv"), sep = "\t")
fwrite(pick("est_reg")[order(group_id, -r2_estimate)], file.path(out, "estimate_regulators.tsv"), sep = "\t")
fwrite(pick("red"), file.path(out, "redundancy.tsv"), sep = "\t")
message("DIAGNOSTICS DONE")
