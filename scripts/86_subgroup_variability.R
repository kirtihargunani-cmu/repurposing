#!/usr/bin/env Rscript
# Step 10 v2: variability of regulator activity within each subgroup (D72; user 2026-10-28: group and patient level, to
# assess intra-group variability). Patient NES = the active build's corrected patient VIPER (61; purity + quality
# corrected), restricted to the subgroup's patients; each parent group computed the same way for comparison.
# Per unit (subgroup or parent):
#   effective_dim     participation ratio of the PCA eigenvalues (sum(l)^2 / sum(l^2)), pc1_var
#   coherence         median Spearman of each patient's NES with the unit's mean profile (leave-one-out mean)
#   own_vs_sibling    subgroups: fraction of patients closer (Spearman) to their own subgroup's leave-one-out mean than to
#                     the sibling subgroup's mean
#   PCs 1-5 x covariates  R^2 with plate (and its chance level: mean over N_PERM permuted plate labels), ESTIMATE score,
#                     RNA-quality index, age, sex
#   top_var           regulators with the largest SD across the unit's patients (top 10)
#   total_var         summed variance of the regulators' NES across the unit's patients (overall spread; little affected by
#                     the number of patients). Added 2026-11-02 AFTER the first run: effective_dim rose in every subgroup,
#                     because splitting removes the parent's main axis and the rest spreads over more axes, so effective_dim
#                     measures how concentrated the variation is, not how much there is.
#   eff_dim_n30, coherence_n30   the same at a fixed size: mean over N_SUB random subsamples of N_FIX (30, the smallest
#                     subgroup) patients. effective_dim and coherence depend on the number of patients, so units of
#                     different size are compared on these (rule fixed 2026-11-02 before looking).
# Units: the 14 D79 subgroups (results/step10_v2/stratum_patients.tsv), their parent restricted to the subgroups'
#   patients (level "group"), and every one of the 41 cohorts with all its patients (level "cohort"; user 2026-11-02:
#   compare subgroup variability with the other disease cohorts).
# Comparison (fixed before looking): each subgroup's eff_dim_n30 and coherence_n30 as a percentile among the 41 cohorts
#   (lower dimension / higher coherence = more homogeneous); "tighter than its cohort" = coherence_n30 above and
#   eff_dim_n30 below the parent cohort's.
# Outputs: results/step10_v2/variability.tsv (one row per unit), variability_pcs.tsv (unit x PC x covariate),
#   variability_vs_cohorts.tsv (one row per subgroup: percentiles among cohorts, vs parent)
suppressPackageStartupMessages({library(data.table); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); source("R/paths.R")
b <- build(); require_build(b, c("patient_nes"))
N_PERM <- 50; N_FIX <- 30; N_SUB <- 20; RNGkind("L'Ecuyer-CMRG"); set.seed(86); out <- "results/step10_v2"  # L'Ecuyer: reproducible across mclapply workers
sp <- fread(file.path(out, "stratum_patients.tsv"))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
qtab <- fread(file.path(p_step6(b), "quality_correction/quality_index.tsv"))[role == "tumor"]
nes_of <- function(gid) { x <- fread(file.path(p_patient_nes(b), paste0(gid, ".patient_nes.tsv.gz")))
  m <- as.matrix(x[, -1]); rownames(m) <- x$regulator; m }
loo_rho <- function(A, B = NULL) {  # each patient vs leave-one-out mean of A (or mean of B)
  s <- rowSums(A); sapply(seq_len(ncol(A)), function(j) {
    ref <- if (is.null(B)) (s - A[, j]) / (ncol(A) - 1) else rowMeans(B)
    cor(A[, j], ref, method = "spearman") })
}
r2 <- function(y, v) { ok <- !is.na(v); if (sum(ok) < 10 || (is.factor(v) && nlevels(droplevels(v[ok])) < 2)) return(NA_real_)
  summary(lm(y[ok] ~ v[ok]))$r.squared }

unit <- function(name, parent, samples, sib = NULL, N, level) {
  A <- N[, intersect(samples, colnames(N)), drop = FALSE]
  pc <- prcomp(t(A), center = TRUE); lam <- pc$sdev^2; ve <- lam / sum(lam)
  meta <- cov[match(colnames(A), sample)]; qi <- qtab[group_id == parent][match(colnames(A), sample), quality_index]
  vars <- list(plate = factor(meta$batch), estimate = meta$estimate_score, quality = qi, age = meta$age, sex = factor(meta$sex))
  pcs <- rbindlist(lapply(seq_len(min(5, ncol(pc$x))), function(k) rbindlist(lapply(names(vars), function(v)
    data.table(unit = name, pc = k, var_explained = round(ve[k], 3), covariate = v, r2 = round(r2(pc$x[, k], vars[[v]]), 3),
               chance = if (v == "plate") round(mean(replicate(N_PERM, r2(pc$x[, k], sample(vars$plate)))), 3) else NA_real_)))))
  own <- loo_rho(A)
  sub <- if (ncol(A) >= N_FIX) sapply(seq_len(N_SUB), function(i) { S <- A[, sample(ncol(A), N_FIX)]
    l <- prcomp(t(S), center = TRUE)$sdev^2; c(sum(l)^2 / sum(l^2), median(loo_rho(S))) }) else matrix(NA_real_, 2, 1)
  sd_r <- sort(apply(A, 1, sd), decreasing = TRUE)
  row <- data.table(unit = name, parent = parent, level = level, n_patients = ncol(A),
    effective_dim = round(sum(lam)^2 / sum(lam^2), 2), pc1_var = round(ve[1], 3),
    coherence = round(median(own), 3), coherence_q10 = round(quantile(own, 0.1), 3),
    total_var = round(sum(apply(A, 1, var)), 1), eff_dim_n30 = round(mean(sub[1, ]), 2), coherence_n30 = round(mean(sub[2, ]), 3),
    own_vs_sibling = if (!is.null(sib)) { B <- N[, intersect(sib, colnames(N)), drop = FALSE]
      round(mean(own > loo_rho(A, B)), 3) } else NA_real_,
    pc1_plate_excess = pcs[pc == 1 & covariate == "plate", round(r2 - chance, 3)],
    pc1_purity_r2 = pcs[pc == 1 & covariate == "estimate", r2], pc1_quality_r2 = pcs[pc == 1 & covariate == "quality", r2],
    top_var = paste(sprintf("%s %.2f", names(sd_r)[1:10], sd_r[1:10]), collapse = ", "))
  list(row = row, pcs = pcs)
}

res <- mclapply(unique(sp$group_id), function(gid) {
  N <- nes_of(gid); st <- sp[group_id == gid]; subs <- unique(st$stratum)
  c(list(unit(gid, gid, st$sample, NULL, N, "group")),
    lapply(subs, function(s) unit(s, gid, st[stratum == s, sample], st[stratum != s, sample], N, "subgroup")))
}, mc.cores = 7)
cohorts <- sub("\\.patient_nes\\.tsv\\.gz$", "", list.files(p_patient_nes(b), "\\.patient_nes\\.tsv\\.gz$"))
resc <- mclapply(cohorts, function(gid) { N <- nes_of(gid)
  list(unit(gid, gid, colnames(N), NULL, N, "cohort")) }, mc.cores = 16)
res <- unlist(c(res, resc), recursive = FALSE)
rows <- rbindlist(lapply(res, `[[`, "row")); pcs <- rbindlist(lapply(res, `[[`, "pcs"))
fwrite(rows, file.path(out, "variability.tsv"), sep = "\t"); fwrite(pcs, file.path(out, "variability_pcs.tsv"), sep = "\t")
co <- rows[level == "cohort"]
cmp <- rows[level == "subgroup", .(subgroup = unit, parent, n_patients, eff_dim_n30, coherence_n30,
  eff_dim_pct_among_cohorts = round(100 * sapply(eff_dim_n30, function(v) mean(co$eff_dim_n30 < v)), 1),
  coherence_pct_among_cohorts = round(100 * sapply(coherence_n30, function(v) mean(co$coherence_n30 < v)), 1),
  total_var, total_var_pct_among_cohorts = round(100 * sapply(total_var, function(v) mean(co$total_var < v)), 1),
  total_var_ratio_to_parent = round(total_var / co$total_var[match(parent, co$unit)], 2),
  parent_eff_dim_n30 = co$eff_dim_n30[match(parent, co$unit)], parent_coherence_n30 = co$coherence_n30[match(parent, co$unit)],
  own_vs_sibling, pc1_plate_excess)]
cmp[, tighter_than_parent := coherence_n30 > parent_coherence_n30 & eff_dim_n30 < parent_eff_dim_n30]  # pre-set rule
cmp[, less_spread_than_parent := total_var_ratio_to_parent < 1]  # added after the first run (see total_var)
fwrite(cmp, file.path(out, "variability_vs_cohorts.tsv"), sep = "\t")
print(cmp)
print(summary(co[, .(n_patients, eff_dim_n30, coherence_n30, effective_dim, coherence)]))
options(width = 220); print(rows[, !"top_var"])
message("SUBGROUP VARIABILITY DONE")
