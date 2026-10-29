#!/usr/bin/env Rscript
# Step 8 subgroup consistency, per group, for the active build (R/paths.R; default v2_tf_xmap): which Step 8 clusters
# are stable and not technical, as input to the decision on which subgroups to keep (user 2026-10-28).
# Main clustering = modules10 at its best K (scripts/80). Per group:
#   significance     M3C decision for modules10 and all (80)
#   ari_all          ARI modules10 vs all at the same K (feature-set agreement)
#   ps_median        prediction strength (Tibshirani & Walther 2005) at best K: N_PS random halves; PAM on each half,
#                    test half assigned to the nearest training medoid; per test cluster, fraction of its pairs kept
#                    together; PS = min over clusters (median over halves). ps_k2 the same at K = 2.
#   plate_v_excess   Cramer's V clusters x plate minus its mean over N_PERM label permutations (many plates inflate V)
#   r2_purity / r2_quality / r2_module   R^2 of ESTIMATE score / RNA-quality index (61) / recurring-module score on
#                    cluster (module score = PC1 of the C8 core members, >= 30% of groups, present in the patient NES)
#   ari_subtype      ARI with PanCanAtlas Subtype_Selected (patients with a label; >= 2 labels)
#   ari_v1           ARI with the v1 modules10 clusters at v1's best K (the 6 groups clustered in v1)
#   concordance      median patient-vs-group MR rho (61), effective_dim (PCs for 80% variance, 61): intra-group spread
# STEP8_SUFFIX=_no_c8: score the clustering without the recurring module (80 EXCLUDE_MODULE=recurring); adds
#   ari_with_c8 = ARI with the original modules10 clusters at the same K
# Output: <step8><STEP8_SUFFIX>/consistency.tsv
suppressPackageStartupMessages({library(data.table); library(cluster); library(mclust); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); source("R/paths.R")
b <- build(); require_build(b, c("step8", "step8_features", "patient_nes"))
SUF <- Sys.getenv("STEP8_SUFFIX"); S8 <- paste0(p_step8(b), SUF); S8F <- paste0(p_step8_features(b), SUF)
orig <- if (nzchar(SUF)) fread(file.path(p_step8(b), "assignments.tsv"))[feature_set == "modules10"] else NULL
N_PS <- 20; N_PERM <- 200; N_CORES <- as.integer(Sys.getenv("N_CORES", "16")); set.seed(84)
dec <- fread(file.path(S8, "decision.tsv"))
asg <- fread(file.path(S8, "assignments.tsv"))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
qtab <- fread(file.path(p_step6(b), "quality_correction/quality_index.tsv"))[role == "tumor"]
conc <- fread(file.path(p_step6(b), "quality_correction/patient_concordance.tsv"))[startsWith(version, "v2")]
edim <- unique(fread(file.path(p_step6(b), "quality_correction/checks.tsv"))[startsWith(version, "v2"), .(group_id, effective_dim)])
core <- fread("results/step6/aracne_project_v2/diagnostics/recurring_module_members_v2.tsv")[frac_groups >= 0.3, regulator]
sub <- fread("data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz"); setnames(sub, 1, "sample")
v1 <- merge(fread("results/step8/assignments.tsv")[feature_set == "modules10"],
            fread("results/step8/decision.tsv")[feature_set == "modules10", .(group_id, K = best_K)], by = c("group_id", "K"))

pred_strength <- function(F, k) median(replicate(N_PS, {
  h <- sample(nrow(F), nrow(F) %/% 2); tr <- F[h, , drop = FALSE]; te <- F[-h, , drop = FALSE]
  med <- pam(tr, k)$medoids
  pred <- apply(te, 1, function(x) which.min(colSums((t(med) - x)^2)))
  cl <- pam(te, k, cluster.only = TRUE)
  min(sapply(unique(cl), function(j) { p <- pred[cl == j]; n <- length(p)
    if (n < 2) return(NA_real_); (sum(table(p)^2) - n) / (n * (n - 1)) }), na.rm = TRUE)
}))
cramer <- function(x, y) { t <- table(x, y); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
  if (min(dim(t)) < 2) return(NA_real_)
  sqrt(suppressWarnings(chisq.test(t, correct = FALSE)$statistic) / (sum(t) * (min(dim(t)) - 1))) }
r2 <- function(y, cl) { ok <- !is.na(y); if (sum(ok) < 10) return(NA_real_); summary(lm(y[ok] ~ factor(cl[ok])))$r.squared }

run <- function(gid) {
  dm <- dec[group_id == gid & feature_set == "modules10"]; k <- dm$best_K
  F <- fread(file.path(S8F, paste0(gid, "_modules10.tsv.gz"))); s <- F$sample; F <- as.matrix(F[, -1])
  cl <- asg[group_id == gid & feature_set == "modules10" & K == k][match(s, sample), cluster]
  ca <- asg[group_id == gid & feature_set == "all" & K == k][match(s, sample), cluster]
  meta <- cov[match(s, sample)]; qi <- qtab[group_id == gid][match(s, sample), quality_index]
  x <- fread(file.path(p_patient_nes(b), paste0(gid, ".patient_nes.tsv.gz")))
  mm <- x[regulator %in% core]; ms <- if (nrow(mm) >= 3) prcomp(t(as.matrix(mm[, s, with = FALSE])), scale. = TRUE)$x[, 1] else rep(NA_real_, length(s))
  v_obs <- cramer(cl, meta$batch); v_null <- mean(replicate(N_PERM, cramer(sample(cl), meta$batch)), na.rm = TRUE)
  st <- sub[match(s, sample), Subtype_Selected]; okst <- !is.na(st)
  w1 <- v1[group_id == gid][match(s, sample), cluster]
  data.table(group_id = gid, n_patients = length(s), n_modules = ncol(F), best_K = k,
    sig_modules10 = dm$decision == "clusters", sig_all = dec[group_id == gid & feature_set == "all", decision] == "clusters",
    RCSI = dm$RCSI, min_cluster = min(table(cl)), cluster_sizes = dm$cluster_sizes,
    ari_all = round(adjustedRandIndex(cl, ca), 3),
    ps_median = round(pred_strength(F, k), 3), ps_k2 = round(pred_strength(F, 2), 3),
    plate_v = round(v_obs, 3), plate_v_excess = round(v_obs - v_null, 3),
    r2_purity = round(r2(meta$estimate_score, cl), 3), r2_quality = round(r2(qi, cl), 3),
    module_members = nrow(mm), r2_module = round(r2(ms, cl), 3),
    n_subtype = sum(okst), n_subtype_labels = uniqueN(st[okst]),
    ari_subtype = if (uniqueN(st[okst]) >= 2) round(adjustedRandIndex(cl[okst], st[okst]), 3) else NA_real_,
    ari_v1 = if (any(!is.na(w1))) round(adjustedRandIndex(cl[!is.na(w1)], w1[!is.na(w1)]), 3) else NA_real_,
    ari_with_c8 = if (!is.null(orig)) round(adjustedRandIndex(cl, orig[group_id == gid & K == k][match(s, sample), cluster]), 3) else NA_real_,
    concordance = conc[group_id == gid, median_rho], effective_dim = edim[group_id == gid, effective_dim])
}
groups <- unique(dec$group_id)
out <- rbindlist(mclapply(groups, function(g) tryCatch(run(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                          mc.cores = N_CORES, mc.preschedule = FALSE))
fwrite(out[order(group_id)], file.path(S8, "consistency.tsv"), sep = "\t", na = "NA")
print(out[order(-sig_modules10, -ps_median), .(group_id, n_patients, best_K, sig_modules10, sig_all, min_cluster, ari_all,
          ps_median, ps_k2, plate_v_excess, r2_purity, r2_quality, r2_module, ari_subtype, ari_v1)], nrows = 100)
message("SUBGROUP CONSISTENCY DONE: ", nrow(out), " groups")
