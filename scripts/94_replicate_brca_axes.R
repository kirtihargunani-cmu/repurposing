#!/usr/bin/env Rscript
# Step 9d: do the 13 TCGA breast axes (82_breast_axes.R) replicate in METABRIC and CPTAC-BRCA? (user 2026-10-12)
# Yardstick as for colorectal (91): each patient scored on the TCGA subtype's modules10 features = mean of member NES,
#   each z-scored within its own cohort-subtype and sign-aligned as in TCGA; modules with < 50% of members present are
#   dropped; features z-scored within cohort-subtype. TCGA loadings are restricted to the modules usable in the cohort.
# Tests per axis x cohort (rule fixed before running):
#   1 emerges    PCA of the cohort's own features; congruence = max |cosine| between the TCGA loading vector and the
#                cohort's PCs 1-10. Null: TCGA loadings shuffled across modules (1,000x), same max -> permutation p.
#   2 captures   share of the cohort's feature variance along the TCGA direction (variance of the projected score / number
#                of features), vs the same null.
#   3 genes      projected axis score (features x TCGA loadings) per patient; per gene Spearman with the cohort's corrected
#                expression; Spearman across genes with the TCGA companion z (results/step10/z_axis_drivers.tsv).
#   verdict      replicated: congruence >= 0.85 ("fair similarity", Lorenzo-Seva & ten Berge 2006) and p < 0.05 and gene
#                rho >= 0.5; partial: p < 0.05 (test 1 or 2) and gene rho >= 0.3; otherwise not replicated.
#                Cohort-subtypes with < 30 patients (CPTAC Basal 27, Her2 14, LumA is 56, LumB 17) are flagged
#                underpowered: their PCA has at most n - 1 components and is unreliable; test 3 is the main evidence there.
# TCGA check: mean-member features vs Step 8 PC1 features; projected vs original axis score.
# Outputs: results/step9/brca_axes_{replication,features_check}.tsv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(94)
N_PERM <- 1000; N_PC <- 10; MIN_N <- 30; out <- "results/step9"
COHORTS <- c("metabric", "cptac_brca")
load_nes <- function(f) { x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator; A }
load_ex <- function(f) { x <- fread(f); E <- as.matrix(x[, -1]); rownames(E) <- x$gene_id; E }
mf <- fread("results/step8/module_features.tsv")[feature_set == "modules10"]
axes <- fread("results/step8/breast_axes_pca.tsv")[keep == TRUE]
L_all <- fread("results/step8/breast_axes_loadings.tsv"); S_all <- fread("results/step8/breast_axes_scores.tsv")
drv <- fread("results/step10/z_axis_drivers.tsv")
cos_ <- function(a, b) abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))

res <- list(); chk <- list()
for (gid in unique(axes$group_id)) {
  tc <- load_nes(sprintf("data/processed/step6/aracne_project_map/%s.patient_nes.tsv.gz", gid))
  m <- mf[group_id == gid]; feats <- unique(m$feature); Zt <- scale(tc)
  sgn <- rbindlist(lapply(feats, function(f) { r <- intersect(m[feature == f, regulator], colnames(Zt))
    data.table(feature = f, regulator = r, sign = sign(cor(Zt[, r, drop = FALSE], rowMeans(Zt[, r, drop = FALSE])))[, 1]) }))
  modscore <- function(A) { Z <- scale(A)
    F <- sapply(feats, function(f) { s <- sgn[feature == f]; r <- s$regulator[s$regulator %in% colnames(Z)]
      if (length(r) < 0.5 * nrow(s)) return(rep(NA_real_, nrow(Z)))
      rowMeans(sweep(Z[, r, drop = FALSE], 2, s[match(r, regulator), sign], `*`)) })
    F <- F[, colSums(is.na(F)) == 0, drop = FALSE]; rownames(F) <- rownames(A); scale(F) }
  Ftc <- modscore(tc)
  s8 <- fread(sprintf("data/processed/step8/features/%s_modules10.tsv.gz", gid)); P <- as.matrix(s8[, -1]); rownames(P) <- s8$sample
  Lg <- L_all[group_id == gid]; Sg <- S_all[group_id == gid]
  for (a in axes[group_id == gid, axis]) {
    l <- setNames(Lg[[a]], Lg$module)
    chk[[paste(gid, a)]] <- data.table(group_id = gid, axis = a, r_features_median = round(median(diag(cor(Ftc[rownames(P), colnames(P)], P))), 3),
      r_projected_vs_original = round(cor(Ftc[Sg$sample, names(l)] %*% l, Sg[[a]])[1], 3)) }
  for (co in COHORTS) {
    A <- load_nes(sprintf("data/processed/step9/%s/%s/patient_nes.tsv.gz", co, gid)); F <- modscore(A); n <- nrow(F)
    E <- load_ex(sprintf("data/processed/step9/%s/%s/expression.tsv.gz", co, gid))[, rownames(F)]
    pcs <- prcomp(F)$rotation; pcs <- pcs[, 1:min(N_PC, ncol(pcs)), drop = FALSE]
    for (a in axes[group_id == gid, axis]) {
      l <- setNames(Lg[[a]], Lg$module); common <- intersect(names(l), colnames(F)); l <- l[common]; l <- l / sqrt(sum(l^2))
      Fc <- F[, common]; pk <- pcs[common, , drop = FALSE]
      cong <- function(v) max(apply(pk, 2, cos_, b = v))
      share <- function(v) var(Fc %*% v)[1] / ncol(Fc)
      obs_c <- cong(l); obs_s <- share(l); best_pc <- which.max(apply(pk, 2, cos_, b = l))
      perm <- replicate(N_PERM, { v <- sample(l); c(cong(v), share(v)) })
      sc <- (Fc %*% l)[, 1]
      g <- intersect(rownames(E), drv$gene_id); tz <- drv[match(g, gene_id)][[paste0(gid, "-", a)]]
      gr <- cor(t(E[g, ]), sc, method = "spearman")[, 1]; okg <- !is.na(tz) & !is.na(gr)
      rho_g <- cor(gr[okg], tz[okg], method = "spearman")
      top <- order(-abs(tz))[1:100]; top <- top[okg[top]]
      p_c <- (1 + sum(perm[1, ] >= obs_c)) / (N_PERM + 1); p_s <- (1 + sum(perm[2, ] >= obs_s)) / (N_PERM + 1)
      verdict <- if (obs_c >= 0.85 && min(p_c, p_s) < 0.05 && rho_g >= 0.5) "replicated" else
                 if (min(p_c, p_s) < 0.05 && rho_g >= 0.3) "partial" else "not replicated"
      res[[paste(gid, a, co)]] <- data.table(group_id = gid, axis = a, cohort = co, n_patients = n,
        n_modules = length(common), n_modules_tcga = nrow(Lg), tcga_var_explained = axes[group_id == gid & axis == a, var_explained],
        congruence = round(obs_c, 3), matched_pc = best_pc, congruence_null95 = round(quantile(perm[1, ], 0.95), 3), congruence_p = signif(p_c, 3),
        variance_share = round(obs_s, 4), variance_share_null95 = round(quantile(perm[2, ], 0.95), 4), variance_share_p = signif(p_s, 3),
        gene_rho = round(rho_g, 3), top100_driver_same_direction = round(mean(sign(gr[top]) == sign(tz[top])), 2),
        underpowered = n < MIN_N, verdict = verdict) }
    message(sprintf("%s %s: %d patients x %d modules", gid, co, n, ncol(F)))
  }
}
rep <- rbindlist(res); fwrite(rep, file.path(out, "brca_axes_replication.tsv"), sep = "\t")
fwrite(rbindlist(chk), file.path(out, "brca_axes_features_check.tsv"), sep = "\t")
print(rep[, .(group_id, axis, cohort, n_patients, congruence, congruence_p, variance_share, variance_share_p, gene_rho, top100_driver_same_direction, verdict)])
print(rbindlist(chk))
message("BRCA AXES REPLICATION DONE")
