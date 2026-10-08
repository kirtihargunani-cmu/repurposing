#!/usr/bin/env Rscript
# Diagnostic: patient heterogeneity in BRCA (ARACNe patient NES, 4 PAM50 subtypes pooled; all scored against
# the same GTEx breast normals with the same BRCA network, so scores are comparable). No pipeline outputs change.
# Three views of the same scores (the corrected views are what-ifs, not a chosen design):
#   as_is            patient NES
#   estimate         ESTIMATE regressed out (subtype kept in the model so subtype differences are not removed)
#   estimate_plate   ESTIMATE + plate regressed out (subtype kept)
# Per view:
#   variance.tsv     share of each regulator's variance explained by subtype / plate / ESTIMATE (each alone),
#                    summarized over regulators
#   pca.tsv          pooled PCs 1-6: variance explained, R^2 on subtype / plate / ESTIMATE / ABSOLUTE
#   separation.tsv   PAM50 silhouette in the top-10-PC space; within-subtype effective dimensionality
#   kmeans.tsv       cluster tendency inside each subtype: k-means (k = 2..8, top 10 PCs) mean silhouette
#                    vs the same on column-shuffled data (no structure)
# Output: results/step6/aracne_project/diagnostics/brca_heterogeneity_*.tsv
suppressPackageStartupMessages({library(data.table); library(cluster)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(55)
SUB <- c("BRCA-Basal", "BRCA-Her2", "BRCA-LumA", "BRCA-LumB"); N_PC <- 10; KS <- 2:8
out <- "results/step6/aracne_project/diagnostics"

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
nes <- lapply(SUB, function(g) { x <- fread(sprintf("data/processed/step6/aracne_project/%s.patient_nes.tsv.gz", g))
  m <- t(as.matrix(x[, -1])); colnames(m) <- x$regulator; m })
regs <- Reduce(intersect, lapply(nes, colnames))
A <- do.call(rbind, lapply(nes, function(m) m[, regs]))
meta <- cov[match(rownames(A), sample)]
sub <- factor(meta$group_id); plate <- factor(meta$batch); est <- meta$estimate_score
message(sprintf("%d patients x %d regulators (shared by all 4 subtypes)", nrow(A), ncol(A)))

remove <- function(M, nuis) {  # regress out nuisance terms, keeping subtype; add column means back
  X <- model.matrix(as.formula(paste("~ sub +", nuis)), data.frame(sub, plate, est))
  q <- qr(X); B <- qr.coef(q, M); B[is.na(B)] <- 0
  keep <- grep("^\\(Intercept\\)|^sub", rownames(B))
  X[, keep, drop = FALSE] %*% B[keep, , drop = FALSE] + qr.resid(q, M)
}
views <- list(as_is = A, estimate = remove(A, "est"), estimate_plate = remove(A, "est + plate"))
r2f <- function(y, f) { X <- model.matrix(~ f); 1 - sum(qr.resid(qr(X), y)^2) / sum((y - mean(y))^2) }
r2c <- function(y, v) { ok <- !is.na(v); cor(y[ok], v[ok])^2 }
r2cols <- function(M, X) { q <- qr(X); 1 - colSums(qr.resid(q, M)^2) / colSums(sweep(M, 2, colMeans(M))^2) }

variance <- pca_t <- sep <- km <- list()
for (v in names(views)) {
  M <- views[[v]]
  vs <- data.table(regulator = colnames(M), subtype = r2cols(M, model.matrix(~ sub)),
                   plate = r2cols(M, model.matrix(~ plate)), estimate = r2cols(M, cbind(1, est)))
  variance[[v]] <- vs[, .(view = v, median_subtype = round(median(subtype), 3), median_plate = round(median(plate), 3),
    median_estimate = round(median(estimate), 3), frac_subtype_gt_0.3 = round(mean(subtype > 0.3), 3),
    frac_plate_gt_0.3 = round(mean(plate > 0.3), 3), frac_estimate_gt_0.3 = round(mean(estimate > 0.3), 3),
    top_subtype = paste(head(regulator[order(-subtype)], 8), collapse = ","))]
  p <- prcomp(M, center = TRUE); ve <- p$sdev^2 / sum(p$sdev^2)
  pca_t[[v]] <- rbindlist(lapply(1:6, function(k) data.table(view = v, pc = k, var_explained = round(ve[k], 3),
    r2_subtype = round(r2f(p$x[, k], sub), 3), r2_plate = round(r2f(p$x[, k], plate), 3),
    r2_estimate = round(r2c(p$x[, k], est), 3), r2_absolute = round(r2c(p$x[, k], meta$absolute_purity), 3),
    top_regulators = paste(names(sort(abs(p$rotation[, k]), decreasing = TRUE))[1:6], collapse = ","))))
  X <- p$x[, 1:N_PC]; d <- dist(X)
  sil <- silhouette(as.integer(sub), d)
  sep[[v]] <- data.table(view = v, subtype = levels(sub), n = as.vector(table(sub)),
    pam50_silhouette = round(tapply(sil[, 3], sub, mean), 3), overall_pam50_silhouette = round(mean(sil[, 3]), 3),
    within_effective_dim = sapply(levels(sub), function(s) { e <- prcomp(M[sub == s, ])$sdev^2; round(sum(e)^2 / sum(e^2), 2) }))
  km[[v]] <- rbindlist(lapply(levels(sub), function(s) {
    Ms <- M[sub == s, ]; Xs <- prcomp(Ms)$x[, 1:N_PC]; ds <- dist(Xs)
    Xn <- apply(Xs, 2, sample); dn <- dist(Xn)  # same per-PC spread, no joint structure
    rbindlist(lapply(KS, function(k) {
      a <- kmeans(Xs, k, nstart = 20, iter.max = 50); b <- kmeans(Xn, k, nstart = 20, iter.max = 50)
      data.table(view = v, subtype = s, k = k, silhouette = round(mean(silhouette(a$cluster, ds)[, 3]), 3),
                 silhouette_null = round(mean(silhouette(b$cluster, dn)[, 3]), 3),
                 smallest_cluster = min(table(a$cluster)))
    }))
  }))
}
for (k in c("variance", "pca_t", "sep", "km")) {
  f <- c(variance = "variance", pca_t = "pca", sep = "separation", km = "kmeans")[[k]]
  fwrite(rbindlist(get(k)), file.path(out, paste0("brca_heterogeneity_", f, ".tsv")), sep = "\t")
}
message("BRCA HETEROGENEITY DONE")
