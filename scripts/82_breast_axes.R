#!/usr/bin/env Rscript
# Step 8c: continuous axes of variation within each breast subtype (BRCA-Basal, Her2, LumA, LumB).
# Step 8 found no robust clusters in breast; the user (2026-10-12) asked for continuous axes instead and delegated the
# choice. Choice: principal components of the Step 8 module features (modules10: one z-scored score per regulator module
# of >= 10 regulators), because modules are many (64-68 per subtype) and overlapping, PCs give a few independent axes
# without hand-picking biology, and fixed loadings can be projected onto other cohorts (METABRIC, CPTAC-BRCA).
# Rule fixed before running: an axis is kept if all of
#   1. it explains more variance than the same-rank PC of 500 column-permuted datasets (95th percentile; parallel
#      analysis, Horn 1965); at most 5 PCs considered;
#   2. it is stable: median |cosine| between its loadings and the best-matching PC of 200 patient bootstraps >= 0.8;
#   3. it is not technical: plate R^2 in excess of chance, (R^2 - chance) / (1 - chance) with chance = (plates - 1) / (n - 1),
#      < 0.2 (a fixed raw-R^2 cutoff was used first and replaced: chance R^2 ranges 0.07-0.31 with group size), and
#      |Spearman| < 0.5 with the RNA-quality index and with ESTIMATE purity.
# Each axis is then described by its top-loading modules (hub + top members) and Spearman correlations with 50 MSigDB
# Hallmark program scores (mean z of program genes in purity + quality corrected tumor logCPM, as in 62), ESTIMATE
# immune, cytolytic activity and the recurring PILRB/CCNL2 module.
# Outputs: results/step8/breast_axes_{pca,annotation,loadings,scores}.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(82)
GROUPS <- c("BRCA-Basal", "BRCA-Her2", "BRCA-LumA", "BRCA-LumB")
MAX_PC <- 5; N_PERM <- 500; N_BOOT <- 200; MIN_COS <- 0.8; MAX_PLATE_EXCESS <- 0.2; MAX_RHO <- 0.5
RECURRING <- c("PILRB", "CCNL2", "ZNF692", "ANKZF1", "HSF4", "MAPK8IP3", "AKAP8L", "GIGYF1")
res_dir <- "results/step8"
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
qtab <- fread("results/step6/aracne_project/quality_correction_map/quality_index.tsv")[role == "tumor"]
memb <- fread(file.path(res_dir, "module_features.tsv"))[feature_set == "modules10"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
hm <- lapply(strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t"), function(l) list(name = sub("^HALLMARK_", "", l[1]), genes = l[-(1:2)]))
hm <- setNames(lapply(hm, function(h) na.omit(unname(sym2ens[h$genes]))), sapply(hm, `[[`, "name"))
sp <- function(a, b) { ok <- !is.na(a) & !is.na(b); cor(a[ok], b[ok], method = "spearman") }
r2f <- function(y, f) { f <- droplevels(factor(f)); 1 - sum(qr.resid(qr(model.matrix(~ f)), y)^2) / sum((y - mean(y))^2) }
zrow <- function(m) (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6)

run <- function(gid) {
  x <- fread(sprintf("data/processed/step8/features/%s_modules10.tsv.gz", gid)); F <- as.matrix(x[, -1]); rownames(F) <- x$sample
  n <- nrow(F); p <- ncol(F); pc <- prcomp(F, center = TRUE, scale. = TRUE); ve <- pc$sdev^2 / sum(pc$sdev^2)
  null <- replicate(N_PERM, { Fp <- apply(F, 2, sample); s <- prcomp(Fp, scale. = TRUE)$sdev[1:MAX_PC]^2; s / p })
  thr <- apply(null, 1, quantile, 0.95)
  L <- pc$rotation[, 1:MAX_PC]
  cosb <- replicate(N_BOOT, { i <- sample(n, replace = TRUE); Lb <- prcomp(F[i, ], scale. = TRUE)$rotation[, 1:(MAX_PC + 2)]
    apply(abs(crossprod(L, Lb)), 1, max) })
  meta <- cov[match(rownames(F), sample)]; qi <- qtab[group_id == gid][match(rownames(F), sample), quality_index]
  # Hallmark scores on purity + quality corrected expression (as in 62_recurring_module.R)
  lc <- cpm(calcNormFactors(DGEList(counts[, rownames(F)])), log = TRUE, prior.count = 1); lc <- lc[apply(lc, 1, median) >= 1, ]
  pur <- meta$estimate_score; ok <- !is.na(pur) & !is.na(qi)
  lcc <- lc; lcc[, ok] <- t(qr.resid(qr(cbind(1, pur[ok] - mean(pur[ok]), qi[ok] - mean(qi[ok]))), t(lc[, ok]))) + rowMeans(lc[, ok])
  zc <- zrow(lcc)
  hs <- sapply(hm, function(g) { g <- intersect(g, rownames(zc)); if (length(g) >= 10) colMeans(zc[g, , drop = FALSE]) else rep(NA, ncol(zc)) })
  cyt <- colMeans(lc[sym2ens[c("GZMA", "PRF1")], ])
  rec_feat <- memb[group_id == gid, .N, by = .(feature, hit = regulator %in% RECURRING)][hit == TRUE][which.max(N), feature]
  rec <- if (length(rec_feat)) F[, rec_feat] else rep(NA_real_, n)
  pca <- rbindlist(lapply(1:MAX_PC, function(k) {
    s <- pc$x[, k]
    plate_r2 <- r2f(s, meta$batch); ch <- (uniqueN(meta$batch) - 1) / (n - 1)
    data.table(group_id = gid, axis = paste0("PC", k), n_patients = n, n_modules = p, var_explained = round(ve[k], 4),
      null_95 = round(thr[k], 4), above_null = ve[k] > thr[k], stability_cos = round(median(cosb[k, ]), 3),
      plate_R2 = round(plate_r2, 3), plate_R2_chance = round(ch, 3), plate_excess = round((plate_r2 - ch) / (1 - ch), 3),
      rho_quality = round(sp(s, qi), 3), rho_purity = round(sp(s, pur), 3)) }))
  # cumulative rule: a PC only counts as above the null if all earlier PCs are
  pca[, above_null := cumprod(above_null) == 1]
  pca[, keep := above_null & stability_cos >= MIN_COS & plate_excess < MAX_PLATE_EXCESS & abs(rho_quality) < MAX_RHO & abs(rho_purity) < MAX_RHO]
  pca[, reason := fifelse(keep, "kept", fifelse(!above_null, "below null", fifelse(stability_cos < MIN_COS, "unstable", "technical")))]
  ann <- rbindlist(lapply(1:MAX_PC, function(k) {
    s <- pc$x[, k]; ld <- L[, k]
    o <- order(-abs(ld))[1:6]
    top <- vapply(o, function(j) { m <- memb[group_id == gid & feature == names(ld)[j], regulator]
      sprintf("%s%s(%s)", if (ld[j] > 0) "+" else "-", names(ld)[j], paste(head(m, 5), collapse = ",")) }, "")
    h <- apply(hs, 2, sp, b = s); ho <- order(-abs(h))[1:5]
    data.table(group_id = gid, axis = paste0("PC", k), top_modules = paste(top, collapse = "; "),
      top_hallmarks = paste(sprintf("%s %+.2f", names(h)[ho], h[ho]), collapse = "; "),
      rho_immune = round(sp(s, meta$estimate_immune), 2), rho_cytolytic = round(sp(s, cyt), 2),
      rho_proliferation = round(sp(s, hs[, "E2F_TARGETS"]), 2), rho_recurring_module = round(sp(s, rec), 2)) }))
  ld <- data.table(group_id = gid, module = rownames(L), round(L, 4))
  sc <- data.table(group_id = gid, sample = rownames(F), round(pc$x[, 1:MAX_PC], 4))
  message(sprintf("%s: %d patients x %d modules; kept %s", gid, n, p, paste(pca[keep == TRUE, axis], collapse = ",")))
  list(pca = pca, ann = ann, ld = ld, sc = sc)
}
res <- lapply(GROUPS, run)
pick <- function(k) rbindlist(lapply(res, `[[`, k))
fwrite(pick("pca"), file.path(res_dir, "breast_axes_pca.tsv"), sep = "\t")
fwrite(pick("ann"), file.path(res_dir, "breast_axes_annotation.tsv"), sep = "\t")
fwrite(pick("ld"), file.path(res_dir, "breast_axes_loadings.tsv"), sep = "\t")
fwrite(pick("sc"), file.path(res_dir, "breast_axes_scores.tsv"), sep = "\t")
message("BREAST AXES DONE")
