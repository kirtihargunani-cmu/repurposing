#!/usr/bin/env Rscript
# Diagnostic: how much of each ARACNe edge is explained by known confounders? No edges are removed.
#   Edge strength proxy = Spearman correlation of regulator and target in the network's own input matrix
#   (cohort tumors, log2(TPM+1)); ARACNe uses MI, which also sees non-linear links, so edges with small |r_raw|
#   are reported separately. Each gene's ranks are residualized on a covariate set and the correlation recomputed
#   (partial Spearman):
#     r_raw        no adjustment
#     r_estimate   ESTIMATE score (purity / stroma / immune)
#     r_plate      TCGA plate (batch); TARGET cohorts: protocol/batch column as recorded
#     r_plate_shuf plate labels shuffled (same plate sizes): attenuation expected from removing that many
#                  degrees of freedom by chance
#     r_both       ESTIMATE + plate
#     r_subtype    subtype / group (only for cohorts pooling > 1 group): edges driven by differences between
#                  subtypes rather than variation inside them
#     r_quality    RNA-quality index (per-sample mean z of long genes minus short genes, top vs bottom genomic-span
#                  quintile of expressed genes; as in scripts/57 and 61)
#     r_all        ESTIMATE + quality + plate (the full technical adjustment)
#   retained_<x> = |r_x| / |r_raw|
# Outputs: data/processed/step5/aracne/<NET>/edge_confounding.tsv.gz
#          results/step5/aracne_edge_confounding.tsv (per network summary; 'collapse' = retains < 50% of |r_raw|,
#          an illustrative cut, not a filter)
# VERSION=v2: networks in aracne_v2/, summary results/step5/aracne_edge_confounding_v2.tsv (all consolidated edges,
#   before the cross-mappable filter of scripts/44)
suppressPackageStartupMessages({library(data.table); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(46)
SFX <- if (identical(Sys.getenv("VERSION", "v1"), "v2")) "_v2" else ""
IN <- paste0("data/processed/step5/aracne_inputs", SFX); NETDIR <- paste0("data/processed/step5/aracne", SFX)
nets <- commandArgs(trailingOnly = TRUE)
if (!length(nets)) nets <- basename(dirname(Sys.glob(file.path(NETDIR, "*", "edge_support.tsv.gz"))))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))
qindex <- function(m) {  # genes x samples log expression -> long-minus-short index per sample
  m <- m[apply(m, 1, median) >= 1 & rownames(m) %in% names(span), , drop = FALSE]
  z <- (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  colMeans(z[span[rownames(z)] >= q[2], , drop = FALSE]) - colMeans(z[span[rownames(z)] <= q[1], , drop = FALSE])
}

zres <- function(R, X) {  # residualize rows of R (genes x samples) on design X, then row-standardize
  E <- if (is.null(X)) R - rowMeans(R) else t(qr.resid(qr(X), t(R)))
  E / sqrt(rowSums(E^2))
}
edge_r <- function(Z, a, b) {
  out <- numeric(length(a))
  for (s in split(seq_along(a), ceiling(seq_along(a) / 50000))) out[s] <- rowSums(Z[a[s], , drop = FALSE] * Z[b[s], , drop = FALSE])
  out
}

run <- function(net) {
  x <- fread(file.path(IN, paste0(net, ".exp.tsv")))
  E0 <- as.matrix(x[, -1]); rownames(E0) <- x$gene; qi <- qindex(E0)
  R <- t(apply(E0, 1, rank)); rm(E0)
  meta <- cov[match(colnames(R), sample)]; stopifnot(!anyNA(meta$sample))
  est <- meta$estimate_score; plate <- factor(fifelse(is.na(meta$batch), "unknown", as.character(meta$batch)))
  grp <- factor(meta$group_id)
  X <- list(raw = NULL, estimate = cbind(1, est))
  if (nlevels(plate) > 1) X <- c(X, list(plate = model.matrix(~ plate), plate_shuf = model.matrix(~ sample(plate)),
                                         both = model.matrix(~ est + plate)))  # single-plate cohorts: no plate terms
  if (nlevels(grp) > 1) X$subtype <- model.matrix(~ grp)
  X$quality <- cbind(1, qi)
  X$all <- if (nlevels(plate) > 1) model.matrix(~ est + qi + plate) else cbind(1, est, qi)
  e <- fread(file.path(NETDIR, net, "edge_support.tsv.gz"), select = c("Regulator", "Target", "MI", "n_boot"))
  ia <- match(e$Regulator, rownames(R)); ib <- match(e$Target, rownames(R))
  for (k in names(X)) { Z <- zres(R, X[[k]]); e[, paste0("r_", k) := round(edge_r(Z, ia, ib), 4)] }
  fwrite(e, file.path(NETDIR, net, "edge_confounding.tsv.gz"), sep = "\t")

  lin <- e[abs(r_raw) >= 0.2]  # edges with a clear monotone link; the rest are reported only by count
  ret <- function(k) abs(lin[[paste0("r_", k)]]) / abs(lin$r_raw)
  flip <- function(k) mean(sign(lin[[paste0("r_", k)]]) != sign(lin$r_raw))
  s <- data.table(net = net, n_samples = ncol(R), n_plates = nlevels(plate), n_groups = nlevels(grp),
    r2_quality_on_plate = if (nlevels(plate) > 1) round(summary(lm(qi ~ plate))$r.squared, 3) else NA_real_,
    edges = nrow(e), frac_weak_linear = round(mean(abs(e$r_raw) < 0.2), 3), median_abs_r_raw = round(median(abs(e$r_raw)), 3))
  for (k in setdiff(names(X), "raw")) {
    s[, paste0("median_retained_", k) := round(median(ret(k)), 3)]
    s[, paste0("collapse_", k) := round(mean(ret(k) < 0.5), 4)]
    s[, paste0("flip_", k) := round(flip(k), 4)]
  }
  # are confounded edges the weakly supported ones? (ESTIMATE + plate, or ESTIMATE alone without plates)
  adj <- if (is.null(X$both)) "estimate" else "both"
  s[, collapse_adj_low_support := round(mean(ret(adj)[lin$n_boot <= median(e$n_boot)] < 0.5), 4)]
  s[, collapse_adj_high_support := round(mean(ret(adj)[lin$n_boot > quantile(e$n_boot, 0.9)] < 0.5), 4)]
  s
}
res <- rbindlist(mclapply(nets, function(n) tryCatch(run(n), error = function(e) { message(n, ": ", conditionMessage(e)); NULL }),
                          mc.cores = 8, mc.preschedule = FALSE), fill = TRUE)
fwrite(res[order(net)], paste0("results/step5/aracne_edge_confounding", SFX, ".tsv"), sep = "\t")
message("EDGE CONFOUNDING DONE: ", nrow(res), " networks")
