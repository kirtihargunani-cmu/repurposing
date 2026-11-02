#!/usr/bin/env Rscript
# Which healthy reference generalizes better: TCGA adjacent normals or GTEx normals? Measures only.
# Per organ unit (groups sharing GTEx donors; adjacent normals pooled within an organ, one per patient):
#  1. Spread at equal n (N_REP subsamples of n = min(n_adjacent, n_gtex, 100) per arm; genes passing
#     filterByExpr in the unit; joint TMM logCPM), raw and after removing age + sex + ESTIMATE within each arm:
#       gene variance ratio      median log2(var_adjacent / var_gtex) over genes, and over genes with mean logCPM >= 3
#       biological CV (edgeR)    common BCV per arm (count-aware; not inflated by low counts)
#       multivariate spread      median distance to the arm centroid in the top 20 PCs of the joint data
#  2. What the extra spread in adjacent normals is: per adjacent sample, distance from the GTEx centroid vs
#     (a) tumor-likeness = projection of (sample - GTEx mean) onto the unit's tumor-vs-GTEx logFC, (b) ESTIMATE,
#     (c) TCGA plate (R^2 and shuffled-plate chance).
#  3. Generalization (the target property: a signature that carries to new patients): split-half reproducibility.
#     Tumors and normals are split by PATIENT into two disjoint halves; signature 1 = tumors A vs normals A,
#     signature 2 = tumors B vs normals B (limma-voom ~ condition + age + ESTIMATE (+ sex)); Spearman(sig1, sig2).
#     Same tumor halves and the same normal n for both references, so only the reference differs.
#     Also: agreement of each half-signature with the full Step 4 signature of the other reference.
# Outputs: results/step4/reference_diagnostics/normal_variance_{spread,adjacent_drivers,splithalf,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); RNGkind("L'Ecuyer-CMRG"); set.seed(33)  # reproducible streams across mclapply workers
N_REP <- as.integer(Sys.getenv("N_REP", "30")); N_MAX <- 100; N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
out <- "results/step4/reference_diagnostics"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
cov <- cov[sample %in% colnames(counts) & !is.na(age) & !is.na(sex) & !is.na(estimate_score)]
UNITS <- list(breast = c("BRCA-Basal", "BRCA-LumA", "BRCA-LumB", "BRCA-Her2"), kidney = c("TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH"),
              lung = c("TCGA-LUAD", "TCGA-LUSC"), colorectum = c("TCGA-COAD", "TCGA-READ"), head_neck = "TCGA-HNSC",
              esophagus = "ESCA-AD", liver = "TCGA-LIHC", prostate = "TCGA-PRAD", stomach = "TCGA-STAD",
              thyroid = "TCGA-THCA", uterus = "UCEC-EEC")

unit_data <- function(u) {
  gs <- UNITS[[u]]; d <- cov[group_id %in% gs]
  adj <- unique(d[role == "adjacent_normal"], by = "patient"); gtx <- unique(d[role == "gtex_normal"], by = "sample")
  tum <- unique(d[role == "tumor" & sample_type_class %in% c("primary", "new_primary")], by = "patient")
  list(adj = adj, gtx = gtx, tum = tum, groups = gs)
}
lcpm <- function(s, keep = NULL) { y <- DGEList(counts[, s]); if (is.null(keep)) keep <- filterByExpr(y, group = rep(1, length(s)))
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); list(e = cpm(y, log = TRUE, prior.count = 2), y = y, keep = keep) }
resid_cov <- function(e, d) {  # remove age + ESTIMATE (+ sex) within one arm, keep the arm mean
  X <- cbind(1, as.numeric(scale(d$age)), as.numeric(scale(d$estimate_score)))
  if (uniqueN(d$sex) > 1) X <- cbind(X, as.numeric(d$sex == d$sex[1]))
  X[is.na(X)] <- 0; f <- lmFit(e, X); e - f$coefficients[, -1, drop = FALSE] %*% t(X[, -1, drop = FALSE]) }

spread_unit <- function(u) {
  U <- unit_data(u); n <- min(nrow(U$adj), nrow(U$gtx), N_MAX); if (n < 10) return(NULL)
  rbindlist(lapply(seq_len(N_REP), function(r) {
    a <- U$adj[sample(.N, n)]; g <- U$gtx[sample(.N, n)]; s <- c(a$sample, g$sample); dd <- rbind(a, g)
    L <- lcpm(s); e <- L$e; ia <- seq_len(n); ig <- n + seq_len(n)
    rows <- list()
    for (adj_cov in c(FALSE, TRUE)) {
      ea <- e[, ia]; eg <- e[, ig]
      if (adj_cov) { ea <- resid_cov(ea, a); eg <- resid_cov(eg, g) }
      va <- apply(ea, 1, var); vg <- apply(eg, 1, var); mu <- rowMeans(e); hi <- mu >= 3
      lr <- log2((va + 1e-4) / (vg + 1e-4))
      p <- prcomp(t(cbind(ea, eg)), rank. = 20); sc <- p$x
      da <- sqrt(rowSums(sweep(sc[ia, ], 2, colMeans(sc[ia, ]))^2)); dg <- sqrt(rowSums(sweep(sc[ig, ], 2, colMeans(sc[ig, ]))^2))
      rows[[length(rows) + 1]] <- data.table(unit = u, rep = r, n_per_arm = n, covariates_removed = adj_cov,
        median_log2_var_ratio = median(lr), median_log2_var_ratio_mean3 = median(lr[hi]), frac_genes_adj_more_var = mean(va > vg),
        pc_spread_ratio = median(da) / median(dg))
    }
    # BCV per arm (raw counts, edgeR common dispersion, intercept-only)
    bcv <- function(ss) { y <- calcNormFactors(DGEList(counts[L$keep, ss])); sqrt(estimateCommonDisp(y)$common.dispersion) }
    out <- rbindlist(rows); out[, `:=`(bcv_adj = bcv(a$sample), bcv_gtex = bcv(g$sample))]; out
  }))
}

drivers_unit <- function(u) {
  U <- unit_data(u); a <- U$adj; g <- U$gtx; L <- lcpm(c(a$sample, g$sample)); e <- L$e
  # within-arm deviations on a common scale (GTEx gene SD), so the lab offset between arms does not count as spread
  ia <- seq_len(nrow(a)); ig <- nrow(a) + seq_len(nrow(g)); gsd <- pmax(apply(e[, ig], 1, sd), 0.1)
  dev <- (e[, ia, drop = FALSE] - rowMeans(e[, ia])) / gsd
  dist_gtex <- sqrt(colMeans(dev^2))                                          # adjacent sample's distance from the adjacent centroid
  dist_gtex_ref <- sqrt(colMeans(((e[, ig] - rowMeans(e[, ig])) / gsd)^2))    # GTEx sample's distance from the GTEx centroid
  # tumor-likeness, same lab only: deviation from the other adjacents aligned with the tumor-vs-ADJACENT signature
  za4 <- fread("results/step4/z_tumor_vs_adjacent.tsv"); ga <- intersect(U$groups, names(za4))
  zt <- rowMeans(as.matrix(za4[match(rownames(e), za4$gene_id), ga, with = FALSE]), na.rm = TRUE); ok <- !is.na(zt)
  tlike <- as.numeric(cor(dev[ok, ], zt[ok]))
  plate_r2 <- function(x) { b <- factor(a$batch); if (nlevels(b) < 2) return(c(NA, NA)); r2 <- summary(lm(x ~ b))$r.squared
    ch <- mean(replicate(200, summary(lm(x ~ sample(b)))$r.squared)); c(r2, ch) }
  pr <- plate_r2(dist_gtex)
  list(per_sample = data.table(unit = u, sample = a$sample, group_id = a$group_id, dist_within_adjacent = dist_gtex, tumor_likeness = tlike,
                               estimate = a$estimate_score, batch = a$batch),
       summary = data.table(unit = u, n_adj = nrow(a), n_gtex = nrow(g),
         median_dist_adj_within = median(dist_gtex), median_dist_gtex_within = median(dist_gtex_ref),
         frac_adj_beyond_gtex_95pct = mean(dist_gtex > quantile(dist_gtex_ref, 0.95)),
         median_tumor_likeness = median(tlike), rho_dist_vs_tumor_likeness = cor(dist_gtex, tlike, method = "spearman"),
         rho_dist_vs_estimate = cor(dist_gtex, a$estimate_score, method = "spearman"), plate_r2 = pr[1], plate_r2_chance = pr[2]))
}

splithalf_unit <- function(u) {
  U <- unit_data(u); n_norm <- floor(min(nrow(U$adj), nrow(U$gtx)) / 2); if (n_norm < 5) return(NULL)
  n_tum <- min(floor(nrow(U$tum) / 2), 150)
  fit_sig <- function(ts, ns) { d <- rbind(ts, ns); d[, cond := factor(ifelse(role == "tumor", "t", "n"), c("n", "t"))]
    d[, `:=`(age_z = as.numeric(scale(age)), est_z = as.numeric(scale(estimate_score)))]
    terms <- c("cond", "age_z", "est_z", if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= 3) "sex")
    X <- model.matrix(reformulate(terms), data = d); y <- DGEList(counts[, d$sample]); k <- filterByExpr(y, X)
    y <- calcNormFactors(y[k, , keep.lib.sizes = FALSE]); fit <- eBayes(lmFit(voom(y, X), X))
    tt <- topTable(fit, coef = "condt", number = Inf, sort.by = "none"); setNames(tt$t, rownames(tt)) }
  rho <- function(x, y) { k <- intersect(names(x), names(y)); cor(x[k], y[k], method = "spearman", use = "complete.obs") }
  zfull <- function(ct) { f <- fread(sprintf("results/step4/z_tumor_vs_%s.tsv", ct)); gs <- intersect(U$groups, names(f))
    setNames(rowMeans(as.matrix(f[, gs, with = FALSE]), na.rm = TRUE), sub("\\..*", "", f$gene_id)) }
  full_gtex <- zfull("gtex"); full_adj <- zfull("adjacent")
  rbindlist(lapply(seq_len(N_REP), function(r) {
    # split patients: those with an adjacent normal are split first, so tumor and normal halves come from disjoint patients
    pa <- sample(U$adj$patient); half <- setNames(rep(1:2, length.out = length(pa)), pa)
    rest <- setdiff(U$tum$patient, pa); half <- c(half, setNames(sample(rep(1:2, length.out = length(rest))), rest))
    tA <- U$tum[half[patient] == 1][sample(.N, min(.N, n_tum))]; tB <- U$tum[half[patient] == 2][sample(.N, min(.N, n_tum))]
    aA <- U$adj[half[patient] == 1][sample(.N, min(.N, n_norm))]; aB <- U$adj[half[patient] == 2][sample(.N, min(.N, n_norm))]
    gg <- U$gtx[sample(.N)]; gA <- gg[seq_len(n_norm)]; gB <- gg[n_norm + seq_len(n_norm)]
    sA_adj <- fit_sig(tA, aA); sB_adj <- fit_sig(tB, aB); sA_gtx <- fit_sig(tA, gA); sB_gtx <- fit_sig(tB, gB)
    data.table(unit = u, rep = r, n_tumor_half = nrow(tA), n_normal_half = n_norm,
               splithalf_adjacent = rho(sA_adj, sB_adj), splithalf_gtex = rho(sA_gtx, sB_gtx),
               adjhalf_vs_full_gtex = rho(sA_adj, full_gtex), gtexhalf_vs_full_adj = rho(sA_gtx, full_adj))
  }))
}

us <- names(UNITS)
sp <- rbindlist(mclapply(us, spread_unit, mc.cores = N_CORES))
dr <- mclapply(us, drivers_unit, mc.cores = N_CORES)
sh <- rbindlist(mclapply(us, splithalf_unit, mc.cores = N_CORES))
fwrite(sp, file.path(out, "normal_variance_spread_reps.tsv"), sep = "\t")
fwrite(rbindlist(lapply(dr, `[[`, "per_sample")), file.path(out, "normal_variance_adjacent_samples.tsv"), sep = "\t")
drs <- rbindlist(lapply(dr, `[[`, "summary")); fwrite(drs, file.path(out, "normal_variance_adjacent_drivers.tsv"), sep = "\t")
fwrite(sh, file.path(out, "normal_variance_splithalf_reps.tsv"), sep = "\t")
spr <- sp[, lapply(.SD, median), by = .(unit, covariates_removed), .SDcols = c("n_per_arm", "median_log2_var_ratio", "median_log2_var_ratio_mean3",
                                                                              "frac_genes_adj_more_var", "pc_spread_ratio", "bcv_adj", "bcv_gtex")]
shs <- sh[, .(n_tumor_half = median(n_tumor_half), n_normal_half = median(n_normal_half),
              splithalf_adjacent = median(splithalf_adjacent), splithalf_gtex = median(splithalf_gtex),
              gtex_better_frac = mean(splithalf_gtex > splithalf_adjacent),
              adjhalf_vs_full_gtex = median(adjhalf_vs_full_gtex), gtexhalf_vs_full_adj = median(gtexhalf_vs_full_adj)), by = unit]
summ <- merge(merge(spr[covariates_removed == FALSE], spr[covariates_removed == TRUE, .(unit, var_ratio_after_cov = median_log2_var_ratio,
              pc_spread_ratio_after_cov = pc_spread_ratio)], by = "unit"), shs, by = "unit", all = TRUE)
summ <- merge(summ, drs, by = "unit", all = TRUE)
num <- names(summ)[sapply(summ, is.numeric)]; summ[, (num) := lapply(.SD, round, 3), .SDcols = num]
fwrite(summ, file.path(out, "normal_variance_summary.tsv"), sep = "\t"); print(summ)
message("NORMAL VARIANCE DONE")
