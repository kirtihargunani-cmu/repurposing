#!/usr/bin/env Rscript
# Step 9b: does the TCGA-COAD 2-cluster split (Step 8; MSI-like vs rest) replicate in GSE39582 and CPTAC-COAD?
# PI decision 2026-10-12: run all three methods + compare with CMS subtypes.
# Common yardstick: every patient (TCGA-COAD, GSE39582, CPTAC-COAD) is scored on the 78 TCGA-COAD module features
#   (results/step8/module_features.tsv, modules10). Module score = mean of the member regulators' NES, each z-scored
#   within its own cohort and sign-aligned as in TCGA (sign of the member's correlation with the TCGA module mean);
#   members missing in a cohort are skipped (module dropped if < 50% present). Each feature is then z-scored within
#   the cohort, so cohort-level shifts (normals vs cohort-average reference, array vs RNA-seq) drop out.
#   MSI labels: TCGA = PanCanAtlas GI.HM-indel (hypermutated-indel) vs other; GSE39582 = dMMR vs pMMR;
#   CPTAC = MSI-H vs MSS.
# 1 transfer     TCGA cluster centroids (mean feature profile per cluster); each replication patient gets the cluster
#                whose centroid it correlates with best (Spearman, as in PAM50 assignment). Tests: Fisher test of
#                assigned MSI-like cluster vs MSI status; in-group proportion (IGP; Kapp & Tibshirani 2007, clusterRepro,
#                1,000 permutations): share of patients whose nearest neighbour is in the same assigned cluster.
# 2 independent  M3C from scratch on each replication cohort's feature matrix (same settings and decision rule as
#                Step 8: K 2-10, PAM, 100 resamples, 100 Monte Carlo references; clusters only if both p < 0.05);
#                compared with the transferred labels (ARI) and with MSI.
# 3 modules      per feature, standardized mean difference (Cohen's d): TCGA MSI-like cluster vs rest, and MSI vs MSS
#                in each cohort; Spearman correlation of the d's across features, sign agreement of TCGA's top-10
#                defining modules, Wilcoxon p (BH) per feature in each cohort.
# CMS            CMScaller 2.0.1 (Eide et al. 2017; nearest-template prediction, 1,000 permutations, FDR < 0.05,
#                unconfident calls left unassigned) on uncorrected data: TCGA COAD/READ counts, GSE39582 array,
#                CPTAC counts. Cross-tabulated with clusters and MSI.
# Outputs: results/step9/crc_{features_check,transfer,igp,independent,module_effects,cms,assignments,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(M3C); library(mclust); library(clusterRepro); library(CMScaller)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(91)
N_CORES <- as.integer(Sys.getenv("N_CORES", "16")); out <- "results/step9"

load_nes <- function(f) { x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator; A }
mf <- fread("results/step8/module_features.tsv")[group_id == "TCGA-COAD" & feature_set == "modules10"]
tc <- load_nes("data/processed/step6/aracne_project_map/TCGA-COAD.patient_nes.tsv.gz")
rep <- list(gse39582 = load_nes("data/processed/step9/gse39582/patient_nes.tsv.gz"),
            cptac_coad = load_nes("data/processed/step9/cptac_coad/patient_nes.tsv.gz"))
gse_s <- fread("results/step9/gse39582_samples.tsv"); cp_s <- fread("results/step9/cptac_coad_samples.tsv")
rep$gse39582 <- rep$gse39582[rownames(rep$gse39582) %in% gse_s[role == "tumor", sample], ]

# signs from TCGA
zc <- function(A) scale(A)
Ztc <- zc(tc); feats <- unique(mf$feature)
sgn <- rbindlist(lapply(feats, function(f) { r <- intersect(mf[feature == f, regulator], colnames(Ztc))
  data.table(feature = f, regulator = r, sign = sign(cor(Ztc[, r, drop = FALSE], rowMeans(Ztc[, r, drop = FALSE])))[, 1]) }))
modscore <- function(A) { Z <- zc(A)
  F <- sapply(feats, function(f) { s <- sgn[feature == f]; r <- s$regulator[s$regulator %in% colnames(Z)]
    if (length(r) < 0.5 * nrow(s)) return(rep(NA_real_, nrow(Z)))
    rowMeans(sweep(Z[, r, drop = FALSE], 2, s[match(r, regulator), sign], `*`)) })
  rownames(F) <- rownames(A); F <- F[, colSums(is.na(F)) == 0, drop = FALSE]; scale(F) }
Ftc <- modscore(tc); Frep <- lapply(rep, modscore)
common <- Reduce(intersect, c(list(colnames(Ftc)), lapply(Frep, colnames)))
Ftc <- Ftc[, common]; Frep <- lapply(Frep, function(F) F[, common])
s8 <- fread(sprintf("data/processed/step8/features/TCGA-COAD_modules10.tsv.gz")); P <- as.matrix(s8[, -1]); rownames(P) <- s8$sample
fchk <- data.table(feature = common, r_mean_vs_pc1 = round(diag(cor(Ftc[rownames(P), common], P[, common])), 3))
fwrite(fchk, file.path(out, "crc_features_check.tsv"), sep = "\t")
message(sprintf("features: %d of %d TCGA modules usable in all cohorts; mean-score vs Step 8 PC1-score r median %.3f",
                length(common), length(feats), median(fchk$r_mean_vs_pc1)))

# TCGA clusters and MSI
asg <- fread("results/step8/assignments.tsv")[group_id == "TCGA-COAD" & feature_set == "modules10" & K == 2]
sub <- fread("data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz"); setnames(sub, 1, "sample")
tl <- data.table(sample = rownames(Ftc)); tl[, cluster := asg$cluster[match(tl$sample, asg$sample)]]
stopifnot(identical(as.vector(sort(table(tl$cluster))), as.vector(sort(table(asg$cluster)))))  # labels follow the patients
tsub <- sub$Subtype_Selected[match(tl$sample, sub$sample)]
tl[, msi := fifelse(tsub == "GI.HM-indel", "MSI", fifelse(is.na(tsub), NA_character_, "MSS"))]
msi_cl <- tl[msi == "MSI", as.integer(names(which.max(table(cluster))))]
tl[, msi_like := cluster == msi_cl]
cent <- sapply(sort(unique(tl$cluster)), function(k) colMeans(Ftc[tl$cluster == k, , drop = FALSE])); colnames(cent) <- sort(unique(tl$cluster))
msi_of <- list(gse39582 = setNames(fifelse(gse_s$mmr.status == "dMMR", "MSI", fifelse(gse_s$mmr.status == "pMMR", "MSS", NA_character_)), gse_s$sample),
               cptac_coad = setNames(fifelse(cp_s$MSI_STATUS == "MSI-H", "MSI", fifelse(cp_s$MSI_STATUS == "MSS", "MSS", NA_character_)), cp_s$sample))

fisher_row <- function(cohort, method, msi_like, msi) { ok <- !is.na(msi) & !is.na(msi_like)
  t <- table(factor(msi_like[ok], c(TRUE, FALSE)), factor(msi[ok], c("MSI", "MSS"))); ft <- fisher.test(t)
  data.table(cohort = cohort, method = method, n = sum(ok), n_msi = sum(msi[ok] == "MSI"),
             msi_in_msi_like = sprintf("%d of %d", t["TRUE", "MSI"], sum(t[, "MSI"])),
             share_msi_captured = round(t["TRUE", "MSI"] / sum(t[, "MSI"]), 2),
             msi_rate_msi_like = round(t["TRUE", "MSI"] / sum(t["TRUE", ]), 2), msi_rate_other = round(t["FALSE", "MSI"] / sum(t["FALSE", ]), 2),
             odds_ratio = round(unname(ft$estimate), 2), p = signif(ft$p.value, 3)) }

# 1 transfer
transfer <- list(fisher_row("TCGA-COAD", "Step 8 clusters (reference)", tl$msi_like, tl$msi)); igp <- list(); asg_all <- list()
igp_tc <- clusterRepro(cent, t(Ftc), Number.of.permutations = 1000)
igp[["TCGA-COAD"]] <- data.table(cohort = "TCGA-COAD", cluster = colnames(cent), msi_like = colnames(cent) == msi_cl,
                                 n = igp_tc$Actual.Size, IGP = round(igp_tc$Actual.IGP, 3), p = igp_tc$p.value)
for (co in names(Frep)) { F <- Frep[[co]]
  r <- cor(t(F), cent, method = "spearman"); lab <- as.integer(colnames(cent)[max.col(r)])
  msi <- msi_of[[co]][rownames(F)]
  transfer[[co]] <- fisher_row(co, "transferred TCGA clusters", lab == msi_cl, msi)
  ig <- clusterRepro(cent, t(F), Number.of.permutations = 1000)
  igp[[co]] <- data.table(cohort = co, cluster = colnames(cent), msi_like = colnames(cent) == msi_cl, n = ig$Actual.Size, IGP = round(ig$Actual.IGP, 3), p = ig$p.value)
  asg_all[[co]] <- data.table(cohort = co, sample = rownames(F), transferred = lab, msi_like = lab == msi_cl, msi = msi)
}
asg_all[["TCGA-COAD"]] <- data.table(cohort = "TCGA-COAD", sample = tl$sample, transferred = tl$cluster, msi_like = tl$msi_like, msi = tl$msi)

# 2 independent M3C
ind <- list(); indsc <- list()
for (co in names(Frep)) { F <- Frep[[co]]
  ids <- sprintf("S%04d", seq_len(nrow(F)))  # M3C breaks on IDs starting with a digit (CPTAC "01CO001")
  back <- setNames(rownames(F), ids); Fm <- F; rownames(Fm) <- ids
  m <- M3C(as.data.frame(t(Fm)), cores = N_CORES, iters = 100, repsref = 100, repsreal = 100, maxK = 10, clusteralg = "pam",
           seed = 123, removeplots = TRUE, silent = TRUE)
  sc <- as.data.table(m$scores)[, cohort := co]; best <- sc[which.max(RCSI)]
  a <- m$realdataresults[[best$K]]$ordered_annotation; cl <- setNames(as.integer(as.character(a$consensuscluster)), back[rownames(a)])
  a2 <- m$realdataresults[[2]]$ordered_annotation; cl2 <- setNames(as.integer(as.character(a2$consensuscluster)), back[rownames(a2)])
  tr <- asg_all[[co]]; tr[, `:=`(independent_bestK = cl[sample], independent_K2 = cl2[sample])]; asg_all[[co]] <- tr
  msi <- tr$msi; k2msi <- names(which.max(table(tr$independent_K2[msi %in% "MSI"])))
  ind[[co]] <- data.table(cohort = co, best_K = best$K, RCSI = round(best$RCSI, 3), MONTECARLO_P = signif(best$MONTECARLO_P, 3),
    NORM_P = signif(best$NORM_P, 3), decision = if (best$MONTECARLO_P < 0.05 & best$NORM_P < 0.05) "clusters" else "no significant structure",
    cluster_sizes = paste(sort(table(cl), decreasing = TRUE), collapse = "/"),
    K2_p = signif(sc[K == 2, NORM_P], 3), ARI_K2_vs_transferred = round(adjustedRandIndex(tr$independent_K2, tr$transferred), 2),
    ARI_bestK_vs_transferred = round(adjustedRandIndex(tr$independent_bestK, tr$transferred), 2),
    K2_msi_fisher_p = signif(fisher.test(table(tr$independent_K2 == as.integer(k2msi), msi)[, c("MSI", "MSS")])$p.value, 3))
  indsc[[co]] <- sc
  message(sprintf("%s independent M3C: best K %d, p %.3g / %.3g -> %s", co, best$K, best$MONTECARLO_P, best$NORM_P, ind[[co]]$decision))
}

# 3 module effects
cohen <- function(x, g) { a <- x[g]; b <- x[!g]; (mean(a) - mean(b)) / sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2)) }
eff <- data.table(feature = common,
  d_tcga_cluster = apply(Ftc, 2, cohen, g = tl$msi_like),
  d_tcga_msi = apply(Ftc[!is.na(tl$msi), ], 2, cohen, g = tl$msi[!is.na(tl$msi)] == "MSI"))
for (co in names(Frep)) { F <- Frep[[co]]; msi <- msi_of[[co]][rownames(F)]; ok <- !is.na(msi)
  eff[, paste0("d_", co, "_msi") := apply(F[ok, ], 2, cohen, g = msi[ok] == "MSI")]
  eff[, paste0("fdr_", co, "_msi") := p.adjust(apply(F[ok, ], 2, function(x) wilcox.test(x ~ msi[ok])$p.value), "BH")] }
eff[, top10_tcga := rank(-abs(d_tcga_cluster)) <= 10]
eff[, members := sapply(feature, function(f) paste(head(mf[feature == f, regulator], 6), collapse = ","))]
num <- names(eff)[sapply(eff, is.numeric)]; eff[, (num) := lapply(.SD, signif, 3), .SDcols = num]
fwrite(eff[order(-abs(d_tcga_cluster))], file.path(out, "crc_module_effects.tsv"), sep = "\t")
modsum <- rbindlist(lapply(names(Frep), function(co) { d <- eff[[paste0("d_", co, "_msi")]]
  data.table(cohort = co, rho_vs_tcga_cluster_d = round(cor(eff$d_tcga_cluster, d, method = "spearman"), 2),
             top10_same_sign = sprintf("%d of 10", sum(sign(d[eff$top10_tcga]) == sign(eff$d_tcga_cluster[eff$top10_tcga]))),
             features_fdr05 = sum(eff[[paste0("fdr_", co, "_msi")]] < 0.05)) }))

# CMS
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
cms_run <- function(m, rnaseq, label) { r <- CMScaller(m, rowNames = "ensg", RNAseq = rnaseq, nPerm = 1000, FDR = 0.05, doPlot = FALSE, verbose = FALSE)
  data.table(cohort = label, sample = colnames(m), CMS = as.character(r$prediction)) }
tcga_cms <- rbindlist(lapply(c("TCGA-COAD", "TCGA-READ"), function(g) { s <- intersect(cov[group_id == g, sample], colnames(counts))
  cms_run(round(counts[, s]), TRUE, g) }))
graw <- fread("data/processed/step9/gse39582/expression_raw.tsv.gz"); G <- as.matrix(graw[, -1]); rownames(G) <- graw$gene_id
craw <- fread("data/processed/step9/cptac_coad/counts_raw.tsv.gz"); C <- as.matrix(craw[, -1]); rownames(C) <- craw$gene_id
cms <- rbind(tcga_cms, cms_run(G[, gse_s[role == "tumor", sample]], FALSE, "gse39582"), cms_run(C, TRUE, "cptac_coad"))
A <- rbindlist(asg_all, fill = TRUE); A <- merge(A, cms, by = c("cohort", "sample"), all.x = TRUE)
r8 <- fread("results/step8/assignments.tsv")[group_id == "TCGA-READ" & feature_set == "modules10" & K == 3]
cms_tab <- rbind(
  A[, .N, by = .(cohort, grouping = fifelse(msi_like, "MSI-like cluster", "other cluster"), CMS)],
  A[!is.na(msi), .N, by = .(cohort, grouping = paste0("MSI status: ", msi), CMS)],
  merge(r8[, .(sample, grouping = paste0("READ cluster ", cluster))], cms[cohort == "TCGA-READ"], by = "sample")[, .N, by = .(cohort, grouping, CMS)])
fwrite(cms_tab[order(cohort, grouping, CMS)], file.path(out, "crc_cms.tsv"), sep = "\t", na = "unassigned")
cms_sum <- A[, .(cms1_share_msi_like = round(mean(CMS[msi_like] %in% "CMS1"), 2), cms1_share_other = round(mean(CMS[!msi_like] %in% "CMS1"), 2),
                 unassigned = round(mean(is.na(CMS)), 2)), by = cohort]

fwrite(A, file.path(out, "crc_assignments.tsv"), sep = "\t", na = "NA")
fwrite(rbindlist(transfer), file.path(out, "crc_transfer.tsv"), sep = "\t")
fwrite(rbindlist(igp), file.path(out, "crc_igp.tsv"), sep = "\t")
fwrite(rbindlist(ind), file.path(out, "crc_independent.tsv"), sep = "\t"); fwrite(rbindlist(indsc), file.path(out, "crc_independent_scores.tsv"), sep = "\t")
fwrite(merge(merge(rbindlist(transfer)[, .(cohort, method, share_msi_captured, odds_ratio, fisher_p = p)], modsum, by = "cohort", all.x = TRUE),
             cms_sum, by = "cohort", all.x = TRUE), file.path(out, "crc_summary.tsv"), sep = "\t", na = "NA")
print(rbindlist(transfer)); print(rbindlist(igp)); print(rbindlist(ind)); print(modsum); print(cms_sum)
message("STEP 9 CRC REPLICATION DONE")
