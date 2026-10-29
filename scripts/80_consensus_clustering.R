#!/usr/bin/env Rscript
# Step 8: consensus clustering of patients on regulator activity, with a null-model test for the number of clusters.
# Input: Step 8 patient NES = mappability filter + gene-level purity/RNA-quality correction + ARACNe VIPER
#        (data/processed/step6/aracne_project_map, script 61). PI decisions 2026-10-09:
#   groups     each of BRCA-Basal, BRCA-Her2, BRCA-LumA, BRCA-LumB, TCGA-COAD, TCGA-READ separately (no pooling)
#   features   main = one score per regulator module of >= 10 regulators (modules: |r| >= 0.7 between patient NES,
#              average linkage, as in script 58; score = PC1 of members, sign = mean member NES; each module z-scored);
#              sensitivity = every module incl. small ones and single regulators ("all")
#   method     M3C (John et al. 2020): consensus clustering (PAM, 100 resamples, 80% of patients each) for K = 2..10,
#              compared with 100 Monte Carlo reference datasets that keep the feature correlation but have no clusters
#              (reverse-PCA). Stability = entropy of the consensus matrix -> RCSI; p-values: Monte Carlo (empirical)
#              and beta-distribution fit (NORM_P).
#   decision   (fixed before running) best K = highest RCSI; clusters accepted only if both MONTECARLO_P and NORM_P of
#              the best K are < 0.05; otherwise "no significant structure" (accepted outcome, PI 2026-10-08).
# Early-warning confounder look for the best K (full checks are Step 9): plate (chi-square), purity covariate and
# RNA-quality index (Kruskal-Wallis), silhouette in feature space.
# Outputs: results/step8/{m3c_scores,decision,assignments,module_features}.tsv;
#          data/processed/step8/features/<group>_<set>.tsv.gz
# FEATURE_SETS=modules10,all (default both); ONLY_GROUPS; N_CORES (default 16); ITERS (default 100)
# EXCLUDE_MODULE=recurring: before building features, drop the group's recurring module (C8 / D26): the |r| >= 0.7
#   module (as in features()) holding the most C8 core members (>= 30% of groups, 62 v2 full set; >= 3 present),
#   plus any other core member; outputs to <step8>_no_c8 and <step8_features>_no_c8 (user 2026-10-28)
suppressPackageStartupMessages({library(data.table); library(edgeR); library(M3C); library(cluster)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); source("R/paths.R")
b <- build(); require_build(b, c("signatures", "patient_nes"))
GROUPS <- p_groups(b)   # every group with an ok Step 4 signature for this build
only <- Sys.getenv("ONLY_GROUPS"); if (nzchar(only)) GROUPS <- intersect(GROUPS, strsplit(only, ",")[[1]])
SETS <- strsplit(Sys.getenv("FEATURE_SETS", "modules10,all"), ",")[[1]]
N_CORES <- as.integer(Sys.getenv("N_CORES", "16")); ITERS <- as.integer(Sys.getenv("ITERS", "100"))
MOD_R <- 0.7; MIN_MOD <- 10; MAXK <- 10; ALPHA <- 0.05
NES_DIR <- p_patient_nes(b)
res_dir <- p_step8(b); feat_dir <- p_step8_features(b)
EXCL <- identical(Sys.getenv("EXCLUDE_MODULE"), "recurring")
if (EXCL) { core <- fread("results/step6/aracne_project_v2/diagnostics/recurring_module_members_v2.tsv")[frac_groups >= 0.3, regulator]
            res_dir <- paste0(res_dir, "_no_c8"); feat_dir <- paste0(feat_dir, "_no_c8") }
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE); dir.create(feat_dir, recursive = TRUE, showWarnings = FALSE)
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
qidx <- fread("results/step6/aracne_project/quality_correction_map/checks.tsv", nrows = 0)  # existence check only
qtab <- fread("results/step6/aracne_project/quality_correction_map/quality_index.tsv")[role == "tumor"]

features <- function(A, set) {
  mod <- cutree(hclust(as.dist(1 - abs(cor(A))), method = "average"), h = 1 - MOD_R)
  sz <- table(mod); keep <- as.integer(names(sz)[if (set == "modules10") sz >= MIN_MOD else sz >= 1])
  F <- sapply(keep, function(m) { S <- A[, mod == m, drop = FALSE]
    if (ncol(S) == 1) return(S[, 1])
    S <- sweep(S, 2, sign(cor(S, rowMeans(S)))[, 1], `*`); s <- prcomp(S, center = TRUE)$x[, 1]; s * sign(cor(s, rowMeans(S))) })
  hub <- vapply(keep, function(m) { S <- A[, mod == m, drop = FALSE]; if (ncol(S) == 1) return(colnames(S))
    colnames(S)[which.max(abs(cor(S, prcomp(S, center = TRUE)$x[, 1])))] }, "")
  colnames(F) <- sprintf("M%03d_%s", seq_along(keep), hub)
  memb <- rbindlist(lapply(seq_along(keep), function(i) data.table(feature = colnames(F)[i], regulator = colnames(A)[mod == keep[i]])))
  list(F = scale(F), members = memb)
}

run <- function(gid, set) {
  x <- fread(file.path(NES_DIR, paste0(gid, ".patient_nes.tsv.gz"))); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator
  n_drop <- 0L
  if (EXCL) {
    mod <- cutree(hclust(as.dist(1 - abs(cor(A))), method = "average"), h = 1 - MOD_R)
    hit <- table(mod[colnames(A) %in% core])
    m <- if (length(hit) && max(hit) >= 3) as.integer(names(hit)[which.max(hit)]) else NA_integer_
    drop <- colnames(A) %in% core | (!is.na(m) & mod == m)
    n_drop <- sum(drop); A <- A[, !drop, drop = FALSE]
    message(sprintf("%s: recurring module removed: %d regulators (%d core members present)", gid, n_drop, sum(x$regulator %in% core)))
  }
  fe <- features(A, set); F <- fe$F
  fwrite(data.table(sample = rownames(F), round(F, 4)), file.path(feat_dir, sprintf("%s_%s.tsv.gz", gid, set)), sep = "\t")
  t0 <- Sys.time()
  r <- M3C(as.data.frame(t(F)), cores = N_CORES, iters = ITERS, repsref = 100, repsreal = 100, maxK = MAXK,
           clusteralg = "pam", seed = 123, removeplots = TRUE, silent = TRUE)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  sc <- as.data.table(r$scores)[, `:=`(group_id = gid, feature_set = set)]
  best <- sc[which.max(RCSI)]
  ok <- best$MONTECARLO_P < ALPHA & best$NORM_P < ALPHA
  back <- setNames(rownames(F), make.names(rownames(F)))  # M3C returns syntactic names (TCGA-AA-0001-01 -> TCGA.AA.0001.01)
  stopifnot(!anyDuplicated(names(back)))
  asg <- rbindlist(lapply(2:MAXK, function(k) { a <- r$realdataresults[[k]]$ordered_annotation
    data.table(group_id = gid, feature_set = set, K = k, sample = unname(back[rownames(a)]),
               cluster = as.integer(as.character(a$consensuscluster))) }))
  stopifnot(!anyNA(asg$sample))
  meta <- cov[match(rownames(F), sample)]; b <- asg[K == best$K][match(rownames(F), sample), cluster]
  qi <- qtab[group_id == gid][match(rownames(F), sample), quality_index]
  pcol <- if (startsWith(gid, "LGG")) "absolute_purity" else "estimate_score"
  sil <- mean(silhouette(b, dist(F))[, 3])
  dec <- data.table(group_id = gid, feature_set = set, n_patients = nrow(F), n_features = ncol(F),
    n_regulators_covered = nrow(fe$members), n_removed_c8 = n_drop, best_K = best$K, RCSI = round(best$RCSI, 3), MONTECARLO_P = signif(best$MONTECARLO_P, 3),
    NORM_P = signif(best$NORM_P, 3), decision = if (ok) "clusters" else "no significant structure",
    cluster_sizes = paste(sort(table(b), decreasing = TRUE), collapse = "/"), silhouette = round(sil, 3),
    plate_p = signif(suppressWarnings(chisq.test(table(b, meta$batch), simulate.p.value = TRUE, B = 2000)$p.value), 3),
    purity_p = signif(kruskal.test(meta[[pcol]], factor(b))$p.value, 3),
    quality_p = signif(kruskal.test(qi, factor(b))$p.value, 3), minutes = mins)
  message(sprintf("%s [%s]: %d patients x %d features; best K = %d, RCSI %.3f, p = %.3g / %.3g -> %s (%.1f min)",
                  gid, set, nrow(F), ncol(F), best$K, best$RCSI, best$MONTECARLO_P, best$NORM_P, dec$decision, mins))
  list(sc = sc, dec = dec, asg = asg, memb = fe$members[, `:=`(group_id = gid, feature_set = set)])
}
res <- list()
for (set in SETS) for (gid in GROUPS) res[[paste(gid, set)]] <- tryCatch(run(gid, set),
  error = function(e) { message(gid, " [", set, "]: ", conditionMessage(e)); NULL })
pick <- function(k) rbindlist(lapply(res, `[[`, k), fill = TRUE)
fwrite(pick("sc"), file.path(res_dir, "m3c_scores.tsv"), sep = "\t")
fwrite(pick("dec"), file.path(res_dir, "decision.tsv"), sep = "\t")
fwrite(pick("asg"), file.path(res_dir, "assignments.tsv"), sep = "\t")
fwrite(pick("memb"), file.path(res_dir, "module_features.tsv"), sep = "\t")
message("STEP 8 DONE: ", length(res), " runs")
