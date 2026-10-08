#!/usr/bin/env Rscript
# Step 6 (+ mechanism heterogeneity): regulator activity per group and per patient.
#
# REGULON=collectri (default) | dorothea_ABC | aracne_project
#   aracne_project: project ARACNe regulons (scripts/44_aracne_regulons.R), one network per cohort - each group
#   uses its cohort's network (subtypes share it; COAD/READ -> COADREAD). Full regulator set, no target cap.
#   Run with ONLY_GROUPS limited to groups whose network is built.
# For each group with a Step 4 signature (tumor vs GTEx; BLCA: tumor vs adjacent):
#  1. Regulon filters (per group, recorded in regulator_filter.tsv):
#       expressed   regulator gene median logCPM >= MIN_LOGCPM in the group's tumors
#                   (complexes AP1 / NFKB: any member expressed)
#       size        >= MIN_TARGETS targets present in the group's gene universe AND
#                   >= MIN_FRAC of the regulon's mapped targets present
#  2. Group mechanism: msviper on the Step 4 z signature; null = 1,000 random splits of the group's normals
#     into two halves (Welch t on logCPM -> z): preserves gene-gene correlation, contains no tumor signal.
#  3. Patient mechanism: viperSignature(tumors vs the group's normals, z-score, 100 permutations) -> viper.
#  4. Concordance per patient over the group's master regulators (nominal p < 0.05, top TOP_MR by |NES|;
#     the normal-split null is conservative, so few regulators reach FDR < 0.05 - FDR is reported, not used):
#       rho       Spearman(patient NES, group NES)
#       frac_hit  fraction of MRs with |patient NES| > 1.96 in the group's direction
#  5. Axes of variation: PCA of patient NES (centered); variance explained, effective dimensionality
#     ((sum lambda)^2 / sum lambda^2), top-loading regulators, and R^2 of PCs 1-5 on Step 3 covariates.
# Outputs: results/step6/<REGULON>/{group_mr,regulator_filter,patient_concordance,axes,
#          axes_covariates,group_summary}.tsv; patient NES matrices in data/processed/step6/<REGULON>/
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper);
                                library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
REGULON <- Sys.getenv("REGULON", "collectri")
MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5; TOP_MR <- 50
N_PERM_GROUP <- 1000; N_PERM_PATIENT <- 100; N_CORES <- as.integer(Sys.getenv("N_CORES", "6"))
COMPLEXES <- list(AP1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"),
                  NFKB = c("NFKB1", "NFKB2", "RELA", "RELB", "REL"))
set.seed(1)

res_dir <- file.path("results/step6", REGULON); proc_dir <- file.path("data/processed/step6", REGULON)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE); dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

if (REGULON == "aracne_project") {
  arac <- readRDS("data/processed/step5/regulons_aracne_project.rds")
  grp <- fread("results/step1/disease_groups_toil.tsv")
  net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
  regs_for <- function(gid) { r <- arac[[net_of[[gid]]]]; if (is.null(r)) stop("no ARACNe network for ", net_of[[gid]]); r }
} else {
  regs_core <- readRDS("results/step5/regulons_core.rds")[[REGULON]]
  regs_for <- function(gid) regs_core
}
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
summ <- fread("results/step4/signature_summary.tsv")[status == "ok"]
contrast_for <- summ[, .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
zmat <- list(gtex = fread("results/step4/z_tumor_vs_gtex.tsv"), adjacent = fread("results/step4/z_tumor_vs_adjacent.tsv"))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")

run_group <- function(gid, contrast) {
  regs_all <- regs_for(gid)
  normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  d <- cov[group_id == gid & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  z <- zmat[[contrast]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]
  sig <- setNames(z[[gid]], z$gene_id)
  y <- calcNormFactors(DGEList(counts[names(sig), d$sample]))
  lc <- cpm(y, log = TRUE, prior.count = 1)
  tum <- lc[, d[role == "tumor", sample], drop = FALSE]; nor <- lc[, d[role == normal_role, sample], drop = FALSE]

  # 1. regulon filters
  med <- apply(tum, 1, median)
  expressed <- function(src) {
    if (src %in% names(COMPLEXES)) { g <- sym2ens[COMPLEXES[[src]]]; return(any(med[g[!is.na(g)]] >= MIN_LOGCPM, na.rm = TRUE)) }
    g <- if (startsWith(src, "ENSG")) src else sym2ens[src]; !is.na(g) && !is.na(med[g]) && med[g] >= MIN_LOGCPM
  }
  filt <- rbindlist(lapply(names(regs_all), function(src) {
    tg <- names(regs_all[[src]]$tfmode); n_in <- sum(tg %in% names(sig))
    data.table(regulator = src, n_targets = length(tg), n_present = n_in, frac_present = round(n_in / length(tg), 3),
               expressed = expressed(src))
  }))
  filt[, keep := expressed & n_present >= MIN_TARGETS & frac_present >= MIN_FRAC]
  filt[, reason := fifelse(keep, "kept", fifelse(!expressed, "not expressed",
                   fifelse(n_present < MIN_TARGETS, "too few targets present", "low fraction of targets present")))]
  reg <- lapply(regs_all[filt[keep == TRUE, regulator]], function(r) {
    k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"

  # 2. group mechanism
  # null = random halves of the normals compared with each other: keeps gene-gene correlation but carries
  # no tumor signal (tumor/normal label permutation leaks a scaled copy of the real signature into the null
  # when the effect is large and groups are unbalanced)
  perm_z <- sapply(seq_len(N_PERM_GROUP), function(i) {
    lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(nor)))
    a <- nor[, lab, drop = FALSE]; b <- nor[, !lab, drop = FALSE]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(a) + ncol(b) - 2)
  })
  rownames(perm_z) <- rownames(lc)
  mr <- msviper(sig, reg, nullmodel = perm_z, minsize = MIN_TARGETS, verbose = FALSE)
  mrt <- data.table(group_id = gid, regulator = names(mr$es$nes), nes = mr$es$nes, p = mr$es$p.value,
                    size = mr$es$size)[, fdr := p.adjust(p, "BH")][order(-abs(nes))]

  # 3. patient mechanism
  vs <- viperSignature(tum, nor, per = N_PERM_PATIENT, method = "zscore", verbose = FALSE)
  pn <- viper(vs, reg, minsize = MIN_TARGETS, verbose = FALSE)  # regulators x patients
  fwrite(data.table(regulator = rownames(pn), round(pn, 3)), file.path(proc_dir, paste0(gid, ".patient_nes.tsv.gz")), sep = "\t")

  # 4. concordance with the group mechanism
  top <- head(mrt[p < 0.05], TOP_MR)[regulator %in% rownames(pn)]
  conc <- if (nrow(top) >= 10) data.table(group_id = gid, sample = colnames(pn),
    rho = apply(pn[top$regulator, , drop = FALSE], 2, cor, y = top$nes, method = "spearman"),
    frac_hit = colMeans(sign(pn[top$regulator, , drop = FALSE]) == sign(top$nes) & abs(pn[top$regulator, , drop = FALSE]) > 1.96),
    n_mr = nrow(top)) else NULL

  # 5. axes of variation
  pc <- prcomp(t(pn), center = TRUE, scale. = FALSE)
  lam <- pc$sdev^2; ve <- lam / sum(lam)
  axes <- data.table(group_id = gid, pc = seq_len(min(10, length(ve))), var_explained = round(ve[seq_len(min(10, length(ve)))], 4),
    top_pos = sapply(seq_len(min(10, length(ve))), function(k) paste(names(sort(pc$rotation[, k], decreasing = TRUE))[1:5], collapse = ",")),
    top_neg = sapply(seq_len(min(10, length(ve))), function(k) paste(names(sort(pc$rotation[, k]))[1:5], collapse = ",")))
  meta <- d[role == "tumor"][match(colnames(pn), sample)]
  covs <- list(estimate_score = meta$estimate_score, age = meta$age, batch = factor(meta$batch),
               site = factor(meta$site), sample_type = factor(meta$sample_type_class))
  r2 <- rbindlist(lapply(1:min(5, ncol(pc$x)), function(k) rbindlist(lapply(names(covs), function(cn) {
    v <- covs[[cn]]; ok <- !is.na(v)
    if (sum(ok) < 10 || (is.factor(v) && nlevels(droplevels(v[ok])) < 2)) return(NULL)
    data.table(group_id = gid, pc = k, covariate = cn, r2 = round(summary(lm(pc$x[ok, k] ~ v[ok]))$r.squared, 3))
  }))))
  gs <- data.table(group_id = gid, contrast = contrast, n_tumor = ncol(tum), n_normal = ncol(nor),
                   n_regulators_kept = nrow(filt[keep == TRUE]), n_mr_fdr05 = mrt[fdr < 0.05, .N], n_mr_p05 = mrt[p < 0.05, .N],
                   effective_dim = round(sum(lam)^2 / sum(lam^2), 2), pc1_var = round(ve[1], 3),
                   n_pc_80pct = which(cumsum(ve) >= 0.8)[1],
                   median_rho = if (!is.null(conc)) round(median(conc$rho), 3) else NA_real_,
                   median_frac_hit = if (!is.null(conc)) round(median(conc$frac_hit), 3) else NA_real_,
                   frac_patients_rho_gt_0.5 = if (!is.null(conc)) round(mean(conc$rho > 0.5), 3) else NA_real_)
  list(filt = filt[, group_id := gid], mr = mrt, conc = conc, axes = axes, r2 = r2, summary = gs)
}

only <- Sys.getenv("ONLY_GROUPS")
jobs <- if (nzchar(only)) contrast_for[group_id %in% strsplit(only, ",")[[1]]] else contrast_for
res <- mclapply(seq_len(nrow(jobs)), function(i)
  tryCatch(run_group(jobs$group_id[i], jobs$contrast[i]),
           error = function(e) list(error = data.table(group_id = jobs$group_id[i], error = conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE)
pick <- function(k) rbindlist(lapply(res, `[[`, k), fill = TRUE)
errs <- pick("error"); if (nrow(errs)) print(errs)
for (k in c("filt", "mr", "conc", "axes", "r2", "summary")) {
  f <- c(filt = "regulator_filter", mr = "group_mr", conc = "patient_concordance", axes = "axes",
         r2 = "axes_covariates", summary = "group_summary")[[k]]
  fwrite(pick(k), file.path(res_dir, paste0(f, ".tsv")), sep = "\t")
}
message("VIPER DONE (", REGULON, "): ", length(res) - nrow(errs), " groups ok, ", nrow(errs), " errors")
