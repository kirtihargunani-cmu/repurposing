#!/usr/bin/env Rscript
# Step 6e: RNA-quality + purity correction at the gene level, then patient-level VIPER (ARACNe regulons).
# Replaces the NES-level purity correction (56) as the Step 8 input (PI decision 2026-10-08).
#   1. Quality index per sample (within group, within role): mean z of long genes minus mean z of short genes
#      (top vs bottom genomic-span quintile of expressed protein-coding genes). Degraded RNA loses long-transcript
#      signal; this index explained the TCGA plate pattern (scripts/57_plate_quality.R). Compared with GTEx RIN and ischemic time in the normals.
#   2. Gene-level linear correction (Gallego Romero et al. 2014): each gene's logCPM is regressed, within the
#      group's tumors, on purity + quality index (purity = ESTIMATE score; ABSOLUTE for LGG) and, within the
#      normals, on the quality index; residual + gene mean kept. Patients without the purity covariate are dropped.
#   3. Patient VIPER exactly as in 50_viper.R (viperSignature vs the group's normals, z-score, 100 permutations;
#      same regulators kept by 50_viper.R's filter), on the corrected matrices.
#   4. Checks vs the purity-corrected NES (56): PCs 1-5 R^2 on quality index, plate (and chance), purity, ABSOLUTE;
#      effective dimensionality; patient concordance with the group's master regulators.
# Outputs: data/processed/step6/aracne_project_qc/<group>.patient_nes.tsv.gz
#          results/step6/aracne_project/quality_correction/{checks,patient_concordance,quality_index}.tsv
# ONLY_GROUPS=a,b (optional); N_CORES (default 16)
suppressPackageStartupMessages({library(data.table); library(edgeR); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(61)
N_PERM_PATIENT <- 100; MIN_TARGETS <- 10; TOP_MR <- 50; N_PC <- 5; N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
ABSOLUTE_GROUPS <- c("LGG-IDHmut-codel", "LGG-IDHmut-noncodel", "LGG-IDHwt")
out_nes <- "data/processed/step6/aracne_project_qc"; res_dir <- "results/step6/aracne_project/quality_correction"
dir.create(out_nes, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

arac <- readRDS("data/processed/step5/regulons_aracne_project.rds")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
summ <- fread("results/step4/signature_summary.tsv")[status == "ok"]
contrast_for <- summ[, .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
zmat <- list(gtex = fread("results/step4/z_tumor_vs_gtex.tsv"), adjacent = fread("results/step4/z_tumor_vs_adjacent.tsv"))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
filt_all <- fread("results/step6/aracne_project/regulator_filter.tsv")[keep == TRUE]
mr_all <- fread("results/step6/aracne_project/group_mr.tsv")
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))

qindex <- function(m) {  # genes x samples logCPM -> long-minus-short index per sample
  m <- m[apply(m, 1, median) >= 1 & rownames(m) %in% names(span), , drop = FALSE]
  z <- (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  colMeans(z[span[rownames(z)] >= q[2], , drop = FALSE]) - colMeans(z[span[rownames(z)] <= q[1], , drop = FALSE])
}
correct <- function(m, X) {  # regress each gene (row) on design X (samples x k, with intercept); keep residual + mean
  t(qr.resid(qr(X), t(m))) + rowMeans(m)
}
r2f <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  1 - sum(qr.resid(qr(model.matrix(~ f)), y)^2) / sum((y - mean(y))^2) }
r2c <- function(y, v) { ok <- !is.na(v); if (sum(ok) < 20) NA_real_ else cor(y[ok], v[ok])^2 }
pc_checks <- function(N, meta, qi, purity, label, gid) {  # N: patients x regulators
  p <- prcomp(N, center = TRUE); ve <- p$sdev^2 / sum(p$sdev^2); k <- min(N_PC, ncol(p$x))
  data.table(group_id = gid, version = label, pc = 1:k, var_explained = round(ve[1:k], 3),
    effective_dim = round(sum(ve)^2 / sum(ve^2), 2),
    r2_quality = round(sapply(1:k, function(i) r2c(p$x[, i], qi)), 3),
    r2_plate = round(sapply(1:k, function(i) r2f(p$x[, i], meta$batch)), 3),
    plate_chance = round((uniqueN(meta$batch) - 1) / (nrow(N) - 1), 3),
    r2_purity = round(sapply(1:k, function(i) r2c(p$x[, i], purity)), 3),
    r2_absolute = round(sapply(1:k, function(i) r2c(p$x[, i], meta$absolute_purity)), 3))
}

run <- function(gid, contrast) {
  normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  d <- cov[group_id == gid & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  pcol <- if (gid %in% ABSOLUTE_GROUPS) "absolute_purity" else "estimate_score"
  d <- d[role != "tumor" | !is.na(get(pcol))]
  z <- zmat[[contrast]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]; sig <- setNames(z[[gid]], z$gene_id)
  lc <- cpm(calcNormFactors(DGEList(counts[names(sig), d$sample])), log = TRUE, prior.count = 1)
  tum <- lc[, d[role == "tumor", sample], drop = FALSE]; nor <- lc[, d[role == normal_role, sample], drop = FALSE]
  mt <- d[role == "tumor"][match(colnames(tum), sample)]; mn <- d[role == normal_role][match(colnames(nor), sample)]
  qt <- qindex(tum); qn <- qindex(nor); pur <- mt[[pcol]]
  tumc <- correct(tum, cbind(1, pur - mean(pur), qt - mean(qt)))
  norc <- correct(nor, cbind(1, qn - mean(qn)))

  regs_all <- arac[[net_of[[gid]]]]
  reg <- lapply(regs_all[filt_all[group_id == gid, regulator]], function(r) {
    k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  vs <- viperSignature(tumc, norc, per = N_PERM_PATIENT, method = "zscore", verbose = FALSE)
  pn <- viper(vs, reg, minsize = MIN_TARGETS, verbose = FALSE)
  fwrite(data.table(regulator = rownames(pn), round(pn, 3)), file.path(out_nes, paste0(gid, ".patient_nes.tsv.gz")), sep = "\t")

  # checks: new (gene-level purity + quality) vs previous (NES-level purity correction, script 56)
  old <- fread(sprintf("data/processed/step6/aracne_project_purity/%s.patient_nes.tsv.gz", gid))
  O <- t(as.matrix(old[, -1])); colnames(O) <- old$regulator
  s <- intersect(rownames(O), colnames(pn)); r <- intersect(colnames(O), rownames(pn))
  meta <- mt[match(s, sample)]; qi <- qt[s]; pu <- meta[[pcol]]
  chk <- rbind(pc_checks(O[s, r], meta, qi, pu, "purity_only (56)", gid),
               pc_checks(t(pn[r, s]), meta, qi, pu, "purity+quality (61)", gid))
  chk[, purity_covariate := pcol]
  top <- head(mr_all[group_id == gid & p < 0.05][order(-abs(nes))], TOP_MR)[regulator %in% r]
  conc <- if (nrow(top) >= 10) rbindlist(lapply(list(`purity_only (56)` = t(O[s, top$regulator]),
                                                     `purity+quality (61)` = pn[top$regulator, s]), function(M)
    data.table(median_rho = round(median(apply(M, 2, cor, y = top$nes, method = "spearman")), 3),
               median_frac_hit = round(median(colMeans(sign(M) == sign(top$nes) & abs(M) > 1.96)), 3))), idcol = "version")[
    , group_id := gid] else NULL
  qtab <- rbind(data.table(group_id = gid, sample = names(qt), role = "tumor", quality_index = round(qt, 4)),
                data.table(group_id = gid, sample = names(qn), role = normal_role, quality_index = round(qn, 4)))
  rin <- if (normal_role == "gtex_normal" && sum(!is.na(mn$gtex_rin)) >= 20) round(cor(qn, mn$gtex_rin, use = "complete"), 3) else NA_real_
  isch <- if (normal_role == "gtex_normal" && sum(!is.na(mn$gtex_ischemic_min)) >= 20)
    round(cor(qn, mn$gtex_ischemic_min, use = "complete", method = "spearman"), 3) else NA_real_
  chk[, `:=`(n_tumor = ncol(tum), n_normal = ncol(nor), r_quality_vs_gtex_rin = rin, rho_quality_vs_gtex_ischemic = isch)]
  list(chk = chk, conc = conc, q = qtab)
}

only <- Sys.getenv("ONLY_GROUPS")
jobs <- if (nzchar(only)) contrast_for[group_id %in% strsplit(only, ",")[[1]]] else contrast_for
res <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(run(jobs$group_id[i], jobs$contrast[i]),
  error = function(e) list(err = data.table(group_id = jobs$group_id[i], error = conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE)
pick <- function(k) rbindlist(lapply(res, `[[`, k), fill = TRUE)
errs <- pick("err"); if (nrow(errs)) print(errs)
fwrite(pick("chk"), file.path(res_dir, "checks.tsv"), sep = "\t")
fwrite(pick("conc"), file.path(res_dir, "patient_concordance.tsv"), sep = "\t")
fwrite(pick("q"), file.path(res_dir, "quality_index.tsv"), sep = "\t")
message("QUALITY CORRECTION DONE: ", length(res) - nrow(errs), " groups ok, ", nrow(errs), " errors")
