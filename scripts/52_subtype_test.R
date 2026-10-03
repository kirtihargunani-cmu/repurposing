#!/usr/bin/env Rscript
# Experiment: does splitting a large group by published molecular subtype recover VIPER signal?
# For each group below, tumors are split by PanCanAtlas Subtype_Selected (independent of VIPER; subtypes with
# >= 30 patients). Each subtype gets its own limma-voom signature vs the group's normals (same genes and
# covariates as Step 4) and msviper with the same normal-split null and the same filtered CollecTRI regulon.
# Control: N_RANDOM random tumor subsets of the same size from the whole group (tests whether any gain comes
# from homogeneity rather than from the smaller sample).
# Output: results/step6/collectri/subtype_test.tsv, subtype_test_mr.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(7)
GROUPS <- c("TCGA-KIRC", "TCGA-HNSC", "TCGA-THCA", "TCGA-PRAD", "TCGA-SKCM", "TCGA-LIHC", "TCGA-LUSC", "TCGA-LUAD", "TCGA-COAD")
MIN_N <- 30; N_PERM <- 1000; N_RANDOM <- 3; CUT <- 0.10

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
zg <- fread("results/step4/z_tumor_vs_gtex.tsv")
regs_all <- readRDS("results/step5/regulons_core.rds")$collectri
filt <- fread("results/step6/collectri/regulator_filter.tsv")[keep == TRUE]
mr_all <- fread("results/step6/collectri/group_mr.tsv")
sub <- fread("data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz")
sub <- sub[as.integer(substr(sampleID, 14, 15)) < 10 & !is.na(Subtype_Selected) & !grepl("\\.NA$|\\.-$", Subtype_Selected)]
sub[, patient := substr(sampleID, 1, 12)]
sub <- unique(sub, by = "patient")

signature <- function(d, genes) {
  d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  d[, age_z := as.numeric(scale(age))][, estimate_z := as.numeric(scale(estimate_score))]
  cls <- d[role == "tumor", .N, by = sample_type_class]
  extra <- setdiff(cls[N >= 3, sample_type_class], c("primary", "new_primary"))
  d <- d[!(role == "tumor" & !sample_type_class %in% c("primary", "new_primary", extra))]
  for (k in extra) d[, (paste0("st_", k)) := as.integer(role == "tumor" & sample_type_class == k)]
  terms <- c("condition", "age_z", "estimate_z", if (length(extra)) paste0("st_", extra))
  if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= 3) terms <- c(terms, "sex")
  design <- model.matrix(reformulate(terms), data = d)
  design <- design[, !colnames(design) %in% nonEstimable(design), drop = FALSE]
  y <- calcNormFactors(DGEList(counts[genes, d$sample]))
  fit <- eBayes(lmFit(voom(y, design), design))
  setNames(zscoreT(fit$t[, "conditiontumor"], df = fit$df.total), genes)
}

run_group <- function(gid) {
  d0 <- cov[group_id == gid & role %in% c("tumor", "gtex_normal") & sample %in% colnames(counts)]
  genes <- zg[!is.na(get(gid)), gene_id]
  nor <- d0[role == "gtex_normal"]; tum <- d0[role == "tumor"]
  tum[, subtype := sub$Subtype_Selected[match(patient, sub$patient)]]
  lc <- cpm(calcNormFactors(DGEList(counts[genes, nor$sample])), log = TRUE, prior.count = 1)
  null <- sapply(seq_len(N_PERM), function(i) {
    lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(lc))); a <- lc[, lab]; b <- lc[, !lab]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(lc) - 2)
  }); rownames(null) <- genes
  reg <- lapply(regs_all[filt[group_id == gid, regulator]], function(r) {
    k <- names(r$tfmode) %in% genes; list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  score <- function(tum_sub, label, kind) {
    sig <- signature(rbind(tum_sub, nor, fill = TRUE), genes)
    m <- msviper(sig, reg, nullmodel = null, minsize = 10, verbose = FALSE)
    r <- data.table(group_id = gid, unit = label, kind = kind, n_tumor = nrow(tum_sub), regulator = names(m$es$nes),
                    nes = m$es$nes, p = m$es$p.value)[, fdr := p.adjust(p, "BH")]
    r
  }
  counts_by <- tum[!is.na(subtype), .N, by = subtype][N >= MIN_N]
  out <- list(score(tum, "all tumors (rerun)", "whole"))
  for (s in counts_by$subtype) {
    ts <- tum[subtype == s]
    out[[length(out) + 1]] <- score(ts, s, "subtype")
    for (k in seq_len(N_RANDOM)) out[[length(out) + 1]] <- score(tum[sample(.N, nrow(ts))], paste0(s, " | random ", k), "random")
  }
  rbindlist(out)
}

res <- rbindlist(mclapply(GROUPS, function(g) tryCatch(run_group(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                          mc.cores = length(GROUPS), mc.preschedule = FALSE))
fwrite(res, "results/step6/collectri/subtype_test_mr.tsv", sep = "\t")
summ <- res[, .(n_tumor = n_tumor[1], pass = sum(fdr < CUT), max_abs_nes = round(max(abs(nes)), 2),
                top = paste(head(regulator[order(-abs(nes))], 5), collapse = ",")), by = .(group_id, unit, kind)]
summ[, subtype := sub(" \\| random.*", "", unit)]
ctrl <- summ[kind == "random", .(random_pass_mean = mean(pass), random_max_nes_mean = round(mean(max_abs_nes), 2)), by = .(group_id, subtype)]
summ <- merge(summ[kind != "random"], ctrl, by = c("group_id", "subtype"), all.x = TRUE)[, subtype := NULL]
fwrite(summ[order(group_id, kind != "whole")], "results/step6/collectri/subtype_test.tsv", sep = "\t")
print(as.data.frame(summ[order(group_id, kind != "whole")]), row.names = FALSE)
message("SUBTYPE TEST DONE")
