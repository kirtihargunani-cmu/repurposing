#!/usr/bin/env Rscript
# Test: does a covariate-adjusted normal-split null sharpen the v2 group master regulators (assumption A3)?
# Group VIPER as 50_viper.R (VERSION=v2, ARACNe full set, filtered regulons, regulators kept by its filter), with
# the null built three ways:
#   A  current 50_viper.R v2: all of the group's GTEx normals in covariates_toil (incl. samples Step 4a excluded:
#      colon muscle wall, kidney medulla, missing collection data), no adjustment
#   B  spec 6a: Step 4a reference normals only (the normals the signature used), no adjustment
#   C  B + each gene's logCPM regressed, within the normals, on the Step 4 model's normal-side covariates before
#      the random splits: age, ESTIMATE, sex (if both sexes >= 3), collection covariates coll_* (GTEx contrast,
#      terms with variation). Rationale: the signature's condition coefficient is conditional on these covariates,
#      so its noise excludes their variation; the null should too. Tissue is NOT regressed (multi-tissue groups):
#      it is not in the Step 4 model, so composition sampling noise is part of the signature's noise as well.
# Calibration (B and C): normals split into disjoint halves H1 / H2; null from 500 random splits of H1, test
#   signature = H2 split in two (same statistic and adjustment); FDR < 0.05 calls are false positives
#   (expected ~0). 3 repeats.
# Decision rule (fixed before running): B replaces A regardless (it is what the spec prescribes; reported).
#   C is adopted over B if, over the test groups, (1) the median FDR < 0.05 count rises, (2) calibration false
#   positives at FDR < 0.05 stay ~0 (mean <= 5 per repeat and not above B's), and (3) known-biology controls
#   (51_compare_mr.R list) gain no wrong-direction call and lose no passes. Otherwise B only.
# Outputs: results/step6/null_adjust_test/{summary,calibration,controls,mr}.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
OUT <- "results/step6/null_adjust_test"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_PERM <- 1000; N_PERM_CAL <- 500; N_CAL <- 3; MIN_TARGETS <- 10
GROUPS <- strsplit(Sys.getenv("GROUPS", paste(
  "ESCA-SCC,TCGA-BLCA,TCGA-STAD,TCGA-UCS,CESC-AD,TCGA-ACC,TCGA-LIHC",   # few FDR < 0.05 MRs in v2
  "TCGA-COAD,TCGA-KIRC,BRCA-LumA,TCGA-LUAD,SARC-UPS", sep = ",")), ",")[[1]]  # strong / most affected by A vs B
COLL <- c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")

arac <- readRDS("data/processed/step5/regulons_aracne_project_v2.rds")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
summ <- fread("results/step4_v2/signature_summary.tsv")[status == "ok"]
contrast_for <- summ[, .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
zmat <- list(gtex = fread("results/step4_v2/z_tumor_vs_gtex.tsv"), adjacent = fread("results/step4_v2/z_tumor_vs_adjacent.tsv"))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
ref <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE, c("sample", COLL), with = FALSE]
filt <- fread("results/step6/aracne_project_v2/regulator_filter.tsv")[keep == TRUE]

ctl <- rbindlist(list(
  data.table(regulator = c("E2F1", "FOXM1", "MYBL2"), expected = 1, why = "proliferation", group_id = list(GROUPS)),
  data.table(regulator = c("ESR1", "FOXA1", "GATA3"), expected = 1, why = "luminal breast identity", group_id = list(c("BRCA-LumA", "BRCA-LumB"))),
  data.table(regulator = c("EPAS1", "HIF1A"), expected = 1, why = "VHL loss / hypoxia", group_id = list("TCGA-KIRC")),
  data.table(regulator = c("HNF4A", "HNF1A"), expected = -1, why = "loss of hepatocyte identity", group_id = list("TCGA-LIHC")),
  data.table(regulator = "MYC", expected = 1, why = "WNT-MYC", group_id = list(c("TCGA-COAD", "TCGA-READ")))
))[, .(group_id = unlist(group_id)), by = .(regulator, expected, why)][group_id %in% GROUPS]

welch_z <- function(m, lab) {
  a <- m[, lab, drop = FALSE]; b <- m[, !lab, drop = FALSE]
  t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
  zscoreT(pmin(pmax(t, -50), 50), df = ncol(a) + ncol(b) - 2)
}
split_null <- function(m, n) sapply(seq_len(n), function(i) welch_z(m, sample(rep(c(TRUE, FALSE), length.out = ncol(m)))))
adjust <- function(m, meta, contrast) {  # residual + gene mean, Step 4 normal-side covariates
  X <- cbind(age = meta$age - mean(meta$age), est = meta$estimate_score - mean(meta$estimate_score))
  if (uniqueN(meta$sex) > 1 && min(table(meta$sex)) >= 3) X <- cbind(X, male = as.numeric(meta$sex == "male"))
  if (contrast == "gtex") { C <- as.matrix(meta[, COLL, with = FALSE]); X <- cbind(X, C[, apply(C, 2, var) > 0, drop = FALSE]) }
  t(qr.resid(qr(cbind(1, X)), t(m))) + rowMeans(m)
}
viper_mr <- function(sig, reg, null) {
  mr <- msviper(sig, reg, nullmodel = null, minsize = MIN_TARGETS, verbose = FALSE)
  data.table(regulator = names(mr$es$nes), nes = mr$es$nes, p = mr$es$p.value)[, fdr := p.adjust(p, "BH")]
}

run <- function(gid) {
  set.seed(64)
  contrast <- contrast_for[group_id == gid, contrast]; normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  z <- zmat[[contrast]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]; sig <- setNames(z[[gid]], z$gene_id)
  regs_all <- arac[[net_of[[gid]]]]
  reg <- lapply(regs_all[filt[group_id == gid, regulator]], function(r) {
    k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  logcpm <- function(d) cpm(calcNormFactors(DGEList(counts[names(sig), d$sample])), log = TRUE, prior.count = 1)
  dA <- cov[group_id == gid & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  dB <- if (contrast == "gtex") merge(dA[role == "tumor" | sample %in% ref$sample], ref, by = "sample", all.x = TRUE) else dA
  lA <- logcpm(dA); lB <- logcpm(dB)
  nA <- lA[, dA[role == normal_role, sample]]; mB <- dB[role == normal_role]; nB <- lB[, mB$sample]
  okC <- !is.na(mB$age) & !is.na(mB$estimate_score) & !is.na(mB$sex)
  mC <- mB[okC]; nC <- adjust(nB[, okC], mC, contrast)
  nulls <- list(A = split_null(nA, N_PERM), B = split_null(nB, N_PERM), C = split_null(nC, N_PERM))
  mr <- rbindlist(lapply(names(nulls), function(a) cbind(arm = a, viper_mr(sig, reg, nulls[[a]]))))
  width <- sapply(nulls, function(n) median(apply(n, 1, sd)))  # per-gene null z spread
  cal <- rbindlist(lapply(seq_len(N_CAL), function(k) {
    h <- sample(rep(c(TRUE, FALSE), length.out = ncol(nB)))
    rbindlist(lapply(c("B", "C"), function(a) {
      m <- if (a == "B") nB else { x <- nB; x[, okC] <- nC; x[, okC, drop = FALSE] }
      hh <- if (a == "B") h else h[okC]
      s <- welch_z(m[, !hh, drop = FALSE], sample(rep(c(TRUE, FALSE), length.out = sum(!hh))))
      r <- viper_mr(s, reg, split_null(m[, hh, drop = FALSE], N_PERM_CAL))
      data.table(group_id = gid, arm = a, rep = k, fp_fdr05 = sum(r$fdr < 0.05), fp_fdr10 = sum(r$fdr < 0.10),
                 max_abs_nes = round(max(abs(r$nes)), 2))
    }))
  }))
  mr[, group_id := gid]
  list(mr = mr, cal = cal, width = data.table(group_id = gid, arm = names(width), null_z_sd = round(width, 3),
                                              n_normals = c(ncol(nA), ncol(nB), ncol(nC))))
}

res <- mclapply(GROUPS, function(g) tryCatch(run(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                mc.cores = length(GROUPS), mc.preschedule = FALSE)
pick <- function(k) rbindlist(lapply(res, `[[`, k))
mr <- pick("mr"); cal <- pick("cal"); width <- pick("width")
fwrite(mr[, .(group_id, arm, regulator, nes = round(nes, 3), p = signif(p, 3), fdr = signif(fdr, 3))], file.path(OUT, "mr.tsv"), sep = "\t")

base <- mr[arm == "A", .(group_id, regulator, nes_A = nes)]
summ_t <- merge(mr, base, by = c("group_id", "regulator"))[, .(
  n_fdr05 = sum(fdr < 0.05), n_fdr10 = sum(fdr < 0.10), max_abs_nes = round(max(abs(nes)), 2),
  rho_vs_A = round(cor(nes, nes_A, method = "spearman"), 4)), by = .(group_id, arm)]
summ_t <- merge(summ_t, width, by = c("group_id", "arm"))
fwrite(summ_t[order(group_id, arm)], file.path(OUT, "summary.tsv"), sep = "\t")
fwrite(cal, file.path(OUT, "calibration.tsv"), sep = "\t")
cv <- merge(ctl, mr[, .(group_id, regulator, arm, nes = round(nes, 2), fdr = signif(fdr, 2))], by = c("group_id", "regulator"))
cv[, verdict := fifelse(sign(nes) != expected, "wrong direction", fifelse(fdr < 0.10, "correct, passes", "correct direction, below cutoff"))]
fwrite(cv, file.path(OUT, "controls.tsv"), sep = "\t")

print(dcast(summ_t, group_id ~ arm, value.var = c("n_fdr05", "max_abs_nes", "null_z_sd")))
print(cal[, .(mean_fp_fdr05 = mean(fp_fdr05), max_fp_fdr05 = max(fp_fdr05), mean_fp_fdr10 = mean(fp_fdr10)), by = arm])
print(dcast(cv[, .N, by = .(arm, verdict)], verdict ~ arm, value.var = "N", fill = 0))
m05 <- summ_t[, .(med = median(n_fdr05)), by = arm]; print(m05)
message("NULL ADJUST TEST DONE: ", length(unique(mr$group_id)), " groups")
