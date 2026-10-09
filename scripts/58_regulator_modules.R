#!/usr/bin/env Rscript
# Diagnostic: what are the modules of co-varying regulators, and does each have a "master" regulator?
# Input: ARACNe patient NES from NES_DIR (default data/processed/step6/aracne_project_qc = gene-level purity + RNA-quality
# corrected, the Step 8 input; data/processed/step6/aracne_project_purity = NES-level purity correction only).
# Modules: average-linkage clustering of regulators on 1 - |r| of patient NES, cut at |r| >= 0.7 (as in script 54).
# Per module (>= 10 regulators):
#   hub                  member most correlated with the module score (first PC of the members): the regulator that
#                        best represents the module
#   frac_tf / frac_kinase  share of members that are Lambert TFs / kinases
#   hub_targets_members  share of the other members that are in the hub's own ARACNe regulon (targets)
#   baseline_targets     same share for all regulators (chance level for that hub): if hub_targets_members is far
#                        above it, the hub sits upstream of the module (candidate master regulator)
#   shared_targets       mean Jaccard overlap of target sets between member pairs (high = members score alike
#                        because they share targets, i.e. redundant regulons)
#   r2_plate / r_length / r_proliferation  what the module score tracks: plate, RNA-quality length index,
#                        proliferation (mean NES of E2F1, FOXM1, MYBL2)
# Output: results/step6/aracne_project/diagnostics/regulator_modules_<qc|purity>.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(Matrix); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
NES_DIR <- Sys.getenv("NES_DIR", "data/processed/step6/aracne_project_map")
TAG <- sub("^aracne_project_?", "", basename(NES_DIR)); MOD_R <- 0.7; MIN_SIZE <- 10; PROLIF <- c("E2F1", "FOXM1", "MYBL2")
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
regs_meta <- fread("results/step5/aracne_regulators.tsv")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
arac <- readRDS("data/processed/step5/regulons_aracne_project.rds")
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))
r2f <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  1 - sum(qr.resid(qr(model.matrix(~ f)), y)^2) / sum((y - mean(y))^2) }

run <- function(f) {
  gid <- sub("\\.patient_nes\\.tsv\\.gz$", "", basename(f))
  x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator
  meta <- cov[match(rownames(A), sample)]; reg <- arac[[net_of[[gid]]]]
  lc <- cpm(calcNormFactors(DGEList(counts[, rownames(A)])), log = TRUE, prior.count = 1)
  lc <- lc[apply(lc, 1, median) >= 1 & rownames(lc) %in% names(span), ]
  z <- (lc - rowMeans(lc)) / apply(lc, 1, sd); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  lidx <- colMeans(z[span[rownames(z)] >= q[2], ]) - colMeans(z[span[rownames(z)] <= q[1], ])
  pro <- intersect(PROLIF, colnames(A)); pscore <- if (length(pro)) rowMeans(scale(A[, pro, drop = FALSE])) else rep(NA, nrow(A))

  R <- cor(A); mod <- cutree(hclust(as.dist(1 - abs(R)), method = "average"), h = 1 - MOD_R)
  # regulator x target incidence (regulators named by symbol; targets Ensembl)
  rn <- intersect(colnames(A), names(reg))
  tg <- lapply(reg[rn], function(r) names(r$tfmode)); alltg <- unique(unlist(tg))
  Mi <- sparseMatrix(i = rep(seq_along(rn), lengths(tg)), j = match(unlist(tg), alltg), x = 1,
                     dims = c(length(rn), length(alltg)), dimnames = list(rn, alltg))
  reg_ens <- sym2ens[colnames(A)]; reg_ens[is.na(reg_ens)] <- colnames(A)[is.na(reg_ens)]
  sizes <- sort(table(mod), decreasing = TRUE); big <- as.integer(names(sizes[sizes >= MIN_SIZE]))
  rbindlist(lapply(big, function(m) {
    mem <- colnames(A)[mod == m]
    S <- A[, mem, drop = FALSE]; S <- sweep(S, 2, sign(cor(S, S[, which.max(colSums(abs(R[mem, mem])))])), `*`)
    score <- prcomp(S, center = TRUE)$x[, 1]; rr <- abs(cor(S, score))[, 1]; hub <- mem[which.max(rr)]
    others <- setdiff(mem, hub); htg <- if (hub %in% rn) tg[[hub]] else character()
    hub_in <- mean(reg_ens[others] %in% htg)
    base <- mean(reg_ens[setdiff(colnames(A), hub)] %in% htg)
    mm <- intersect(mem, rn); jac <- NA_real_
    if (length(mm) >= 2) { sub <- Mi[mm, , drop = FALSE]; inter <- as.matrix(tcrossprod(sub)); sz <- rowSums(sub)
      un <- outer(sz, sz, "+") - inter; J <- inter / un; jac <- mean(J[upper.tri(J)]) }
    data.table(group_id = gid, module = m, size = length(mem),
      frac_tf = round(mean(mem %in% regs_meta[is_lambert_tf == TRUE, symbol]), 2),
      frac_kinase = round(mean(mem %in% regs_meta[is_kinase == TRUE, symbol]), 2),
      hub = hub, hub_is_tf = hub %in% regs_meta[is_lambert_tf == TRUE, symbol], hub_r = round(max(rr), 2),
      hub_targets_members = round(hub_in, 3), baseline_targets = round(base, 3), shared_targets = round(jac, 3),
      r2_plate = round(r2f(score, meta$batch), 2), r_length = round(cor(score, lidx), 2),
      r_proliferation = round(cor(score, pscore), 2),
      members_top = paste(head(mem[order(-rr)], 10), collapse = ","))
  }))
}
files <- Sys.glob(file.path(NES_DIR, "*.patient_nes.tsv.gz"))
res <- rbindlist(mclapply(files, function(f) tryCatch(run(f), error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL }),
                          mc.cores = 8, mc.preschedule = FALSE))
fwrite(res[order(group_id, -size)], sprintf("results/step6/aracne_project/diagnostics/regulator_modules_%s.tsv", TAG), sep = "\t")
message("REGULATOR MODULES DONE: ", uniqueN(res$group_id), " groups, ", nrow(res), " modules")
