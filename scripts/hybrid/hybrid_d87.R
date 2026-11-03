#!/usr/bin/env Rscript
# D87 (rule fixed 2026-11-02 before running, docs/DECISIONS.md): hybrid drug-side regulons = CollecTRI targets with
# directions learned from LINCS drug x line profiles (scripts/79c).
# Learning: for each CollecTRI regulator scorable on the profile genes (>= 10 targets, >= 50% present, as 79d) and each
#   of its targets, the Pearson correlation across drugs between the target and the regulator's leave-one-out activity
#   (mean of the CollecTRI-signed z of its other targets). Own directions per line with >= 300 learning profiles,
#   pooled (weighted mean Fisher z, weights n - 3) over lines with >= 100 for every other line. Drugs in the D84 ligand
#   list are never used for learning.
# Variants: H1 flip CollecTRI's sign where the correlation is opposite with p < 0.01 (primary); H2 sign of the
#   correlation; H3 sign of the correlation as direction, |r| as weight. Unlearned links keep CollecTRI's sign.
# Check: D84 main-analysis signatures (knockdown + overexpression, own_ok), scored as 79n; ligand profiles as 79v.
# ARM=specific (default, the rule) | raw (report; feeds D85). N_CORES (default 20; <= 24 by the user's thread limit).
# Outputs: results/lincs/hybrid/{edges_<arm>.tsv.gz, signatures_<arm>.tsv.gz, ligands_<arm>.tsv}; regulons in
#          data/processed/lincs/hybrid/regulons_<arm>.rds (not in git). Decisions: scripts/hybrid/hybrid_d87_decide.R
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../.."))
setwd(root)
MIN_TARGETS <- 10; MIN_FRAC <- 0.5; OWN_N <- 300; POOL_N <- 100; P_FLIP <- 0.01
ARM <- Sys.getenv("ARM", "specific"); stopifnot(ARM %in% c("specific", "raw"))
N_CORES <- min(as.integer(Sys.getenv("N_CORES", "20")), 24L)
proc <- "data/processed/lincs"; res <- "results/lincs/hybrid"
dir.create(res, showWarnings = FALSE, recursive = TRUE); dir.create(file.path(proc, "hybrid"), showWarnings = FALSE)
t0 <- Sys.time(); tick <- function(...) message(sprintf("[%4.1f min] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), sprintf(...))

## drug profiles (genes x profiles)
pg <- fread(file.path(proc, "drug_line_profiles_genes.tsv")); idx <- fread(file.path(proc, "drug_line_profiles_index.tsv"))
pf <- if (ARM == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32"
P <- matrix(readBin(file.path(proc, pf), "numeric", size = 4, n = nrow(pg) * nrow(idx)), nrow = nrow(pg))
dimnames(P) <- list(pg$gene_id, paste(idx$pert_id, idx$ModelID, sep = "|"))
lig <- fread("results/lincs/d84_ligands.tsv")[profile == ARM]
cmp <- fread("results/lincs/compounds.tsv")
lig_ids <- cmp[cmap_name %in% lig$drug, unique(pert_id)]
learn <- idx[, .(col = .I, pert_id, ModelID)][!pert_id %in% lig_ids]
np <- learn[, .N, by = ModelID]
own <- np[N >= OWN_N, ModelID]; pooled <- np[N >= POOL_N, ModelID]
tick("profiles %s: %d genes x %d; learning on %d (ligand drugs left out: %d compounds); own lines %d, pooled over %d",
     ARM, nrow(P), ncol(P), nrow(learn), length(lig_ids), length(own), length(pooled))

## CollecTRI, links on the profile genes
cl <- readRDS("results/step5/regulons_core.rds")$collectri
lr <- Filter(function(x) { k <- names(x$tfmode) %in% rownames(P); sum(k) >= MIN_TARGETS && mean(k) >= MIN_FRAC }, cl)
lr <- lapply(lr, function(x) { k <- names(x$tfmode) %in% rownames(P); x$tfmode[k] })
tick("CollecTRI: %d regulators, %d scorable on profile genes (%d links)", length(cl), length(lr), sum(lengths(lr)))

edge_cor <- function(line) {
  Z <- P[, learn[ModelID == line, col], drop = FALSE]; Z <- Z - rowMeans(Z)
  rbindlist(lapply(names(lr), function(tf) {
    s <- lr[[tf]]; M <- Z[names(s), , drop = FALSE] * s            # signed targets
    A <- matrix(colSums(M), nrow(M), ncol(M), byrow = TRUE) - M      # leave-one-out signed sum
    A <- A - rowMeans(A); Zt <- M * s                                # Zt: unsigned target rows (already centred)
    r <- rowSums(Zt * A) / sqrt(rowSums(Zt^2) * rowSums(A^2))
    data.table(regulator = tf, target = names(s), prior = as.integer(s), r = r)
  }))[, `:=`(line = line, n = ncol(Z))]
}
E <- rbindlist(mclapply(pooled, edge_cor, mc.cores = min(N_CORES, length(pooled))))
tick("link correlations: %d rows", nrow(E))
E[!is.finite(r), r := NA_real_]
pool <- E[!is.na(r), .(r = tanh(sum((n - 3) * atanh(pmin(pmax(r, -0.999999), 0.999999))) / sum(n - 3)), n = sum(n - 3) + 3L, lines = .N),
          by = .(regulator, target, prior)][, line := "POOLED"]
E <- rbind(E[line %in% own, .(regulator, target, prior, r, n, line, lines = 1L)], pool, use.names = TRUE)
E[, p := 2 * pnorm(-abs(atanh(pmin(pmax(r, -0.999999), 0.999999)) * sqrt(n - 3)))]
E[, flip := !is.na(r) & sign(r) != prior & p < P_FLIP]
fwrite(E[, .(line, regulator, target, prior, r = round(r, 4), n, p = signif(p, 3), flip)], file.path(res, sprintf("edges_%s.tsv.gz", ARM)), sep = "\t")

## regulons per line key x variant: full CollecTRI with learned links re-signed
build <- function(key, variant) {
  e <- E[line == key]; setkey(e, regulator, target)
  out <- lapply(names(cl), function(tf) {
    x <- cl[[tf]]; if (!tf %in% names(lr)) return(x)
    m <- e[.(tf, names(x$tfmode)), on = .(regulator, target)]
    has <- !is.na(m$r)
    if (variant == "H1") x$tfmode[m$flip %in% TRUE] <- -x$tfmode[m$flip %in% TRUE]
    if (variant %in% c("H2", "H3")) x$tfmode[has & m$r != 0] <- sign(m$r[has & m$r != 0])
    if (variant == "H3") x$likelihood <- ifelse(has, abs(m$r), if (any(has)) mean(abs(m$r[has])) else 1)
    x })
  names(out) <- names(cl); out
}
keys <- c(own, "POOLED"); VARS <- if (ARM == "specific") c("H1", "H2", "H3") else "H1"
REG <- setNames(lapply(keys, function(k) setNames(lapply(VARS, function(v) build(k, v)), VARS)), keys)
REG$CollecTRI <- cl
saveRDS(REG, file.path(proc, "hybrid", sprintf("regulons_%s.rds", ARM)))
tick("regulons built: %d keys x %d variants", length(keys), length(VARS))

## benchmark: D84 main-analysis signatures (as 79n)
genes <- fread(file.path(proc, "level5_bing_genes.tsv"))
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])
rd <- function(f, n) { X <- matrix(readBin(file.path(proc, f), "numeric", size = 4, n = nrow(genes) * n), nrow = nrow(genes)); X[ok_g, , drop = FALSE] }
mk <- fread(file.path(proc, "knockdown_bing_sigs.tsv.gz")); mo <- fread(file.path(proc, "overexpression_bing_sigs.tsv.gz"))
X <- cbind(rd("knockdown_bing.f32", nrow(mk)), rd("overexpression_bing.f32", nrow(mo)))
rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- c(mk$sig_id, mo$sig_id)
d84 <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")
sg <- d84[set == "collectri" & own_ok == TRUE, .(sig_id, technology, direction, cmap_name, cell_iname, nes_d84 = nes)]
sg <- merge(sg, rbind(mk[, .(sig_id, ModelID)], mo[, .(sig_id, ModelID)]), by = "sig_id")
sg[, key := fifelse(ModelID %in% own, ModelID, "POOLED")]
X <- X[, sg$sig_id]

prep <- function(rg, genes, noself) {
  r <- lapply(names(rg), function(n) { tm <- rg[[n]]$tfmode; lk <- rg[[n]]$likelihood; e <- sym2ens[n]
    k <- names(tm) %in% genes & (!noself | is.na(e) | names(tm) != e)
    if (sum(k) < MIN_TARGETS || sum(k) / length(tm) < MIN_FRAC) return(NULL)
    list(tfmode = tm[k], likelihood = lk[k]) })
  names(r) <- names(rg); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r
}
nes_at <- function(M, reg, tfs) {   # NES of each column's own TF
  v <- viper(M, reg, method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  if (ncol(M) == 1) v <- matrix(v, ncol = 1, dimnames = list(names(v), colnames(M)))
  i <- match(tfs, rownames(v)); out <- rep(NA_real_, length(tfs)); ok <- !is.na(i); out[ok] <- v[cbind(i[ok], which(ok))]; out
}
jobs <- CJ(key = unique(sg$key), variant = c("CollecTRI", VARS), sorted = FALSE)
B <- rbindlist(mclapply(seq_len(nrow(jobs)), function(j) {
  k <- jobs$key[j]; v <- jobs$variant[j]; s <- sg[key == k]
  rg <- if (v == "CollecTRI") cl else REG[[k]][[v]]
  M <- X[, s$sig_id, drop = FALSE]; if (ncol(M) == 1) M <- cbind(M, M)    # viper needs 2 columns
  data.table(sig_id = s$sig_id, variant = v, nes = nes_at(M, prep(rg, rownames(X), TRUE), c(s$cmap_name, if (nrow(s) == 1) s$cmap_name))[seq_len(nrow(s))])
}, mc.cores = N_CORES, mc.preschedule = FALSE))
B <- merge(sg, B, by = "sig_id")
chk <- B[variant == "CollecTRI", cor(nes, nes_d84, use = "complete.obs")]
tick("benchmark: %d signatures; CollecTRI NES vs D84 run r = %.4f", nrow(sg), chk)
B[, variant := if (ARM == "raw") fifelse(variant == "CollecTRI", variant, paste0(variant, "_raw")) else variant]
fwrite(B[, .(sig_id, technology, direction, cmap_name, cell_iname, ModelID, key, variant, nes = round(nes, 4))],
       file.path(res, sprintf("signatures_%s.tsv.gz", ARM)), sep = "\t")

## ligand profiles (as 79v: drug profile NES, curated sets pruned as 79d, no self-removal)
lg <- unique(lig[, .(drug, tf, direction, ModelID)])
lg <- merge(lg, cmp[, .(drug = cmap_name, pert_id)], by = "drug", allow.cartesian = TRUE)
lg[, colname := paste(pert_id, ModelID, sep = "|")]; lg <- lg[colname %in% colnames(P)]
lg[, key := fifelse(ModelID %in% own, ModelID, "POOLED")]
LJ <- CJ(key = unique(lg$key), variant = c("CollecTRI", VARS), sorted = FALSE)
L <- rbindlist(mclapply(seq_len(nrow(LJ)), function(j) {
  k <- LJ$key[j]; v <- LJ$variant[j]; s <- lg[key == k]
  rg <- if (v == "CollecTRI") cl else REG[[k]][[v]]
  M <- P[, unique(s$colname), drop = FALSE]; if (ncol(M) == 1) M <- cbind(M, M)
  vv <- viper(M, prep(rg, rownames(P), FALSE), method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  i <- match(s$tf, rownames(vv)); c2 <- match(s$colname, colnames(vv))
  cbind(s, variant = v, nes = ifelse(is.na(i), NA_real_, vv[cbind(ifelse(is.na(i), 1L, i), c2)]))
}, mc.cores = N_CORES, mc.preschedule = FALSE))
L <- L[, .(nes = mean(nes)), by = .(drug, tf, direction, ModelID, key, variant)]
L[, variant := if (ARM == "raw") fifelse(variant == "CollecTRI", variant, paste0(variant, "_raw")) else variant]
fwrite(L[, .(arm = ARM, drug, tf, direction, ModelID, key, variant, nes = round(nes, 4))], file.path(res, sprintf("ligands_%s.tsv", ARM)), sep = "\t")
tick("ligands: %d drug x TF x line rows. DONE", uniqueN(L[, .(drug, tf, ModelID)]))
