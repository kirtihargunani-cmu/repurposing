#!/usr/bin/env Rscript
# Diagnostic: is the weak CollecTRI signal a regulon problem? Rerun group msviper with the precomputed,
# tissue-specific TCGA ARACNe networks (aracne.networks) on the same Step 4 signatures, same normal-split
# null and same filters (regulator expressed, >= 10 targets present, >= 50% of targets present).
# Full regulator set (TFs, co-factors, signaling). Output: results/step6/aracne_precomputed/{group_mr,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(11)
NET <- c("TCGA-PRAD" = "prad", "TCGA-LUSC" = "lusc", "BRCA-LumA" = "brca", "TCGA-KIRC" = "kirc", "TCGA-COAD" = "coad",
         "TCGA-THCA" = "thca", "TCGA-LUAD" = "luad", "TCGA-HNSC" = "hnsc", "TCGA-LIHC" = "lihc", "TCGA-STAD" = "stad",
         "TCGA-BLCA" = "blca", "TCGA-PAAD" = "paad")
N_PERM <- 1000; MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5
out_dir <- "results/step6/aracne_precomputed"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

nets <- readRDS("data/processed/step5/regulons_aracne.rds")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
tfs <- fread("results/step5/tf_universe.tsv")$gene_id
regs_meta <- fread("results/step5/aracne_regulators.tsv")
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
zg <- fread("results/step4/z_tumor_vs_gtex.tsv")

run <- function(gid) {
  z <- zg[!is.na(get(gid)), .(gene_id, v = get(gid))]; sig <- setNames(z$v, z$gene_id)
  d <- cov[group_id == gid & role %in% c("tumor", "gtex_normal") & sample %in% colnames(counts)]
  lc <- cpm(calcNormFactors(DGEList(counts[names(sig), d$sample])), log = TRUE, prior.count = 1)
  tum <- lc[, d[role == "tumor", sample]]; nor <- lc[, d[role == "gtex_normal", sample]]
  med <- apply(tum, 1, median)
  null <- sapply(seq_len(N_PERM), function(i) {
    lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(nor))); a <- nor[, lab]; b <- nor[, !lab]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(nor) - 2)
  }); rownames(null) <- names(sig)
  net <- nets[[paste0("aracne_", NET[[gid]])]]
  src_gene <- sym2ens[names(net)]
  keep <- vapply(seq_along(net), function(i) {
    g <- src_gene[i]; tg <- names(net[[i]]$tfmode); n_in <- sum(tg %in% names(sig))
    !is.na(g) && !is.na(med[g]) && med[g] >= MIN_LOGCPM && n_in >= MIN_TARGETS && n_in / length(tg) >= MIN_FRAC
  }, logical(1))
  reg <- lapply(net[keep], function(r) { k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  m <- msviper(sig, reg, nullmodel = null, minsize = MIN_TARGETS, verbose = FALSE)
  data.table(group_id = gid, network = NET[[gid]], regulator = names(m$es$nes), nes = m$es$nes, p = m$es$p.value,
             size = m$es$size)[, fdr := p.adjust(p, "BH")][, is_lambert_tf := sym2ens[regulator] %in% tfs][
             , is_kinase := regulator %in% regs_meta[is_kinase == TRUE, symbol]]
}
res <- rbindlist(mclapply(names(NET), function(g) tryCatch(run(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                          mc.cores = 6, mc.preschedule = FALSE))
fwrite(res, file.path(out_dir, "group_mr.tsv"), sep = "\t")
ct <- fread("results/step6/collectri/group_mr.tsv")[, .(ct_fdr10 = sum(fdr < 0.10), ct_fdr25 = sum(fdr < 0.25)), by = group_id]
s <- res[, .(tested = .N, fdr10 = sum(fdr < 0.10), fdr25 = sum(fdr < 0.25), fdr10_tf = sum(fdr < 0.10 & is_lambert_tf),
             fdr10_kinase = sum(fdr < 0.10 & is_kinase), max_abs_nes = round(max(abs(nes)), 2),
             top_up = paste(head(regulator[order(-nes)], 6), collapse = ","),
             top_down = paste(head(regulator[order(nes)], 4), collapse = ",")), by = group_id]
s <- merge(s, ct, by = "group_id")
fwrite(s, file.path(out_dir, "summary.tsv"), sep = "\t")
print(as.data.frame(s), row.names = FALSE)
message("ARACNE PRECOMPUTED TEST DONE")
