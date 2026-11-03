#!/usr/bin/env Rscript
# Exploratory (2026-11-03): evaluate the lagged ARACNe networks (scripts/explore/lagged_aracne_inputs.py): regulator at
# 6 h -> each landmark gene's late change (24 h beyond its own 6 h value), across the drugs of one line.
# Link sign = Spearman of regulator row and target row of the input. Compared with the same line's contemporaneous drug
# network (D86), the lineage tumor network and DepMap baseline (overlap over chance), and checked against the line's D84
# knockdown / overexpression signatures (link agrees if the target moves as sign x perturbation direction; for the
# activity variant every regulator with signatures in the line, measured or not).
# Outputs results/lincs/pnet/lagged_{overlap,agreement,esr1_mcf7}.tsv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
res <- "results/lincs/pnet"; pl <- "data/processed/lincs"
LINES <- data.table(line = c("MCF7", "PC3", "A375"), tumor = c("BRCA", "PRAD", "SKCM"))
idm <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; e2s <- setNames(idm$symbol, idm$gene_id); s2e <- setNames(idm$gene_id, idm$symbol)
genes <- fread(file.path(pl, "pnet/MCF7_specific/all.exp.tsv"), select = 1)$gene; regs81 <- readLines(file.path(pl, "pnet/regulators.txt"))
signed <- function(netf, expf, pfx) {
  n <- fread(netf); setnames(n, c("reg", "target", "mi", "p")); n <- n[!grepl("^(R6|A6)_", target)]
  x <- fread(expf); m <- as.matrix(x[, -1]); rownames(m) <- x$gene
  rk <- t(apply(m[unique(c(n$reg, n$target)), , drop = FALSE], 1, rank))
  n[, rho := mapply(function(a, b) cor(rk[a, ], rk[b, ]), reg, target)]
  n[, regulator := sub(paste0("^", pfx, "_"), "", reg)]; n[, .(regulator, target, mi, mode = sign(rho), rho)] }
cont <- function(line) {   # contemporaneous drug network with signs
  signed(file.path(pl, sprintf("pnet/%s_specific/net_all/network.txt", line)), file.path(pl, sprintf("pnet/%s_specific/all.exp.tsv", line)), "NONE") }
tum <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
tumor_edges <- function(net, regset) rbindlist(lapply(names(tum[[net]]), function(r) { e <- if (startsWith(r, "ENSG")) r else s2e[r]
  if (is.na(e) || !e %in% regset) return(NULL); tm <- tum[[net]][[r]]$tfmode; k <- names(tm) %in% genes & names(tm) != e
  if (!any(k)) return(NULL); data.table(regulator = e, target = names(tm)[k], mode = sign(tm[k])) }))
dep <- fread(file.path(pl, "pnet_depmap/ALL/net_all/network.txt")); setnames(dep, c("regulator", "target", "mi", "p"))
key <- function(x) paste(x$regulator, x$target)

## benchmark signatures
gl <- fread(file.path(pl, "level5_bing_genes.tsv")); ok_g <- !is.na(gl$ensembl_id) & gl$ensembl_id != "" & !duplicated(gl$ensembl_id)
rdX <- function(f, n) { X <- matrix(readBin(file.path(pl, f), "numeric", size = 4, n = nrow(gl) * n), nrow = nrow(gl)); X <- X[ok_g, , drop = FALSE]; rownames(X) <- gl$ensembl_id[ok_g]; X }
mk <- fread(file.path(pl, "knockdown_bing_sigs.tsv.gz")); mo <- fread(file.path(pl, "overexpression_bing_sigs.tsv.gz"))
X <- cbind(rdX("knockdown_bing.f32", nrow(mk)), rdX("overexpression_bing.f32", nrow(mo))); colnames(X) <- c(mk$sig_id, mo$sig_id)
sg <- unique(fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[own_ok == TRUE, .(sig_id, cmap_name, direction, cell_iname)])
sg[, regulator := s2e[cmap_name]]
agree <- function(ed, line) { s <- sg[cell_iname == line & sig_id %in% colnames(X)]
  m <- merge(ed[, .(regulator, target, mode)], s, by = "regulator", allow.cartesian = TRUE)[target %in% rownames(X)]
  if (!nrow(m)) return(NULL); m[, z := X[cbind(match(target, rownames(X)), match(sig_id, colnames(X)))]]
  m[, .(rows = .N, signatures = uniqueN(sig_id), regulators = uniqueN(regulator), agree = round(mean(sign(z) == mode * direction), 3))] }

OV <- list(); AG <- list(); ES <- NULL
for (i in seq_len(nrow(LINES))) {
  ln <- LINES$line[i]; tn <- LINES$tumor[i]
  lm <- signed(file.path(pl, sprintf("pnet_lag_mrna/%s/net_all/network.txt", ln)), file.path(pl, sprintf("pnet_lag_mrna/%s/all.exp.tsv", ln)), "R6")
  la <- signed(file.path(pl, sprintf("pnet_lag_act/%s/net_all/network.txt", ln)), file.path(pl, sprintf("pnet_lag_act/%s/all.exp.tsv", ln)), "A6")
  ct <- cont(ln); tu <- tumor_edges(tn, regs81); tua <- tumor_edges(tn, unique(la$regulator))
  N <- length(regs81) * (length(genes) - 1)
  nets <- list(lagged_mrna = lm, contemporaneous_drug = ct, tumor = tu, depmap = dep)
  for (a in names(nets)) for (b in names(nets)) if (a < b) { A <- key(nets[[a]]); B <- key(nets[[b]]); k <- length(intersect(A, B))
    OV[[length(OV) + 1]] <- data.table(line = ln, a = a, b = b, edges_a = length(A), edges_b = length(B), shared = k, fold = round(k / (length(A) * length(B) / N), 2)) }
  Na <- uniqueN(la$regulator) * (length(genes) - 1)
  k <- length(intersect(key(la), key(tua)))
  OV[[length(OV) + 1]] <- data.table(line = ln, a = "lagged_act", b = "tumor (same regulators)", edges_a = nrow(la), edges_b = nrow(tua), shared = k, fold = round(k / (nrow(la) * nrow(tua) / Na), 2))
  # agreement, on the regulators each pair of networks shares with the benchmark
  for (nm in c("lagged_mrna", "contemporaneous_drug", "tumor")) { x <- agree(nets[[nm]], ln); if (!is.null(x)) AG[[length(AG) + 1]] <- cbind(line = ln, network = nm, x) }
  for (nm in c("lagged_act", "tumor_act_regs")) { x <- agree(if (nm == "lagged_act") la else tua[regulator %in% la$regulator], ln); if (!is.null(x)) AG[[length(AG) + 1]] <- cbind(line = ln, network = nm, x) }
  if (ln == "MCF7") ES <- la[regulator == s2e["ESR1"]][order(-mi)][, .(target = e2s[target], mi = round(mi, 3), rho = round(rho, 3))]
}
OV <- rbindlist(OV); AG <- rbindlist(AG)
fwrite(OV, file.path(res, "lagged_overlap.tsv"), sep = "\t"); fwrite(AG, file.path(res, "lagged_agreement.tsv"), sep = "\t")
if (!is.null(ES)) fwrite(ES, file.path(res, "lagged_esr1_mcf7.tsv"), sep = "\t")
options(width = 200); print(OV); print(AG)
cat("\nESR1 (MCF7, early activity -> late change), links:", nrow(ES), " positive", sum(ES$rho > 0), " negative", sum(ES$rho < 0), "\n"); print(head(ES, 30))
lum <- c("GATA3", "FOXA1", "ERBB4", "ANXA9", "HEBP1", "MINDY1", "SCCPDH", "SPDEF", "TBX3", "AR", "XBP1", "MYB", "GREB1", "CCND1", "NRIP1", "IGFBP4")
cat("luminal / estrogen-response genes among ESR1 lagged targets:\n"); print(ES[target %in% lum])
