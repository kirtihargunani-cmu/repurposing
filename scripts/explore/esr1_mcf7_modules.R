## Exploratory (2026-11-03): do ESR1's BRCA ARACNe targets split into identity vs estrogen-response modules across MCF7 drugs (ligand drugs left out)?
suppressPackageStartupMessages(library(data.table))
proc <- "data/processed/lincs"
pg <- fread(file.path(proc, "drug_line_profiles_genes.tsv")); idx <- fread(file.path(proc, "drug_line_profiles_index.tsv"))
lig <- fread("results/lincs/d84_ligands.tsv"); cmp <- fread("results/lincs/compounds.tsv")
for (arm in c("specific", "raw")) {
P <- matrix(readBin(file.path(proc, if (arm == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32"), "numeric", size = 4, n = nrow(pg) * nrow(idx)), nrow = nrow(pg))
rownames(P) <- pg$gene_id
cols <- which(idx$ModelID == "ACH-000019" & !idx$pert_id %in% cmp[cmap_name %in% lig$drug, pert_id])
Z <- P[, cols]; Z <- Z - rowMeans(Z)
r <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")$BRCA$ESR1
tg <- intersect(names(r$tfmode), rownames(Z)); s <- sign(r$tfmode[tg]); M <- Z[tg, ] * s
A <- matrix(colSums(M), nrow(M), ncol(M), byrow = TRUE) - M; A <- A - rowMeans(A)
loo <- rowSums(M * A) / sqrt(rowSums(M^2) * rowSums(A^2))
sym <- setNames(pg$symbol, pg$gene_id)
pc <- prcomp(t(Z[tg, ]))$rotation[, 1]
cat(sprintf("\n== %s: MCF7 %d drugs; ESR1 BRCA ARACNe targets on profiles %d (%d positive)\n", arm, ncol(Z), length(tg), sum(s > 0)))
cat(sprintf("targets agreeing with ARACNe sign vs rest-of-regulon activity: %.0f%% (|r|>0.1: agree %d, contradict %d)\n", 100 * mean(loo > 0), sum(loo > 0.1), sum(loo < -0.1)))
canon <- c("GREB1", "TFF1", "PGR", "MYB", "XBP1", "PDZK1", "CA12", "NRIP1", "IGFBP4", "SIAH2", "AREG", "CCND1")
lum <- c("GATA3", "FOXA1", "SPDEF", "AR", "TBX3", "ERBB4", "ANXA9", "SCCPDH", "MINDY1", "HEBP1")
show <- function(g) { e <- names(sym)[sym %in% g]; e <- intersect(e, rownames(Z)); data.table(gene = sym[e], in_aracne = e %in% tg, aracne_sign = ifelse(e %in% tg, s[e], NA), loo_r = ifelse(e %in% tg, round(loo[e], 2), NA), pc1 = ifelse(e %in% tg, round(pc[e] * sqrt(length(tg)), 2), NA)) }
cat("canonical estrogen-response genes:\n"); print(show(canon)); cat("luminal-identity genes:\n"); print(show(lum))
# correlation between canonical and luminal module scores across drugs
ce <- intersect(names(sym)[sym %in% canon], rownames(Z)); le <- intersect(names(sym)[sym %in% lum], rownames(Z))
cat(sprintf("across MCF7 drugs, corr(mean canonical, mean luminal) = %.2f\n", cor(colMeans(Z[ce, ]), colMeans(Z[le, ]))))
}
