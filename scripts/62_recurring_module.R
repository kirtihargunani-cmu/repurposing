#!/usr/bin/env Rscript
# Diagnostic: what is the regulator module that recurs across groups (PILRB, HSF4, CCNL2, AKAP8L, ...)?
# Biology (immune signaling, apoptosis, DNA repair, ...) or a technical axis (plate, RNA quality, GC content)?
# Input: ARACNe patient NES from NES_DIR (default aracne_project_qc; aracne_project_cqn -> outputs with suffix _cqn).
#   1. Per group: rebuild modules (as 58: |r| >= 0.7, average linkage); take the module with the most seed members
#      (>= 3). Core = regulators found in that module in >= 30% of the groups where it occurs.
#   2. Per group, Spearman of the module score (PC1 of members, sign = mean member NES) with
#        - 50 MSigDB Hallmark program scores (mean z of the program's genes in the purity + quality corrected tumor
#          expression, i.e. the data the scores come from), incl. interferon, inflammatory, TNFa/NF-kB, IL6/JAK/STAT3,
#          apoptosis, DNA repair, p53, UV response, unfolded protein response
#        - immune: ESTIMATE immune score; cytolytic activity (mean log GZMA, PRF1; Rooney et al. 2015)
#        - technical: plate (R^2 vs chance), RNA-quality index, GC index (mean z of high-GC minus low-GC genes, top vs
#          bottom quintile of longest-transcript GC; GC bias is a known RNA-seq artifact, Risso 2011, Hansen 2012),
#          both on uncorrected expression
#   3. Member features vs all regulators: transcript GC, length, expression, chromosome.
#   4. Targets: ARACNe targets of the core members vs all targets in each network: GC, length and Hallmark enrichment.
# Outputs: results/step6/aracne_project/diagnostics/recurring_module_{members,associations,features,targets}.tsv
# NES_DIR=data/processed/step6/aracne_project_v2_patient (v2, check C8): v2 regulons / regulator list / gene set,
#   ESTIMATE purity in every group; outputs results/step6/aracne_project_v2/diagnostics/recurring_module_*_v2.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(62)
SEED <- c("PILRB", "HSF4", "CCNL2", "AKAP8L", "ANKZF1", "MAPK8IP3", "ZNF692", "CENPT", "AGAP6", "ZNF276", "STK36", "GIGYF1")
MOD_R <- 0.7; ABSOLUTE_GROUPS <- c("LGG-IDHmut-codel", "LGG-IDHmut-noncodel", "LGG-IDHwt")
NES_DIR <- Sys.getenv("NES_DIR", "data/processed/step6/aracne_project_map"); out <- "results/step6/aracne_project/diagnostics"
TAG <- sub("^aracne_project_?", "", basename(NES_DIR)); SFX <- paste0("_", TAG)
V2 <- TAG == "v2_patient"
if (V2) { SFX <- "_v2"; out <- "results/step6/aracne_project_v2/diagnostics"; ABSOLUTE_GROUPS <- character(0) }
dir.create(out, recursive = TRUE, showWarnings = FALSE)

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE & role == "tumor"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
if (V2) counts <- counts[rownames(counts) %in% fread("results/step3/genes_protein_coding_v23_v2.tsv")$gene_id, ]
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol); ens2sym <- setNames(idmap$symbol, idmap$gene_id)
gcl <- fread("data/interim/gene_gc_v23.tsv"); gc <- setNames(gcl$gc_longest_tx, gcl$gene_id); txl <- setNames(gcl$longest_tx_len, gcl$gene_id)
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id)); chrom <- setNames(pm$chrom, sub("\\..*", "", pm$id))
hm <- lapply(strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t"), function(l) list(name = sub("^HALLMARK_", "", l[1]), genes = l[-(1:2)]))
hm <- setNames(lapply(hm, function(h) na.omit(unname(sym2ens[h$genes]))), sapply(hm, `[[`, "name"))
regs_meta <- fread(if (V2) "results/step5/aracne_regulators_v2.tsv" else "results/step5/aracne_regulators.tsv")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
arac <- readRDS(if (V2) "data/processed/step5/regulons_aracne_project_v2.rds" else "data/processed/step5/regulons_aracne_project.rds")

zrow <- function(m) (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6)
index <- function(z, val) { v <- val[rownames(z)]; ok <- !is.na(v); q <- quantile(v[ok], c(0.2, 0.8))
  colMeans(z[ok & v >= q[2], , drop = FALSE]) - colMeans(z[ok & v <= q[1], , drop = FALSE]) }
r2f <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  1 - sum(qr.resid(qr(model.matrix(~ f)), y)^2) / sum((y - mean(y))^2) }
sp <- function(a, b) { ok <- !is.na(a) & !is.na(b); if (sum(ok) < 20) NA_real_ else cor(a[ok], b[ok], method = "spearman") }

run <- function(f) {
  gid <- sub("\\.patient_nes\\.tsv\\.gz$", "", basename(f))
  x <- fread(f); A <- t(as.matrix(x[, -1])); colnames(A) <- x$regulator
  R <- cor(A); mod <- cutree(hclust(as.dist(1 - abs(R)), method = "average"), h = 1 - MOD_R)
  hits <- table(mod[intersect(SEED, colnames(A))]); if (!length(hits) || max(hits) < 3) return(NULL)
  m <- as.integer(names(hits)[which.max(hits)]); mem <- colnames(A)[mod == m]
  S <- A[, mem, drop = FALSE]; S <- sweep(S, 2, sign(cor(S, rowMeans(S)))[, 1], `*`)
  score <- prcomp(S, center = TRUE)$x[, 1]; score <- score * sign(cor(score, rowMeans(S)))
  meta <- cov[match(rownames(A), sample)]
  lc <- cpm(calcNormFactors(DGEList(counts[, rownames(A)])), log = TRUE, prior.count = 1)
  lc <- lc[apply(lc, 1, median) >= 1, ]
  zr <- zrow(lc)
  pcol <- if (gid %in% ABSOLUTE_GROUPS) "absolute_purity" else "estimate_score"
  qi <- index(zr[rownames(zr) %in% names(span), ], span); gci <- index(zr, gc)
  pur <- meta[[pcol]]; ok <- !is.na(pur)
  lcc <- lc; lcc[, ok] <- t(qr.resid(qr(cbind(1, pur[ok] - mean(pur[ok]), qi[ok] - mean(qi[ok]))), t(lc[, ok]))) + rowMeans(lc[, ok])
  zc <- zrow(lcc)
  hs <- sapply(hm, function(g) { g <- intersect(g, rownames(zc)); if (length(g) >= 10) colMeans(zc[g, , drop = FALSE]) else rep(NA, ncol(zc)) })
  cyt <- if (all(sym2ens[c("GZMA", "PRF1")] %in% rownames(lc))) colMeans(lc[sym2ens[c("GZMA", "PRF1")], ]) else rep(NA, ncol(lc))
  assoc <- rbind(
    data.table(feature = colnames(hs), kind = "hallmark", value = apply(hs, 2, sp, b = score)),
    data.table(feature = c("estimate_immune", "cytolytic_GZMA_PRF1", "quality_index", "gc_index", "purity"),
               kind = c("immune", "immune", "technical", "technical", "purity"),
               value = c(sp(score, meta$estimate_immune), sp(score, cyt), sp(score, qi), sp(score, gci), sp(score, pur))),
    data.table(feature = "plate_R2", kind = "technical", value = r2f(score, meta$batch)),
    data.table(feature = "plate_R2_chance", kind = "technical", value = (uniqueN(meta$batch) - 1) / (nrow(A) - 1)))
  assoc[, `:=`(group_id = gid, module_size = length(mem), value = round(value, 3))]
  list(members = data.table(group_id = gid, regulator = mem, module_size = length(mem)), assoc = assoc)
}
files <- Sys.glob(file.path(NES_DIR, "*.patient_nes.tsv.gz"))
res <- mclapply(files, function(f) tryCatch(run(f), error = function(e) { message(basename(f), ": ", conditionMessage(e)); NULL }),
                mc.cores = 12, mc.preschedule = FALSE)
members <- rbindlist(lapply(res, `[[`, "members")); assoc <- rbindlist(lapply(res, `[[`, "assoc"))
ng <- uniqueN(members$group_id)
freq <- members[, .(n_groups = .N), by = regulator][, frac_groups := round(n_groups / ng, 2)][order(-n_groups)]
core <- freq[frac_groups >= 0.3, regulator]
fwrite(freq, file.path(out, paste0("recurring_module_members", SFX, ".tsv")), sep = "\t")
fwrite(assoc[order(kind, feature, group_id)], file.path(out, paste0("recurring_module_associations", SFX, ".tsv")), sep = "\t")

# member features vs all ARACNe regulators
allreg <- regs_meta$symbol; feat <- function(syms, label) { e <- na.omit(sym2ens[syms])
  data.table(set = label, n = length(e), median_gc = round(median(gc[e], na.rm = TRUE), 3),
             median_tx_len = as.integer(median(txl[e], na.rm = TRUE)), median_span = as.integer(median(span[e], na.rm = TRUE)),
             frac_tf = round(mean(syms %in% regs_meta[is_lambert_tf == TRUE, symbol]), 2),
             frac_chr16_17_19_22 = round(mean(chrom[e] %in% c("chr16", "chr17", "chr19", "chr22"), na.rm = TRUE), 2),
             frac_chr19 = round(mean(chrom[e] == "chr19", na.rm = TRUE), 2)) }
features <- rbind(feat(core, "core members"), feat(allreg, "all ARACNe regulators"))
fwrite(features, file.path(out, paste0("recurring_module_features", SFX, ".tsv")), sep = "\t")

# targets of the core members in each network vs all targets: GC, length, Hallmark enrichment (hypergeometric)
tg <- rbindlist(lapply(names(arac), function(net) { reg <- arac[[net]]
  alltg <- unique(unlist(lapply(reg, function(r) names(r$tfmode))))
  ct <- unique(unlist(lapply(reg[intersect(core, names(reg))], function(r) names(r$tfmode))))
  if (length(ct) < 20) return(NULL)
  enr <- rbindlist(lapply(names(hm), function(h) { g <- intersect(hm[[h]], alltg); k <- sum(ct %in% g)
    data.table(net = net, hallmark = h, in_targets = k, expected = round(length(ct) * length(g) / length(alltg), 1),
               p = phyper(k - 1, length(g), length(alltg) - length(g), length(ct), lower.tail = FALSE)) }))
  enr[, `:=`(n_core_targets = length(ct), gc_core_targets = round(median(gc[ct], na.rm = TRUE), 3),
             gc_all_targets = round(median(gc[alltg], na.rm = TRUE), 3),
             len_core_targets = as.integer(median(txl[ct], na.rm = TRUE)), len_all_targets = as.integer(median(txl[alltg], na.rm = TRUE)))]
}))
tg[, fdr := p.adjust(p, "BH"), by = net]
fwrite(tg[order(net, p)], file.path(out, paste0("recurring_module_targets", SFX, ".tsv")), sep = "\t")
message(sprintf("RECURRING MODULE DONE: found in %d of %d groups; %d core members", ng, length(files), length(core)))
