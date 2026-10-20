#!/usr/bin/env Rscript
# Step 4 v2 acceptance checks (SPEC_rebuild_v2.md, "Acceptance checks"; reported whatever they show, D19).
#   C1  GC tilt: Spearman of (tumor-vs-GTEx z minus tumor-vs-adjacent z) with gene GC (longest transcript), per group
#       with both contrasts
#   C3  agreement of tumor-vs-GTEx with tumor-vs-adjacent z (Spearman), from the signature summaries
#   C4  collection signal left: Spearman of tumor-vs-GTEx z with the within-GTEx ischemic-time z of the same group
#       (results/step4/reference_diagnostics/z_ischemic.tsv, scripts/31)
#   C6  proliferation retained: Hallmark E2F targets / G2M checkpoint in tumor-vs-GTEx z (limma::cameraPR)
#   C11 colon / kidney: split-half reproducibility (N_SPLIT random halves of the patients; equal normal n per half for
#       both references) of tumor-vs-GTEx (v2: 4a reference, collection covariates, voomByGroup) and tumor-vs-adjacent
#       (v1 model, voomByGroup) signatures
#   C12 voomByGroup vs voom: whole-signature z Spearman and top-500 overlap per group (from the v2 summary; the
#       regulator part follows in Step 6)
# Every check is computed for v1 (results/step4/) and v2 (results/step4_v2/), on the genes both share.
# Outputs: results/step4_v2/checks/{c1_gc,c3_agreement,c4_ischemic,c6_proliferation,c11_splithalf,c12_voom,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
source("scripts/vendor/voomByGroup.R")
OUT <- "results/step4_v2/checks"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_SPLIT <- 10; N_CORES <- as.integer(Sys.getenv("N_CORES", "12")); MIN_CLASS <- 3
RNGkind("L'Ecuyer-CMRG"); set.seed(1)
med <- function(x) round(median(x, na.rm = TRUE), 3)

V <- list(v1 = "results/step4", v2 = "results/step4_v2")
Z <- lapply(V, function(d) list(gtex = fread(file.path(d, "z_tumor_vs_gtex.tsv")), adjacent = fread(file.path(d, "z_tumor_vs_adjacent.tsv"))))
S <- lapply(V, function(d) fread(file.path(d, "signature_summary.tsv"))[status == "ok"])
vec <- function(m, g) { x <- m[[g]]; names(x) <- m$gene_id; x[!is.na(x)] }
common <- function(...) Reduce(intersect, lapply(list(...), names))
sp <- function(a, b) { g <- common(a, b); round(cor(a[g], b[g], method = "spearman"), 3) }

# C1 GC tilt
gcl <- fread("data/interim/gene_gc_v23.tsv")[, gene_id := sub("\\..*", "", gene_id)]
gc <- setNames(gcl$gc_longest_tx, gcl$gene_id)
both <- intersect(intersect(names(Z$v1$adjacent), names(Z$v2$adjacent)), names(Z$v2$gtex))[-(1:2)]
c1 <- rbindlist(lapply(both, function(g) data.table(group_id = g, t(sapply(c("v1", "v2"), function(v) {
  d <- vec(Z[[v]]$gtex, g); a <- vec(Z[[v]]$adjacent, g); k <- common(d, a); sp(d[k] - a[k], gc) })))))
fwrite(c1, file.path(OUT, "c1_gc.tsv"), sep = "\t")

# C3 agreement of the two contrasts
c3 <- merge(S$v1[contrast == "gtex", .(group_id, v1 = spearman_gtex_vs_adjacent)],
            S$v2[contrast == "gtex", .(group_id, v2 = spearman_gtex_vs_adjacent)], by = "group_id")[!is.na(v1) | !is.na(v2)]
fwrite(c3, file.path(OUT, "c3_agreement.tsv"), sep = "\t")

# C4 ischemic-time signal left
zi <- fread("results/step4/reference_diagnostics/z_ischemic.tsv")
gg <- intersect(intersect(names(zi), names(Z$v2$gtex)), names(Z$v1$gtex))[-1]
gg <- setdiff(gg, "symbol")
c4 <- rbindlist(lapply(gg, function(g) { i <- vec(zi, g); if (length(i) < 1000) return(NULL)
  data.table(group_id = g, v1 = sp(vec(Z$v1$gtex, g), i), v2 = sp(vec(Z$v2$gtex, g), i)) }))
fwrite(c4, file.path(OUT, "c4_ischemic.tsv"), sep = "\t")

# C6 proliferation (Hallmark E2F / G2M)
hm <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
hm <- setNames(lapply(hm, `[`, -(1:2)), sapply(hm, `[`, 1))[c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")]
c6 <- rbindlist(lapply(setdiff(names(Z$v2$gtex), c("gene_id", "symbol")), function(g) rbindlist(lapply(c("v1", "v2"), function(v) {
  m <- Z[[v]]$gtex; if (!g %in% names(m)) return(NULL); z <- m[[g]]; names(z) <- m$symbol; z <- z[!is.na(z) & !is.na(names(z))]
  r <- as.data.table(cameraPR(z, ids2indices(hm, names(z))), keep.rownames = "hallmark")
  r[, .(group_id = g, version = v, hallmark, direction = Direction, fdr = signif(FDR, 3))] }))))
fwrite(c6, file.path(OUT, "c6_proliferation.tsv"), sep = "\t")

# C11 split-half, colon / kidney (same samples and model as 30_signatures.R VERSION=v2)
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]
ref <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE, .(sample, coll_isch, coll_fast, coll_slow, coll_nonbp)]
COLL <- c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")
pc <- fread("results/step3/genes_protein_coding_v23_v2.tsv")
counts <- readRDS("data/processed/step4/counts_protein_coding.rds"); counts <- counts[rownames(counts) %in% pc$gene_id, ]
fit_z <- function(d, coll) {
  d <- copy(d); d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  cls <- d[role == "tumor", .N, by = sample_type_class]
  d <- d[!(role == "tumor" & sample_type_class %in% cls[sample_type_class != "primary" & N < MIN_CLASS, sample_type_class])]
  extra <- setdiff(cls[N >= MIN_CLASS, sample_type_class], "primary")
  for (k in extra) d[, (paste0("st_", k)) := as.integer(role == "tumor" & sample_type_class == k)]
  d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  d[, age_z := as.numeric(scale(age))]; d[, estimate_z := as.numeric(scale(estimate_score))]
  terms <- c("condition", "age_z", "estimate_z", if (length(extra)) paste0("st_", extra))
  if (coll) terms <- c(terms, COLL[sapply(COLL, function(k) var(d[[k]]) > 0)])
  if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= MIN_CLASS) terms <- c(terms, "sex")
  X <- model.matrix(reformulate(terms), data = d); X <- X[, !colnames(X) %in% nonEstimable(X), drop = FALSE]
  y <- DGEList(counts[, d$sample]); y <- calcNormFactors(y[filterByExpr(y, X), , keep.lib.sizes = FALSE])
  f <- eBayes(lmFit(voomByGroup(y, group = d$condition, design = X, print = FALSE), X))
  setNames(zscoreT(f$t[, "conditiontumor"], df = f$df.total), rownames(f))
}
c11 <- rbindlist(mclapply(c("TCGA-COAD", "TCGA-READ", "TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH"), function(gid) {
  d <- cov[group_id == gid & sample %in% colnames(counts)]
  tum <- d[role == "tumor"]; adj <- d[role == "adjacent_normal"]
  gtx <- merge(d[role == "gtex_normal" & sample %in% ref$sample], ref, by = "sample")
  for (k in COLL) tum[, (k) := 0]
  n_norm <- floor(min(nrow(adj), nrow(gtx)) / 2)
  rbindlist(lapply(seq_len(N_SPLIT), function(r) {
    ht <- sample(rep(1:2, length.out = nrow(tum))); ha <- sample(nrow(adj)); hg <- sample(nrow(gtx))
    A <- list(adj[ha[1:n_norm]], adj[ha[(n_norm + 1):(2 * n_norm)]]); G <- list(gtx[hg[1:n_norm]], gtx[hg[(n_norm + 1):(2 * n_norm)]])
    za <- lapply(1:2, function(h) fit_z(rbind(tum[ht == h], A[[h]], fill = TRUE), FALSE))
    zg <- lapply(1:2, function(h) fit_z(rbind(tum[ht == h], G[[h]], fill = TRUE), TRUE))
    data.table(group_id = gid, rep = r, n_tumor_half = floor(nrow(tum) / 2), n_normal_half = n_norm,
               splithalf_adjacent = sp(za[[1]], za[[2]]), splithalf_gtex = sp(zg[[1]], zg[[2]]))
  }))
}, mc.cores = min(N_CORES, 5)))
fwrite(c11, file.path(OUT, "c11_splithalf.tsv"), sep = "\t")
c11s <- c11[, .(n_tumor_half = n_tumor_half[1], n_normal_half = n_normal_half[1], splithalf_gtex = med(splithalf_gtex),
                splithalf_adjacent = med(splithalf_adjacent), gtex_better_frac = mean(splithalf_gtex > splithalf_adjacent)), by = group_id]

# C12 voomByGroup vs voom
c12 <- S$v2[, .(group_id, contrast, rho_vs_voom, top500_vs_voom)]
fwrite(c12, file.path(OUT, "c12_voom.tsv"), sep = "\t")

summ <- data.table(
  check = c("C1 GC tilt, median (v1 -> v2)", "C3 GTEx vs adjacent agreement, median (v1 -> v2)",
            "C4 ischemic-time signal, median (v1 -> v2)", "C6 E2F / G2M up at FDR < 0.05, groups (v1 -> v2)",
            "C11 split-half GTEx vs adjacent, medians", "C12 z rho with voom, median (range)", "C12 top-500 overlap with voom, median (range)"),
  value = c(sprintf("%s -> %s (%d groups)", med(c1$v1), med(c1$v2), nrow(c1)),
            sprintf("%s -> %s (%d groups)", med(c3$v1), med(c3$v2), c3[!is.na(v2), .N]),
            sprintf("%s -> %s (%d groups)", med(c4$v1), med(c4$v2), nrow(c4)),
            paste(sapply(c("v1", "v2"), function(v) sprintf("%d / %d", c6[version == v & hallmark == "HALLMARK_E2F_TARGETS" & direction == "Up" & fdr < 0.05, .N],
                                                         c6[version == v & hallmark == "HALLMARK_G2M_CHECKPOINT" & direction == "Up" & fdr < 0.05, .N])), collapse = " -> "),
            c11s[, paste(sprintf("%s %s vs %s", sub("TCGA-", "", group_id), splithalf_gtex, splithalf_adjacent), collapse = "; ")],
            sprintf("%s (%s-%s)", med(c12$rho_vs_voom), min(c12$rho_vs_voom), max(c12$rho_vs_voom)),
            sprintf("%s (%s-%s)", med(c12$top500_vs_voom), min(c12$top500_vs_voom), max(c12$top500_vs_voom))))
fwrite(summ, file.path(OUT, "summary.tsv"), sep = "\t")
print(summ); print(c1); print(c4); print(c11s)
print(dcast(c6, group_id ~ version + hallmark, value.var = "fdr")[order(group_id)])
message("STEP 4 V2 CHECKS DONE")
