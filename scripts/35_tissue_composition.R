#!/usr/bin/env Rscript
# Tissue composition of tumors vs normal tissue, per disease cohort (docs/SPEC_tissue_composition.md, D61).
# Samples per group (40 groups with a tumor-vs-GTEx signature): primary_per_patient tumors; GTEx normals of the v2
# reference (results/step4_v2/gtex_reference.tsv, use == TRUE; colon / kidney sub-region filter); TCGA adjacent
# normals (one per patient). Genes: v2 gene set. Expression: TMM log-CPM (prior 2) per group over all its samples.
#  1. Scores per sample
#     1a MCP-counter (Becht 2016; package + gene lists pinned, scripts/04b) on gene symbols: 10 immune / stromal populations
#     1b HPA v24 single-cell types, closely related types grouped into families (D67): marker of a unit = max nTPM over
#        its member types >= 10 and >= 4x the highest type outside the unit, gene in the v2 set; panel = up to 20
#        markers with the largest fold (>= 8 needed); score = mean log-CPM
#     + proliferation (Hallmark G2M mean log-CPM; not a cell type)
#  2. Differences per group x cell type: tumor vs GTEx, tumor vs adjacent, adjacent vs GTEx (adjacent arms need >= 10
#     samples). Age and sex removed within each arm first (linear; arm mean kept). Median difference, Cliff's delta,
#     Wilcoxon p, BH FDR over all tests. Also tumor vs GTEx after the v2 collection correction (GTEx scores adjusted to
#     the anchor: coll_* covariates of scripts/32 set to 0), and, for colon / kidney, tumor vs the unfiltered GTEx set.
#  3. Composition share of the signature: limma-voom tumor vs GTEx (v2 genes, Step 4 samples incl. the v2 GTEx
#     reference) with and without the ESTIMATE term (~ condition + age + sex [+ ESTIMATE] + tumor sample-type
#     indicators); z regressed on the HPA gene x cell-type specificity matrix (share of each gene's total nTPM per
#     type) by ridge regression; R^2 out-of-fold (10-fold CV, lambda chosen by inner CV); top cell types.
#  4. Positive controls fixed in the spec.
# Outputs: results/step4/composition/composition_{panels,differences,share,controls}.tsv, composition_scores.tsv.gz
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(MCPcounter); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); RNGkind("L'Ecuyer-CMRG"); set.seed(35)
OUT <- "results/step4/composition"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_CORES <- as.integer(Sys.getenv("N_CORES", "8")); MIN_ADJ <- 10; MIN_CLASS <- 3

genes <- fread("results/step3/genes_protein_coding_v23_v2.tsv"); sym <- setNames(genes$symbol, genes$gene_id)
counts <- readRDS("data/processed/step4/counts_protein_coding.rds"); counts <- counts[rownames(counts) %in% genes$gene_id, ]
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]
ref <- fread("results/step4_v2/gtex_reference.tsv")
summ <- fread("results/step4/signature_summary.tsv")[status == "ok" & contrast == "gtex"]
gids <- summ$group_id

# ---- 1b HPA panels
hpa <- fread("data/raw/composition/rna_single_cell_type.tsv", col.names = c("gene_id", "symbol", "cell_type", "ntpm"))
H <- dcast(hpa[gene_id %in% genes$gene_id], gene_id ~ cell_type, value.var = "ntpm", fun.aggregate = max, fill = 0)
Hm <- as.matrix(H[, -1]); rownames(Hm) <- H$gene_id
# Cell-type families (user 2026-10-18, D67): closely related HPA types share markers, so a single type rarely has
# genes >= 4x every other type. Each family is scored as one unit; its markers are compared with types OUTSIDE it.
FAMILIES <- list(
  "Enterocytes and Paneth cells" = c("Proximal enterocytes", "Distal enterocytes", "Paneth cells"),
  "Squamous epithelium and keratinocytes" = c("Basal keratinocytes", "Suprabasal keratinocytes", "Basal squamous epithelial cells", "Squamous epithelial cells"),
  "Alveolar cells" = c("Alveolar cells type 1", "Alveolar cells type 2"),
  "Airway basal, club and ionocytes" = c("Basal respiratory cells", "Club cells", "Ionocytes"),
  "Breast epithelium" = c("Breast glandular cells", "Breast myoepithelial cells"),
  "Distal nephron" = c("Distal tubular cells", "Collecting duct cells"),
  "Trophoblasts" = c("Cytotrophoblasts", "Syncytiotrophoblasts", "Extravillous trophoblasts"),
  "Endothelial cells (blood + lymphatic)" = c("Endothelial cells", "Lymphatic endothelial cells"),
  "Tissue macrophages" = c("Macrophages", "Kupffer cells", "Hofbauer cells", "Microglial cells"),
  "Smooth muscle and myoid cells" = c("Smooth muscle cells", "Peritubular cells"),
  "Fibroblasts and stromal cells" = c("Fibroblasts", "Endometrial stromal cells", "Ovarian stromal cells"),
  "Gonadal somatic cells" = c("Leydig cells", "Sertoli cells", "Granulosa cells"),
  "Gland secretory cells" = c("Mucus glandular cells", "Serous glandular cells", "Salivary duct cells"),
  "Endocrine cells (pancreatic + enteroendocrine)" = c("Pancreatic endocrine cells", "Enteroendocrine cells"),
  "Ductal cells and cholangiocytes" = c("Ductal cells", "Cholangiocytes"))
stopifnot(all(unlist(FAMILIES) %in% colnames(Hm)))
units <- c(FAMILIES, setNames(as.list(setdiff(colnames(Hm), unlist(FAMILIES))), setdiff(colnames(Hm), unlist(FAMILIES))))
cand <- rbindlist(lapply(names(units), function(u) { m <- units[[u]]
  inside <- apply(Hm[, m, drop = FALSE], 1, max); outside <- apply(Hm[, setdiff(colnames(Hm), m), drop = FALSE], 1, max)
  k <- inside >= 10 & inside >= 4 * outside
  data.table(gene_id = rownames(Hm)[k], cell_type = u, ntpm = inside[k], next_ntpm = outside[k], n_members = length(m)) }))
cand[, fold := ntpm / pmax(next_ntpm, 0.1)]
cand <- cand[order(cell_type, -fold)]
panels <- cand[, head(.SD, 20), by = cell_type][, if (.N >= 8) .SD, by = cell_type]
unit_tab <- data.table(cell_type = names(units), members = sapply(units, paste, collapse = "; "),
                       n_candidate_markers = sapply(names(units), function(u) cand[cell_type == u, .N]))
fwrite(merge(unit_tab, panels[, .(n_markers = .N, markers = paste(sym[gene_id], collapse = ","), min_fold = round(min(fold), 1)), by = cell_type],
             by = "cell_type", all.x = TRUE)[, `:=`(scored = !is.na(n_markers), hpa_release = "HPA v24")],
       file.path(OUT, "composition_panels.tsv"), sep = "\t")
spec_mat <- Hm / pmax(rowSums(Hm), 1e-9)  # gene x cell type specificity (share of total nTPM)
hm <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
g2m <- hm[[which(sapply(hm, `[`, 1) == "HALLMARK_G2M_CHECKPOINT")]][-(1:2)]
mcp_genes <- read.table("data/raw/composition/mcpcounter_genes.txt", sep = "\t", header = TRUE, colClasses = "character", check.names = FALSE)

samples_for <- function(gid, unfiltered = FALSE) {
  t <- cov[group_id == gid & role == "tumor" & sample %in% colnames(counts)]
  r <- ref[grepl(paste0("(^|;)", gid, "($|;)"), groups) & sample %in% colnames(counts)]
  n <- cov[group_id == gid & role == "gtex_normal" & sample %in% r[if (unfiltered) excluded %in% c(NA, "") | grepl("^sub-region", excluded) else use == TRUE, sample]]
  a <- cov[group_id == gid & role == "adjacent_normal" & sample %in% colnames(counts)]
  list(t = t, n = n, a = a, ref = r)
}
lcpm <- function(s) { y <- calcNormFactors(DGEList(counts[, s])); cpm(y, log = TRUE, prior.count = 2) }
score_samples <- function(e) {
  es <- e; rownames(es) <- sym[rownames(e)]; es <- es[!is.na(rownames(es)) & !duplicated(rownames(es)), ]
  mcp <- MCPcounter.estimate(es, featuresType = "HUGO_symbols", genes = mcp_genes)
  rownames(mcp) <- paste0("MCP: ", rownames(mcp))
  hp <- t(sapply(split(panels$gene_id, panels$cell_type), function(gs) colMeans(e[intersect(gs, rownames(e)), , drop = FALSE])))
  rownames(hp) <- paste0("HPA: ", rownames(hp))
  pr <- colMeans(es[intersect(g2m, rownames(es)), , drop = FALSE])
  rbind(mcp, hp, "Proliferation (G2M)" = pr)
}
resid_arm <- function(S, d) {  # remove age (+ sex) within one arm; keep the arm mean
  X <- cbind(1, as.numeric(scale(d$age))); if (uniqueN(d$sex) > 1) X <- cbind(X, as.numeric(d$sex == d$sex[1]))
  X[is.na(X)] <- 0; b <- lmFit(S, X)$coefficients; b[is.na(b)] <- 0; S - b[, -1, drop = FALSE] %*% t(X[, -1, drop = FALSE]) }
cliff <- function(x, y) { r <- rank(c(x, y)); nx <- length(x); ny <- length(y)
  U <- sum(r[seq_len(nx)]) - nx * (nx + 1) / 2; 2 * U / (nx * ny) - 1 }
compare <- function(S1, S2, lab) rbindlist(lapply(rownames(S1), function(k) { x <- S1[k, ]; y <- S2[k, ]
  data.table(cell_type = k, comparison = lab, n1 = length(x), n2 = length(y), median_1 = median(x), median_2 = median(y),
             median_diff = median(x) - median(y), cliff_delta = cliff(x, y), p = suppressWarnings(wilcox.test(x, y)$p.value)) }))

per_group <- function(gid) {
  S <- samples_for(gid); d <- rbind(S$t, S$n, S$a, fill = TRUE)
  extra <- NULL
  if (gid %in% c("TCGA-COAD", "TCGA-READ", "TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH")) {
    U <- samples_for(gid, unfiltered = TRUE); extra <- U$n[!sample %in% S$n$sample]; d <- rbind(d, extra, fill = TRUE) }
  e <- lcpm(d$sample); sc <- score_samples(e)
  arm <- function(dd) resid_arm(sc[, dd$sample, drop = FALSE], dd)
  T <- arm(S$t); N <- arm(S$n)
  res <- list(compare(T, N, "tumor_vs_gtex"))
  # collection-corrected GTEx (anchor = coll_* 0)
  rn <- S$ref[match(S$n$sample, sample)]
  Xc <- as.matrix(rn[, .(coll_isch, coll_fast, coll_slow, coll_nonbp)]); Xc <- Xc[, apply(Xc, 2, function(v) sd(v, na.rm = TRUE) > 0), drop = FALSE]
  if (ncol(Xc)) { X <- cbind(1, scale(S$n$age), Xc); X[is.na(X)] <- 0; b <- lmFit(sc[, S$n$sample, drop = FALSE], X)$coefficients; b[is.na(b)] <- 0
    Nc <- N - b[, colnames(Xc), drop = FALSE] %*% t(Xc); res[[2]] <- compare(T, Nc, "tumor_vs_gtex_collection_corrected") }
  if (nrow(S$a) >= MIN_ADJ) { A <- arm(S$a); res <- c(res, list(compare(T, A, "tumor_vs_adjacent"), compare(A, N, "adjacent_vs_gtex"))) }
  if (!is.null(extra) && nrow(extra)) { Nu <- arm(rbind(S$n, extra)); res <- c(res, list(compare(T, Nu, "tumor_vs_gtex_unfiltered"))) }
  role <- setNames(c(rep("tumor", nrow(S$t)), rep("gtex", nrow(S$n)), rep("adjacent", nrow(S$a)), rep("gtex_excluded", NROW(extra))), d$sample)
  list(diff = rbindlist(res)[, group_id := gid],
       scores = data.table(group_id = gid, sample = colnames(sc), role = role[colnames(sc)], t(round(sc, 4)), check.names = FALSE))
}
res <- mclapply(gids, function(g) tryCatch(per_group(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }), mc.cores = N_CORES)
res <- Filter(Negate(is.null), res)
D <- rbindlist(lapply(res, `[[`, "diff")); D[, fdr := p.adjust(p, "BH")]
num <- c("median_1", "median_2", "median_diff", "cliff_delta"); D[, (num) := lapply(.SD, round, 4), .SDcols = num]
setcolorder(D, c("group_id", "cell_type", "comparison")); fwrite(D, file.path(OUT, "composition_differences.tsv"), sep = "\t")
fwrite(rbindlist(lapply(res, `[[`, "scores"), fill = TRUE), file.path(OUT, "composition_scores.tsv.gz"), sep = "\t")

# ---- 3. composition share of the signature (with / without the ESTIMATE term)
ridge_r2 <- function(z, X, k = 10) {
  X <- scale(X); X[is.na(X)] <- 0; z <- z - mean(z); fold <- sample(rep(seq_len(k), length.out = length(z)))
  lam <- 10^seq(-1, 4, length.out = 11)
  fitb <- function(Xa, za, l) solve(crossprod(Xa) + l * diag(ncol(Xa)), crossprod(Xa, za))
  pred <- numeric(length(z))
  for (f in seq_len(k)) { tr <- fold != f
    inner <- sample(rep(1:5, length.out = sum(tr))); Xt <- X[tr, ]; zt <- z[tr]
    cvl <- sapply(lam, function(l) sum(sapply(1:5, function(i) { b <- fitb(Xt[inner != i, ], zt[inner != i], l); sum((zt[inner == i] - Xt[inner == i, ] %*% b)^2) })))
    b <- fitb(Xt, zt, lam[which.min(cvl)]); pred[!tr] <- X[!tr, ] %*% b }
  bfull <- fitb(X, z, lam[6])
  list(r2 = 1 - sum((z - pred)^2) / sum(z^2), coef = setNames(as.numeric(bfull), colnames(X)))
}
share_group <- function(gid) {
  S <- samples_for(gid); t <- S$t; cls <- t[, .N, by = sample_type_class]
  t <- t[sample_type_class == "primary" | sample_type_class %in% cls[N >= MIN_CLASS, sample_type_class]]
  d <- rbind(t, S$n, fill = TRUE)[!is.na(age) & !is.na(sex) & !is.na(estimate_score)]
  d[, `:=`(cond = factor(ifelse(role == "tumor", "t", "n"), c("n", "t")), age_z = as.numeric(scale(age)), est_z = as.numeric(scale(estimate_score)))]
  extra <- setdiff(unique(d[role == "tumor", sample_type_class]), "primary")
  for (k in extra) d[, (paste0("st_", k)) := as.integer(role == "tumor" & sample_type_class == k)]
  base <- c("cond", "age_z", if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= MIN_CLASS) "sex", if (length(extra)) paste0("st_", extra))
  rbindlist(lapply(c(TRUE, FALSE), function(withest) {
    X <- model.matrix(reformulate(c(base, if (withest) "est_z")), data = d)
    y <- DGEList(counts[, d$sample]); k <- filterByExpr(y, X); y <- calcNormFactors(y[k, , keep.lib.sizes = FALSE])
    fit <- eBayes(lmFit(voom(y, X), X)); tt <- topTable(fit, coef = "condt", number = Inf, sort.by = "none")
    z <- setNames(zscoreT(tt$t, df = fit$df.total), rownames(tt)); gs <- intersect(names(z), rownames(spec_mat))
    r <- ridge_r2(z[gs], spec_mat[gs, ]); o <- order(r$coef)
    data.table(group_id = gid, estimate_term = withest, n_tumor = sum(d$role == "tumor"), n_gtex = sum(d$role == "gtex_normal"),
               genes = length(gs), composition_share_r2 = round(r$r2, 4),
               top_up_in_tumor = paste(head(names(r$coef)[rev(o)], 5), collapse = "; "),
               top_up_in_normal = paste(head(names(r$coef)[o], 5), collapse = "; "))
  }))
}
SH <- rbindlist(Filter(Negate(is.null), mclapply(gids, function(g) tryCatch(share_group(g), error = function(e) { message(g, " share: ", conditionMessage(e)); NULL }), mc.cores = N_CORES)))
fwrite(SH, file.path(OUT, "composition_share.tsv"), sep = "\t")

# ---- 4. positive controls (fixed in the spec)
ctrl <- rbindlist(list(
  data.table(group_id = c("BRCA-Basal", "BRCA-Her2", "BRCA-LumA", "BRCA-LumB"), cell_type = "HPA: Adipocytes", comparison = "tumor_vs_gtex", expect = "lower"),
  data.table(group_id = c("UCEC-EEC", "UCEC-SER", "TCGA-UCS"), cell_type = "HPA: Smooth muscle and myoid cells", comparison = "tumor_vs_gtex", expect = "lower"),
  data.table(group_id = "TCGA-LUAD", cell_type = "HPA: Alveolar cells", comparison = "tumor_vs_gtex", expect = "lower"),
  data.table(group_id = c("TCGA-COAD", "TCGA-READ"), cell_type = "HPA: Enterocytes and Paneth cells", comparison = "tumor_vs_gtex", expect = "similar"),
  data.table(group_id = c("TCGA-COAD", "TCGA-READ"), cell_type = "HPA: Smooth muscle and myoid cells", comparison = "tumor_vs_gtex", expect = "smaller gap than unfiltered"),
  data.table(group_id = c("TCGA-KIRC", "TCGA-KIRP"), cell_type = "HPA: Proximal tubular cells", comparison = "tumor_vs_gtex", expect = "lower"),
  data.table(group_id = "TGCT-SEM", cell_type = c("MCP: T cells", "MCP: B lineage"), comparison = "tumor_vs_gtex", expect = "higher"),
  data.table(group_id = "TCGA-SKCM", cell_type = "HPA: Melanocytes", comparison = "tumor_vs_gtex", expect = "higher"),
  data.table(group_id = "TCGA-SKCM", cell_type = "HPA: Squamous epithelium and keratinocytes", comparison = "tumor_vs_gtex", expect = "lower")))
C <- merge(ctrl, D, by = c("group_id", "cell_type", "comparison"), all.x = TRUE)
unf <- D[comparison == "tumor_vs_gtex_unfiltered", .(group_id, cell_type, delta_unfiltered = cliff_delta)]
C <- merge(C, unf, by = c("group_id", "cell_type"), all.x = TRUE)
C[, pass := fcase(is.na(cliff_delta), NA,
                  expect == "lower", cliff_delta < 0 & fdr < 0.05, expect == "higher", cliff_delta > 0 & fdr < 0.05,
                  expect == "similar", abs(cliff_delta) < 0.33,
                  expect == "smaller gap than unfiltered", abs(cliff_delta) < abs(delta_unfiltered))]
fwrite(C, file.path(OUT, "composition_controls.tsv"), sep = "\t")
print(C[, .(group_id, cell_type, expect, cliff_delta, fdr = signif(fdr, 2), delta_unfiltered, pass)])
print(SH[, .(group_id, estimate_term, composition_share_r2)][order(group_id, estimate_term)])
message("COMPOSITION DONE: ", uniqueN(D$group_id), " groups, ", uniqueN(D$cell_type), " cell types; controls pass ", sum(C$pass, na.rm = TRUE), "/", nrow(C))
