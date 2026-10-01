#!/usr/bin/env Rscript
# Step 4: limma-voom tumor-vs-normal signatures for every included Step 1 group.
#
# Contrasts (per group):
#   gtex      tumor vs GTEx normal of the group's tissue(s)        - primary signature
#   adjacent  tumor vs TCGA/TARGET adjacent normal (>= MIN_NORMAL) - check on the GTEx contrast
# Samples: primary_per_patient rows only (one tumor per patient, one normal per donor).
# Model:   ~ condition + sex + age + ESTIMATE score + tumor sample-type indicators
#   - sex dropped when a group has one sex; age and ESTIMATE score standardized
#   - sample type: indicator per non-primary tumor class (metastatic / recurrent / post_treatment)
#     with >= MIN_CLASS tumors, so `condition` is primary tumor vs normal; tumors in rarer classes
#     are dropped. new_primary counts as primary.
#   - data source / plate are NOT in the model: in the gtex contrast they are collinear with condition.
# Genes: GENCODE v23 protein-coding, then edgeR::filterByExpr per contrast; TMM; voom.
# Signature value: moderated t -> signed z (limma::zscoreT, df.total).
#
# Outputs:
#   data/processed/step4/<group>__<contrast>.tsv.gz   full topTable + z (not in git)
#   results/step4/signature_summary.tsv               one row per group x contrast
#   results/step4/z_tumor_vs_gtex.tsv, z_tumor_vs_adjacent.tsv   gene x group z matrices
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_NORMAL <- 10
MIN_CLASS <- 3
N_CORES <- 24

cov <- fread("results/step3/covariates_toil.tsv")
groups <- fread("results/step1/disease_groups_toil.tsv")[status == "include"]
pc <- fread("results/step3/genes_protein_coding_v23.tsv")
cov <- cov[primary_per_patient == TRUE & group_id %in% groups$group_id]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]

cnt_f <- "data/raw/solid_toil/TcgaTargetGtex_gene_expected_count.gz"
cache <- "data/processed/step4/counts_protein_coding.rds"  # rebuilt when the sample set changes
dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
counts <- if (file.exists(cache)) readRDS(cache) else NULL
if (is.null(counts) || !all(intersect(unique(cov$sample), colnames(counts)) %in% colnames(counts)) ||
    length(setdiff(unique(cov$sample), colnames(counts))) > 13) {
  hdr <- names(fread(cmd = paste("zcat", cnt_f, "| head -1")))
  use <- intersect(unique(cov$sample), hdr)
  message(length(use), " of ", uniqueN(cov$sample), " samples in count matrix; reading...")
  x <- fread(cmd = paste("zcat", cnt_f), select = c(hdr[1], use))
  ids <- sub("\\..*", "", x[[1]])
  keep <- ids %in% pc$gene_id
  counts <- as.matrix(x[keep, -1])
  rownames(counts) <- ids[keep]
  rm(x); invisible(gc())
  counts <- round(pmax(2^counts - 1, 0))  # file is log2(expected_count + 1)
  saveRDS(counts, cache)
}
message(nrow(counts), " protein-coding genes x ", ncol(counts), " samples")
symbol <- setNames(pc$symbol, pc$gene_id)

dir.create("data/processed/step4", recursive = TRUE, showWarnings = FALSE)
dir.create("results/step4", recursive = TRUE, showWarnings = FALSE)

run <- function(gid, contrast) {
  normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  d <- cov[group_id == gid & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  out <- list(group_id = gid, contrast = contrast, n_tumor = sum(d$role == "tumor"),
              n_normal = sum(d$role == normal_role), status = "ok", note = "")
  if (out$n_normal < MIN_NORMAL) return(modifyList(out, list(status = "skipped",
                                                            note = sprintf("< %d normals", MIN_NORMAL))))
  # tumor sample-type classes
  cls <- d[role == "tumor", .N, by = sample_type_class]
  rare <- cls[sample_type_class != "primary" & N < MIN_CLASS, sample_type_class]
  dropped <- d[role == "tumor" & sample_type_class %in% rare, .N]
  d <- d[!(role == "tumor" & sample_type_class %in% rare)]
  extra <- setdiff(cls[N >= MIN_CLASS, sample_type_class], "primary")
  for (k in extra) d[, (paste0("st_", k)) := as.integer(role == "tumor" & sample_type_class == k)]
  # covariates; drop samples missing any used covariate
  n0 <- nrow(d)
  d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  terms <- c("condition", "age_z", "estimate_z", if (length(extra)) paste0("st_", extra))
  if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= MIN_CLASS) terms <- c(terms, "sex")
  d[, age_z := as.numeric(scale(age))]
  d[, estimate_z := as.numeric(scale(estimate_score))]
  design <- model.matrix(reformulate(terms), data = d)
  if (!is.fullrank(design)) {
    ne <- nonEstimable(design)
    design <- design[, !colnames(design) %in% ne, drop = FALSE]
    out$note <- paste("non-estimable dropped:", paste(ne, collapse = ","))
  }
  y <- DGEList(counts[, d$sample])
  keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  v <- voom(y, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "conditiontumor", number = Inf, sort.by = "none")
  tt <- data.table(gene_id = rownames(tt), symbol = symbol[rownames(tt)], tt)
  tt[, z := zscoreT(t, df = fit$df.total)]
  fwrite(tt, sprintf("data/processed/step4/%s__%s.tsv.gz", gid, contrast), sep = "\t")
  modifyList(out, list(
    n_tumor = sum(d$condition == "tumor"), n_normal = sum(d$condition == "normal"),
    n_dropped_missing_cov = n0 - nrow(d), n_dropped_rare_type = dropped,
    covariates = paste(setdiff(colnames(design), c("(Intercept)", "conditiontumor")), collapse = ","),
    n_genes = nrow(tt), n_de_fdr05_lfc1 = tt[adj.P.Val < 0.05 & abs(logFC) > 1, .N],
    n_up = tt[adj.P.Val < 0.05 & logFC > 1, .N], n_down = tt[adj.P.Val < 0.05 & logFC < -1, .N],
    z = list(setNames(tt$z, tt$gene_id))))
}

jobs <- CJ(group_id = groups$group_id, contrast = c("gtex", "adjacent"))
only <- Sys.getenv("ONLY_GROUPS")  # e.g. ONLY_GROUPS=TCGA-LUAD,BRCA-Basal for a quick test
if (nzchar(only)) jobs <- jobs[group_id %in% strsplit(only, ",")[[1]]]
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  tryCatch(run(jobs$group_id[i], jobs$contrast[i]),
           error = function(e) list(group_id = jobs$group_id[i], contrast = jobs$contrast[i],
                                    status = "error", note = conditionMessage(e)))
}, mc.cores = N_CORES)

summ <- rbindlist(lapply(res, function(r) { r$z <- NULL; as.data.table(r) }), fill = TRUE)
if (summ[status == "error", .N]) print(summ[status == "error", .(group_id, contrast, note)])
for (ct in c("gtex", "adjacent")) {
  zs <- Filter(function(r) identical(r$contrast, ct) && identical(r$status, "ok"), res)
  if (!length(zs)) next
  genes <- sort(unique(unlist(lapply(zs, function(r) names(r$z[[1]])))))
  zm <- data.table(gene_id = genes, symbol = symbol[genes])
  for (r in zs) zm[, (r$group_id) := round(r$z[[1]][genes], 3)]
  fwrite(zm, sprintf("results/step4/z_tumor_vs_%s.tsv", ct), sep = "\t", na = "NA")
}
# agreement between the two contrasts, where both exist
both <- Filter(function(r) identical(r$status, "ok"), res)
agree <- rbindlist(lapply(unique(sapply(both, `[[`, "group_id")), function(g) {
  a <- Filter(function(r) r$group_id == g, both)
  if (length(a) < 2) return(NULL)
  za <- a[[1]]$z[[1]]; zb <- a[[2]]$z[[1]]; gs <- intersect(names(za), names(zb))
  data.table(group_id = g, spearman_gtex_vs_adjacent = round(cor(za[gs], zb[gs], method = "spearman"), 3))
}))
if (nrow(agree)) summ <- merge(summ, agree, by = "group_id", all.x = TRUE)
setorder(summ, group_id, contrast)
fwrite(summ, "results/step4/signature_summary.tsv", sep = "\t")
message("SIGNATURES DONE: ", summ[status == "ok", .N], " ok, ", summ[status == "skipped", .N],
        " skipped, ", summ[status == "error", .N], " errors")
