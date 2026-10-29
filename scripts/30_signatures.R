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
# Genes: GENESET=protein_coding (default; used for drug matching) or GENESET=extended
#   (protein-coding + lincRNA + antisense, for mechanism work), then edgeR::filterByExpr per
#   contrast; TMM; voom. Extended outputs go to *_extended directories.
# Normals: >= MIN_NORMAL per contrast, unless the group's min_normals (config/subtype_splits.tsv) is lower.
# Signature value: moderated t -> signed z (limma::zscoreT, df.total).
#
# Outputs:
#   data/processed/step4/<group>__<contrast>.tsv.gz   full topTable + z (not in git)
#   results/step4/signature_summary.tsv               one row per group x contrast
#   results/step4/z_tumor_vs_gtex.tsv, z_tumor_vs_adjacent.tsv   gene x group z matrices
#
# VERSION=v2 (rebuild spec SPEC_rebuild_v2.md, Step 4b; protein-coding only):
#   genes      v2 gene set (results/step3/genes_protein_coding_v23_v2.tsv; no mappability < 0.8, no readthrough; D46)
#   GTEx       only the Step 4a reference samples (results/step4_v2/gtex_reference.tsv, use == TRUE: sub-region
#              filter for colon / kidney, collection data present; D59)
#   model      gtex contrast adds the collection covariates coll_isch, coll_fast, coll_slow, coll_nonbp (centered per
#              GTEx tissue at its anchor in 4a; 0 for tumors; terms without variation in the group are left out), so
#              the condition coefficient is tumor vs a GTEx donor at the anchor (D57). Adjacent contrast: v1 model.
#   variance   voomByGroup (scripts/vendor/voomByGroup.R, group = tumor / normal; D69) in both contrasts; voom
#              is also fitted and kept for check C12 (columns logFC_voom / t_voom / z_voom per group file;
#              z_tumor_vs_<contrast>_voom.tsv; agreement columns in the summary)
#   outputs    data/processed/step4_v2/, results/step4_v2/
# STRATA=<file> (VERSION=v2 only; Step 10 v2, D72): one signature per subgroup in the file (columns stratum, group_id,
#   sample; results/step10_v2/stratum_patients.tsv): the subgroup's tumors vs the parent group's normals, same model
#   and normals as the parent; outputs data/processed/step10_v2/, results/step10_v2/ (group_id column = stratum)
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_NORMAL <- 10
MIN_CLASS <- 3
N_CORES <- as.integer(Sys.getenv("N_CORES", "24"))

cov <- fread("results/step3/covariates_toil.tsv")
groups <- fread("results/step1/disease_groups_toil.tsv")[status == "include"]
GENESET <- Sys.getenv("GENESET", "protein_coding")
stopifnot(GENESET %in% c("protein_coding", "extended"))
VERSION <- Sys.getenv("VERSION", "v1")
stopifnot(VERSION %in% c("v1", "v2"), VERSION == "v1" || GENESET == "protein_coding")
V2 <- VERSION == "v2"
STRATA <- Sys.getenv("STRATA"); ST <- nzchar(STRATA); stopifnot(!ST || V2)
if (ST) { sp <- fread(STRATA); parent_of <- setNames(unique(sp[, .(stratum, group_id)])$group_id, unique(sp[, .(stratum, group_id)])$stratum) }
sfx <- if (GENESET == "extended") "_extended" else ""
pc <- fread(sprintf("results/step3/genes_%s_v23.tsv", GENESET))
if (V2) {
  source("scripts/vendor/voomByGroup.R")
  pc <- fread("results/step3/genes_protein_coding_v23_v2.tsv")
  ref <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE, .(sample, anchor_class, coll_isch, coll_fast, coll_slow, coll_nonbp)]
  COLL <- c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")
}
cov <- cov[primary_per_patient == TRUE & group_id %in% groups$group_id]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]

cnt_f <- "data/raw/solid_toil/TcgaTargetGtex_gene_expected_count.gz"
cache <- sprintf("data/processed/step4%s/counts_%s.rds", sfx, GENESET)  # rebuilt when the sample set changes
dir.create(dirname(cache), recursive = TRUE, showWarnings = FALSE)
hdr <- names(fread(cmd = paste("zcat", cnt_f, "| head -1")))
use <- intersect(unique(cov$sample), hdr)
counts <- if (file.exists(cache)) readRDS(cache) else NULL
if (is.null(counts) || !all(use %in% colnames(counts))) {
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
if (V2) counts <- counts[rownames(counts) %in% pc$gene_id, ]
message(nrow(counts), " ", GENESET, " genes x ", ncol(counts), " samples (", VERSION, ")")
symbol <- setNames(pc$symbol, pc$gene_id)

vsfx <- if (V2) "_v2" else ""
proc_dir <- paste0("data/processed/step4", sfx, vsfx); res_dir <- paste0("results/step4", sfx, vsfx)
if (ST) { proc_dir <- "data/processed/step10_v2"; res_dir <- "results/step10_v2" }
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
min_normals <- setNames(suppressWarnings(as.integer(groups$min_normals)), groups$group_id)

run <- function(gid, contrast) {
  normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  par <- if (ST) parent_of[[gid]] else gid
  d <- cov[group_id == par & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  if (ST) d <- d[role != "tumor" | sample %in% sp[stratum == gid, sample]]
  if (V2 && contrast == "gtex") {  # Step 4a reference samples + collection covariates (0 for tumors)
    d <- merge(d[role == "tumor" | sample %in% ref$sample], ref, by = "sample", all.x = TRUE)
    for (k in COLL) d[role == "tumor", (k) := 0]
  }
  d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  out <- list(group_id = gid, contrast = contrast, n_tumor = sum(d$role == "tumor"),
              n_normal = sum(d$role == normal_role), status = "ok", note = "")
  min_n <- if (is.na(min_normals[par])) MIN_NORMAL else min(min_normals[par], MIN_NORMAL)
  if (out$n_normal < min_n) return(modifyList(out, list(status = "skipped",
                                                       note = sprintf("< %d normals", min_n))))
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
  if (V2 && contrast == "gtex") terms <- c(terms, COLL[sapply(COLL, function(k) var(d[[k]]) > 0)])
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
  v <- if (V2) voomByGroup(y, group = d$condition, design = design, print = FALSE) else voom(y, design)
  fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = "conditiontumor", number = Inf, sort.by = "none")
  tt <- data.table(gene_id = rownames(tt), symbol = symbol[rownames(tt)], tt)
  if ("biotype" %in% names(pc)) tt[, biotype := pc$biotype[match(gene_id, pc$gene_id)]]
  tt[, z := zscoreT(t, df = fit$df.total)]
  if (V2) {  # voom fit kept for check C12
    f0 <- eBayes(lmFit(voom(y, design), design))
    tt[, `:=`(logFC_voom = f0$coefficients[gene_id, "conditiontumor"], t_voom = f0$t[gene_id, "conditiontumor"],
              z_voom = zscoreT(f0$t[gene_id, "conditiontumor"], df = f0$df.total[match(gene_id, rownames(f0))]))]
  }
  fwrite(tt, sprintf("%s/%s__%s.tsv.gz", proc_dir, gid, contrast), sep = "\t")
  modifyList(out, list(
    n_tumor = sum(d$condition == "tumor"), n_normal = sum(d$condition == "normal"),
    n_dropped_missing_cov = n0 - nrow(d), n_dropped_rare_type = dropped,
    covariates = paste(setdiff(colnames(design), c("(Intercept)", "conditiontumor")), collapse = ","),
    n_genes = nrow(tt), n_de_fdr05_lfc1 = tt[adj.P.Val < 0.05 & abs(logFC) > 1, .N],
    n_up = tt[adj.P.Val < 0.05 & logFC > 1, .N], n_down = tt[adj.P.Val < 0.05 & logFC < -1, .N],
    normals_by_anchor = if (V2 && contrast == "gtex") d[role == normal_role, .N, keyby = anchor_class][, paste0(anchor_class, ":", N, collapse = ";")] else NA_character_,
    rho_vs_voom = if (V2) round(cor(tt$z, tt$z_voom, method = "spearman"), 4) else NA_real_,
    top500_vs_voom = if (V2) length(intersect(tt[order(-abs(z))][1:500, gene_id], tt[order(-abs(z_voom))][1:500, gene_id])) / 500 else NA_real_,
    z = list(setNames(tt$z, tt$gene_id)),
    z_voom = if (V2) list(setNames(tt$z_voom, tt$gene_id)) else NULL))
}

jobs <- CJ(group_id = if (ST) unique(sp$stratum) else groups$group_id, contrast = c("gtex", "adjacent"))
only <- Sys.getenv("ONLY_GROUPS")  # e.g. ONLY_GROUPS=TCGA-LUAD,BRCA-Basal for a quick test
if (nzchar(only)) jobs <- jobs[group_id %in% strsplit(only, ",")[[1]]]
res <- mclapply(seq_len(nrow(jobs)), function(i) {
  tryCatch(run(jobs$group_id[i], jobs$contrast[i]),
           error = function(e) list(group_id = jobs$group_id[i], contrast = jobs$contrast[i],
                                    status = "error", note = conditionMessage(e)))
}, mc.cores = N_CORES)

summ <- rbindlist(lapply(res, function(r) { r$z <- NULL; r$z_voom <- NULL; as.data.table(r) }), fill = TRUE)
if (summ[status == "error", .N]) print(summ[status == "error", .(group_id, contrast, note)])
for (ct in c("gtex", "adjacent")) for (zk in c("z", if (V2) "z_voom")) {
  zs <- Filter(function(r) identical(r$contrast, ct) && identical(r$status, "ok"), res)
  if (!length(zs)) next
  genes <- sort(unique(unlist(lapply(zs, function(r) names(r[[zk]][[1]])))))
  zm <- data.table(gene_id = genes, symbol = symbol[genes])
  for (r in zs) zm[, (r$group_id) := round(r[[zk]][[1]][genes], 3)]
  fwrite(zm, sprintf("%s/z_tumor_vs_%s%s.tsv", res_dir, ct, if (zk == "z_voom") "_voom" else ""), sep = "\t", na = "NA")
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
fwrite(summ, file.path(res_dir, "signature_summary.tsv"), sep = "\t")
message("SIGNATURES DONE: ", summ[status == "ok", .N], " ok, ", summ[status == "skipped", .N],
        " skipped, ", summ[status == "error", .N], " errors")
