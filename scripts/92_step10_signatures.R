#!/usr/bin/env Rscript
# Step 10: disease signatures for the Step 8/9 strata (results/step8/step10_strata.tsv), for drug matching.
# User decisions 2026-10-12:
#   subgroups  COAD-CMS1/3-like, COAD-CMS2/4-like, READ-C1/C2/C3: subgroup tumors vs healthy GTEx tissue of the organ
#              (the parent group's GTEx normals), exactly as Step 4 (30_signatures.R)
#   breast     each kept axis (82_breast_axes.R): patients split into thirds by axis score (within subtype);
#              top third ("high") and bottom third ("low") each vs GTEx breast, as Step 4 -> two signatures per axis
#              (option C, used for drug matching). "high" = positive PC score; what that means: breast_axes_annotation.tsv.
#   companion  per breast axis, which genes follow the axis (option A; explanation only, not for matching):
#              limma-voom across all of the subtype's tumors, ~ axis score (z) + age + ESTIMATE (+ sex), moderated t -> z.
# Step 4 model: ~ condition + age_z + estimate_z (+ sex if both sexes >= 3); protein-coding genes, filterByExpr, TMM, voom,
#   eBayes, moderated t -> signed z (zscoreT). All tumors here are primary, one per patient (Step 8 patients).
# Checks: Spearman of each stratum signature with its parent group's Step 4 signature; high vs low end of each axis;
#   (high - low) difference vs the companion axis z.
# Outputs: data/processed/step10/<stratum>.tsv.gz (full tables, not in git)
#          results/step10/{signature_summary,z_strata_vs_gtex,z_axis_drivers,axis_driver_top}.tsv
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
N_CORES <- as.integer(Sys.getenv("N_CORES", "16")); MIN_CLASS <- 3; N_TOP <- 25
proc_dir <- "data/processed/step10"; res_dir <- "results/step10"
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
pc <- fread("results/step3/genes_protein_coding_v23.tsv"); symbol <- setNames(pc$symbol, pc$gene_id)
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
strata <- fread("results/step8/step10_strata.tsv")
pat <- fread("results/step8/step10_patient_strata.tsv")
ax <- fread("results/step8/breast_axes_scores.tsv")
step4 <- fread("results/step4/z_tumor_vs_gtex.tsv")
fname <- function(s) gsub("/", "_", s)

# job list: stratum name, parent group, tumor samples
jobs <- list()
for (s in strata[kind == "subgroup", stratum]) {
  g <- strata[stratum == s, group_id]; jobs[[s]] <- list(stratum = s, group_id = g, end = "", tumors = pat[stratum == s, sample]) }
for (s in strata[kind == "continuous axis", stratum]) {
  g <- strata[stratum == s, group_id]; a <- sub(".*-", "", s); sc <- ax[group_id == g, c("sample", a), with = FALSE]
  q <- quantile(sc[[a]], c(1/3, 2/3))
  jobs[[paste0(s, "-high")]] <- list(stratum = paste0(s, "-high"), group_id = g, end = "high", tumors = sc[get(a) >= q[2], sample])
  jobs[[paste0(s, "-low")]] <- list(stratum = paste0(s, "-low"), group_id = g, end = "low", tumors = sc[get(a) <= q[1], sample]) }

fit_z <- function(d, terms, coef) {
  if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= MIN_CLASS) terms <- c(terms, "sex")
  design <- model.matrix(reformulate(terms), data = d)
  y <- DGEList(counts[, d$sample]); keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); v <- voom(y, design); fit <- eBayes(lmFit(v, design))
  tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  tt <- data.table(gene_id = rownames(tt), symbol = symbol[rownames(tt)], tt)
  tt[, z := zscoreT(t, df = fit$df.total)]
  list(tt = tt, covariates = paste(setdiff(colnames(design), c("(Intercept)", coef)), collapse = ","))
}

run <- function(j) {
  d <- cov[(group_id == j$group_id & role == "gtex_normal") | (role == "tumor" & sample %in% j$tumors)]
  d <- unique(d, by = "sample")[sample %in% colnames(counts)]
  n0 <- nrow(d); d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  d[, `:=`(age_z = as.numeric(scale(age)), estimate_z = as.numeric(scale(estimate_score)))]
  stopifnot(all(d[role == "tumor", sample_type_class] %in% c("primary", "new_primary")))
  f <- fit_z(d, c("condition", "age_z", "estimate_z"), "conditiontumor"); tt <- f$tt
  fwrite(tt, sprintf("%s/%s.tsv.gz", proc_dir, fname(j$stratum)), sep = "\t")
  p4 <- setNames(step4[[j$group_id]], step4$gene_id); gs <- intersect(tt$gene_id, names(p4)[!is.na(p4)])
  list(summary = data.table(stratum = j$stratum, group_id = j$group_id, end = j$end,
    n_tumor = sum(d$condition == "tumor"), n_normal = sum(d$condition == "normal"), n_dropped_missing_cov = n0 - nrow(d),
    covariates = f$covariates, n_genes = nrow(tt), n_up = tt[adj.P.Val < 0.05 & logFC > 1, .N], n_down = tt[adj.P.Val < 0.05 & logFC < -1, .N],
    rho_vs_parent_step4 = round(cor(tt[match(gs, gene_id), z], p4[gs], method = "spearman"), 3)),
    z = setNames(tt$z, tt$gene_id))
}
res <- mclapply(jobs, function(j) tryCatch(run(j), error = function(e) { message(j$stratum, ": ", conditionMessage(e)); NULL }),
                mc.cores = N_CORES, mc.preschedule = FALSE)
stopifnot(!any(vapply(res, is.null, TRUE)))
summ <- rbindlist(lapply(res, `[[`, "summary"))
genes <- sort(unique(unlist(lapply(res, function(r) names(r$z)))))
zm <- data.table(gene_id = genes, symbol = symbol[genes])
for (r in res) zm[, (r$summary$stratum) := round(r$z[genes], 3)]
fwrite(zm, file.path(res_dir, "z_strata_vs_gtex.tsv"), sep = "\t", na = "NA")

# companion: genes that follow each breast axis, across all of the subtype's tumors
axes <- strata[kind == "continuous axis"]
comp <- mclapply(seq_len(nrow(axes)), function(i) {
  g <- axes$group_id[i]; a <- sub(".*-", "", axes$stratum[i])
  sc <- ax[group_id == g]; d <- cov[role == "tumor" & sample %in% sc$sample][sample %in% colnames(counts)]
  d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  d[, `:=`(axis_z = as.numeric(scale(sc[match(d$sample, sample)][[a]])), age_z = as.numeric(scale(age)),
           estimate_z = as.numeric(scale(estimate_score)))]
  f <- fit_z(d, c("axis_z", "age_z", "estimate_z"), "axis_z")
  list(axis = axes$stratum[i], n = nrow(d), z = setNames(f$tt$z, f$tt$gene_id))
}, mc.cores = N_CORES)
cg <- sort(unique(unlist(lapply(comp, function(r) names(r$z)))))
cz <- data.table(gene_id = cg, symbol = symbol[cg]); for (r in comp) cz[, (r$axis) := round(r$z[cg], 3)]
fwrite(cz, file.path(res_dir, "z_axis_drivers.tsv"), sep = "\t", na = "NA")
top <- rbindlist(lapply(comp, function(r) { z <- sort(r$z)
  data.table(axis = r$axis, n_tumors = r$n,
             genes_up_toward_high = paste(symbol[names(rev(tail(z, N_TOP)))], collapse = ","),
             genes_up_toward_low = paste(symbol[names(head(z, N_TOP))], collapse = ",")) }))
fwrite(top, file.path(res_dir, "axis_driver_top.tsv"), sep = "\t")

# checks per axis: high vs low end; does (high - low) follow the companion axis z
chk <- rbindlist(lapply(axes$stratum, function(s) {
  h <- zm[[paste0(s, "-high")]]; l <- zm[[paste0(s, "-low")]]; cmp <- cz[match(zm$gene_id, gene_id)][[s]]
  ok <- !is.na(h) & !is.na(l) & !is.na(cmp)
  data.table(stratum = paste0(s, c("-high", "-low")), rho_high_vs_low = round(cor(h[ok], l[ok], method = "spearman"), 3),
             rho_diff_vs_axis_drivers = round(cor(h[ok] - l[ok], cmp[ok], method = "spearman"), 3)) }))
summ <- merge(summ, chk, by = "stratum", all.x = TRUE)
summ <- summ[order(match(group_id, unique(strata$group_id)), stratum)]
fwrite(summ, file.path(res_dir, "signature_summary.tsv"), sep = "\t")
print(summ[, .(stratum, n_tumor, n_normal, n_up, n_down, rho_vs_parent_step4, rho_high_vs_low, rho_diff_vs_axis_drivers)])
message("STEP 10 SIGNATURES DONE: ", nrow(summ), " signatures, ", length(comp), " axis driver lists")
