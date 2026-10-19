#!/usr/bin/env Rscript
# Can GTEx collection effects (post-mortem, ischemia, autolysis, handling) be corrected, and how? Measures only.
# Per GTEx tissue used as a Step 4 reference (one sample per donor; genes = protein-coding minus mappability < 0.8
# and readthrough loci, as v2; filterByExpr + TMM logCPM):
#
# Q1. Variance: does correcting for collection variables lower donor-to-donor variance in donors NOT used to fit?
#     5-fold cross-validation over donors; per gene, held-out residual variance (MSE) of each model vs M0.
#       M0  age + sex + ESTIMATE                                   (what Step 4 already adjusts)
#       M1  M0 + log ischemic time + death circumstance (ventilator / fast / slow; Hardy 0 / 1-2 / 3-4)
#       M2  M1 + autolysis score (SMATSSCR) + time in fixative (SMTSPAX)
#       M3  M2 + RIN
#       M4  M3 + nucleic-acid isolation batch (SMNABTCHT; unseen batches get no batch effect)
#     Shifting all samples by the same amount does not change variance, so Q1 cannot choose the target (L1 vs L2).
# Q2. Target: which reference looks most like real best-preserved donors?
#     Best-preserved (BP) = ventilator death (Hardy 0) + ischemic time <= the tissue's 33rd percentile
#     (+ autolysis <= 1 where scored). BP donors are split into 5 folds; each fold is hidden, the model is fit on
#     everyone else, and three references are compared with the hidden BP donors' average profile (all at the
#     GTEx-average age / sex / ESTIMATE, using the training fit):
#       L1  model prediction at the average donor's collection values (= correction level 1)
#       L2  model prediction at the training BP donors' median collection values (= level 2)
#       L3  plain average of the training BP donors (= level 3, no model)
#     Score: RMSE over genes to the hidden BP average; gap closed by L2 = 1 - RMSE(L2) / RMSE(L1).
# Q3. Check only (never subtracted): does moving GTEx from L1 to L2 move it toward TCGA adjacent normals of the
#     same organ? Spearman over genes of (L2 - L1) with (adjacent mean - L1).
# Outputs: results/step4/reference_diagnostics/gtex_collection_{variance,anchor,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); RNGkind("L'Ecuyer-CMRG"); set.seed(34)
K <- 5; MIN_DONORS <- 40; MIN_BP <- 10; N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
out <- "results/step4/reference_diagnostics"

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
sa <- fread("data/raw/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
            select = c("SAMPID", "SMTSD", "SMTSISCH", "SMATSSCR", "SMTSPAX", "SMRIN", "SMNABTCHT"))
sj <- fread("data/raw/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
# v2 gene set
gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
lowmap <- sub("\\..*", "", gm[!is.na(map) & map < 0.8, id])
rt <- fread("data/interim/readthrough_genes_v23.tsv")[drop == TRUE, gene_id]
genes <- setdiff(rownames(counts), c(lowmap, rt))

g <- unique(cov[role == "gtex_normal" & sample %in% colnames(counts), .(sample, age, sex, estimate_score)], by = "sample")
g <- merge(g, sa, by.x = "sample", by.y = "SAMPID")
g[, donor := sub("^(GTEX-[^-]+).*", "\\1", sample)]
g[sj, on = c(donor = "SUBJID"), hardy := i.DTHHRDY]
g[, death := factor(fcase(hardy == 0, "ventilator", hardy %in% 1:2, "fast", hardy %in% 3:4, "slow"), c("ventilator", "fast", "slow"))]
g <- g[!is.na(age) & !is.na(sex) & !is.na(estimate_score) & !is.na(SMTSISCH) & !is.na(death)]
g <- unique(g, by = c("SMTSD", "donor"))
adj_of <- list("Breast - Mammary Tissue" = c("BRCA-Basal", "BRCA-LumA", "BRCA-LumB", "BRCA-Her2"),
               "Colon - Transverse" = c("TCGA-COAD", "TCGA-READ"), "Kidney - Cortex" = c("TCGA-KIRC", "TCGA-KIRP", "TCGA-KICH"),
               "Liver" = "TCGA-LIHC", "Lung" = c("TCGA-LUAD", "TCGA-LUSC"), "Prostate" = "TCGA-PRAD", "Stomach" = "TCGA-STAD",
               "Thyroid" = "TCGA-THCA", "Uterus" = "UCEC-EEC", "Esophagus - Mucosa" = "ESCA-AD")

design_for <- function(d, model, ref = NULL) {  # ref: data used to standardize / define levels (training set)
  if (is.null(ref)) ref <- d
  z <- function(x, r) (x - mean(r, na.rm = TRUE)) / sd(r, na.rm = TRUE)
  X <- cbind(age = z(d$age, ref$age), est = z(d$estimate_score, ref$estimate_score),
             sex = if (uniqueN(ref$sex) > 1) as.numeric(d$sex == ref$sex[1]) else NULL)
  if (model >= 1) X <- cbind(X, isch = z(log(d$SMTSISCH), log(ref$SMTSISCH)),
                             fast = as.numeric(d$death == "fast"), slow = as.numeric(d$death == "slow"))
  if (model >= 2) { if (all(!is.na(ref$SMATSSCR))) X <- cbind(X, autol = z(d$SMATSSCR, ref$SMATSSCR))
                    if (all(!is.na(ref$SMTSPAX))) X <- cbind(X, pax = z(d$SMTSPAX, ref$SMTSPAX)) }
  if (model >= 3 && all(!is.na(ref$SMRIN))) X <- cbind(X, rin = z(d$SMRIN, ref$SMRIN))
  if (model >= 4) { lv <- names(which(table(ref$SMNABTCHT) >= 3)); lv <- lv[-1]
    if (length(lv)) X <- cbind(X, sapply(lv, function(l) as.numeric(d$SMNABTCHT == l))) }
  X[is.na(X)] <- 0
  keep <- apply(cbind(X), 2, function(v) sd(v) > 0 || nrow(X) < nrow(ref)); cbind(int = 1, X[, keep, drop = FALSE])
}
fit_beta <- function(E, X) { b <- lmFit(E, X)$coefficients; b[is.na(b)] <- 0; b }
align <- function(X, cols) { M <- matrix(0, nrow(X), length(cols), dimnames = list(NULL, cols)); k <- intersect(cols, colnames(X)); M[, k] <- X[, k]; M }

per_tissue <- function(tis) {
  d <- g[SMTSD == tis]; if (nrow(d) < MIN_DONORS && !tis %in% names(adj_of)) return(NULL)
  if (nrow(d) < 20) return(NULL)
  adj <- if (tis %in% names(adj_of)) unique(cov[role == "adjacent_normal" & group_id %in% adj_of[[tis]] & sample %in% colnames(counts)], by = "patient") else NULL
  s_all <- c(d$sample, adj$sample)
  y <- DGEList(counts[genes, s_all]); k <- filterByExpr(y, group = c(rep("g", nrow(d)), rep("a", length(adj$sample))))
  y <- calcNormFactors(y[k, , keep.lib.sizes = FALSE]); E_all <- cpm(y, log = TRUE, prior.count = 2)
  E <- E_all[, d$sample]
  # Q1: cross-validated variance
  fold <- sample(rep(seq_len(K), length.out = nrow(d)))
  mse <- sapply(0:4, function(m) {
    r <- matrix(NA, nrow(E), nrow(d))
    for (f in seq_len(K)) { tr <- fold != f; Xtr <- design_for(d[tr], m); b <- fit_beta(E[, tr], Xtr)
      Xte <- align(design_for(d[!tr], m, d[tr]), colnames(Xtr)); r[, !tr] <- E[, !tr] - b %*% t(Xte) }
    rowMeans(r^2) })
  var_row <- data.table(tissue = tis, n_donors = nrow(d), genes = nrow(E),
                        t(setNames(round(apply(mse[, -1], 2, function(v) median(v / mse[, 1])), 4), paste0("mse_ratio_M", 1:4, "_vs_M0"))),
                        t(setNames(round(apply(mse[, -1], 2, function(v) mean(v < mse[, 1])), 3), paste0("frac_genes_lower_M", 1:4))))
  # Q2: anchor test with each collection model (M1-M4)
  q33 <- quantile(d$SMTSISCH, 1 / 3)
  bp <- d[death == "ventilator" & SMTSISCH <= q33 & (is.na(SMATSSCR) | SMATSSCR <= 1), which = TRUE]
  anc <- NULL; toward <- NULL
  if (length(bp) >= MIN_BP) {
    bfold <- sample(rep(seq_len(K), length.out = length(bp)))
    anc <- rbindlist(lapply(1:4, function(m) rbindlist(lapply(seq_len(K), function(f) {
      hid <- bp[bfold == f]; tr <- setdiff(seq_len(nrow(d)), hid); dtr <- d[tr]
      Xtr <- design_for(dtr, m); b <- fit_beta(E[, tr], Xtr); cn <- colnames(Xtr)
      demo <- intersect(cn, c("age", "est", "sex"))
      # hidden BP donors, demographics moved to the training average with the training fit
      Xh <- align(design_for(d[hid], m, dtr), cn)
      T <- rowMeans(E[, hid, drop = FALSE] - b[, demo, drop = FALSE] %*% t(Xh[, demo, drop = FALSE]))
      avg <- colMeans(Xtr); L1 <- as.numeric(b %*% avg)
      bpt <- intersect(bp, tr); Xbp <- Xtr[match(bpt, tr), , drop = FALSE]
      a2 <- avg; coll <- setdiff(cn, c("int", demo)); a2[coll] <- apply(Xbp[, coll, drop = FALSE], 2, median)
      L2 <- as.numeric(b %*% a2)
      L3 <- rowMeans(E[, bpt, drop = FALSE] - b[, demo, drop = FALSE] %*% t(Xbp[, demo, drop = FALSE]))
      rmse <- function(x) sqrt(mean((x - T)^2))
      data.table(tissue = tis, model = paste0("M", m), fold = f, n_hidden = length(hid), n_bp_train = length(bpt),
                 rmse_L1 = rmse(L1), rmse_L2 = rmse(L2), rmse_L3 = rmse(L3),
                 rho_shift_L2_vs_true = cor(L2 - L1, T - L1, method = "spearman"))
    }))))
    if (!is.null(adj)) {  # Q3: check only
      Xf <- design_for(d, 2); b <- fit_beta(E, Xf); cn <- colnames(Xf); demo <- intersect(cn, c("age", "est", "sex"))
      avg <- colMeans(Xf); a2 <- avg; coll <- setdiff(cn, c("int", demo)); a2[coll] <- apply(Xf[bp, coll, drop = FALSE], 2, median)
      L1 <- as.numeric(b %*% avg); L2 <- as.numeric(b %*% a2); A <- rowMeans(E_all[, adj$sample, drop = FALSE])
      toward <- data.table(tissue = tis, n_adjacent = length(adj$sample), rho_L2shift_vs_adjacent_gap = round(cor(L2 - L1, A - L1, method = "spearman"), 3),
                           gap_closed_toward_adjacent = round(1 - sqrt(mean((A - L2)^2)) / sqrt(mean((A - L1)^2)), 4))
    }
  }
  list(var = var_row, anc = anc, toward = toward, n_bp = length(bp))
}

tissues <- g[, .N, by = SMTSD][N >= 20, SMTSD]
res <- mclapply(tissues, function(t) tryCatch(per_tissue(t), error = function(e) { message(t, ": ", conditionMessage(e)); NULL }), mc.cores = N_CORES)
res <- Filter(Negate(is.null), res)
V <- rbindlist(lapply(res, `[[`, "var")); V[, n_best_preserved := sapply(res, `[[`, "n_bp")]
A <- rbindlist(lapply(res, `[[`, "anc")); Tw <- rbindlist(lapply(res, `[[`, "toward"))
fwrite(V, file.path(out, "gtex_collection_variance.tsv"), sep = "\t")
fwrite(A, file.path(out, "gtex_collection_anchor_folds.tsv"), sep = "\t")
As <- A[, .(rmse_L1 = mean(rmse_L1), rmse_L2 = mean(rmse_L2), rmse_L3 = mean(rmse_L3), rho_shift_L2_vs_true = mean(rho_shift_L2_vs_true)), by = .(tissue, model)]
As[, `:=`(gap_closed_L2 = 1 - rmse_L2 / rmse_L1, gap_closed_L3 = 1 - rmse_L3 / rmse_L1)]
num <- names(As)[sapply(As, is.numeric)]; As[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(As, file.path(out, "gtex_collection_anchor.tsv"), sep = "\t")
if (nrow(Tw)) fwrite(Tw, file.path(out, "gtex_collection_toward_adjacent.tsv"), sep = "\t")
print(V); print(As[order(tissue, model)]); print(Tw)
message("Overall: median MSE ratio vs M0: ", paste(sapply(paste0("mse_ratio_M", 1:4, "_vs_M0"), function(c) round(median(V[[c]]), 4)), collapse = " / "))
message("Overall: median gap closed by L2 per model: ", paste(As[, round(median(gap_closed_L2), 3), by = model]$V1, collapse = " / "),
        "; by L3: ", paste(As[, round(median(gap_closed_L3), 3), by = model]$V1, collapse = " / "))
message("GTEX COLLECTION DONE")
