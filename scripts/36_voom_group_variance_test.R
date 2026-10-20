#!/usr/bin/env Rscript
# Test of assumption A2 (limma treats tumors and normals as equally variable): does modelling the variance of
# tumors and normals separately change or improve the Step 4 signatures and the Step 6 master regulators?
# Methods (everything else identical; v2 Step 4 model: v2 gene set, v2 GTEx reference (Step 4a, use == TRUE),
# ~ condition + age + sex + ESTIMATE + tumor sample type + coll_isch + coll_fast + coll_slow (+ coll_nonbp)):
#   voom   limma::voom (v1 / spec default; one mean-variance trend for all samples)
#   vbg    voomByGroup (Law, Dong, You; github YOU-k/voomByGroup 38c8abf, vendored in scripts/vendor/): separate
#          mean-variance trends, hence separate precision weights, for tumors and normals
#   vqw    limma::voomWithQualityWeights(var.group = condition): one trend, plus one variance factor per arm
# Groups: those with known-biology controls (Step 7) plus large / small contrasts of other tissues.
# Measures (per group x method):
#   bcv_tumor, bcv_normal   edgeR common BCV per arm (size of the heteroscedasticity)
#   rho_vs_voom, top500     Spearman of z with voom; overlap of the 500 genes with largest |z|
#   splithalf_rho           tumors and normals each split into random halves (N_SPLIT times); signature from
#                           each half; Spearman between the two half-signatures (reproducibility)
#   MR (msviper, v1 project ARACNe regulons - v2 networks not built yet; null and filters as scripts/50_viper.R):
#                           NES Spearman with voom, FDR < 0.05 overlap (Jaccard), known-biology controls,
#                           split-half NES Spearman (N_SPLIT_MR halves)
# Decision rule (fixed before running): vbg (or vqw) is adopted if (1) the arms differ in BCV by >= 20% in most
# groups, (2) median split-half reproducibility of z and of NES is not lower than voom (difference >= -0.005) and
# higher in most groups, (3) known-biology controls do not get worse. If (2)-(3) hold but the change is
# negligible (median rho_vs_voom >= 0.99 and NES rho >= 0.98), it is reported as "no practical effect".
# Outputs: results/step4/voom_group_variance/{summary,splithalf,controls,mr_overlap}.tsv
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
source("scripts/vendor/voomByGroup.R")
OUT <- "results/step4/voom_group_variance"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_CORES <- as.integer(Sys.getenv("N_CORES", "12")); N_SPLIT <- as.integer(Sys.getenv("N_SPLIT", "10")); N_SPLIT_MR <- min(3L, N_SPLIT); N_PERM <- 1000
MIN_CLASS <- 3; MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5; TOP <- 500
METHODS <- c("voom", "vbg", "vqw")
GROUPS <- if (nzchar(Sys.getenv("ONLY_GROUPS"))) strsplit(Sys.getenv("ONLY_GROUPS"), ",")[[1]] else c("BRCA-LumA", "BRCA-LumB", "BRCA-Basal", "TCGA-KIRC", "TCGA-LIHC", "TCGA-PRAD", "TCGA-COAD", "TCGA-READ",
            "TCGA-LUAD", "TCGA-THCA", "TCGA-STAD", "ESCA-SCC", "PCPG-PGL")
CONTROLS <- rbindlist(list(
  data.table(regulator = c("E2F1", "FOXM1", "MYBL2"), expected = 1, group_id = list(GROUPS)),
  data.table(regulator = c("ESR1", "FOXA1", "GATA3"), expected = 1, group_id = list(c("BRCA-LumA", "BRCA-LumB"))),
  data.table(regulator = c("MYC", "FOXM1"), expected = 1, group_id = list("BRCA-Basal")),
  data.table(regulator = c("EPAS1", "HIF1A"), expected = 1, group_id = list("TCGA-KIRC")),
  data.table(regulator = c("HNF4A", "HNF1A"), expected = -1, group_id = list("TCGA-LIHC")),
  data.table(regulator = c("AR", "ERG"), expected = 1, group_id = list("TCGA-PRAD")),
  data.table(regulator = "MYC", expected = 1, group_id = list(c("TCGA-COAD", "TCGA-READ")))))
CONTROLS <- unique(CONTROLS[, .(group_id = unlist(group_id)), by = .(regulator, expected)])
COMPLEXES <- list(AP1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"),
                  NFKB = c("NFKB1", "NFKB2", "RELA", "RELB", "REL"))
RNGkind("L'Ecuyer-CMRG"); set.seed(1)

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]
ref <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE]
genes <- fread("results/step3/genes_protein_coding_v23_v2.tsv")$gene_id
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
counts <- counts[rownames(counts) %in% genes, ]
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
arac <- readRDS("data/processed/step5/regulons_aracne_project.rds")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)

# samples + design for one group (v2 Step 4 model)
setup <- function(gid) {
  d <- cov[group_id == gid & (role == "tumor" | (role == "gtex_normal" & sample %in% ref$sample)) & sample %in% colnames(counts)]
  d <- merge(d, ref[, .(sample, coll_isch, coll_fast, coll_slow, coll_nonbp)], by = "sample", all.x = TRUE)
  for (k in c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")) d[role == "tumor", (k) := 0]
  d[, condition := factor(ifelse(role == "tumor", "tumor", "normal"), levels = c("normal", "tumor"))]
  cls <- d[role == "tumor", .N, by = sample_type_class]
  rare <- cls[sample_type_class != "primary" & N < MIN_CLASS, sample_type_class]
  d <- d[!(role == "tumor" & sample_type_class %in% rare)]
  extra <- setdiff(cls[N >= MIN_CLASS, sample_type_class], "primary")
  for (k in extra) d[, (paste0("st_", k)) := as.integer(role == "tumor" & sample_type_class == k)]
  d <- d[!is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  coll <- c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")
  coll <- coll[sapply(coll, function(k) var(d[[k]]) > 0)]
  terms <- c("condition", "age_z", "estimate_z", if (length(extra)) paste0("st_", extra), coll)
  if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= MIN_CLASS) terms <- c(terms, "sex")
  list(d = d, terms = terms)
}
design_of <- function(d, terms) {
  d <- copy(d); d[, age_z := as.numeric(scale(age))]; d[, estimate_z := as.numeric(scale(estimate_score))]
  X <- model.matrix(reformulate(terms), data = d)
  X[, !colnames(X) %in% nonEstimable(X), drop = FALSE]
}
# z signature for one method on a sample subset
signature <- function(d, terms, method) {
  X <- design_of(d, terms)
  y <- DGEList(counts[, d$sample])
  y <- calcNormFactors(y[filterByExpr(y, X), , keep.lib.sizes = FALSE])
  v <- switch(method,
    voom = voom(y, X),
    vbg  = voomByGroup(y, group = d$condition, design = X, print = FALSE),
    vqw  = voomWithQualityWeights(y, X, var.group = d$condition))
  fit <- eBayes(lmFit(v, X))
  setNames(zscoreT(fit$t[, "conditiontumor"], df = fit$df.total), rownames(fit))
}
half <- function(d) {  # random halves of tumors and of normals
  a <- d[, .(sample, h = sample(rep(1:2, length.out = .N))), by = condition]
  list(d[sample %in% a[h == 1, sample]], d[sample %in% a[h == 2, sample]])
}
rho <- function(a, b) { g <- intersect(names(a), names(b)); cor(a[g], b[g], method = "spearman") }
top_overlap <- function(a, b) { g <- intersect(names(a), names(b)); a <- a[g]; b <- b[g]
  length(intersect(names(sort(-abs(a)))[1:TOP], names(sort(-abs(b)))[1:TOP])) / TOP }

run_group <- function(gid) {
  s <- setup(gid); d <- s$d
  y0 <- DGEList(counts[, d$sample]); y0 <- calcNormFactors(y0[filterByExpr(y0, design_of(d, s$terms)), , keep.lib.sizes = FALSE])
  bcv <- sapply(c("tumor", "normal"), function(k) sqrt(estimateCommonDisp(y0[, d$condition == k])$common))
  lc <- cpm(y0, log = TRUE, prior.count = 1)
  tum <- lc[, d$condition == "tumor"]; nor <- lc[, d$condition == "normal"]
  # regulon filters + null (as scripts/50_viper.R; the null is method-independent)
  regs_all <- arac[[net_of[[gid]]]]; med <- apply(tum, 1, median)
  expressed <- function(src) {
    g <- if (src %in% names(COMPLEXES)) sym2ens[COMPLEXES[[src]]] else sym2ens[src]
    any(med[g[!is.na(g)]] >= MIN_LOGCPM, na.rm = TRUE) }
  perm_z <- sapply(seq_len(N_PERM), function(i) {
    lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(nor))); a <- nor[, lab]; b <- nor[, !lab]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(a) + ncol(b) - 2) })
  rownames(perm_z) <- rownames(lc)
  mr_of <- function(sig) {
    g <- intersect(names(sig), rownames(perm_z)); sig <- sig[g]
    keep <- names(regs_all)[sapply(names(regs_all), function(src) {
      tg <- names(regs_all[[src]]$tfmode); n <- sum(tg %in% g)
      n >= MIN_TARGETS && n / length(tg) >= MIN_FRAC && expressed(src) })]
    reg <- lapply(regs_all[keep], function(r) { k <- names(r$tfmode) %in% g
      list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) }); class(reg) <- "regulon"
    m <- msviper(sig, reg, nullmodel = perm_z[g, ], minsize = MIN_TARGETS, verbose = FALSE)
    data.table(regulator = names(m$es$nes), nes = m$es$nes, p = m$es$p.value)[, fdr := p.adjust(p, "BH")] }
  z <- lapply(setNames(METHODS, METHODS), function(m) signature(d, s$terms, m))
  mr <- lapply(z, mr_of)
  halves <- lapply(seq_len(N_SPLIT), function(i) half(d))
  sh <- rbindlist(lapply(seq_len(N_SPLIT), function(i) rbindlist(lapply(METHODS, function(m) {
    za <- signature(halves[[i]][[1]], s$terms, m); zb <- signature(halves[[i]][[2]], s$terms, m)
    nes_rho <- if (i <= N_SPLIT_MR) { a <- mr_of(za); b <- merge(a, mr_of(zb), by = "regulator")
      cor(b$nes.x, b$nes.y, method = "spearman") } else NA_real_
    data.table(group_id = gid, rep = i, method = m, z_rho = rho(za, zb), nes_rho = nes_rho) }))))
  summ <- rbindlist(lapply(METHODS, function(m) {
    b <- merge(mr$voom, mr[[m]], by = "regulator")
    hv <- mr$voom[fdr < 0.05, regulator]; hm <- mr[[m]][fdr < 0.05, regulator]
    data.table(group_id = gid, method = m, n_tumor = sum(d$condition == "tumor"), n_normal = sum(d$condition == "normal"),
               bcv_tumor = bcv[["tumor"]], bcv_normal = bcv[["normal"]], n_genes = length(z[[m]]),
               rho_vs_voom = rho(z$voom, z[[m]]), top500_vs_voom = top_overlap(z$voom, z[[m]]),
               nes_rho_vs_voom = cor(b$nes.x, b$nes.y, method = "spearman"), n_mr_fdr05 = length(hm),
               mr_jaccard_vs_voom = length(intersect(hv, hm)) / max(1, length(union(hv, hm))),
               z_splithalf = median(sh[method == m, z_rho]), nes_splithalf = median(sh[method == m, nes_rho], na.rm = TRUE)) }))
  ctl <- rbindlist(lapply(METHODS, function(m) merge(CONTROLS[group_id == gid], mr[[m]][, .(regulator, nes, fdr)],
                                                     by = "regulator", all.x = TRUE)[, method := m]))
  message(gid, " done")
  list(summ = summ, sh = sh, ctl = ctl)
}

res <- mclapply(GROUPS, function(g) tryCatch(run_group(g), error = function(e) { message(g, ": ", conditionMessage(e)); NULL }),
                mc.cores = min(N_CORES, length(GROUPS)), mc.set.seed = TRUE)
res <- Filter(Negate(is.null), res)
summ <- rbindlist(lapply(res, `[[`, "summ")); sh <- rbindlist(lapply(res, `[[`, "sh")); ctl <- rbindlist(lapply(res, `[[`, "ctl"))
ctl[, verdict := fifelse(is.na(nes), "not tested", fifelse(sign(nes) != expected, "wrong direction",
                 fifelse(fdr < 0.05, "passes", "right direction, FDR >= 0.05")))]
num <- names(summ)[sapply(summ, is.double)]; summ[, (num) := lapply(.SD, round, 4), .SDcols = num]
fwrite(summ, file.path(OUT, "summary.tsv"), sep = "\t")
fwrite(sh, file.path(OUT, "splithalf.tsv"), sep = "\t", na = "NA")
fwrite(ctl, file.path(OUT, "controls.tsv"), sep = "\t", na = "NA")
w <- dcast(summ, group_id + n_tumor + n_normal + bcv_tumor + bcv_normal ~ method,
           value.var = c("rho_vs_voom", "nes_rho_vs_voom", "z_splithalf", "nes_splithalf"))
print(w)
print(ctl[, .N, by = .(method, verdict)][order(method, verdict)])
message("VOOM GROUP VARIANCE TEST DONE: ", nrow(summ) / length(METHODS), " groups")
