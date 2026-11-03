#!/usr/bin/env Rscript
# D86 perturbation-network pilot: regulons from the line networks (scripts/pnet/pnet_run_aracne.sh) and the pre-set comparison.
# Rule fixed 2026-11-02 before running (docs/DECISIONS.md D86).
#   Regulons: viper::aracne2regulon on each network.txt and the matrix it was built from (as scripts/44; no
#     cross-mappability filter: that corrects RNA-seq read assignment, not L1000 beads); regulators -> HGNC symbols.
#   Sets per line, on landmark genes, regulator removed from its own targets, >= 10 targets (no 50% rule):
#     P  perturbation network of the line built on all compounds (arm raw or specific)
#     T  tumor ARACNe v2 TF, lineage network (MCF7 -> BRCA, PC3 -> PRAD, A375 -> SKCM)
#     C  CollecTRI
#   Primary: D84 knockdown / overexpression signatures (main analysis: own_ok) in the 3 lines on regulators scored by
#     P, T and C; correct = NES sign as expected (down for knockdown, up for overexpression), paired by signature.
#   Scale up if, in either arm, P - T >= 5 points with sign-test p < 0.05 over signatures where P and T differ, and
#     P - C >= -2 points; both qualify -> the higher P goes forward. Else stop and report.
#   Reported: ligand profiles in the 3 lines scored held out (the fold network without the drug's fold; T and C on the
#     same profiles); edge stability (share of the all-compound network's edges found in >= 4 of 5 fold networks);
#     regulators with >= 10 targets per set.
# Outputs: data/processed/lincs/pnet/regulons.rds; results/lincs/pnet/{signatures,summary,ligands,stability}.tsv
suppressPackageStartupMessages({library(data.table); library(viper)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../.."))
setwd(root)
MIN_T <- 10; D <- "data/processed/lincs/pnet"; res <- "results/lincs/pnet"
LINES <- c(MCF7 = "ACH-000019", PC3 = "ACH-000090", A375 = "ACH-000219"); TNET <- c(MCF7 = "BRCA", PC3 = "PRAD", A375 = "SKCM")
ARMS <- c("raw", "specific"); SETS <- c("all", paste0("f", 1:5))
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; ens2sym <- setNames(idmap$symbol, idmap$gene_id)

## 1. regulons
rf <- file.path(D, "regulons.rds")
if (!file.exists(rf)) {
  R <- list()
  for (l in names(LINES)) for (a in ARMS) for (s in SETS) {
    x <- fread(file.path(D, sprintf("%s_%s/%s.exp.tsv", l, a, s))); m <- as.matrix(x[, -1]); rownames(m) <- x$gene
    f <- file.path(D, sprintf("%s_%s/net_%s/network.txt", l, a, s)); stopifnot(file.exists(f))
    tmp <- tempfile(fileext = ".txt"); ed <- fread(f)[, 1:3]; fwrite(ed, tmp, sep = "\t", col.names = FALSE)
    rg <- aracne2regulon(tmp, m, format = "3col", verbose = FALSE); unlink(tmp)
    names(rg) <- ifelse(is.na(ens2sym[names(rg)]), names(rg), ens2sym[names(rg)])
    R[[sprintf("%s_%s/%s", l, a, s)]] <- rg
    message(sprintf("%s_%s/%s: %d edges, %d regulators", l, a, s, nrow(ed), length(rg)))
  }
  saveRDS(R, rf)
}
R <- readRDS(rf)

## 2. benchmark signatures (landmark genes)
genes <- fread("data/processed/lincs/level5_bing_genes.tsv")
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
ld <- function(f, ids) { X <- matrix(readBin(file.path("data/processed/lincs", f), "numeric", size = 4, n = nrow(genes) * length(ids)), nrow = nrow(genes))
  X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- ids; X }
mk <- fread("data/processed/lincs/knockdown_bing_sigs.tsv.gz"); mo <- fread("data/processed/lincs/overexpression_bing_sigs.tsv.gz")
X <- cbind(ld("knockdown_bing.f32", mk$sig_id), ld("overexpression_bing.f32", mo$sig_id))
lm <- genes$ensembl_id[ok_g & genes$feature_space == "landmark"]; X <- X[lm, ]
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])
sg <- unique(fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[own_ok == TRUE & cell_iname %in% names(LINES),
             .(sig_id, cmap_name, cell_iname, technology, direction)])

prep <- function(rg, genes_in) {
  r <- lapply(names(rg), function(n) { tm <- rg[[n]]$tfmode; lk <- rg[[n]]$likelihood; e <- sym2ens[n]
    k <- names(tm) %in% genes_in & (is.na(e) | names(tm) != e)
    if (sum(k) < MIN_T) return(NULL); list(tfmode = tm[k], likelihood = lk[k]) })
  names(r) <- names(rg); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r
}
Tn <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds"); Cr <- readRDS("results/step5/regulons_core.rds")$collectri
nes_of <- function(reg, M, tf) { v <- viper(M, reg, method = "none", minsize = MIN_T, verbose = FALSE)
  v[cbind(match(tf, rownames(v)), seq_len(ncol(M)))] }

out <- list(); nreg <- list()
for (l in names(LINES)) {
  T_ <- prep(Tn[[TNET[l]]], lm); C_ <- prep(Cr, lm)
  for (a in ARMS) {
    P_ <- prep(R[[sprintf("%s_%s/all", l, a)]], lm)
    common <- Reduce(intersect, list(names(P_), names(T_), names(C_)))
    nreg[[paste(l, a)]] <- data.table(line = l, arm = a, P = length(P_), T = length(T_), C = length(C_), common = length(common))
    s <- sg[cell_iname == l & cmap_name %in% common]
    if (!nrow(s)) next
    M <- X[, s$sig_id, drop = FALSE]
    out[[paste(l, a)]] <- data.table(line = l, arm = a, s, nes_P = nes_of(P_, M, s$cmap_name), nes_T = nes_of(T_, M, s$cmap_name), nes_C = nes_of(C_, M, s$cmap_name))
  }
}
nreg <- rbindlist(nreg); print(nreg)
S <- rbindlist(out)
for (k in c("P", "T", "C")) S[[paste0("ok_", k)]] <- sign(S[[paste0("nes_", k)]]) == S$direction
fwrite(S, file.path(res, "signatures.tsv"), sep = "\t")

dec <- rbindlist(lapply(ARMS, function(a) { x <- S[arm == a]
  b <- sum(x$ok_P & !x$ok_T); w <- sum(!x$ok_P & x$ok_T)
  data.table(arm = a, signatures = nrow(x), tfs = uniqueN(x$cmap_name), correct_P = mean(x$ok_P), correct_T = mean(x$ok_T), correct_C = mean(x$ok_C),
             P_minus_T = mean(x$ok_P) - mean(x$ok_T), P_minus_C = mean(x$ok_P) - mean(x$ok_C), P_only = b, T_only = w,
             sign_p = if (b + w) binom.test(b, b + w)$p.value else 1) }))
dec[, qualifies := P_minus_T >= 0.05 & sign_p < 0.05 & P_minus_C >= -0.02]
dec[, decision := if (any(qualifies)) fifelse(qualifies & correct_P == max(correct_P[qualifies]), "scale_up", "-") else "stop"]
fwrite(dec, file.path(res, "summary.tsv"), sep = "\t")
print(dec[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 3) else v)])
print(S[, .(sigs = .N, P = round(mean(ok_P), 3), T = round(mean(ok_T), 3), C = round(mean(ok_C), 3)), by = .(arm, line)])
print(S[, .(sigs = .N, P = round(mean(ok_P), 3), T = round(mean(ok_T), 3), C = round(mean(ok_C), 3)), by = .(arm, technology)])

## 3. reported: edge stability
st <- rbindlist(lapply(names(LINES), function(l) rbindlist(lapply(ARMS, function(a) {
  e <- function(s) { f <- fread(file.path(D, sprintf("%s_%s/net_%s/network.txt", l, a, s)))[, 1:2]; paste(f[[1]], f[[2]]) }
  full <- e("all"); hits <- Reduce(`+`, lapply(paste0("f", 1:5), function(s) full %in% e(s)))
  data.table(line = l, arm = a, edges_all = length(full), share_in_4of5 = round(mean(hits >= 4), 3)) }))))
fwrite(st, file.path(res, "stability.tsv"), sep = "\t"); print(st)

## 4. reported: ligands, held out
idx <- fread("data/processed/lincs/drug_line_profiles_index.tsv"); pg <- fread("data/processed/lincs/drug_line_profiles_genes.tsv")
folds <- fread(file.path(res, "folds.tsv"))
ann <- fread("results/lincs/drug_annotations.tsv", colClasses = "character")
tg <- function(i) unique(unlist(strsplit(paste(ann$lincs_target[i], ann$hub_target[i], sep = "|"), "\\|")))
moa <- tolower(paste(ann$hub_moa, ann$lincs_moa))
ago <- grepl("agonist|activator", gsub("antagonist", "", moa)); ant <- grepl("antagonist|inhibitor|degrader|blocker", moa)
lig <- rbindlist(lapply(which(ago != ant), function(i) data.table(pert_id = ann$pert_id[i], drug = ann$cmap_name[i], tf = tg(i), direction = fifelse(ago[i], 1L, -1L))))
lig <- unique(lig[tf != "" & tf != "NA"])
lig <- merge(lig, folds[, .(pert_id, ModelID, fold)], by = "pert_id")
lig <- lig[ModelID %in% LINES]; lig[, line := names(LINES)[match(ModelID, LINES)]]
LG <- list()
for (a in ARMS) {
  pf <- if (a == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32"
  P <- matrix(readBin(file.path("data/processed/lincs", pf), "numeric", size = 4, n = nrow(pg) * nrow(idx)), nrow = nrow(pg))
  rownames(P) <- pg$gene_id; P <- P[lm, ]
  for (i in seq_len(nrow(lig))) { r <- lig[i]
    j <- which(idx$pert_id == r$pert_id & idx$ModelID == r$ModelID); M <- P[, j, drop = FALSE]
    Pf <- prep(R[[sprintf("%s_%s/f%d", r$line, a, r$fold)]], lm); T_ <- prep(Tn[[TNET[r$line]]], lm); C_ <- prep(Cr, lm)
    if (!all(r$tf %in% names(Pf), r$tf %in% names(T_), r$tf %in% names(C_))) next
    LG[[length(LG) + 1]] <- data.table(arm = a, r, nes_P = nes_of(Pf, M, r$tf), nes_T = nes_of(T_, M, r$tf), nes_C = nes_of(C_, M, r$tf)) }
}
LG <- rbindlist(LG)
if (nrow(LG)) {
  for (k in c("P", "T", "C")) LG[[paste0("ok_", k)]] <- sign(LG[[paste0("nes_", k)]]) == LG$direction
  fwrite(LG, file.path(res, "ligands.tsv"), sep = "\t")
  print(LG[, .(profiles = .N, drugs = uniqueN(drug), P = round(mean(ok_P), 2), T = round(mean(ok_T), 2), C = round(mean(ok_C), 2)), by = .(arm, tf, direction)][order(arm, tf)])
} else message("no ligand profiles on regulators scored by all three sets")
message("PNET EVAL DONE")
