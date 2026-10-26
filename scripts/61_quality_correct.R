#!/usr/bin/env Rscript
# Step 6e: RNA-quality + purity correction at the gene level, then patient-level VIPER (ARACNe regulons).
# Replaces the NES-level purity correction (56) as the Step 8 input (PI decision 2026-10-08).
#   1. Quality index per sample (within group, within role): mean z of long genes minus mean z of short genes
#      (top vs bottom genomic-span quintile of expressed protein-coding genes). Degraded RNA loses long-transcript
#      signal; this index explained the TCGA plate pattern (scripts/57_plate_quality.R). Compared with GTEx RIN and ischemic time in the normals.
#   2. Gene-level linear correction (Gallego Romero et al. 2014): each gene's logCPM is regressed, within the
#      group's tumors, on purity + quality index (purity = ESTIMATE score; ABSOLUTE for LGG) and, within the
#      normals, on the quality index; residual + gene mean kept. Patients without the purity covariate are dropped.
#   3. Patient VIPER exactly as in 50_viper.R (viperSignature vs the group's normals, z-score, 100 permutations;
#      same regulators kept by 50_viper.R's filter), on the corrected matrices.
#   4. Checks vs the purity-corrected NES (56): PCs 1-5 R^2 on quality index, plate (and chance), purity, ABSOLUTE;
#      effective dimensionality; patient concordance with the group's master regulators.
# Outputs: data/processed/step6/aracne_project_qc/<group>.patient_nes.tsv.gz
#          results/step6/aracne_project/quality_correction/{checks,patient_concordance,quality_index}.tsv
# NORM=tmm (default) | cqn. cqn (PI decision 2026-10-08): counts are normalized with conditional quantile normalization
#   (Hansen, Irizarry & Wu 2012; per-sample smooth GC-content and transcript-length effects; GC and length of the longest
#   GENCODE v23 transcript, data/interim/gene_gc_v23.tsv) instead of TMM logCPM, and readthrough genes
#   (data/interim/readthrough_genes_v23.tsv, drop == TRUE: fused or clone-named loci whose transcripts are all tagged
#   readthrough) are removed from the expression data and from the regulons (as regulators and targets). Purity +
#   quality correction and patient VIPER then run as above. Outputs go to aracne_project_cqn / quality_correction_cqn
#   and are compared with the TMM version (aracne_project_qc).
# MAPPABILITY=1 (default with NORM=tmm since 2026-10-09; MAPPABILITY=0 reproduces aracne_project_qc): mappability filter of Saha & Battle 2018 (as used for GTEx v8
#   trans-eQTLs; resources from figshare collection 4297352, GRCh38/GENCODE v26, in data/raw/crossmap/): genes with
#   gene mappability < 0.8 are removed from the expression data and as regulators/targets, and regulon edges between
#   cross-mappable gene pairs (symmetric cross-mappability > 0) are removed. Outputs: aracne_project_map /
#   quality_correction_map, compared with aracne_project_qc.
# VERSION=v2 (SPEC_rebuild_v2.md, Step 6b; D45, D57): v2 Step 4 signatures (results/step4_v2), v2 regulons (already
#   without cross-mappable edges and low-mappability genes, scripts/44), v2 group filters / MRs (results/step6/aracne_project_v2).
#   Tumors: purity = ESTIMATE in every group (LGG included). GTEx normals: the Step 4a reference samples only, regressed
#   on quality index + the 4a collection covariates (coll_*, centered at each tissue's anchor) and set to their value
#   at the anchor (residual + intercept: coll_* = 0, quality index at its mean), as in Step 4. Adjacent normals as v1.
#   Outputs: data/processed/step6/aracne_project_v2_patient/, results/step6/aracne_project_v2/quality_correction/
#   (checks compare with v1 aracne_project_map).
# ONLY_GROUPS=a,b (optional); N_CORES (default 16)
suppressPackageStartupMessages({library(data.table); library(edgeR); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(61)
N_PERM_PATIENT <- 100; MIN_TARGETS <- 10; TOP_MR <- 50; N_PC <- 5; N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
ABSOLUTE_GROUPS <- c("LGG-IDHmut-codel", "LGG-IDHmut-noncodel", "LGG-IDHwt")
V2 <- identical(Sys.getenv("VERSION", "v1"), "v2")
NORM <- Sys.getenv("NORM", "tmm"); stopifnot(NORM %in% c("tmm", "cqn"), !V2 || NORM == "tmm")
MAPF <- !V2 && Sys.getenv("MAPPABILITY", if (NORM == "cqn") "0" else "1") == "1"; stopifnot(!(MAPF && NORM == "cqn"))  # default on (PI 2026-10-09); v2: in the regulons
if (V2) ABSOLUTE_GROUPS <- character(0)  # ESTIMATE everywhere (D45)
TAG <- if (V2) "v2" else if (NORM == "cqn") "cqn" else if (MAPF) "map" else "qc"
out_nes <- if (V2) "data/processed/step6/aracne_project_v2_patient" else paste0("data/processed/step6/aracne_project_", TAG)
res_dir <- if (V2) "results/step6/aracne_project_v2/quality_correction" else
  paste0("results/step6/aracne_project/quality_correction", if (TAG == "qc") "" else paste0("_", TAG))
BASE <- if (V2) c(dir = "aracne_project_map", label = "v1: purity+quality, TMM + mappability filter (61)") else
  if (TAG != "qc") c(dir = "aracne_project_qc", label = "purity+quality, TMM (61)") else
  c(dir = "aracne_project_purity", label = "purity_only (56)")
NEW_LABEL <- switch(TAG, v2 = "v2: purity (ESTIMATE) + quality, GTEx at anchor (61)",
                    cqn = "purity+quality, CQN + no readthrough (61)", map = "purity+quality, TMM + mappability filter (61)",
                    qc = "purity+quality (61)")
VS <- if (V2) "_v2" else ""
if (MAPF) {
  gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
  lowmap <- sub("\\..*", "", gm[!is.na(map) & map < 0.8 & !grepl("_PAR_Y", id), id])  # Y copies of PAR genes (map 0) skipped
  cm <- fread(cmd = "zcat data/raw/crossmap/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz", header = FALSE,
              col.names = c("a", "b", "s"))[s > 0]
  cm[, `:=`(a = sub("\\..*", "", a), b = sub("\\..*", "", b), s = NULL)]
  cm <- unique(rbind(cm, cm[, .(a = b, b = a)])); setkey(cm, a, b)
  # filter each network's regulons once: drop edges between cross-mappable genes (one keyed join per network)
  filter_net <- function(reg) {
    ens <- vapply(names(reg), function(rn) if (startsWith(rn, "ENSG")) rn else { e <- sym2ens[rn]; if (is.na(e)) rn else e }, "")
    n <- lengths(lapply(reg, `[[`, "tfmode"))
    ed <- data.table(rn = rep(names(reg), n), a = rep(ens, n), b = unlist(lapply(reg, function(r) names(r$tfmode)), use.names = FALSE))
    ed[, cross := !is.na(cm[ed, on = .(a, b), mult = "first", which = TRUE])]
    message(sprintf("mappability filter: %d of %d edges cross-mappable", sum(ed$cross), nrow(ed)))
    drop <- split(ed$b[ed$cross], ed$rn[ed$cross])
    for (rn in names(drop)) { k <- !(names(reg[[rn]]$tfmode) %in% drop[[rn]])
      reg[[rn]] <- list(tfmode = reg[[rn]]$tfmode[k], likelihood = reg[[rn]]$likelihood[k]) }
    reg
  }
}
if (NORM == "cqn") {
  suppressPackageStartupMessages(library(cqn))
  gcl <- fread("data/interim/gene_gc_v23.tsv"); gcv <- setNames(gcl$gc_longest_tx, gcl$gene_id); txl <- setNames(gcl$longest_tx_len, gcl$gene_id)
  rt <- fread("data/interim/readthrough_genes_v23.tsv")[drop == TRUE, gene_id]
}
dir.create(out_nes, recursive = TRUE, showWarnings = FALSE); dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

arac <- readRDS(paste0("data/processed/step5/regulons_aracne_project", VS, ".rds"))
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
summ <- fread(paste0("results/step4", VS, "/signature_summary.tsv"))[status == "ok"]
contrast_for <- summ[, .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
zmat <- list(gtex = fread(paste0("results/step4", VS, "/z_tumor_vs_gtex.tsv")), adjacent = fread(paste0("results/step4", VS, "/z_tumor_vs_adjacent.tsv")))
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
filt_all <- fread(paste0("results/step6/aracne_project", VS, "/regulator_filter.tsv"))[keep == TRUE]
mr_all <- fread(paste0("results/step6/aracne_project", VS, "/group_mr.tsv"))
COLL <- c("coll_isch", "coll_fast", "coll_slow", "coll_nonbp")
ref4a <- if (V2) fread("results/step4_v2/gtex_reference.tsv")[use == TRUE, c("sample", COLL), with = FALSE] else NULL
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
if (MAPF) { arac <- lapply(arac, filter_net); rm(cm); invisible(gc()) }
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))

qindex <- function(m) {  # genes x samples logCPM -> long-minus-short index per sample
  m <- m[apply(m, 1, median) >= 1 & rownames(m) %in% names(span), , drop = FALSE]
  z <- (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  colMeans(z[span[rownames(z)] >= q[2], , drop = FALSE]) - colMeans(z[span[rownames(z)] <= q[1], , drop = FALSE])
}
correct <- function(m, X) {  # regress each gene (row) on design X (samples x k, with intercept); keep residual + mean
  t(qr.resid(qr(X), t(m))) + rowMeans(m)
}
correct_at0 <- function(m, X) {  # as correct(), but keep residual + intercept = value where the centered covariates are 0
  q <- qr(X); t(qr.resid(q, t(m))) + qr.coef(q, t(m))[1, ]
}
r2f <- function(y, f) { f <- droplevels(factor(f)); if (nlevels(f) < 2) return(NA_real_)
  1 - sum(qr.resid(qr(model.matrix(~ f)), y)^2) / sum((y - mean(y))^2) }
r2c <- function(y, v) { ok <- !is.na(v); if (sum(ok) < 20) NA_real_ else cor(y[ok], v[ok])^2 }
pc_checks <- function(N, meta, qi, purity, label, gid) {  # N: patients x regulators
  p <- prcomp(N, center = TRUE); ve <- p$sdev^2 / sum(p$sdev^2); k <- min(N_PC, ncol(p$x))
  data.table(group_id = gid, version = label, pc = 1:k, var_explained = round(ve[1:k], 3),
    effective_dim = round(sum(ve)^2 / sum(ve^2), 2),
    r2_quality = round(sapply(1:k, function(i) r2c(p$x[, i], qi)), 3),
    r2_plate = round(sapply(1:k, function(i) r2f(p$x[, i], meta$batch)), 3),
    plate_chance = round((uniqueN(meta$batch) - 1) / (nrow(N) - 1), 3),
    r2_purity = round(sapply(1:k, function(i) r2c(p$x[, i], purity)), 3),
    r2_absolute = round(sapply(1:k, function(i) r2c(p$x[, i], meta$absolute_purity)), 3))
}

run <- function(gid, contrast) {
  normal_role <- if (contrast == "gtex") "gtex_normal" else "adjacent_normal"
  d <- cov[group_id == gid & role %in% c("tumor", normal_role) & sample %in% colnames(counts)]
  if (V2 && contrast == "gtex") d <- merge(d[role == "tumor" | sample %in% ref4a$sample], ref4a, by = "sample", all.x = TRUE)
  pcol <- if (gid %in% ABSOLUTE_GROUPS) "absolute_purity" else "estimate_score"
  d <- d[role != "tumor" | !is.na(get(pcol))]
  z <- zmat[[contrast]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]; sig <- setNames(z[[gid]], z$gene_id)
  if (NORM == "cqn") {
    genes <- names(sig)[!names(sig) %in% rt & names(sig) %in% names(gcv)]
    cq <- cqn(counts[genes, d$sample], x = gcv[genes], lengths = txl[genes], sizeFactors = colSums(counts[, d$sample]),
              verbose = FALSE)
    lc <- cq$y + cq$offset  # log2 expression with per-sample GC and length effects removed
  } else {
    genes <- if (MAPF) setdiff(names(sig), lowmap) else names(sig)
    lc <- cpm(calcNormFactors(DGEList(counts[genes, d$sample])), log = TRUE, prior.count = 1)
  }
  tum <- lc[, d[role == "tumor", sample], drop = FALSE]; nor <- lc[, d[role == normal_role, sample], drop = FALSE]
  mt <- d[role == "tumor"][match(colnames(tum), sample)]; mn <- d[role == normal_role][match(colnames(nor), sample)]
  qt <- qindex(tum); qn <- qindex(nor); pur <- mt[[pcol]]
  tumc <- correct(tum, cbind(1, pur - mean(pur), qt - mean(qt)))
  if (V2 && contrast == "gtex") {  # GTEx normals at the anchor (Step 4a collection covariates)
    C <- as.matrix(mn[, COLL, with = FALSE]); C <- C[, apply(C, 2, var) > 0, drop = FALSE]
    norc <- correct_at0(nor, cbind(1, qn - mean(qn), C))
  } else norc <- correct(nor, cbind(1, qn - mean(qn)))

  regs_all <- arac[[net_of[[gid]]]]
  keep_regs <- filt_all[group_id == gid, regulator]
  if (NORM == "cqn") keep_regs <- keep_regs[!(keep_regs %in% rt) & !(sym2ens[keep_regs] %in% rt)]
  if (MAPF) keep_regs <- keep_regs[!(sym2ens[keep_regs] %in% lowmap)]
  reg <- setNames(lapply(keep_regs, function(rn) { r <- regs_all[[rn]]
    k <- names(r$tfmode) %in% rownames(lc)
    list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) }), keep_regs)
  reg <- reg[vapply(reg, function(r) length(r$tfmode) >= MIN_TARGETS, logical(1))]
  class(reg) <- "regulon"
  vs <- viperSignature(tumc, norc, per = N_PERM_PATIENT, method = "zscore", verbose = FALSE)
  pn <- viper(vs, reg, minsize = MIN_TARGETS, verbose = FALSE)
  fwrite(data.table(regulator = rownames(pn), round(pn, 3)), file.path(out_nes, paste0(gid, ".patient_nes.tsv.gz")), sep = "\t")

  # checks: new (gene-level purity + quality) vs previous (NES-level purity correction, script 56)
  old <- fread(sprintf("data/processed/step6/%s/%s.patient_nes.tsv.gz", BASE[["dir"]], gid))
  O <- t(as.matrix(old[, -1])); colnames(O) <- old$regulator
  s <- intersect(rownames(O), colnames(pn)); r <- intersect(colnames(O), rownames(pn))
  meta <- mt[match(s, sample)]; qi <- qt[s]; pu <- meta[[pcol]]
  chk <- rbind(pc_checks(O[s, r], meta, qi, pu, BASE[["label"]], gid),
               pc_checks(t(pn[r, s]), meta, qi, pu, NEW_LABEL, gid))
  chk[, purity_covariate := pcol]
  top <- head(mr_all[group_id == gid & p < 0.05][order(-abs(nes))], TOP_MR)[regulator %in% r]
  conc <- if (nrow(top) >= 10) rbindlist(lapply(setNames(list(t(O[s, top$regulator]), pn[top$regulator, s]),
                                                         c(BASE[["label"]], NEW_LABEL)), function(M)
    data.table(median_rho = round(median(apply(M, 2, cor, y = top$nes, method = "spearman")), 3),
               median_frac_hit = round(median(colMeans(sign(M) == sign(top$nes) & abs(M) > 1.96)), 3))), idcol = "version")[
    , group_id := gid] else NULL
  qtab <- rbind(data.table(group_id = gid, sample = names(qt), role = "tumor", quality_index = round(qt, 4)),
                data.table(group_id = gid, sample = names(qn), role = normal_role, quality_index = round(qn, 4)))
  rin <- if (normal_role == "gtex_normal" && sum(!is.na(mn$gtex_rin)) >= 20) round(cor(qn, mn$gtex_rin, use = "complete"), 3) else NA_real_
  isch <- if (normal_role == "gtex_normal" && sum(!is.na(mn$gtex_ischemic_min)) >= 20)
    round(cor(qn, mn$gtex_ischemic_min, use = "complete", method = "spearman"), 3) else NA_real_
  chk[, `:=`(n_tumor = ncol(tum), n_normal = ncol(nor), r_quality_vs_gtex_rin = rin, rho_quality_vs_gtex_ischemic = isch)]
  list(chk = chk, conc = conc, q = qtab)
}

only <- Sys.getenv("ONLY_GROUPS")
jobs <- if (nzchar(only)) contrast_for[group_id %in% strsplit(only, ",")[[1]]] else contrast_for
res <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(run(jobs$group_id[i], jobs$contrast[i]),
  error = function(e) list(err = data.table(group_id = jobs$group_id[i], error = conditionMessage(e)))),
  mc.cores = N_CORES, mc.preschedule = FALSE)
pick <- function(k) rbindlist(lapply(res, `[[`, k), fill = TRUE)
errs <- pick("err"); if (nrow(errs)) print(errs)
fwrite(pick("chk"), file.path(res_dir, "checks.tsv"), sep = "\t")
fwrite(pick("conc"), file.path(res_dir, "patient_concordance.tsv"), sep = "\t")
fwrite(pick("q"), file.path(res_dir, "quality_index.tsv"), sep = "\t")
message("QUALITY CORRECTION DONE: ", length(res) - nrow(errs), " groups ok, ", nrow(errs), " errors")
