#!/usr/bin/env Rscript
# Step 7 (D65 / D68 / D81): BRCA_LUM test - do luminal breast master regulators change materially on a luminal-only network?
#
# Groups BRCA-LumA, BRCA-LumB. Networks (TF-only build, cross-mappable edges removed):
#   BRCA       pooled breast network (1,040 tumors; what LumA / LumB use now)
#   BRCA_LUM   luminal-only network (769 LumA + LumB tumors)
#   BRCA_SUB1-3 size-matched controls: 769 random BRCA tumors of all subtypes (47_brca_size_controls.py; regulons in
#              regulons_aracne_project_v2_tf_brca_controls.rds). They separate "fewer tumors" from "luminal-only".
# Group NES as in 50_viper.R steps 1-2 (v2: Step 4 v2 z signature, Step 4a GTEx normals, regulon filters, msviper with
# a normal-split null), except that one set of 1,000 normal splits per group is shared by all five networks, so NES
# differences come from the network alone.
# Cutoff (D81, fixed 2026-11-01 before the comparison). Material if, in LumA or LumB, either fires:
#   (a) rho_w(BRCA_LUM, BRCA) < every rho_w(BRCA_SUB_k, BRCA) and <= median_k rho_w(BRCA_SUB_k, BRCA) - 0.10
#       rho_w = weighted Spearman: Pearson correlation of NES ranks over regulators scored on both networks, each
#       regulator weighted by its mean |NES| on the two
#   (c) ESR1, FOXA1, GATA3 or XBP1: sign(NES) differs between BRCA_LUM and BRCA, |NES| >= 1.96 on at least one of the
#       two, and the same flip happens on no control
#   Context only (not in the rule): unweighted Spearman, top-100 overlap by |NES|.
# Sanity: the pooled-BRCA NES here vs results/step6/aracne_project_v2_tf/group_mr.tsv (same inputs, other random splits).
# Outputs: results/step7/brca_lum_test/{group_mr,agreement,known_biology,verdict}.tsv
# DRYRUN=1: every network replaced by pooled BRCA (tests the code without showing the result; all rho = 1, no flips).
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5; N_PERM_GROUP <- 1000   # as 50_viper.R
RHO_MARGIN <- 0.10; FLIP_NES <- 1.96; TOP_N <- 100                          # D81
KNOWN <- c("ESR1", "FOXA1", "GATA3", "XBP1")
GROUPS <- c("BRCA-LumA", "BRCA-LumB")
DRY <- identical(Sys.getenv("DRYRUN"), "1")
out <- if (DRY) "results/step7/brca_lum_test_dryrun" else "results/step7/brca_lum_test"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

prod <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
cf <- "data/processed/step5/regulons_aracne_project_v2_tf_brca_controls.rds"
ctrl <- if (file.exists(cf)) readRDS(cf) else list()
CTRLS <- paste0("BRCA_SUB", 1:3)
nets <- c(list(BRCA = prod$BRCA, BRCA_LUM = prod$BRCA_LUM), ctrl[CTRLS])
if (DRY) nets <- setNames(rep(list(prod$BRCA), 5), c("BRCA", "BRCA_LUM", CTRLS))
miss <- setdiff(c("BRCA", "BRCA_LUM", CTRLS), names(Filter(Negate(is.null), nets)))
if (length(miss)) stop("networks missing: ", paste(miss, collapse = ", "),
                       " (controls: 43_run_aracne.sh then VERSION=v2 REGSET=tf OUT_TAG=_brca_controls 44_aracne_regulons.R BRCA_SUB1 BRCA_SUB2 BRCA_SUB3)")

idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
z <- fread("results/step4_v2/z_tumor_vs_gtex.tsv")
ref <- fread("results/step4_v2/gtex_reference.tsv")[use == TRUE]
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")

score_group <- function(gid) {
  d <- cov[group_id == gid & role %in% c("tumor", "gtex_normal") & sample %in% colnames(counts)]
  d <- d[role == "tumor" | sample %in% ref$sample]
  zz <- z[, c("gene_id", gid), with = FALSE][!is.na(get(gid))]
  sig <- setNames(zz[[gid]], zz$gene_id)
  lc <- cpm(calcNormFactors(DGEList(counts[names(sig), d$sample])), log = TRUE, prior.count = 1)
  tum <- lc[, d[role == "tumor", sample], drop = FALSE]; nor <- lc[, d[role == "gtex_normal", sample], drop = FALSE]
  med <- apply(tum, 1, median)
  perm_z <- sapply(seq_len(N_PERM_GROUP), function(i) {   # one null for all networks
    lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(nor)))
    a <- nor[, lab, drop = FALSE]; b <- nor[, !lab, drop = FALSE]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(a) + ncol(b) - 2)
  })
  rownames(perm_z) <- rownames(lc)
  rbindlist(mclapply(names(nets), function(net) {
    ra <- nets[[net]]
    keep <- vapply(names(ra), function(src) {
      g <- if (startsWith(src, "ENSG")) src else sym2ens[src]
      tg <- names(ra[[src]]$tfmode); n_in <- sum(tg %in% names(sig))
      !is.na(g) && !is.na(med[g]) && med[g] >= MIN_LOGCPM && n_in >= MIN_TARGETS && n_in / length(tg) >= MIN_FRAC
    }, logical(1))
    reg <- lapply(ra[keep], function(r) { k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
    class(reg) <- "regulon"
    mr <- msviper(sig, reg, nullmodel = perm_z, minsize = MIN_TARGETS, verbose = FALSE)
    data.table(group_id = gid, network = net, regulator = names(mr$es$nes), nes = mr$es$nes, p = mr$es$p.value,
               size = mr$es$size)[, fdr := p.adjust(p, "BH")]
  }, mc.cores = length(nets)))
}
mr <- rbindlist(lapply(GROUPS, score_group))
fwrite(mr, file.path(out, "group_mr.tsv"), sep = "\t")

wcor <- function(a, b) {  # weighted Spearman: weighted Pearson of ranks, weight = mean |NES|
  ra <- rank(a); rb <- rank(b); w <- (abs(a) + abs(b)) / 2; w <- w / sum(w)
  ma <- sum(w * ra); mb <- sum(w * rb)
  sum(w * (ra - ma) * (rb - mb)) / sqrt(sum(w * (ra - ma)^2) * sum(w * (rb - mb)^2))
}
agree <- rbindlist(lapply(GROUPS, function(gid) rbindlist(lapply(setdiff(names(nets), "BRCA"), function(net) {
  m <- merge(mr[group_id == gid & network == "BRCA", .(regulator, a = nes)],
             mr[group_id == gid & network == net, .(regulator, b = nes)], by = "regulator")
  ta <- head(m[order(-abs(a)), regulator], TOP_N); tb <- head(m[order(-abs(b)), regulator], TOP_N)
  data.table(group_id = gid, network = net, n_common = nrow(m), rho_w = round(wcor(m$a, m$b), 4),
             rho = round(cor(m$a, m$b, method = "spearman"), 4), top100_overlap = length(intersect(ta, tb)) / TOP_N)
}))))
fwrite(agree, file.path(out, "agreement.tsv"), sep = "\t")

kb <- dcast(mr[regulator %in% KNOWN], group_id + regulator ~ network, value.var = "nes")
flip_on <- function(net) !is.na(kb[[net]]) & !is.na(kb$BRCA) & sign(kb[[net]]) != sign(kb$BRCA) &
  pmax(abs(kb[[net]]), abs(kb$BRCA)) >= FLIP_NES
kb[, flip_lum := flip_on("BRCA_LUM")]
kb[, n_ctrl_flip := Reduce(`+`, lapply(CTRLS, function(n) as.integer(flip_on(n))))]
kb[, fires := flip_lum & n_ctrl_flip == 0]
fwrite(kb, file.path(out, "known_biology.tsv"), sep = "\t")

verdict <- rbindlist(lapply(GROUPS, function(gid) {
  ag <- agree[group_id == gid]; lum <- ag[network == "BRCA_LUM", rho_w]; ct <- ag[network %in% CTRLS, rho_w]
  a_fires <- lum < min(ct) && lum <= median(ct) - RHO_MARGIN
  c_fires <- any(kb[group_id == gid, fires])
  data.table(group_id = gid, rho_w_lum = lum, rho_w_ctrl_min = min(ct), rho_w_ctrl_median = median(ct),
             a_fires = a_fires, c_fires = c_fires,
             c_genes = paste(kb[group_id == gid & fires == TRUE, regulator], collapse = ","), material = a_fires || c_fires)
}))
fwrite(verdict, file.path(out, "verdict.tsv"), sep = "\t")

old <- fread("results/step6/aracne_project_v2_tf/group_mr.tsv")[group_id %in% GROUPS]
chk <- merge(old[, .(group_id, regulator, nes_50 = nes)], mr[network == "BRCA", .(group_id, regulator, nes)], by = c("group_id", "regulator"))
message("sanity, pooled BRCA vs 50_viper.R: ", paste(chk[, .(r = round(cor(nes_50, nes), 4)), by = group_id][, paste(group_id, r)], collapse = "; "))
print(agree); print(kb); print(verdict)
message("BRCA_LUM TEST (D81): ", if (any(verdict$material)) "MATERIAL CHANGE" else "no material change", if (DRY) " [DRYRUN]")
