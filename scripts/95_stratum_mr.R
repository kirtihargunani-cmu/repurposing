#!/usr/bin/env Rscript
# Step 10b: master regulators (MRs) of each Step 10 stratum, in the format of results/step6/aracne_project/group_mr.tsv,
# for drug reversal / specificity scoring on the drug side (user 2026-10-12). Method exactly as 50_viper.R
# (REGULON=aracne_project, the version the drug side's group_mr.tsv comes from; no mappability filter there, so none here):
#   signature  the stratum's Step 10 z (results/step10/z_strata_vs_gtex.tsv; stratum tumors vs the parent's GTEx normals)
#   network    the parent's project ARACNe network (BRCA; COADREAD), full regulator set, no target cap
#   filter     regulator expressed (median logCPM >= 1) in the STRATUM's tumors (user 2026-10-12: judged on the
#              subgroup's own patients); >= 10 targets present and >= 50% of mapped targets present
#   null       1,000 random halves of the parent's GTEx normals (Welch t -> z), msviper; BH FDR
# Check: the same code on each parent group's Step 4 signature and all its tumors should reproduce group_mr.tsv
#   (random null differs; NES Spearman and FDR < 0.10 overlap reported).
# Also writes the patient membership of every stratum (incl. breast axis ends) for per-stratum Celligner weights and
# target gate on the drug side.
# Outputs: results/step10/{stratum_mr,stratum_regulator_filter,stratum_mr_summary,stratum_patients,stratum_mr_parent_check}.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(95)
MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5; N_PERM <- 1000; MR_FDR <- 0.10
N_CORES <- as.integer(Sys.getenv("N_CORES", "12"))
COMPLEXES <- list(AP1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"), NFKB = c("NFKB1", "NFKB2", "RELA", "RELB", "REL"))
res_dir <- "results/step10"
arac <- readRDS("data/processed/step5/regulons_aracne_project.rds")
net_of <- function(gid) if (gid %in% c("TCGA-COAD", "TCGA-READ")) "COADREAD" else "BRCA"
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
zs <- fread(file.path(res_dir, "z_strata_vs_gtex.tsv")); z4 <- fread("results/step4/z_tumor_vs_gtex.tsv")
summ10 <- fread(file.path(res_dir, "signature_summary.tsv"))
gmr <- fread("results/step6/aracne_project/group_mr.tsv")

# stratum membership: COAD/READ subgroups + breast axis thirds (as 92_step10_signatures.R), restricted to patients with
# the Step 4 covariates (those that entered the signature)
pat <- fread("results/step8/step10_patient_strata.tsv")
ax <- fread("results/step8/breast_axes_scores.tsv"); strata <- fread("results/step8/step10_strata.tsv")
brs <- rbindlist(lapply(strata[kind == "continuous axis", stratum], function(s) {
  g <- strata[stratum == s, group_id]; a <- sub(".*-", "", s); v <- ax[group_id == g][[a]]; smp <- ax[group_id == g, sample]
  q <- quantile(v, c(1/3, 2/3))
  rbind(data.table(group_id = g, stratum = paste0(s, "-high"), sample = smp[v >= q[2]]),
        data.table(group_id = g, stratum = paste0(s, "-low"), sample = smp[v <= q[1]])) }))
mem <- rbind(pat, brs)
okcov <- cov[role == "tumor" & !is.na(age) & !is.na(estimate_score) & !is.na(sex), sample]
mem <- mem[sample %in% okcov & sample %in% colnames(counts)]
stopifnot(all(mem[, .N, by = stratum][match(summ10$stratum, stratum), N] == summ10$n_tumor))
mem[, weak := stratum %in% c("BRCA-Her2-PC4-high", "BRCA-Her2-PC4-low")]
fwrite(mem, file.path(res_dir, "stratum_patients.tsv"), sep = "\t")

perm_null <- function(nor) {  # normals: genes x samples logCPM
  pz <- sapply(seq_len(N_PERM), function(i) { lab <- sample(rep(c(TRUE, FALSE), length.out = ncol(nor)))
    a <- nor[, lab, drop = FALSE]; b <- nor[, !lab, drop = FALSE]
    t <- (rowMeans(a) - rowMeans(b)) / sqrt(apply(a, 1, var) / ncol(a) + apply(b, 1, var) / ncol(b) + 1e-6)
    zscoreT(pmin(pmax(t, -50), 50), df = ncol(a) + ncol(b) - 2) })
  rownames(pz) <- rownames(nor); pz }

run <- function(name, gid, sig, tumors) {
  regs_all <- arac[[net_of(gid)]]
  nor_s <- cov[group_id == gid & role == "gtex_normal" & sample %in% colnames(counts), sample]
  y <- calcNormFactors(DGEList(counts[names(sig), c(tumors, nor_s)])); lc <- cpm(y, log = TRUE, prior.count = 1)
  med <- apply(lc[, tumors, drop = FALSE], 1, median)
  expressed <- function(src) {
    if (src %in% names(COMPLEXES)) { g <- sym2ens[COMPLEXES[[src]]]; return(any(med[g[!is.na(g)]] >= MIN_LOGCPM, na.rm = TRUE)) }
    g <- if (startsWith(src, "ENSG")) src else sym2ens[src]; !is.na(g) && !is.na(med[g]) && med[g] >= MIN_LOGCPM }
  filt <- rbindlist(lapply(names(regs_all), function(src) { tg <- names(regs_all[[src]]$tfmode); n_in <- sum(tg %in% names(sig))
    data.table(regulator = src, n_targets = length(tg), n_present = n_in, frac_present = round(n_in / length(tg), 3), expressed = expressed(src)) }))
  filt[, keep := expressed & n_present >= MIN_TARGETS & frac_present >= MIN_FRAC]
  filt[, reason := fifelse(keep, "kept", fifelse(!expressed, "not expressed",
                   fifelse(n_present < MIN_TARGETS, "too few targets present", "low fraction of targets present")))]
  reg <- lapply(regs_all[filt[keep == TRUE, regulator]], function(r) { k <- names(r$tfmode) %in% names(sig)
    list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) }); class(reg) <- "regulon"
  mr <- msviper(sig, reg, nullmodel = perm_null(lc[, nor_s]), minsize = MIN_TARGETS, verbose = FALSE)
  mrt <- data.table(group_id = name, parent_group = gid, regulator = names(mr$es$nes), nes = mr$es$nes, p = mr$es$p.value,
                    size = mr$es$size)[, fdr := p.adjust(p, "BH")][order(-abs(nes))]
  list(mr = mrt, filt = filt[, `:=`(group_id = name, parent_group = gid)])
}

jobs <- c(lapply(summ10$stratum, function(s) { gid <- summ10[stratum == s, group_id]; z <- zs[!is.na(get(s))]
            list(name = s, gid = gid, sig = setNames(z[[s]], z$gene_id), tumors = mem[stratum == s, sample], check = FALSE) }),
          lapply(unique(summ10$group_id), function(g) { z <- z4[!is.na(get(g))]
            list(name = g, gid = g, sig = setNames(z[[g]], z$gene_id), check = TRUE,
                 tumors = intersect(cov[group_id == g & role == "tumor", sample], colnames(counts))) }))
res <- mclapply(jobs, function(j) tryCatch(c(run(j$name, j$gid, j$sig, j$tumors), check = j$check),
                                           error = function(e) { message(j$name, ": ", conditionMessage(e)); NULL }),
                mc.cores = N_CORES, mc.preschedule = FALSE)
stopifnot(!any(vapply(res, is.null, TRUE)))
st <- Filter(function(r) !r$check, res); ck <- Filter(function(r) r$check, res)
mr <- rbindlist(lapply(st, `[[`, "mr")); fwrite(mr, file.path(res_dir, "stratum_mr.tsv"), sep = "\t")
fwrite(rbindlist(lapply(st, `[[`, "filt")), file.path(res_dir, "stratum_regulator_filter.tsv"), sep = "\t")

# parent reproduction check
chk <- rbindlist(lapply(ck, function(r) { a <- r$mr; b <- gmr[group_id == a$group_id[1]]; m <- merge(a, b, by = "regulator")
  sa <- a[fdr < MR_FDR, regulator]; sb <- b[fdr < MR_FDR, regulator]
  data.table(group_id = a$group_id[1], n_regulators_now = nrow(a), n_regulators_step6 = nrow(b),
             nes_spearman = round(cor(m$nes.x, m$nes.y, method = "spearman"), 4),
             mr_fdr10_now = length(sa), mr_fdr10_step6 = length(sb), jaccard_fdr10 = round(length(intersect(sa, sb)) / length(union(sa, sb)), 3)) }))
fwrite(chk, file.path(res_dir, "stratum_mr_parent_check.tsv"), sep = "\t"); print(chk)

# summary: MRs per stratum, similarity to parent, regulators specific vs siblings
sib_of <- function(s) { if (startsWith(s, "COAD-")) return(setdiff(c("COAD-CMS1/3-like", "COAD-CMS2/4-like"), s))
  if (startsWith(s, "READ-")) return(setdiff(paste0("READ-C", 1:3), s))
  sub("-(high|low)$", ifelse(endsWith(s, "-high"), "-low", "-high"), s) }
sm <- rbindlist(lapply(unique(mr$group_id), function(s) { a <- mr[group_id == s]; g <- a$parent_group[1]; p <- gmr[group_id == g]
  m <- merge(a, p, by = "regulator"); sig_s <- a[fdr < MR_FDR, regulator]; sig_p <- p[fdr < MR_FDR, regulator]
  sibs <- mr[group_id %in% sib_of(s)]; sib_hit <- sibs[fdr < MR_FDR, unique(regulator)]
  sib_nes <- sibs[, .(sib_nes = mean(nes)), by = regulator]; d <- merge(a, sib_nes, by = "regulator")
  top_diff <- d[fdr < MR_FDR][order(-abs(nes - sib_nes))][1:min(.N, 8)]
  data.table(stratum = s, parent_group = g, n_tumor = summ10[stratum == s, n_tumor], n_regulators = nrow(a),
             n_mr_fdr10 = length(sig_s), n_mr_fdr05 = a[fdr < 0.05, .N], parent_mr_fdr10 = length(sig_p),
             nes_spearman_vs_parent = round(cor(m$nes.x, m$nes.y, method = "spearman"), 3),
             mr_shared_with_parent = length(intersect(sig_s, sig_p)), mr_not_in_siblings = length(setdiff(sig_s, sib_hit)),
             top_differing_from_siblings = paste(sprintf("%s %+.1f vs %+.1f", top_diff$regulator, top_diff$nes, top_diff$sib_nes), collapse = "; "),
             weak = s %in% c("BRCA-Her2-PC4-high", "BRCA-Her2-PC4-low")) }))
fwrite(sm, file.path(res_dir, "stratum_mr_summary.tsv"), sep = "\t")
print(sm[, .(stratum, n_tumor, n_mr_fdr10, parent_mr_fdr10, nes_spearman_vs_parent, mr_not_in_siblings)])
message("STRATUM MR DONE: ", uniqueN(mr$group_id), " strata")
