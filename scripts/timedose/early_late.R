#!/usr/bin/env Rscript
# D89 (definitions fixed 2026-11-03 before running, docs/DECISIONS.md): early and late-only drug signatures, scored with
# tumor ARACNe v2 TF, and the early-vs-late test per regulator.
# Inputs: scripts/timedose/td_profiles.py time arm, line-specific version (time_T6_specific = E, time_T24_specific = L,
#   time_A0_specific = A0; time_index.tsv).
# Signatures: E, L, R = L - proj_E(L) (late-only), A0 (baseline). Activity: viper (method "none"), every ARACNe v2 TF
#   network, regulators / pruning as 79d. Saved as data/processed/lincs/timedose/early_late/<sig>_<net>.f32 (regulators x
#   pairs, column-major) + <net>.regulators.txt; the late-only profiles as time_R_specific.f32.
# Test: d = a(L) - a(E) per network x line x regulator; robust z over the line's pairs (lines with < 30 pairs pooled as
#   "other"), two-sided p, BH within network x line; FDR < 0.05. Onset fraction |a(R)| / (|a(E)| + |a(R)|).
# Outputs: results/lincs/early_late/{tests_matched.tsv.gz (lineage-matched networks), line_summary.tsv,
#          regulator_summary.tsv, baseline.tsv, antier_mcf7.tsv, drug_summary.tsv}
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
MIN_TARGETS <- 10; MIN_FRAC <- 0.5; MIN_LINE <- 30; FDR <- 0.05
N_CORES <- min(as.integer(Sys.getenv("N_CORES", "20")), 24L)
proc <- "data/processed/lincs"; td <- file.path(proc, "timedose"); od <- file.path(td, "early_late"); res <- "results/lincs/early_late"
dir.create(od, FALSE, TRUE); dir.create(res, FALSE, TRUE)
genes <- fread(file.path(proc, "drug_line_profiles_genes.tsv"))$gene_id
I <- fread(file.path(td, "time_index.tsv")); n <- nrow(I)
rd <- function(a) { M <- matrix(readBin(file.path(td, sprintf("time_%s_specific.f32", a)), "numeric", size = 4, n = length(genes) * n), nrow = length(genes))
  dimnames(M) <- list(genes, paste(I$pert_id, I$ModelID, sep = "|")); M }
E <- rd("T6"); L <- rd("T24"); A0 <- rd("A0")
u <- colSums(L * E) / pmax(colSums(E * E), 1e-12)
R <- L - sweep(E, 2, u, `*`)
writeBin(as.vector(R), file.path(td, "time_R_specific.f32"), size = 4)   # profile-major, as the other arms
cat(sprintf("pairs %d; share of late signal persisting from early (median R^2 of L on E) %.2f\n", n,
            median(colSums(L * E)^2 / (colSums(L * L) * colSums(E * E)))))
SIG <- list(E = E, L = L, R = R, A0 = A0); rm(E, L, R, A0)

## regulons as 79d
arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
groups <- fread("results/step4_v2/signature_summary.tsv")[status == "ok", unique(group_id)]
filt <- fread("results/step6/aracne_project_v2_tf/regulator_filter.tsv")[keep == TRUE & group_id %in% groups][, net := net_of[group_id]]
prune <- function(rg, cand) { r <- lapply(rg[intersect(cand, names(rg))], function(x) { k <- names(x$tfmode) %in% genes
  if (sum(k) < MIN_TARGETS || mean(k) < MIN_FRAC) return(NULL); list(tfmode = x$tfmode[k], likelihood = x$likelihood[k]) })
  r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r }
nets <- sort(unique(filt$net))
jobs <- CJ(sig = names(SIG), net = nets)
done <- mclapply(seq_len(nrow(jobs)), function(j) {
  f <- file.path(od, sprintf("%s_%s.f32", jobs$sig[j], jobs$net[j])); rf <- file.path(od, paste0(jobs$net[j], ".regulators.txt"))
  if (file.exists(f) && file.exists(rf)) return(TRUE)
  v <- viper(SIG[[jobs$sig[j]]], prune(arac[[jobs$net[j]]], unique(filt[net == jobs$net[j], regulator])), method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  writeBin(as.vector(v), f, size = 4); if (jobs$sig[j] == "E") writeLines(rownames(v), rf); TRUE
}, mc.cores = N_CORES, mc.preschedule = FALSE)
if (!all(vapply(done, isTRUE, logical(1)))) stop("scoring failed: ", paste(unlist(done[!vapply(done, isTRUE, logical(1))]), collapse = "; "))
rm(SIG); invisible(gc())
act <- function(s, net) { r <- readLines(file.path(od, paste0(net, ".regulators.txt")))
  matrix(readBin(file.path(od, sprintf("%s_%s.f32", s, net)), "numeric", size = 4, n = length(r) * n), nrow = length(r), dimnames = list(r, NULL)) }

## lineage matching (as 79n)
keys <- fread(cmd = "grep -v '^#' config/celligner_tissue_keys.tsv")
ln <- unique(keys[group_id %in% names(net_of), .(lineage = unlist(strsplit(lineages, ";")), net = net_of[group_id]), by = group_id][, .(lineage, net)])
mdl <- fread("data/raw/depmap_24q4/Model.csv"); lin <- setNames(mdl$OncotreeLineage, mdl$ModelID); cname <- setNames(mdl$CellLineName, mdl$ModelID)
nl <- I[, .N, by = ModelID]; I[, unit := fifelse(ModelID %in% nl[N >= MIN_LINE, ModelID], ModelID, "other")]

RES <- mclapply(nets, function(net) {
  a <- lapply(setNames(c("E", "L", "R", "A0"), c("E", "L", "R", "A0")), act, net = net)
  d <- a$L - a$E
  mcol <- which(paste(lin[I$ModelID], net) %in% ln[, paste(lineage, net)])
  cnt <- list(); det <- list()
  for (u in unique(I$unit)) { k <- which(I$unit == u)
    x <- d[, k, drop = FALSE]; med <- apply(x, 1, median); md <- apply(x, 1, mad); z <- (x - med) / pmax(md, 1e-9)
    q <- matrix(p.adjust(2 * pnorm(-abs(z)), "BH"), nrow(z))
    cnt[[u]] <- data.table(net = net, unit = u, tests = length(z), sig = sum(q < FDR), gain = sum(q < FDR & z > 0), loss = sum(q < FDR & z < 0))
    km <- intersect(k, mcol)
    if (length(km)) { j <- match(km, k)
      det[[u]] <- data.table(net = net, col = rep(km, each = nrow(x)), regulator = rep(rownames(x), length(km)),
                             aE = as.vector(a$E[, km]), aL = as.vector(a$L[, km]), aR = as.vector(a$R[, km]), aA0 = as.vector(a$A0[, km]),
                             z = as.vector(z[, j]), fdr = as.vector(q[, j])) } }
  bl <- rbindlist(lapply(mcol, function(i) data.table(net = net, col = i,
    rho_A0_E = cor(a$A0[, i], a$E[, i], method = "spearman"), rho_A0_L = cor(a$A0[, i], a$L[, i], method = "spearman"),
    rho_A0_R = cor(a$A0[, i], a$R[, i], method = "spearman"), rho_E_L = cor(a$E[, i], a$L[, i], method = "spearman"))))
  list(cnt = rbindlist(cnt), det = rbindlist(det), bl = bl)
}, mc.cores = min(N_CORES, length(nets)), mc.preschedule = FALSE)
bad <- vapply(RES, inherits, logical(1), "try-error"); if (any(bad)) stop(RES[[which(bad)[1]]])
CNT <- rbindlist(lapply(RES, `[[`, "cnt")); D <- rbindlist(lapply(RES, `[[`, "det")); B <- rbindlist(lapply(RES, `[[`, "bl")); rm(RES)
D[, `:=`(pert_id = I$pert_id[col], ModelID = I$ModelID[col], drug = I$cmap_name[col], line = cname[I$ModelID[col]])]
D[, onset := abs(aR) / (abs(aE) + abs(aR) + 1e-12)]
D[, change := fifelse(fdr < FDR, fifelse(z > 0, "late gain", "late loss"), "no change")]
fwrite(D[, .(line, drug, pert_id, net, regulator, aE = round(aE, 3), aL = round(aL, 3), aR = round(aR, 3), aA0 = round(aA0, 3),
             z = round(z, 2), fdr = signif(fdr, 3), onset = round(onset, 3), change)], file.path(res, "tests_matched.tsv.gz"), sep = "\t")

## summaries
ls_all <- CNT[, .(tests = sum(tests), sig = sum(sig), gain = sum(gain), loss = sum(loss)), by = unit][, unit_name := fifelse(unit == "other", "other (< 30 pairs)", cname[unit])]
ls_m <- D[, .(pairs = uniqueN(col), tests = .N, share_sig = round(mean(fdr < FDR), 4), late_gain = round(mean(change == "late gain"), 4),
              late_loss = round(mean(change == "late loss"), 4), median_onset = round(median(onset), 3)), by = .(line, net)]
fwrite(merge(ls_all[, .(unit_name, tests_all_nets = tests, share_sig_all_nets = round(sig / tests, 4))], ls_m, by.x = "unit_name", by.y = "line", all = TRUE),
       file.path(res, "line_summary.tsv"), sep = "\t")
rs <- D[, .(pairs = .N, late_gain = round(mean(change == "late gain"), 3), late_loss = round(mean(change == "late loss"), 3),
            median_onset = round(median(onset), 3), median_aE = round(median(aE), 2), median_aL = round(median(aL), 2)), by = .(net, regulator)]
fwrite(rs[order(-(late_gain + late_loss))], file.path(res, "regulator_summary.tsv"), sep = "\t")
D[, top_A0 := frank(-abs(aA0), ties.method = "first") <= 25, by = .(net, col)]
bs <- rbind(D[, .(regulators = "top 25 by |A0|", share_sig = round(mean(fdr < FDR), 4), median_onset = round(median(onset), 3),
                  share_A0_sign_eq_E = round(mean(sign(aA0) == sign(aE)), 3), share_A0_sign_eq_L = round(mean(sign(aA0) == sign(aL)), 3)), by = top_A0][top_A0 == TRUE],
            D[, .(top_A0 = NA, regulators = "all", share_sig = round(mean(fdr < FDR), 4), median_onset = round(median(onset), 3),
                  share_A0_sign_eq_E = round(mean(sign(aA0) == sign(aE)), 3), share_A0_sign_eq_L = round(mean(sign(aA0) == sign(aL)), 3))])
bb <- B[, .(median_rho_A0_E = round(median(rho_A0_E), 3), median_rho_A0_L = round(median(rho_A0_L), 3), median_rho_A0_R = round(median(rho_A0_R), 3),
            median_rho_E_L = round(median(rho_E_L), 3), pairs = .N)]
fwrite(cbind(bs[, -"top_A0"], bb), file.path(res, "baseline.tsv"), sep = "\t")
tc <- fread("results/lincs/time_dose/time_classes.tsv")[, .(pert_id, ModelID, time_class)]
ds <- merge(D[, .(sig = sum(fdr < FDR), gain = sum(change == "late gain"), loss = sum(change == "late loss"), regs = .N), by = .(col, pert_id, ModelID, drug, line, net)],
            tc, by = c("pert_id", "ModelID"), all.x = TRUE)
fwrite(ds, file.path(res, "drug_summary.tsv"), sep = "\t")
ae <- D[ModelID == "ACH-000019" & net == "BRCA" & regulator %in% c("ESR1", "GATA3", "FOXA1", "XBP1", "PGR", "MYB") &
        drug %in% c("fulvestrant", "tamoxifen", "raloxifene", "toremifene", "bazedoxifene", "clomifene", "ZK-164015", "ICI-182780")]
fwrite(ae[, .(drug, regulator, aE = round(aE, 2), aL = round(aL, 2), aR = round(aR, 2), aA0 = round(aA0, 2), z = round(z, 2), fdr = signif(fdr, 3), change)][order(regulator, drug)],
       file.path(res, "antier_mcf7.tsv"), sep = "\t")

options(width = 200)
cat("\n== per line (lineage-matched networks)\n"); print(ls_m[order(-pairs)])
cat("\n== all networks, per line unit\n"); print(ls_all[, .(unit_name, tests, share_sig = round(sig / tests, 4), gain_share = round(gain / pmax(sig, 1), 3))])
cat("\n== regulators most often changed late (>= 50 pairs)\n"); print(head(rs[pairs >= 50][order(-(late_gain + late_loss))], 25))
cat("\n== baseline (A0) vs early / late\n"); print(bs); print(bb)
cat("\n== pairs: significant changes by time class\n"); print(ds[, .(pairs = .N, median_sig = as.numeric(median(sig)), mean_share = round(mean(sig / regs), 3)), by = time_class])
cat("\n== anti-estrogens in MCF7 (BRCA network)\n"); print(ae[regulator %in% c("ESR1", "GATA3", "FOXA1", "XBP1"), .(drug, regulator, aE = round(aE, 2), aL = round(aL, 2), aR = round(aR, 2), aA0 = round(aA0, 2), z = round(z, 2), change)][order(regulator, drug)])
