#!/usr/bin/env Rscript
# Step 6d: cross-reference group-level VIPER results from the three regulon sets (CollecTRI, DoRothEA A-C,
# project ARACNe). Same signatures, null and filters in all three; only the regulons differ.
#   per_group.tsv   per group: regulators tested / passing FDR < FDR_CUT per set; pairwise Spearman of NES on
#                   regulators tested in both sets; consensus = passing in >= 2 sets with the same sign
#   consensus.tsv   per group x regulator passing in >= 2 sets (same sign): NES and FDR in each set
#   controls.tsv    known-biology controls side by side (verdict per set)
# FDR_CUT=0.10 (default). Output: results/step6/comparison/
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
CUT <- as.numeric(Sys.getenv("FDR_CUT", "0.10")); SETS <- c("collectri", "dorothea_ABC", "aracne_project")
out <- "results/step6/comparison"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

mr <- rbindlist(lapply(SETS, function(s) fread(file.path("results/step6", s, "group_mr.tsv"))[
  , .(group_id, regulator, nes, fdr, set = s)]))
w <- dcast(mr, group_id + regulator ~ set, value.var = c("nes", "fdr"))
pass <- function(s) !is.na(w[[paste0("fdr_", s)]]) & w[[paste0("fdr_", s)]] < CUT
sgn <- function(s) sign(w[[paste0("nes_", s)]])
P <- sapply(SETS, pass); S <- sapply(SETS, sgn)
# consensus: passes in >= 2 sets and every passing set agrees on direction
npass <- rowSums(P)
agree <- vapply(seq_len(nrow(w)), function(i) { s <- S[i, P[i, ]]; length(unique(s)) == 1 }, logical(1))
w[, `:=`(n_sets_pass = npass, consensus = npass >= 2 & agree)]

rho <- function(a, b) w[, { x <- get(paste0("nes_", a)); y <- get(paste0("nes_", b)); ok <- !is.na(x) & !is.na(y)
  .(n = sum(ok), rho = if (sum(ok) >= 10) round(cor(x[ok], y[ok], method = "spearman"), 2) else NA_real_) }, by = group_id]
pr <- Reduce(function(a, b) merge(a, b, by = "group_id"), list(
  setnames(rho("collectri", "dorothea_ABC"), c("n", "rho"), c("n_shared_ct_do", "rho_ct_do")),
  setnames(rho("collectri", "aracne_project"), c("n", "rho"), c("n_shared_ct_ar", "rho_ct_ar")),
  setnames(rho("dorothea_ABC", "aracne_project"), c("n", "rho"), c("n_shared_do_ar", "rho_do_ar"))))
cnt <- w[, .(tested_ct = sum(!is.na(nes_collectri)), tested_do = sum(!is.na(nes_dorothea_ABC)),
             tested_ar = sum(!is.na(nes_aracne_project)),
             pass_ct = sum(fdr_collectri < CUT, na.rm = TRUE), pass_do = sum(fdr_dorothea_ABC < CUT, na.rm = TRUE),
             pass_ar = sum(fdr_aracne_project < CUT, na.rm = TRUE), consensus = sum(consensus),
             consensus_all3 = sum(n_sets_pass == 3 & consensus),
             top_consensus = paste(head(regulator[consensus][order(-abs(nes_aracne_project[consensus]))], 8), collapse = ",")),
         by = group_id]
fwrite(merge(cnt, pr, by = "group_id")[order(-consensus)], file.path(out, "per_group.tsv"), sep = "\t")
fwrite(w[consensus == TRUE][order(group_id, -n_sets_pass, -abs(nes_aracne_project))], file.path(out, "consensus.tsv"), sep = "\t", na = "NA")

ctl <- rbindlist(lapply(SETS, function(s) fread(file.path("results/step6", s, "controls.tsv"))[, set := s]))
ctl <- dcast(unique(ctl[, .(group_id, regulator, why, set, verdict)]), group_id + regulator + why ~ set, value.var = "verdict")
fwrite(ctl, file.path(out, "controls.tsv"), sep = "\t", na = "NA")
message(sprintf("COMPARISON DONE: %d groups; consensus regulators (>= 2 sets, FDR < %.2f): %d", nrow(cnt), CUT, sum(w$consensus)))
