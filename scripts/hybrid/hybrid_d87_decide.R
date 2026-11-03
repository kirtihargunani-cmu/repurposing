#!/usr/bin/env Rscript
# D87 decision (rule fixed 2026-11-02 before running, docs/DECISIONS.md), from scripts/hybrid/hybrid_d87.R outputs.
# Adopt H1 as the curated half of the drug-side consensus if, vs CollecTRI on the D84 main-analysis signatures:
#   (1) correct share +>= 2 points; (2) sign test over signatures where they differ p < 0.05; (3) median per-TF
#   difference >= 0 (TFs with >= 3 signatures); (4) ESR1 negative in most MCF7 anti-estrogen profiles and the ligand
#   check no more than 5 points below CollecTRI. Everything else is reported only.
# Outputs: results/lincs/hybrid/{summary.tsv, tf.tsv, decision.tsv, flips.tsv, ligand_summary.tsv}
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../.."))
setwd(root); res <- "results/lincs/hybrid"
MIN_GAIN <- 0.02; P_MAX <- 0.05; MIN_SIGS <- 3; LIG_TOL <- 0.05; MCF7 <- "ACH-000019"

S <- rbind(fread(file.path(res, "signatures_specific.tsv.gz")),
           fread(file.path(res, "signatures_raw.tsv.gz"))[variant != "CollecTRI"])
d84 <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[set == "aracne_v2_tf", .(sig_id, nes_a = nes)]
lm <- fread("data/processed/lincs/level5_landmark_genes.tsv")
lmsym <- lm[[intersect(c("gene_symbol", "symbol", "pr_gene_symbol"), names(lm))[1]]]
W <- dcast(S, sig_id + technology + direction + cmap_name + cell_iname + key ~ variant, value.var = "nes")
W <- merge(W, d84, by = "sig_id", all.x = TRUE)
VARS <- setdiff(unique(S$variant), "CollecTRI")
W <- W[!is.na(CollecTRI)]
ok <- function(x, d) sign(x) == d
cons <- function(a, b) fifelse(is.na(a), b, fifelse(is.na(b), a, (a + b) / 2))

cmp <- function(w, v, base = "CollecTRI", label = "") {
  w <- w[!is.na(get(v)) & !is.na(get(base))]
  h <- ok(w[[v]], w$direction); c <- ok(w[[base]], w$direction)
  win <- sum(h & !c); los <- sum(!h & c)
  t <- w[, .(h = mean(ok(get(v), direction)), c = mean(ok(get(base), direction)), n = .N), by = cmap_name][n >= MIN_SIGS]
  data.table(comparison = label, variant = v, base = base, signatures = nrow(w), tfs = uniqueN(w$cmap_name),
             correct_variant = mean(h), correct_base = mean(c), diff = mean(h) - mean(c), wins = win, losses = los,
             p_sign = if (win + los) binom.test(win, win + los)$p.value else NA_real_,
             tfs_3plus = nrow(t), median_tf_diff = if (nrow(t)) median(t$h - t$c) else NA_real_,
             tf_better = sum(t$h > t$c), tf_worse = sum(t$h < t$c))
}
sm <- rbindlist(c(
  lapply(VARS, function(v) cmp(W, v, label = "all")),
  lapply(c("knockdown", "overexpression"), function(k) cmp(W[(direction < 0) == (k == "knockdown")], "H1", label = k)),
  list(cmp(W[cmap_name %in% lmsym], "H1", label = "TF measured on L1000"),
       cmp(W[!cmap_name %in% lmsym], "H1", label = "TF not measured on L1000"),
       cmp(W[key != "POOLED"], "H1", label = "lines with own directions"),
       cmp(W[key == "POOLED"], "H1", label = "lines with pooled directions"))))
W[, `:=`(cons_C = cons(nes_a, CollecTRI), cons_H1 = cons(nes_a, H1))]
sm <- rbind(sm, cmp(W, "cons_H1", "cons_C", "consensus with ARACNe"))
sm[, c("correct_variant", "correct_base", "diff", "median_tf_diff") := lapply(.SD, round, 4), .SDcols = c("correct_variant", "correct_base", "diff", "median_tf_diff")]
sm[, p_sign := signif(p_sign, 3)]
fwrite(sm, file.path(res, "summary.tsv"), sep = "\t")

tf <- melt(W, id.vars = c("sig_id", "cmap_name", "direction"), measure.vars = c("CollecTRI", VARS), variable.name = "variant", value.name = "nes")[!is.na(nes)]
tf <- tf[, .(n_sigs = .N, correct = round(mean(ok(nes, direction)), 3), median_nes = round(median(nes * direction), 3)), by = .(cmap_name, variant)]
fwrite(dcast(tf, cmap_name ~ variant, value.var = c("n_sigs", "correct")), file.path(res, "tf.tsv"), sep = "\t")

## ligands and ESR1
L <- rbind(fread(file.path(res, "ligands_specific.tsv")), fread(file.path(res, "ligands_raw.tsv")))
LW <- dcast(L, arm + drug + tf + direction + ModelID ~ variant, value.var = "nes")
lsum <- rbindlist(lapply(intersect(c(VARS), names(LW)), function(v) {
  x <- LW[!is.na(get(v)) & !is.na(CollecTRI)]
  x[, .(variant = v, profiles = .N, correct_variant = round(mean(ok(get(v), direction)), 4), correct_collectri = round(mean(ok(CollecTRI, direction)), 4)), by = arm] }))
fwrite(lsum, file.path(res, "ligand_summary.tsv"), sep = "\t")
esr <- LW[tf == "ESR1" & direction < 0 & ModelID == MCF7]
esr_s <- esr[arm == "specific"]

## flips
E <- rbind(cbind(arm = "specific", fread(file.path(res, "edges_specific.tsv.gz"))), cbind(arm = "raw", fread(file.path(res, "edges_raw.tsv.gz"))))
fl <- E[!is.na(r), .(links = .N, flipped_H1 = round(mean(flip), 4), sign_disagrees = round(mean(sign(r) != prior), 4),
                     median_abs_r = round(median(abs(r)), 3)), by = .(arm, line)]
foc <- E[arm == "specific" & regulator %in% c("ESR1", "AR", "FOXA1", "GATA3") & (line == MCF7 | line == "POOLED") & !is.na(r),
         .(links = .N, flipped_H1 = sum(flip), sign_disagrees = sum(sign(r) != prior), median_r_signed = round(median(r * prior), 3)), by = .(arm, regulator, line)]
fwrite(rbind(fl, foc, fill = TRUE), file.path(res, "flips.tsv"), sep = "\t")

## decision
a <- sm[comparison == "all" & variant == "H1"]; lh <- lsum[arm == "specific" & variant == "H1"]
crit <- data.table(
  criterion = c("1 gain >= 2 points", "2 sign test p < 0.05", "3 median per-TF diff >= 0", "4a ESR1 negative in most MCF7 anti-estrogen profiles", "4b ligands no more than 5 points below CollecTRI"),
  value = c(sprintf("%+.1f points (%.1f%% vs %.1f%%, %d signatures)", 100 * a$diff, 100 * a$correct_variant, 100 * a$correct_base, a$signatures),
            sprintf("%d wins, %d losses, p = %.3g", a$wins, a$losses, a$p_sign),
            sprintf("%+.3f over %d TFs (%d better, %d worse)", a$median_tf_diff, a$tfs_3plus, a$tf_better, a$tf_worse),
            sprintf("%d of %d negative (median NES %.2f; CollecTRI %d of %d, %.2f)", sum(esr_s$H1 < 0, na.rm = TRUE), sum(!is.na(esr_s$H1)), median(esr_s$H1, na.rm = TRUE),
                    sum(esr_s$CollecTRI < 0, na.rm = TRUE), sum(!is.na(esr_s$CollecTRI)), median(esr_s$CollecTRI, na.rm = TRUE)),
            sprintf("%.1f%% vs %.1f%% (%d profiles)", 100 * lh$correct_variant, 100 * lh$correct_collectri, lh$profiles)),
  met = c(a$diff >= MIN_GAIN, isTRUE(a$p_sign < P_MAX) && a$wins > a$losses, isTRUE(a$median_tf_diff >= 0),
          mean(esr_s$H1 < 0, na.rm = TRUE) > 0.5, lh$correct_variant >= lh$correct_collectri - LIG_TOL))
crit <- rbind(crit, data.table(criterion = "DECISION", value = if (all(crit$met)) "adopt H1 as the curated half of the consensus" else "keep CollecTRI", met = all(crit$met)))
fwrite(crit, file.path(res, "decision.tsv"), sep = "\t")

options(width = 200)
cat("\n== benchmark (vs CollecTRI unless noted)\n"); print(sm)
cat("\n== ligands\n"); print(lsum)
cat("\n== ESR1, anti-estrogens in MCF7\n"); print(esr[, .(arm, drug, CollecTRI, H1, H2 = if ("H2" %in% names(esr)) H2 else NA, H1_raw = if ("H1_raw" %in% names(esr)) H1_raw else NA)])
cat("\n== flips\n"); print(fl); print(foc)
cat("\n== DECISION\n"); print(crit)
