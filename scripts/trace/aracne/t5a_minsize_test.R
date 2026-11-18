#!/usr/bin/env Rscript
# D138.3: is a 10-target floor the right place? MEASURED, not asserted.
#
# D16 set ">= 10 targets and >= 50% of targets measured" with the rationale "avoid scoring
# regulators that are absent or barely measured". No empirical basis for 10 was recorded, and
# VIPER's own default is minsize = 25, so our floor is MORE permissive than the package default.
#
# Test: take regulators whose BRCA_LUMA regulon has plenty of LINCS-space targets, subsample the
# regulon down to k targets, run the real VIPER on the real MCF7 drug profiles, and ask how well
# the k-target activity reproduces the full-regulon activity. Reliability is the correlation across
# the 1,221 MCF7 compounds, per regulator, median over 10 independent subsamples.
suppressPackageStartupMessages({library(data.table); library(viper)})
setwd("/home/kirti/projects/repurposing")
SCR <- "/tmp/claude-1014/-home-kirti/ff915099-39db-42da-910a-c733ac0b89dd/scratchpad"
set.seed(20261117)
KS <- c(3, 5, 8, 10, 15, 25, 50, 100); NREP <- 10; NET <- "BRCA_LUMA"; N_FULL <- 150

P <- as.matrix(fread(file.path(SCR, "mcf7_profiles.tsv")), rownames = "gene")
R <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")[[NET]]
R <- lapply(R, function(r) { k <- names(r$tfmode) %in% rownames(P)
                             list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
sz <- vapply(R, function(r) length(r$tfmode), integer(1))
big <- names(sz)[sz >= N_FULL]
cat(sprintf("network %s: %d regulators, %d with >= %d LINCS-space targets\n",
            NET, length(R), length(big), N_FULL))

run <- function(reg) { class(reg) <- "regulon"
                       viper(P, reg, method = "none", minsize = 1, verbose = FALSE) }
full <- run(R[big])
out <- rbindlist(lapply(KS, function(k) {
  rbindlist(lapply(seq_len(NREP), function(i) {
    sub <- lapply(R[big], function(r) { j <- sample.int(length(r$tfmode), k)
                                        list(tfmode = r$tfmode[j], likelihood = r$likelihood[j]) })
    nes <- run(sub)
    g <- intersect(rownames(full), rownames(nes))
    data.table(k = k, rep = i, regulator = g,
               r = vapply(g, function(x) cor(full[x, ], nes[x, ], method = "spearman"), numeric(1)))
  }))
}))
s <- out[, .(median_r = median(r), q10 = quantile(r, .1), frac_above_0.8 = mean(r > 0.8)),
         by = k][order(k)]
cat("\nreliability of a k-target regulon vs its own full regulon (>= ", N_FULL,
    " targets), across 1,221 MCF7 compounds:\n", sep = "")
print(s)
fwrite(out, "results/lincs/trace_aracne/minsize_reliability_raw.tsv.gz", sep = "\t")
fwrite(s, "results/lincs/trace_aracne/minsize_reliability.tsv", sep = "\t")
cat("\nregulators tested:", length(big), " reps per k:", NREP, "\n")
