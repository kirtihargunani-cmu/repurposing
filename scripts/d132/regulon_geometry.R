#!/usr/bin/env Rscript
# D132 diagnostic (D132.1 item 5): is the inflation a simple function of regulon SIZE and internal COHERENCE?
# For every regulator actually scored by the null run, computes from the SAME regulons and the SAME null profiles:
#   n_targets      targets present on the LINCS gene space (what 79d scores on)
#   eff_n_weights  (sum|lik|)^2 / sum lik^2, the weight-concentration count (NOT a correlation measure)
#   mean_abs_corr  mean |pairwise correlation| among the regulon's targets across the null profiles
#   quad_ratio     sqrt(w' R w / w' w), w = tfmode * likelihood: the ANALYTIC prediction of the inflation from the
#                  target-target correlation alone. If this tracks the measured inflation, an analytic correction is
#                  available and no empirical NES null is needed in production.
# It is an approximation of aREA (which rank-transforms and combines a one- and a two-tail component), so it is a
# diagnostic for the measured inflation, never a substitute for it.
suppressPackageStartupMessages({library(data.table)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../.."))
setwd(root)
NULLNAME <- Sys.getenv("NULLNAME", "replicate"); TAGN <- substr(NULLNAME, 1, 3)
proc <- "data/processed/lincs"
genes <- fread(file.path(proc, "drug_line_profiles_genes.tsv"))
idx <- fread(file.path(proc, sprintf("d132/%s_index.tsv", NULLNAME)))
P <- matrix(readBin(file.path(proc, sprintf("d132/%s_profiles.f32", NULLNAME)), "numeric", size = 4,
                    n = nrow(genes) * nrow(idx)), nrow = nrow(genes))
rownames(P) <- genes$gene_id
cat(sprintf("null %s: %d genes x %d profiles\n", NULLNAME, nrow(P), ncol(P)))
Z <- t(scale(t(P)))                       # genes standardised across null profiles
Z[!is.finite(Z)] <- 0
m <- ncol(Z)

out <- list()
for (SET in c("tf", "collectri")) {
  if (SET == "tf") {
    arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
    dir <- file.path(proc, paste0("drug_viper_specific_v2_tf_d132", TAGN))
  } else {
    arac <- list(ALL = readRDS("results/step5/regulons_core.rds")[["collectri"]])
    dir <- file.path(proc, paste0("drug_viper_specific_collectri_d132", TAGN))
  }
  for (f in list.files(dir, pattern = "\\.regulators\\.txt$", full.names = TRUE)) {
    net <- sub("\\.regulators\\.txt$", "", basename(f))
    regs <- readLines(f)
    ra <- arac[[net]]
    rows <- rbindlist(lapply(regs, function(r) {
      o <- ra[[r]]; k <- names(o$tfmode) %in% rownames(Z)
      tg <- names(o$tfmode)[k]; w <- as.numeric(o$tfmode[k]) * as.numeric(o$likelihood[k])
      n <- length(tg)
      if (n < 2) return(NULL)
      Zt <- Z[tg, , drop = FALSE]
      v <- as.numeric(w %*% Zt) / sqrt(m - 1)          # w' Z' , so var = w' R w
      quad <- sqrt(sum(v^2) / sum(w^2))
      # mean |corr| among targets, exact for small regulons and on a seeded 300-target subsample for large ones
      sel <- if (n > 300) { set.seed(132); sort(sample.int(n, 300)) } else seq_len(n)
      C <- tcrossprod(Zt[sel, , drop = FALSE]) / (m - 1)
      p <- length(sel)
      data.table(set = SET, net = net, regulator = r, n_targets = n,
                 eff_n_weights = sum(abs(w))^2 / sum(w^2),
                 mean_abs_corr = (sum(abs(C)) - sum(abs(diag(C)))) / (p * (p - 1)),
                 quad_ratio = quad)
    }))
    out[[paste(SET, net)]] <- rows
  }
  cat(sprintf("%s: %d networks done\n", SET, length(list.files(dir, pattern = "\\.regulators\\.txt$"))))
}
D <- rbindlist(out)
dir.create("results/lincs/d132", showWarnings = FALSE, recursive = TRUE)
fwrite(D, sprintf("results/lincs/d132/regulon_geometry_%s.tsv.gz", NULLNAME), sep = "\t")
cat(sprintf("wrote %d (set, net, regulator) rows\n", nrow(D)))
