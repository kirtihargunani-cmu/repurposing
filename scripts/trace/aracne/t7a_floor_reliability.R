#!/usr/bin/env Rscript
# D139-S1: does raising the regulon-size floor buy RELIABILITY, and where does it stop paying?
#
# D138.3 measured reliability by SUBSAMPLING a large regulon to k targets. That answers "how good
# is a k-target regulon", but raising the floor does not shrink anyone's regulon - it removes the
# small ones. So the question here is different: as the floor rises, does the DRUG SCORE become
# more reproducible, and what does it cost in coverage?
#
# Method: split-half. Each regulon is split into two random halves; VIPER is run on each half, so
# NES_A and NES_B are two independent estimates of the same regulator's activity on the same
# profiles. Per regulator, reliability = cor(NES_A, NES_B) across profiles, Spearman-Brown
# corrected to the full regulon. Per SCORE, rev_A and rev_B are built from the two halves and
# correlated - that is the reproducibility of the drug score itself at a given floor.
# Split-half halves the regulon, so these are conservative (lower bounds) in absolute terms; the
# comparison ACROSS floors is the point.
suppressPackageStartupMessages({library(data.table); library(viper)})
setwd("/home/kirti/projects/repurposing")
SCR <- "/tmp/claude-1014/-home-kirti/ff915099-39db-42da-910a-c733ac0b89dd/scratchpad"
set.seed(20261117)
NETS <- c("BRCA_LUMA", "SKCM", "COADREAD"); NREP <- 5
FLOORS <- c(10, 25, 50, 75, 100, 150, 200)

P <- as.matrix(fread(file.path(SCR, "profiles2k.tsv")), rownames = "gene")
MR <- fread("results/step6/aracne_project_v2_tf/group_mr.tsv")
G4N <- c(BRCA_LUMA = "BRCA-LumA", SKCM = "TCGA-SKCM", COADREAD = "TCGA-COAD")
ALL <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
res <- list(); rel <- list()

for (NET in NETS) {
  R <- lapply(ALL[[NET]], function(r) { k <- names(r$tfmode) %in% rownames(P)
                                        list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  sz <- vapply(R, function(r) length(r$tfmode), integer(1))
  R <- R[sz >= 10]; sz <- sz[sz >= 10]
  m <- MR[group_id == G4N[[NET]]]
  w <- setNames(m$nes, m$regulator)
  cat(sprintf("\n%s: %d regulators >=10 targets, %d with a disease weight\n",
              NET, length(R), sum(names(R) %in% names(w))))
  for (rep in seq_len(NREP)) {
    halves <- lapply(R, function(r) { n <- length(r$tfmode); i <- sample.int(n)
      a <- i[seq_len(floor(n/2))]; b <- i[(floor(n/2)+1):n]
      list(A = list(tfmode = r$tfmode[a], likelihood = r$likelihood[a]),
           B = list(tfmode = r$tfmode[b], likelihood = r$likelihood[b])) })
    run <- function(w2) { x <- lapply(halves, `[[`, w2); class(x) <- "regulon"
                          viper(P, x, method = "none", minsize = 1, eset.filter = FALSE, verbose = FALSE) }
    A <- run("A"); B <- run("B")
    g <- intersect(rownames(A), rownames(B))
    rr <- vapply(g, function(x) cor(A[x, ], B[x, ], method = "spearman"), numeric(1))
    rel[[length(rel)+1]] <- data.table(net = NET, rep = rep, regulator = g, size = sz[g],
                                       r_half = rr, r_full_sb = 2*rr/(1+rr))
    for (f in FLOORS) {
      k <- g[sz[g] >= f & g %in% names(w)]
      if (length(k) < 2) next
      s <- w[k]; nrm <- sqrt(sum(s^2))
      revA <- -as.vector(t(A[k, , drop=FALSE]) %*% s) / nrm
      revB <- -as.vector(t(B[k, , drop=FALSE]) %*% s) / nrm
      res[[length(res)+1]] <- data.table(net = NET, rep = rep, floor = f, n_regulators = length(k),
                                         weight_mass = sum(s^2)/sum(w[intersect(g, names(w))]^2),
                                         score_reliability = cor(revA, revB, method = "spearman"))
    }
  }
}
RES <- rbindlist(res); REL <- rbindlist(rel)
fwrite(RES, "results/lincs/trace_aracne/floor_score_reliability.tsv", sep = "\t")
fwrite(REL, "results/lincs/trace_aracne/floor_regulator_reliability.tsv.gz", sep = "\t")

cat("\n=== SCORE reproducibility vs floor (split-half, median over", NREP, "reps) ===\n")
print(RES[, .(n_regulators = round(median(n_regulators)),
              pct_weight_kept = round(100*median(weight_mass), 1),
              score_reliability = round(median(score_reliability), 4)), by = .(net, floor)][order(net, floor)])
cat("\n=== individual REGULATOR reliability by its own size (Spearman-Brown corrected) ===\n")
REL[, bin := cut(size, c(9,24,49,99,199,1e6), labels=c("10-24","25-49","50-99","100-199","200+"))]
print(REL[, .(regulators = uniqueN(regulator), median_r_full = round(median(r_full_sb), 3),
              pct_above_0.8 = round(100*mean(r_full_sb > 0.8), 1)), by = bin][order(bin)])
