#!/usr/bin/env Rscript
# Step 8b: characterize the Step 8 clusterings (scripts/80_consensus_clustering.R).
#   agreement.tsv          ARI between the main (modules10) and sensitivity (all) assignments, each at its best K and
#                          both at the main best K (robustness to the feature set)
#   known_subtypes.tsv     best-K clusters (main) x PanCanAtlas Subtype_Selected (TCGA molecular subtypes)
#   defining_modules.tsv   per cluster: module features with the largest mean difference vs the other patients, with
#                          the module's top member regulators
# Outputs in results/step8/
suppressPackageStartupMessages({library(data.table); library(mclust)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
d <- fread("results/step8/decision.tsv"); a <- fread("results/step8/assignments.tsv"); mf <- fread("results/step8/module_features.tsv")
ab <- merge(a, d[, .(group_id, feature_set, best_K)], by.x = c("group_id", "feature_set", "K"), by.y = c("group_id", "feature_set", "best_K"))
w <- dcast(ab, group_id + sample ~ feature_set, value.var = "cluster")
k2 <- merge(a, d[feature_set == "modules10", .(group_id, Kmain = best_K)], by = "group_id")[K == Kmain]
w2 <- dcast(k2, group_id + sample ~ feature_set, value.var = "cluster")
agr <- merge(w[, .(ARI_best_K = round(adjustedRandIndex(modules10, all), 2)), by = group_id],
             w2[, .(ARI_main_K = round(adjustedRandIndex(modules10, all), 2)), by = group_id], by = "group_id")
fwrite(agr, "results/step8/agreement.tsv", sep = "\t")
sub <- fread("data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz"); setnames(sub, 1, "sample")
m <- merge(ab[feature_set == "modules10"], sub[, .(sample, Subtype_Selected)], by = "sample", all.x = TRUE)
fwrite(m[, .N, by = .(group_id, cluster, Subtype_Selected)][order(group_id, cluster, -N)], "results/step8/known_subtypes.tsv", sep = "\t", na = "NA")
dm <- rbindlist(lapply(unique(d$group_id), function(g) rbindlist(lapply(c("modules10", "all"), function(set) {
  F <- fread(sprintf("data/processed/step8/features/%s_%s.tsv.gz", g, set)); M <- as.matrix(F[, -1]); rownames(M) <- F$sample
  cl <- ab[group_id == g & feature_set == set][match(rownames(M), sample), cluster]
  rbindlist(lapply(sort(unique(cl)), function(k) { dlt <- colMeans(M[cl == k, , drop = FALSE]) - colMeans(M[cl != k, , drop = FALSE])
    top <- head(order(-abs(dlt)), 8)
    data.table(group_id = g, feature_set = set, cluster = k, n = sum(cl == k), feature = colnames(M)[top], mean_diff = round(dlt[top], 2),
               top_members = sapply(colnames(M)[top], function(f) paste(head(mf[group_id == g & feature_set == set & feature == f, regulator], 8), collapse = ","))) }))
}))))
fwrite(dm, "results/step8/defining_modules.tsv", sep = "\t")
message("STEP 8 CHARACTERIZATION DONE")
