#!/usr/bin/env Rscript
# Step 10 v2: what defines each kept subgroup (D79), from the Step 8 features it was clustered on (v2 of scripts/81).
# Per subgroup, both feature sets of its source clustering (results/step10_v2/strata.tsv: step8_v2_tf_no_c8, or
# step8_v2_tf for COAD), K = 2: every feature's mean difference subgroup minus sibling (features are z-scored, so the
# difference is in SD units) and a Wilcoxon p (BH within subgroup x feature set).
#   size_class   "module10" (module of >= 10 regulators: also a modules10 feature), "small module" (2-9 regulators) or
#                "single" (one regulator): small modules and singles exist only in the all set, the finer programs
#                between patients (user 2026-10-28)
#   members      the feature's regulators (first 10), hub = the member most correlated with the module score
# Outputs: results/step10_v2/defining_features.tsv (all features), defining_top.tsv (per subgroup x feature set: top 10 by
#          |difference| with FDR < 0.05, plus the top 5 small modules / singles of the all set)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
out <- "results/step10_v2"
st <- fread(file.path(out, "strata.tsv")); sp <- fread(file.path(out, "stratum_patients.tsv"))
feat_dir <- function(src) sub("^results/(step8_v2_tf)(.*)$", "data/processed/\\1/features\\2", src)

res <- rbindlist(lapply(unique(st$group_id), function(g) {
  s <- st[group_id == g]; src <- s$source[1]
  mf <- fread(file.path(src, "module_features.tsv"))[group_id == g]
  rbindlist(lapply(c("modules10", "all"), function(set) {
    F <- fread(file.path(feat_dir(src), sprintf("%s_%s.tsv.gz", g, set))); M <- as.matrix(F[, -1]); rownames(M) <- F$sample
    lab <- sp[match(rownames(M), sample), stratum]
    sz <- mf[feature_set == set, .N, by = feature]; mem <- mf[feature_set == set, .(members = paste(head(regulator, 10), collapse = ",")), by = feature]
    rbindlist(lapply(s$stratum, function(k) {
      a <- lab == k; b <- !is.na(lab) & lab != k
      p <- apply(M, 2, function(x) suppressWarnings(wilcox.test(x[a], x[b])$p.value))
      data.table(group_id = g, stratum = k, feature_set = set, feature = colnames(M),
                 mean_diff = round(colMeans(M[a, , drop = FALSE]) - colMeans(M[b, , drop = FALSE]), 3), p = signif(p, 3))[
        , fdr := signif(p.adjust(p, "BH"), 3)][
        , n_regulators := sz[match(feature, sz$feature), N]][
        , size_class := fifelse(n_regulators >= 10, "module10", fifelse(n_regulators >= 2, "small module", "single"))][
        , members := mem[match(feature, mem$feature), members]]
    }))
  }))
}))
fwrite(res[order(stratum, feature_set, -abs(mean_diff))], file.path(out, "defining_features.tsv"), sep = "\t")
top <- rbind(res[fdr < 0.05][order(-abs(mean_diff)), head(.SD, 10), by = .(stratum, feature_set)][, pick := "top"],
             res[feature_set == "all" & size_class != "module10" & fdr < 0.05][order(-abs(mean_diff)), head(.SD, 5), by = stratum][, pick := "small program"])
fwrite(top[order(stratum, pick, feature_set, -abs(mean_diff))], file.path(out, "defining_top.tsv"), sep = "\t")
print(res[, .(features = .N, fdr05 = sum(fdr < 0.05), small_fdr05 = sum(fdr < 0.05 & size_class != "module10")), by = .(stratum, feature_set)])
message("SUBGROUP CHARACTERIZATION DONE")
