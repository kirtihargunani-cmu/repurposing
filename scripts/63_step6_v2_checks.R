#!/usr/bin/env Rscript
# Step 6 v2 acceptance checks (SPEC_rebuild_v2.md, "Acceptance checks"; reported whatever they show, D19).
#   C2  GC-only signature vs group MR NES: each group's kept regulons (regulator_filter.tsv) scored with aREA on a
#       signature = minus z-scored gene GC (longest transcript), Spearman with the group's msviper NES. v1 recomputed the
#       same way for comparison. High = regulator ranking predictable from GC alone.
#   C5  known-biology controls: scripts/51_compare_mr.R with REGULON=aracne_project_v2 (controls.tsv); summarised here.
#   C7  master regulators at FDR < 0.10 (and < 0.05) per group, v1 vs v2: counts, Jaccard, NES Spearman on
#       regulators scored in both.
#   C9  network size: edges and regulators per network before / after the cross-mappable filter, v1 vs v2
#       (bootstrap support needs the bootstrap files, kept on the build machine: scripts/45 there).
#   C10 mappability: regulators / targets with gene mappability < 0.8, and cross-mappable edges left in the v2
#       regulons (both expected exactly 0).
# SET=full (default): reference = v1 (aracne_project), test = v2 full set (aracne_project_v2) -> results/step6/checks_v2/
# SET=tf: reference = v2 full set, test = v2 TF-only set (aracne_project_v2_tf) -> results/step6/checks_v2_tf/
#   (v1 had no TF-only set; columns *_ref / *_test then mean full / TF-only)
# Output: {c2_gc,c5_controls,c5_controls_changed,c7_mr,c9_networks,c10_mappability,summary}.tsv
suppressPackageStartupMessages({library(data.table); library(viper)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
SET <- Sys.getenv("SET", "full"); stopifnot(SET %in% c("full", "tf"))
OUT <- if (SET == "tf") "results/step6/checks_v2_tf" else "results/step6/checks_v2"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
V <- if (SET == "tf") c(ref = "aracne_project_v2", test = "aracne_project_v2_tf") else c(ref = "aracne_project", test = "aracne_project_v2")
regs <- lapply(V, function(v) readRDS(paste0("data/processed/step5/regulons_", v, ".rds")))
mr <- lapply(V, function(v) fread(file.path("results/step6", v, "group_mr.tsv")))
filt <- lapply(V, function(v) fread(file.path("results/step6", v, "regulator_filter.tsv"))[keep == TRUE])
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
groups <- intersect(unique(mr$ref$group_id), unique(mr$test$group_id))

# C2: GC-only signature
gc <- fread("data/interim/gene_gc_v23.tsv")[!is.na(gc_longest_tx)]
# sign as the Step 4 tilt (lower GC = higher in tumor-vs-GTEx z): rho > 0 = MR ranking follows GC
gsig <- matrix(-as.numeric(scale(gc$gc_longest_tx)), dimnames = list(gc$gene_id, "gc"))
c2 <- rbindlist(lapply(names(V), function(v) rbindlist(lapply(groups, function(gid) {
  reg <- regs[[v]][[net_of[[gid]]]][filt[[v]][group_id == gid, regulator]]
  reg <- lapply(reg, function(r) { k <- names(r$tfmode) %in% rownames(gsig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  reg <- reg[lengths(lapply(reg, `[[`, "tfmode")) >= 10]; class(reg) <- "regulon"
  es <- viper::aREA(gsig, reg, minsize = 10)$nes[, 1]
  m <- mr[[v]][group_id == gid & regulator %in% names(es)]
  data.table(version = v, group_id = gid, n = nrow(m), rho_gc = round(cor(es[m$regulator], m$nes, method = "spearman"), 3),
             gc_sign_agree_fdr10 = round(mean(sign(es[m[fdr < 0.10, regulator]]) == sign(m[fdr < 0.10, nes])), 3))
}))))
c2w <- dcast(c2, group_id ~ version, value.var = c("rho_gc", "gc_sign_agree_fdr10"))
fwrite(c2w, file.path(OUT, "c2_gc.tsv"), sep = "\t", na = "NA")

# C5: controls from 51_compare_mr.R
ctl <- rbindlist(lapply(names(V), function(v) {
  f <- file.path("results/step6", V[[v]], "controls.tsv"); if (!file.exists(f)) return(NULL)
  cbind(version = v, fread(f))
}), fill = TRUE)
if (nrow(ctl)) {
  ctl[, kind := fifelse(why == "proliferation", "proliferation", "tissue-specific")]
  c5 <- dcast(ctl[, .N, by = .(version, kind, verdict)], kind + verdict ~ version, value.var = "N", fill = 0)
  fwrite(c5, file.path(OUT, "c5_controls.tsv"), sep = "\t")
  ch <- dcast(ctl, group_id + regulator + why ~ version, value.var = "verdict")[ref != test]
  fwrite(ch, file.path(OUT, "c5_controls_changed.tsv"), sep = "\t", na = "NA")
}

# C7: master regulators v1 vs v2
c7 <- rbindlist(lapply(groups, function(gid) {
  a <- mr$ref[group_id == gid]; b <- mr$test[group_id == gid]; com <- intersect(a$regulator, b$regulator)
  row <- data.table(group_id = gid, n_scored_ref = nrow(a), n_scored_test = nrow(b),
                    rho_nes = round(cor(a[match(com, regulator), nes], b[match(com, regulator), nes], method = "spearman"), 3))
  for (cut in c(0.10, 0.05)) {
    s1 <- a[fdr < cut, regulator]; s2 <- b[fdr < cut, regulator]; tg <- sprintf("fdr%02d", round(cut * 100))
    row[, (paste0("n_mr_ref_", tg)) := length(s1)][, (paste0("n_mr_test_", tg)) := length(s2)]
    row[, (paste0("jaccard_", tg)) := if (length(union(s1, s2))) round(length(intersect(s1, s2)) / length(union(s1, s2)), 3) else NA_real_]
  }
  row
}))
fwrite(c7, file.path(OUT, "c7_mr.tsv"), sep = "\t", na = "NA")

# C9: network size
size <- function(reg) data.table(regulators = length(reg), edges = sum(lengths(lapply(reg, `[[`, "tfmode"))))
rm_tab <- fread(paste0("results/step5/crossmap_edges_removed_v2", if (SET == "tf") "_tf" else "", ".tsv"))
c9 <- merge(rbindlist(lapply(names(regs$ref), function(n) cbind(network = n, size(regs$ref[[n]])))),
            rbindlist(lapply(names(regs$test), function(n) cbind(network = n, size(regs$test[[n]])))),
            by = "network", all = TRUE, suffixes = c("_ref", "_test"))
c9 <- merge(c9, rm_tab[, .(network, edges_test_before_filter = edges, cross_mappable_removed = cross_mappable, frac_removed = frac)],
            by = "network", all = TRUE)
fwrite(c9, file.path(OUT, "c9_networks.tsv"), sep = "\t", na = "NA")

# C10: mappability in the test regulons
gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
# _PAR_Y rows (Y copies of pseudoautosomal genes, mappability 0) skipped: the X copy is the gene (as scripts/22)
lowmap <- sub("\\..*", "", gm[!is.na(map) & map < 0.8 & !grepl("_PAR_Y", id), id])
cm <- fread(cmd = "zcat data/raw/crossmap/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz", header = FALSE,
            col.names = c("a", "b", "s"))[s > 0]
cm[, `:=`(a = sub("\\..*", "", a), b = sub("\\..*", "", b), s = NULL)]
cm <- unique(rbind(cm, cm[, .(a = b, b = a)])); setkey(cm, a, b)
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
c10 <- rbindlist(lapply(names(regs$test), function(n) {
  reg <- regs$test[[n]]; k <- lengths(lapply(reg, `[[`, "tfmode"))
  ens <- ifelse(startsWith(names(reg), "ENSG"), names(reg), sym2ens[names(reg)])
  ed <- data.table(a = rep(ens, k), b = unlist(lapply(reg, function(r) names(r$tfmode)), use.names = FALSE))
  data.table(network = n, regulators_unmapped_id = sum(is.na(ens)), regulators_lowmap = sum(ens %in% lowmap),
             targets_lowmap = uniqueN(ed$b[ed$b %in% lowmap]),
             cross_mappable_edges = sum(!is.na(cm[ed, on = .(a, b), mult = "first", which = TRUE])))
}))
fwrite(c10, file.path(OUT, "c10_mappability.tsv"), sep = "\t")

med <- function(x) round(median(x, na.rm = TRUE), 3)
summ <- data.table(
  check = c("C2 GC-only vs MR NES, median rho", "C2 groups with |rho| > 0.4",
            "C7 MRs FDR < 0.10, median per group", "C7 Jaccard ref-test FDR < 0.10, median", "C7 NES rho ref-test, median",
            "C9 edges per network, median", "C9 cross-mappable edges removed, median fraction",
            "C10 low-mappability regulators + targets", "C10 cross-mappable edges left"),
  ref = c(med(c2w$rho_gc_ref), sum(abs(c2w$rho_gc_ref) > 0.4), med(c7$n_mr_ref_fdr10), NA, NA, med(c9$edges_ref), NA, NA, NA),
  test = c(med(c2w$rho_gc_test), sum(abs(c2w$rho_gc_test) > 0.4), med(c7$n_mr_test_fdr10), med(c7$jaccard_fdr10), med(c7$rho_nes),
         med(c9$edges_test), med(c9$frac_removed), sum(c10$regulators_lowmap + c10$targets_lowmap), sum(c10$cross_mappable_edges)))
fwrite(summ, file.path(OUT, "summary.tsv"), sep = "\t", na = "NA")
print(summ)
message("STEP 6 V2 CHECKS DONE")
