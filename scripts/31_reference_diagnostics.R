#!/usr/bin/env Rscript
# Step 4 pre-rebuild diagnostics: how much do lab/source, post-mortem and RNA-quality differences between
# TCGA tumors and GTEx normals put into the tumor-vs-GTEx signatures? Measures only; no pipeline output changes.
#
# A. Source effect (healthy vs healthy): per group with >= MIN_NORMAL TCGA adjacent normals and GTEx normals,
#    limma-voom adjacent vs GTEx ~ source + age + sex + ESTIMATE. Both arms are normal tissue of the same organ,
#    so the difference = lab/protocol + post-mortem/surgical handling + field effect + sampling differences.
#    Reported: genes changed, GC / length / mappability tilt, agreement across tissues (a shared component
#    would be lab/protocol and transferable to groups without adjacent normals), overlap with the Step 4
#    tumor-vs-GTEx signature.
# B. Post-mortem within GTEx: per group, the group's GTEx donors only, ~ ischemic time + ventilator death
#    (Hardy scale 0) + age + sex. Reported per term: genes changed, tilt, cross-tissue agreement, and
#    alignment with the source effect and the tumor signature.
# C. RNA-quality index computed jointly on tumors + GTEx + adjacent normals (the Step 6e index was
#    computed within role, so tumor-vs-normal differences were never measured): long-minus-short gene
#    index on the joint TMM logCPM; differences by role; within GTEx vs RIN and ischemic time.
# D. Mappability (Saha & Battle 2018, gene mappability < 0.8) and readthrough loci among Step 4 genes.
# Outputs: results/step4/reference_diagnostics/*.tsv; gene x group z matrices z_source / z_ischemic / z_ventilator.
suppressPackageStartupMessages({library(data.table); library(limma); library(edgeR); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(31)
MIN_NORMAL <- 10; N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
out <- "results/step4/reference_diagnostics"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient == TRUE]
groups <- fread("results/step1/disease_groups_toil.tsv")[status == "include"]
summ4 <- fread("results/step4/signature_summary.tsv")[status == "ok"]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
cov <- cov[sample %in% colnames(counts) & group_id %in% groups$group_id]
cov[sample_type_class == "new_primary", sample_type_class := "primary"]
# GTEx donor death circumstance (Hardy scale; 0 = ventilator case)
subj <- fread("data/raw/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
cov[role == "gtex_normal", donor := sub("^(GTEX-[^-]+).*", "\\1", sample)]
cov[subj, on = c(donor = "SUBJID"), hardy := i.DTHHRDY]

# gene annotations: GC, transcript length, genomic span, mappability, readthrough
gcl <- fread("data/interim/gene_gc_v23.tsv")[, gene_id := sub("\\..*", "", gene_id)]
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))
gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
mapv <- setNames(gm$map, sub("\\..*", "", gm$id))
rt <- fread("data/interim/readthrough_genes_v23.tsv")
ann <- data.table(gene_id = rownames(counts))
ann[, `:=`(gc = gcl$gc_longest_tx[match(gene_id, gcl$gene_id)], txlen = log10(gcl$longest_tx_len[match(gene_id, gcl$gene_id)]),
           span = log10(span[gene_id]), map = mapv[gene_id], readthrough = gene_id %in% rt[drop == TRUE, gene_id])]

lcpm_of <- function(s) { y <- calcNormFactors(DGEList(counts[, s])); cpm(y, log = TRUE, prior.count = 2) }
qindex <- function(m) {  # same definition as scripts/61_quality_correct.R, on whatever samples are passed
  m <- m[apply(m, 1, median) >= 1 & rownames(m) %in% names(span), , drop = FALSE]
  z <- (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  colMeans(z[span[rownames(z)] >= q[2], , drop = FALSE]) - colMeans(z[span[rownames(z)] <= q[1], , drop = FALSE])
}
fit_z <- function(s, d, terms, coefs) {  # voom fit on samples s with design data d; z per coefficient
  design <- model.matrix(reformulate(terms), data = d)
  design <- design[, !colnames(design) %in% nonEstimable(design), drop = FALSE]
  y <- DGEList(counts[, s]); keep <- filterByExpr(y, design)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE]); fit <- eBayes(lmFit(voom(y, design), design))
  lapply(setNames(intersect(coefs, colnames(design)), intersect(coefs, colnames(design))), function(k) {
    tt <- topTable(fit, coef = k, number = Inf, sort.by = "none")
    list(z = setNames(zscoreT(tt$t, df = fit$df.total), rownames(tt)),
         n_fdr05_lfc1 = sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1), n_fdr05 = sum(tt$adj.P.Val < 0.05), n_genes = nrow(tt))
  })
}
tilt <- function(z) { a <- ann[match(names(z), gene_id)]
  sapply(c("gc", "txlen", "span", "map"), function(v) round(cor(z, a[[v]], method = "spearman", use = "complete.obs"), 3)) }
sex_term <- function(d) if (uniqueN(d$sex) > 1 && min(table(d$sex)) >= 3) "sex" else NULL

per_group <- function(gid) {
  d <- cov[group_id == gid]; res <- list(group_id = gid)
  # C. joint quality index (tumor + GTEx + adjacent), one TMM over all of the group's samples
  qi <- qindex(lcpm_of(d$sample)); d[, qi := qi[sample]]
  res$quality <- d[, .(group_id = gid, sample, role, source, qi, gtex_rin, gtex_ischemic_min, hardy)]
  # A. source effect: adjacent vs GTEx (both healthy)
  a <- d[role %in% c("adjacent_normal", "gtex_normal") & !is.na(age) & !is.na(estimate_score) & !is.na(sex)]
  if (a[role == "adjacent_normal", .N] >= MIN_NORMAL && a[role == "gtex_normal", .N] >= MIN_NORMAL) {
    a[, `:=`(src = factor(ifelse(role == "adjacent_normal", "tcga", "gtex"), c("gtex", "tcga")),
             age_z = as.numeric(scale(age)), estimate_z = as.numeric(scale(estimate_score)))]
    res$source <- fit_z(a$sample, a, c("src", "age_z", "estimate_z", sex_term(a)), "srctcga")[[1]]
    res$source$n <- a[, .N, by = role][, setNames(N, role)]
  }
  # A2. split-half control: tumor-vs-GTEx and adjacent-vs-GTEx share the GTEx arm, so their agreement is partly
  #     shared sampling noise. Split GTEx donors in half: source effect vs half 1, tumor signature vs half 2.
  if (!is.null(res$source)) {
    gs <- a[role == "gtex_normal", sample]; h1 <- sample(gs, floor(length(gs) / 2)); h2 <- setdiff(gs, h1)
    t <- d[role == "tumor" & !is.na(age) & !is.na(estimate_score) & !is.na(sex) & sample_type_class == "primary", sample]
    mk <- function(s, lvl) { x <- d[match(s, sample)]
      x[, `:=`(cond = factor(ifelse(role %in% lvl, "b", "a"), c("a", "b")), age_z = as.numeric(scale(age)), estimate_z = as.numeric(scale(estimate_score)))]; x }
    x1 <- mk(c(a[role == "adjacent_normal", sample], h1), "adjacent_normal"); x2 <- mk(c(t, h2), "tumor")
    s1 <- fit_z(x1$sample, x1, c("cond", "age_z", "estimate_z", sex_term(x1)), "condb")[[1]]$z
    s2 <- fit_z(x2$sample, x2, c("cond", "age_z", "estimate_z", sex_term(x2)), "condb")[[1]]$z
    k <- intersect(names(s1), names(s2)); res$split_rho <- round(cor(s1[k], s2[k], method = "spearman"), 3)
  }
  # B. post-mortem within GTEx donors
  g <- d[role == "gtex_normal" & !is.na(age) & !is.na(sex) & !is.na(gtex_ischemic_min) & !is.na(hardy)]
  if (nrow(g) >= 20) {
    g[, `:=`(isch_z = as.numeric(scale(gtex_ischemic_min)), vent = as.integer(hardy == 0), age_z = as.numeric(scale(age)))]
    terms <- c("isch_z", if (g[, sum(vent) >= 5 && sum(!vent) >= 5]) "vent", "age_z", sex_term(g))
    res$pm <- fit_z(g$sample, g, terms, c("isch_z", "vent")); res$pm_n <- c(n = nrow(g), n_vent = g[, sum(vent)])
  }
  res
}
gids <- summ4[contrast == "gtex", unique(group_id)]
res <- mclapply(gids, function(g) tryCatch(per_group(g), error = function(e) list(group_id = g, error = conditionMessage(e))),
                mc.cores = N_CORES)
errs <- Filter(function(r) !is.null(r$error), res); if (length(errs)) print(rbindlist(lapply(errs, as.data.table)))
res <- Filter(function(r) is.null(r$error), res); names(res) <- sapply(res, `[[`, "group_id")

zt <- fread("results/step4/z_tumor_vs_gtex.tsv"); ztv <- function(g) setNames(zt[[g]], zt$gene_id)[!is.na(zt[[g]])]
za <- fread("results/step4/z_tumor_vs_adjacent.tsv")
zmat <- function(lst) { genes <- sort(unique(unlist(lapply(lst, function(z) names(z)))))
  m <- data.table(gene_id = genes); for (g in names(lst)) m[, (g) := round(lst[[g]][genes], 3)]; m }
rho <- function(x, y) { k <- intersect(names(x), names(y)); round(cor(x[k], y[k], method = "spearman"), 3) }
cross <- function(lst) { n <- names(lst); if (length(n) < 2) return(NA)
  round(median(unlist(lapply(combn(n, 2, simplify = FALSE), function(p) rho(lst[[p[1]]], lst[[p[2]]])))), 3) }
# leave-one-TISSUE-out: drop every group that shares a GTEx donor sample with g (breast subtypes, colon/rectum, kidney, ...)
gtex_of <- cov[role == "gtex_normal", .(s = list(sample)), by = group_id]; gtex_of <- setNames(gtex_of$s, gtex_of$group_id)
shares <- function(g, h) length(intersect(gtex_of[[g]], gtex_of[[h]])) > 0
loo_consensus <- function(lst, g) { o <- Filter(function(h) h != g && !shares(g, h), names(lst)); zm <- zmat(lst[o])
  m <- as.matrix(zm[, -1]); rownames(m) <- zm$gene_id; apply(m, 1, median, na.rm = TRUE) }

# A summary
src <- Filter(Negate(is.null), lapply(res, function(r) r$source)); zsrc <- lapply(src, `[[`, "z")
fwrite(zmat(zsrc), file.path(out, "z_source.tsv"), sep = "\t", na = "NA")
src_cons <- { m <- zmat(zsrc); setNames(apply(as.matrix(m[, -1]), 1, median, na.rm = TRUE), m$gene_id) }
sA <- rbindlist(lapply(names(src), function(g) {
  s <- src[[g]]; tg <- ztv(g); zad <- setNames(za[[g]], za$gene_id)[!is.na(za[[g]])]
  # how much of tumor-vs-GTEx is explained by the source effect (linear, shared genes)
  k <- intersect(names(tg), names(s$z)); r2 <- summary(lm(tg[k] ~ s$z[k]))$r.squared
  data.table(group_id = g, n_adjacent = s$n[["adjacent_normal"]], n_gtex = s$n[["gtex_normal"]], genes = s$n_genes,
             src_fdr05 = s$n_fdr05, src_fdr05_lfc1 = s$n_fdr05_lfc1, t(tilt(s$z)),
             rho_src_vs_tumor_gtex = rho(s$z, tg), r2_tumor_gtex_on_src = round(r2, 3),
             rho_src_vs_tumor_adj = if (length(zad)) rho(s$z, zad) else NA,
             rho_src_vs_loo_consensus = rho(s$z, loo_consensus(zsrc, g)), split_half_rho_src_vs_tumor = res[[g]]$split_rho,
             gtex_minus_adj_vs_src = if (length(zad)) rho(tg[intersect(names(tg), names(zad))] - zad[intersect(names(tg), names(zad))], s$z) else NA)
}))
fwrite(sA, file.path(out, "source_effect_summary.tsv"), sep = "\t")
# shared source component vs every group's tumor-vs-GTEx signature (groups with adjacent normals use leave-one-out)
sAll <- rbindlist(lapply(gids, function(g) {
  cons <- if (g %in% names(zsrc)) loo_consensus(zsrc, g) else src_cons
  data.table(group_id = g, has_adjacent = g %in% names(zsrc), rho_tumor_gtex_vs_src_consensus = rho(ztv(g), cons))
}))
fwrite(sAll, file.path(out, "source_consensus_vs_signatures.tsv"), sep = "\t")
# what the shared source component is: top genes and Hallmark enrichment (rank-based, limma::cameraPR)
sym <- fread("results/step3/genes_protein_coding_v23.tsv"); symv <- setNames(sym$symbol, sym$gene_id)
cz <- src_cons[!is.na(src_cons)]; names(cz) <- symv[names(cz)]; cz <- cz[!is.na(names(cz)) & !duplicated(names(cz))]
hm <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
hm <- setNames(lapply(hm, function(x) x[-(1:2)]), sapply(hm, `[`, 1))
cam <- as.data.table(cameraPR(cz, ids2indices(hm, names(cz))), keep.rownames = "hallmark")
fwrite(cam, file.path(out, "source_consensus_hallmarks.tsv"), sep = "\t")
o <- order(cz); fwrite(data.table(rank = 1:40, up_in_tcga_normals = names(cz)[rev(o)][1:40], up_in_gtex = names(cz)[o][1:40]),
                     file.path(out, "source_consensus_top_genes.tsv"), sep = "\t")

# B summary
pmz <- function(term) { l <- lapply(res, function(r) r$pm[[term]]$z); Filter(Negate(is.null), l) }
zis <- pmz("isch_z"); zve <- pmz("vent")
fwrite(zmat(zis), file.path(out, "z_ischemic.tsv"), sep = "\t", na = "NA"); fwrite(zmat(zve), file.path(out, "z_ventilator.tsv"), sep = "\t", na = "NA")
sB <- rbindlist(lapply(names(zis), function(g) {
  r <- res[[g]]; tg <- ztv(g)
  ti <- tilt(r$pm$isch_z$z); names(ti) <- paste0("isch_", names(ti))
  data.table(group_id = g, n_gtex = r$pm_n[["n"]], n_vent = r$pm_n[["n_vent"]],
             isch_fdr05 = r$pm$isch_z$n_fdr05, isch_fdr05_lfc1 = r$pm$isch_z$n_fdr05_lfc1, t(ti),
             rho_isch_vs_tumor_gtex = rho(r$pm$isch_z$z, tg),
             rho_isch_vs_src = if (g %in% names(zsrc)) rho(r$pm$isch_z$z, zsrc[[g]]) else NA,
             vent_fdr05 = if (!is.null(r$pm$vent)) r$pm$vent$n_fdr05 else NA,
             rho_vent_vs_tumor_gtex = if (!is.null(r$pm$vent)) rho(r$pm$vent$z, tg) else NA,
             rho_vent_vs_src = if (!is.null(r$pm$vent) && g %in% names(zsrc)) rho(r$pm$vent$z, zsrc[[g]]) else NA)
}))
fwrite(sB, file.path(out, "postmortem_summary.tsv"), sep = "\t")

# C summary
q <- rbindlist(lapply(res, `[[`, "quality")); fwrite(q, file.path(out, "quality_index_joint.tsv"), sep = "\t")
dco <- function(x, y) round((mean(x) - mean(y)) / sqrt((var(x) + var(y)) / 2), 2)
sC <- q[, {
  t <- qi[role == "tumor"]; n <- qi[role == "gtex_normal"]; a <- qi[role == "adjacent_normal"]; gg <- .SD[role == "gtex_normal"]
  .(n_tumor = length(t), n_gtex = length(n), n_adj = length(a),
    mean_tumor = round(mean(t), 3), mean_gtex = round(mean(n), 3), mean_adj = if (length(a)) round(mean(a), 3) else NA_real_,
    d_tumor_vs_gtex = dco(t, n), d_adj_vs_gtex = if (length(a) >= MIN_NORMAL) dco(a, n) else NA_real_,
    d_tumor_vs_adj = if (length(a) >= MIN_NORMAL) dco(t, a) else NA_real_,
    rho_gtex_qi_rin = round(gg[, cor(qi, gtex_rin, method = "spearman", use = "complete.obs")], 3),
    rho_gtex_qi_ischemic = round(gg[, cor(qi, gtex_ischemic_min, method = "spearman", use = "complete.obs")], 3))
}, by = group_id]
fwrite(sC, file.path(out, "quality_by_role.tsv"), sep = "\t")

# D. mappability / readthrough among genes in any Step 4 signature
g4 <- sub("\\..*", "", zt$gene_id); a4 <- ann[gene_id %in% g4]
sD <- data.table(step4_genes = nrow(a4), no_mappability_value = a4[is.na(map), .N], map_lt_0.8 = a4[map < 0.8, .N],
                 readthrough_drop = a4[readthrough == TRUE, .N], either = a4[map < 0.8 | readthrough == TRUE, .N])
fwrite(sD, file.path(out, "mappability_genes.tsv"), sep = "\t")

overview <- data.table(
  item = c("source effect: groups", "source: median genes FDR<0.05 & |lfc|>1", "source: median GC tilt", "source: median length (span) tilt",
           "source: median cross-tissue agreement (rho)", "source: median rho with tumor-vs-GTEx", "source: median R2 of tumor-vs-GTEx on source",
           "shared source component vs tumor-vs-GTEx (all groups): median rho",
           "ischemic time: median genes FDR<0.05", "ischemic: median cross-tissue agreement", "ischemic: median rho with tumor-vs-GTEx",
           "ischemic: median rho with source effect", "ventilator death: median rho with tumor-vs-GTEx",
           "quality index: median d tumor vs GTEx", "quality index: median d adjacent vs GTEx", "quality index: median rho with GTEx RIN",
           "genes mappability < 0.8 or readthrough"),
  value = c(nrow(sA), median(sA$src_fdr05_lfc1), median(sA$gc), median(sA$span), cross(zsrc), median(sA$rho_src_vs_tumor_gtex),
            median(sA$r2_tumor_gtex_on_src), median(sAll$rho_tumor_gtex_vs_src_consensus),
            median(sB$isch_fdr05), cross(zis), median(sB$rho_isch_vs_tumor_gtex), median(sB$rho_isch_vs_src, na.rm = TRUE),
            median(sB$rho_vent_vs_tumor_gtex, na.rm = TRUE),
            median(sC$d_tumor_vs_gtex), median(sC$d_adj_vs_gtex, na.rm = TRUE), median(sC$rho_gtex_qi_rin, na.rm = TRUE), sD$either))
fwrite(overview, file.path(out, "overview.tsv"), sep = "\t"); print(overview)
message("REFERENCE DIAGNOSTICS DONE")
