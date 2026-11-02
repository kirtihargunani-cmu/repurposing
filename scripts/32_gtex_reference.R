#!/usr/bin/env Rscript
# Step 4a (v2 rebuild, SPEC_rebuild_v2.md; D57, D59): prepare the GTEx healthy reference.
# One row per GTEx sample used by any group (primary_per_patient GTEx rows of results/step3/covariates_toil.tsv:
# one sample per donor per tissue).
#  1. Drop samples absent from the Toil expected-count matrix or without ischemic time / death circumstance.
#  2. Sub-region filter (D59): Colon - Transverse keeps samples with mean log-CPM of 8 epithelial markers >= 7
#     (others = muscle wall); Kidney - Cortex keeps samples with mean log-CPM of 8 proximal-tubule markers >= 5
#     (others = medulla-dominant or near-empty). Log-CPM: TMM over the tissue's GTEx samples, prior count 2,
#     protein-coding genes.
#  3. Best-preserved donors per tissue (after the filter): ventilator death (Hardy 0), ischemic time <= the
#     tissue's 33rd percentile, autolysis score <= 1 where scored.
#  4. Anchor class per tissue: >= 45 best-preserved -> "observed"; 10-44 -> "model"; < 10 -> "average" (flagged).
#  5. Collection covariates (GTEx only; Step 4 sets them to 0 for tumors), centered at the tissue's anchor:
#       coll_isch   (log ischemic min - anchor) / SD of log ischemic min in the tissue; anchor = median over the
#                   best-preserved donors ("observed", "model") or the tissue mean ("average")
#       coll_fast, coll_slow   death-circumstance indicators minus the anchor value (0 for "observed" / "model",
#                   whose best-preserved donors are all ventilator cases; tissue mean for "average")
#       coll_nonbp  1 for donors that are not best-preserved, in "observed" tissues only; else 0
# Outputs: results/step4_v2/gtex_reference.tsv (sample-level; `use` = TRUE for samples kept),
#          results/step4_v2/gtex_reference_tissues.tsv (per tissue: counts, anchor class, anchor values)
suppressPackageStartupMessages({library(data.table); library(edgeR)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
OUT <- "results/step4_v2"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_OBSERVED <- 45; N_MODEL <- 10; ISCH_Q <- 1 / 3; MAX_AUTOL <- 1
SUBREGION <- list(
  "Colon - Transverse" = list(name = "epithelium", min = 7,
                              genes = c("CDX1", "CDX2", "VIL1", "KRT20", "CDH17", "LGALS4", "CEACAM5", "EPCAM")),
  "Kidney - Cortex" = list(name = "proximal_tubule", min = 5,
                           genes = c("SLC34A1", "LRP2", "CUBN", "SLC22A6", "SLC22A8", "SLC5A12", "ALDOB", "GATM")))

cov <- fread("results/step3/covariates_toil.tsv")[role == "gtex_normal" & primary_per_patient == TRUE]
g <- cov[, .(groups = paste(sort(unique(group_id)), collapse = ";")), by = .(sample, donor = patient, tissue = sub("^GTEx ", "", sample_type),
                                                                            ischemic_min = gtex_ischemic_min, death = gtex_death,
                                                                            autolysis = gtex_autolysis)]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
g[, excluded := fcase(!sample %in% colnames(counts), "not in Toil expression data",
                      is.na(ischemic_min) | is.na(death) | death == "", "missing collection data", default = NA_character_)]

# 2. sub-region filter
pc <- fread("results/step3/genes_protein_coding_v23.tsv"); s2e <- setNames(sub("\\..*", "", pc$gene_id), pc$symbol)
g[, `:=`(subregion = NA_character_, subregion_score = NA_real_)]
for (t in names(SUBREGION)) {
  R <- SUBREGION[[t]]; s <- g[tissue == t & is.na(excluded), sample]
  y <- calcNormFactors(DGEList(counts[, s])); e <- cpm(y, log = TRUE, prior.count = 2)
  mk <- intersect(s2e[R$genes], rownames(e)); stopifnot(length(mk) == length(R$genes))
  sc <- colMeans(e[mk, , drop = FALSE])
  g[match(s, sample), `:=`(subregion = R$name, subregion_score = round(sc, 3))]
  g[sample %in% s & subregion_score < R$min, excluded := sprintf("sub-region: %s score < %g", R$name, R$min)]
}
g[, use := is.na(excluded)]

# 3-5. best-preserved donors, anchor class, collection covariates (per tissue, kept samples only)
g[, `:=`(best_preserved = FALSE, anchor_class = NA_character_, coll_isch = NA_real_, coll_fast = NA_real_,
         coll_slow = NA_real_, coll_nonbp = NA_real_)]
tis <- rbindlist(lapply(sort(unique(g$tissue)), function(t) {
  k <- g$tissue == t & g$use; d <- g[k]
  li <- log(d$ischemic_min); q <- quantile(d$ischemic_min, ISCH_Q)
  bp <- d$death == "ventilator" & d$ischemic_min <= q & (is.na(d$autolysis) | d$autolysis <= MAX_AUTOL)
  n_bp <- sum(bp); cls <- if (n_bp >= N_OBSERVED) "observed" else if (n_bp >= N_MODEL) "model" else "average"
  fast <- as.numeric(d$death == "fast"); slow <- as.numeric(d$death == "slow")
  a_isch <- if (cls == "average") mean(li) else median(li[bp])
  a_fast <- if (cls == "average") mean(fast) else 0; a_slow <- if (cls == "average") mean(slow) else 0
  sdl <- sd(li); if (!is.finite(sdl) || sdl == 0) sdl <- 1
  g[k, `:=`(best_preserved = bp, anchor_class = cls, coll_isch = round((li - a_isch) / sdl, 5),
            coll_fast = round(fast - a_fast, 5), coll_slow = round(slow - a_slow, 5),
            coll_nonbp = if (cls == "observed") as.numeric(!bp) else 0)]
  data.table(tissue = t, n_samples = sum(g$tissue == t), n_excluded = sum(g$tissue == t & !g$use), n_used = nrow(d),
             ischemic_q33_min = round(q, 1), n_ventilator = sum(d$death == "ventilator"), n_best_preserved = n_bp,
             anchor_class = cls, anchor_ischemic_min = round(exp(a_isch), 1), anchor_fast = round(a_fast, 3),
             anchor_slow = round(a_slow, 3), groups = paste(sort(unique(unlist(strsplit(d$groups, ";")))), collapse = ";"))
}))
setorder(g, tissue, sample)
fwrite(g, file.path(OUT, "gtex_reference.tsv"), sep = "\t", na = "NA")
fwrite(tis, file.path(OUT, "gtex_reference_tissues.tsv"), sep = "\t")
print(tis[, .(tissue, n_samples, n_excluded, n_used, n_best_preserved, anchor_class, anchor_ischemic_min)])
message(sprintf("GTEX REFERENCE DONE: %d samples, %d used; anchor classes: %s", nrow(g), sum(g$use),
                paste(tis[, .N, by = anchor_class][, paste0(anchor_class, " ", N)], collapse = ", ")))
