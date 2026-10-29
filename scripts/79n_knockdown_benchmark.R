#!/usr/bin/env Rscript
# Which regulon set reads a TF's inhibition correctly on perturbation data: LINCS TF knockdowns (scripts/79m) scored with
# the v2 TF ARACNe networks, CollecTRI and DoRothEA ABC; per TF, the set the drug side uses (user 2026-10-28).
# Why: the tumor-network ESR1 regulon encodes luminal identity across patients and reads anti-estrogens in MCF7 as ESR1
# activation; CollecTRI reads them as inhibition. Knockdowns test every TF with data the same way.
# Rules (fixed 2026-10-28, before any result):
#   signatures  shRNA consensus (trt_sh.cgs) and CRISPR (trt_xpr, is_hiq) in cancer lines; kept only where the TF is
#               expressed in the line (DepMap 24Q4 log2(TPM + 1) >= 1) and, if the TF's own mRNA is on the profile,
#               it is down (z <= -1.5: the knockdown worked)
#   scoring     viper (method = "none"; profiles are MODZ z) per regulon set, regulons pruned to the profile genes
#               (>= 10 targets, >= 50% present, as 79d); the knocked-down TF is removed from its own targets.
#               Per signature: percentile of the TF's NES among all regulators scored (0 = most inhibited) and its sign.
#               ARACNe: networks of the line's lineage (config/celligner_tissue_keys.tsv) holding the TF, mean NES and
#               percentile over them; no matching network -> all networks holding the TF (aracne_matched = FALSE)
#   pass        per TF x set: >= 3 signatures, median percentile <= 0.25 and >= 60% of signatures with NES < 0
#   choice      ARACNe if it passes; else CollecTRI, else DoRothEA (first that passes); none passes -> "unreliable"
#               (left out of the drug-side reversal score); no usable knockdown -> ARACNe, flagged "untested"
# Outputs: results/lincs/knockdown_benchmark_signatures.tsv.gz (signature x set), knockdown_benchmark_tf.tsv (TF x set),
#          knockdown_regulon_choice.tsv (TF -> set), knockdown_benchmark_summary.tsv
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_TARGETS <- 10; MIN_FRAC <- 0.5; MIN_EXPR <- 1; KD_Z <- -1.5; MIN_SIGS <- 3; MAX_PCT <- 0.25; MIN_NEG <- 0.6
N_CORES <- as.integer(Sys.getenv("N_CORES", "16"))
proc <- "data/processed/lincs"; res <- "results/lincs"

genes <- fread(file.path(proc, "level5_bing_genes.tsv"))
meta <- fread(file.path(proc, "knockdown_bing_sigs.tsv.gz"))
X <- matrix(readBin(file.path(proc, "knockdown_bing.f32"), "numeric", size = 4, n = nrow(genes) * nrow(meta)), nrow = nrow(genes))
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- meta$sig_id
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])

# expression of the TF in the line; knockdown effect on its own mRNA
ex <- fread("data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
setnames(ex, 1, "ModelID"); setnames(ex, sub(" \\(.*", "", names(ex)))
exl <- melt(ex[ModelID %in% meta$ModelID, c("ModelID", intersect(unique(meta$cmap_name), names(ex))), with = FALSE],
            id.vars = "ModelID", variable.name = "cmap_name", value.name = "expr", variable.factor = FALSE)
meta <- merge(meta, exl, by = c("ModelID", "cmap_name"), all.x = TRUE)
meta[, own_z := { e <- sym2ens[cmap_name]; ifelse(is.na(e), NA_real_, X[cbind(match(e, rownames(X)), match(sig_id, colnames(X)))]) }]
meta[, keep := !is.na(expr) & expr >= MIN_EXPR & (is.na(own_z) | own_z <= KD_Z)]
message(sprintf("signatures: %d; TF expressed and knockdown confirmed: %d (%d TFs)", nrow(meta), sum(meta$keep),
                uniqueN(meta[keep == TRUE, cmap_name])))
meta <- meta[keep == TRUE]; X <- X[, meta$sig_id]

prep <- function(rg) {
  r <- lapply(names(rg), function(n) { tm <- rg[[n]]$tfmode; lk <- rg[[n]]$likelihood; e <- sym2ens[n]
    k <- names(tm) %in% rownames(X) & (is.na(e) | names(tm) != e)          # no self-target
    if (sum(k) < MIN_TARGETS || sum(k) / length(tm) < MIN_FRAC) return(NULL)
    list(tfmode = tm[k], likelihood = lk[k]) })
  names(r) <- names(rg); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r
}
score <- function(reg, cols) {  # TF NES and percentile in each signature
  v <- viper(X[, cols, drop = FALSE], reg, method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  m <- meta[match(cols, sig_id)]; i <- match(m$cmap_name, rownames(v)); ok <- !is.na(i)
  pct <- apply(v, 2, function(x) (rank(x) - 1) / (length(x) - 1))
  data.table(sig_id = cols[ok], nes = v[cbind(i[ok], which(ok))], pct = pct[cbind(i[ok], which(ok))], n_reg = nrow(v))
}

core <- readRDS("results/step5/regulons_core.rds")
cur <- rbindlist(lapply(c("collectri", "dorothea_ABC"), function(s) cbind(set = s, score(prep(core[[s]]), meta$sig_id))))

arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
keys <- fread(cmd = "grep -v '^#' config/celligner_tissue_keys.tsv")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
ln <- unique(keys[group_id %in% names(net_of), .(lineage = unlist(strsplit(lineages, ";")), net = net_of[group_id]), by = group_id][, .(lineage, net)])
ln <- unique(ln[net %in% names(arac)])
lin <- fread("data/raw/depmap_24q4/Model.csv")[, .(ModelID, lineage = OncotreeLineage)]
meta <- merge(meta, lin, by = "ModelID", all.x = TRUE)
ar <- rbindlist(mclapply(names(arac), function(n) cbind(net = n, score(prep(arac[[n]]), meta$sig_id)), mc.cores = N_CORES))
ar <- merge(ar, meta[, .(sig_id, lineage)], by = "sig_id")
ar[, matched := paste(lineage, net) %in% ln[, paste(lineage, net)]]
ar <- ar[, if (any(matched)) .(nes = mean(nes[matched]), pct = mean(pct[matched]), n_reg = round(mean(n_reg[matched])), aracne_matched = TRUE, nets = paste(net[matched], collapse = ";"))
           else .(nes = mean(nes), pct = mean(pct), n_reg = round(mean(n_reg)), aracne_matched = FALSE, nets = paste(net, collapse = ";")), by = sig_id]
sg <- rbind(cbind(set = "aracne_v2_tf", ar), cur, fill = TRUE)
sg <- merge(meta[, .(sig_id, technology, cmap_name, cell_iname, lineage, pert_time, expr, own_z)], sg, by = "sig_id")
fwrite(sg, file.path(res, "knockdown_benchmark_signatures.tsv.gz"), sep = "\t")

tf <- sg[, .(n_sigs = .N, n_lines = uniqueN(cell_iname), median_pct = round(median(pct), 3), frac_neg = round(mean(nes < 0), 3),
             median_nes = round(median(nes), 2), frac_aracne_matched = if (set[1] == "aracne_v2_tf") round(mean(aracne_matched), 2) else NA_real_,
             shRNA = sum(technology == "shRNA"), CRISPR = sum(technology == "CRISPR")), by = .(cmap_name, set)]
tf[, pass := n_sigs >= MIN_SIGS & median_pct <= MAX_PCT & frac_neg >= MIN_NEG]
fwrite(tf[order(cmap_name, set)], file.path(res, "knockdown_benchmark_tf.tsv"), sep = "\t")

allreg <- unique(c(unlist(lapply(arac, names)), names(core$collectri), names(core$dorothea_ABC)))
allreg <- allreg[!startsWith(allreg, "ENSG")]
pick <- function(x) { for (s in c("aracne_v2_tf", "collectri", "dorothea_ABC")) if (isTRUE(x[set == s, pass])) return(s)
  if (any(x$n_sigs >= MIN_SIGS)) "unreliable" else "untested" }
ch <- tf[, .(choice = pick(.SD)), by = cmap_name]
ch <- merge(data.table(cmap_name = allreg), ch, by = "cmap_name", all.x = TRUE)[is.na(choice), choice := "untested"]
ch[, regulon_used := fifelse(choice == "untested", "aracne_v2_tf", fifelse(choice == "unreliable", "none", choice))]
fwrite(ch, file.path(res, "knockdown_regulon_choice.tsv"), sep = "\t")

sm <- sg[, .(signatures = .N, tfs = uniqueN(cmap_name), auc = round(1 - mean(pct), 3), frac_neg = round(mean(nes < 0), 3)), by = .(set, technology)]
sm2 <- tf[n_sigs >= MIN_SIGS, .(tfs_tested = .N, tfs_pass = sum(pass), frac_pass = round(mean(pass), 3)), by = set]
fwrite(merge(sm, sm2, by = "set"), file.path(res, "knockdown_benchmark_summary.tsv"), sep = "\t")
print(merge(sm, sm2, by = "set")); print(ch[, .N, by = choice])
print(tf[cmap_name %in% c("ESR1", "AR", "GATA3", "FOXA1", "HNF4A", "MYC", "E2F1", "FOXM1", "TP53", "SOX2"),
         .(cmap_name, set, n_sigs, median_pct, frac_neg, pass)][order(cmap_name, set)], nrows = 60)
message("KNOCKDOWN BENCHMARK DONE")
