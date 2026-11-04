#!/usr/bin/env Rscript
# D91 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): identity-linked regulators per tumor ARACNe v2 TF network.
#   I1 lineage: regulons of the network applied with VIPER (method "none") to DepMap 24Q4 baseline expression of all
#      cancer lines (log2 TPM + 1, per-gene z across lines); eta^2 of activity by OncotreeLineage (lineages >= 10 lines)
#   I2 subtype: networks used by >= 2 disease groups: eta^2 of per-tumor activity (Step 6 patient_nes) by disease group;
#      BRCA also: eta^2 over breast lines by DepMap ModelSubtypeFeatures class (luminal / HER2+ / basal, >= 5 lines each;
#      "luminal" in the label -> luminal, else "HER2+" -> HER2+, else "basal" or "TNBC" -> basal)
#   identity = any eta^2 >= 0.5. Regulators: those scored on the drug side (79d, _d85 run).
# Output: results/lincs/d91/identity_regulators.tsv
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
ETA <- 0.5; MIN_LIN <- 10; MIN_SUB <- 5; MIN_TARGETS <- 10
N_CORES <- min(as.integer(Sys.getenv("N_CORES", "20")), 24L)
res <- "results/lincs/d91"; dir.create(res, FALSE, TRUE)
vd <- Sys.getenv("VD", "data/processed/lincs/drug_viper_specific_v2_tf_d85")   # D93: VD / NETS / EXTRA_RDS / OUT_FILE
nets <- sub("\\.regulators\\.txt$", "", list.files(vd, "\\.regulators\\.txt$"))
if (nzchar(Sys.getenv("NETS"))) nets <- strsplit(Sys.getenv("NETS"), ",")[[1]]
idm <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE & !is.na(entrez_id)]
e2g <- setNames(idm$gene_id, as.character(idm$entrez_id))
ex <- fread("data/raw/depmap_24q4/OmicsExpressionProteinCodingGenesTPMLogp1.csv")
ids <- ex[[1]]; ex <- as.matrix(ex[, -1]); rownames(ex) <- ids
g <- e2g[sub(".*\\((\\d+)\\)$", "\\1", colnames(ex))]; keep <- !is.na(g) & !duplicated(g)
ex <- ex[, keep]; colnames(ex) <- g[keep]
mdl <- fread("data/raw/depmap_24q4/Model.csv")[OncotreePrimaryDisease != "Non-Cancerous"]
ex <- ex[rownames(ex) %in% mdl$ModelID, ]
Z <- t(scale(ex)); Z <- Z[is.finite(rowSums(Z)), ]
mdl <- mdl[match(colnames(Z), ModelID)]
lin <- mdl$OncotreeLineage; okl <- lin %in% names(which(table(lin) >= MIN_LIN))
sub_lab <- function(s) fifelse(grepl("luminal", s), "luminal", fifelse(grepl("HER2\\+", s), "HER2+", fifelse(grepl("basal|TNBC", s), "basal", NA_character_)))
bl <- mdl$OncotreeLineage == "Breast"; bsub <- sub_lab(mdl$ModelSubtypeFeatures); bsub[!bl] <- NA
bsub[bsub %in% names(which(table(bsub) < MIN_SUB))] <- NA
eta2 <- function(v, f) { k <- !is.na(f) & is.finite(v); v <- v[k]; f <- f[k]; if (length(unique(f)) < 2) return(NA_real_)
  m <- ave(v, f); sum((m - mean(v))^2) / sum((v - mean(v))^2) }
arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
if (nzchar(Sys.getenv("EXTRA_RDS"))) { ex <- readRDS(Sys.getenv("EXTRA_RDS")); arac <- c(arac, ex[setdiff(names(ex), names(arac))]) }
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
pdir <- "data/processed/step6/aracne_project_v2_tf"
groups41 <- fread("results/step4_v2/signature_summary.tsv")[status == "ok", unique(group_id)]   # the 41 disease groups
cat(sprintf("DepMap: %d lines (%d in lineages with >= %d lines), %d genes; breast subtype lines: %s\n", ncol(Z), sum(okl), MIN_LIN, nrow(Z),
            paste(names(table(bsub)), table(bsub), collapse = ", ")))
out <- rbindlist(mclapply(nets, function(net) {
  regs <- readLines(file.path(vd, paste0(net, ".regulators.txt")))
  rg <- lapply(arac[[net]][intersect(regs, names(arac[[net]]))], function(x) { k <- names(x$tfmode) %in% rownames(Z)
    if (sum(k) < MIN_TARGETS) return(NULL); list(tfmode = x$tfmode[k], likelihood = x$likelihood[k]) })
  rg <- Filter(Negate(is.null), rg); class(rg) <- "regulon"
  v <- viper(Z, rg, method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  d <- data.table(net = net, regulator = rownames(v),
                  eta2_lineage = apply(v, 1, function(x) eta2(x[okl], lin[okl])),
                  eta2_breast_lines = if (net == "BRCA") apply(v, 1, eta2, f = bsub) else NA_real_)
  gs <- intersect(names(net_of)[net_of == net], groups41)
  if (length(gs) >= 2) {
    P <- lapply(gs, function(gg) { f <- file.path(pdir, paste0(gg, ".patient_nes.tsv.gz")); if (!file.exists(f)) return(NULL)
      x <- fread(f); m <- as.matrix(x[, -1]); rownames(m) <- x$regulator; list(m = m, g = gg) })
    P <- Filter(Negate(is.null), P); common <- Reduce(intersect, lapply(P, function(p) rownames(p$m)))
    M <- do.call(cbind, lapply(P, function(p) p$m[common, , drop = FALSE])); f <- unlist(lapply(P, function(p) rep(p$g, ncol(p$m))))
    et <- setNames(apply(M, 1, eta2, f = f), common)
    d[, eta2_subtype_tumors := et[regulator]]; d[, groups := paste(gs, collapse = ";")]
  } else d[, `:=`(eta2_subtype_tumors = NA_real_, groups = gs[1])]
  d
}, mc.cores = min(N_CORES, length(nets)), mc.preschedule = FALSE), fill = TRUE)
out[, identity := pmax(eta2_lineage, eta2_subtype_tumors, eta2_breast_lines, na.rm = TRUE) >= ETA]
out[is.na(identity), identity := FALSE]
fwrite(out[, .(net, regulator, eta2_lineage = round(eta2_lineage, 3), eta2_subtype_tumors = round(eta2_subtype_tumors, 3),
               eta2_breast_lines = round(eta2_breast_lines, 3), identity, groups)], Sys.getenv("OUT_FILE", file.path(res, "identity_regulators.tsv")), sep = "\t")
options(width = 200); print(out[, .(regulators = .N, identity = sum(identity), by_lineage = sum(eta2_lineage >= ETA, na.rm = TRUE),
              by_subtype = sum(eta2_subtype_tumors >= ETA, na.rm = TRUE), by_breast_lines = sum(eta2_breast_lines >= ETA, na.rm = TRUE)), by = net])
print(out[net == "BRCA" & regulator %in% c("ESR1", "FOXA1", "GATA3", "XBP1", "SPDEF", "RARA", "SOX13", "CREB3L4", "FOXC1", "TFAP2C", "FOXM1", "E2F1", "MYC")])
