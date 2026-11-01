#!/usr/bin/env Rscript
# Step 5 diagnostic (D82): which network edges are subtype identity, and which are signaling?
#
# A pooled network can call a gene a regulator's target only because both genes differ between subtypes (both high in
# luminal breast tumors, say). Such "identity" edges describe what a cell is, not what a drug can change within hours.
# For every edge of the TF-only networks (cross-mappable edges removed: regulons_aracne_project_v2_tf.rds), per label set:
#   histology   our subtype groups (Step 1 group_id): BRCA PAM50, CESC, ESCA, LGG, PCPG, SARC, TGCT, UCEC, COAD vs READ
#   molecular   TCGA PanCanAtlas Subtype_mRNA; where a cohort has none, Subtype_Selected (source column says which)
# Tumors whose label has >= MIN_LABEL tumors; >= 2 labels needed. On the labelled tumors of the network's ARACNe input
# (log2 TPM + 1):
#   rho_p   pooled Spearman, regulator vs target
#   rho_w   within-subtype: Pearson of ranks after removing each subtype's mean rank
#   class   identity: rho_w keeps < 50% of rho_p or flips sign; signaling: keeps >= 50%, same sign; weak: |rho_p| < 0.1
#   subtype_specific (flag; A9): |rho| >= 0.3 in one subtype of >= 50 tumors and < 0.1 or opposite sign in another of >= 50
# Outputs: results/step5/identity_edges_v2_tf/
#   edges_<NET>_<set>.tsv.gz  every edge: regulator, target, mode (ARACNe sign), rho_p, rho_w, retention, class, flag,
#                             per-subtype rho
#   regulon_summary.tsv       per network x set x regulator: targets, shares identity / signaling / weak / subtype-specific
#   network_summary.tsv       per network x set: tumors, labels, edge shares; most identity-driven regulators
#   label_sets.tsv            labels and tumor counts used
# ONLY_NETS=ESCA,CESC: run only these networks (summaries then cover only them).
# Diagnostic only; whether to use signaling-only regulons is a separate decision (D82).
suppressPackageStartupMessages({library(data.table); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_LABEL <- 10; RETAIN <- 0.5; WEAK <- 0.1; SPEC_N <- 50; SPEC_HI <- 0.3; SPEC_LO <- 0.1   # D82
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
out <- "results/step5/identity_edges_v2_tf"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

reg <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
reg <- reg[setdiff(names(reg), "BRCA_LUM")]          # BRCA_LUM is one subtype set; its test is 66 / 67
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
cov <- unique(fread("results/step3/covariates_toil.tsv", select = c("sample", "group_id", "role"))[role == "tumor"], by = "sample")
pc <- fread("data/raw/tcga_pancan/TCGASubtype.20170308.tsv.gz"); setnames(pc, 1, "sid")
HIST <- c("BRCA", "CESC", "ESCA", "LGG", "PCPG", "SARC", "TGCT", "UCEC", "COADREAD")
clean <- function(x) { x <- trimws(as.character(x)); x[x %in% c("", "NA", "-")] <- NA; x }

labels_for <- function(net, samples) {
  sets <- list()
  if (net %in% HIST) sets$histology <- list(source = "group_id", lab = clean(cov$group_id[match(samples, cov$sample)]))
  bc15 <- substr(samples, 1, 15)
  m <- clean(pc$Subtype_mRNA[match(bc15, pc$sid)]); s <- clean(pc$Subtype_Selected[match(bc15, pc$sid)])
  mol <- if (sum(!is.na(m)) > 0) list(source = "Subtype_mRNA", lab = m) else list(source = "Subtype_Selected", lab = s)
  if (!(net == "BRCA")) sets$molecular <- mol   # BRCA: PAM50 is both; kept once as histology
  lapply(sets, function(st) {
    tab <- table(st$lab); keep <- names(tab)[tab >= MIN_LABEL]
    st$lab[!st$lab %in% keep] <- NA; st
  })
}

rank_rows <- function(m) t(apply(m, 1, rank))
cor_rows <- function(a, B) {   # Pearson of a (vector) with each row of B, rows already centred
  as.vector(B %*% a) / (sqrt(sum(a^2)) * sqrt(rowSums(B^2)))
}

run_net <- function(net) {
  x <- fread(file.path("data/processed/step5/aracne_inputs_v2", paste0(net, ".exp.tsv")))
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene
  sets <- labels_for(net, colnames(m))
  rg <- reg[[net]]
  res <- list()
  for (sn in names(sets)) {
    lab <- sets[[sn]]$lab; ok <- !is.na(lab)
    if (length(unique(lab[ok])) < 2) next
    lab <- lab[ok]; mm <- m[, ok, drop = FALSE]
    R <- rank_rows(mm)
    Rp <- R - rowMeans(R)
    Rw <- R; for (l in unique(lab)) { j <- lab == l; Rw[, j] <- R[, j] - rowMeans(R[, j, drop = FALSE]) }
    big <- names(which(table(lab) >= SPEC_N))
    Rk <- lapply(setNames(big, big), function(l) { j <- lab == l; Rl <- t(apply(mm[, j, drop = FALSE], 1, rank)); Rl - rowMeans(Rl) })
    ed <- rbindlist(lapply(names(rg), function(r) {
      ge <- if (startsWith(r, "ENSG")) r else sym2ens[[r]]
      tg <- intersect(names(rg[[r]]$tfmode), rownames(mm))
      if (is.null(ge) || is.na(ge) || !ge %in% rownames(mm) || !length(tg)) return(NULL)
      d <- data.table(regulator = r, target = tg, mode = round(rg[[r]]$tfmode[tg], 3),
                      rho_p = cor_rows(Rp[ge, ], Rp[tg, , drop = FALSE]), rho_w = cor_rows(Rw[ge, ], Rw[tg, , drop = FALSE]))
      for (l in big) d[[paste0("rho_", make.names(l))]] <- round(cor_rows(Rk[[l]][ge, ], Rk[[l]][tg, , drop = FALSE]), 3)
      d
    }))
    ed[, retention := round(rho_w / rho_p, 3)]
    ed[, class := fifelse(abs(rho_p) < WEAK, "weak", fifelse(sign(rho_w) == sign(rho_p) & abs(rho_w) >= RETAIN * abs(rho_p), "signaling", "identity"))]
    kc <- grep("^rho_(?!p$|w$)", names(ed), perl = TRUE, value = TRUE)
    if (length(kc) >= 2) {
      K <- as.matrix(ed[, ..kc])
      hi <- abs(K) >= SPEC_HI; mx <- apply(K, 1, function(v) v[which.max(abs(v))])
      lo <- abs(K) < SPEC_LO | sign(K) != sign(mx)
      ed[, subtype_specific := rowSums(hi) >= 1 & rowSums(lo & !hi) >= 1]
    } else ed[, subtype_specific := NA]
    ed[, `:=`(rho_p = round(rho_p, 3), rho_w = round(rho_w, 3))]
    fwrite(ed, file.path(out, sprintf("edges_%s_%s.tsv.gz", net, sn)), sep = "\t")
    rs <- ed[, .(n_targets = .N, frac_identity = round(mean(class == "identity"), 3), frac_signaling = round(mean(class == "signaling"), 3),
                 frac_weak = round(mean(class == "weak"), 3), frac_subtype_specific = round(mean(subtype_specific), 3)), by = regulator]
    rs[, `:=`(network = net, label_set = sn)]
    tab <- table(lab)
    ns <- data.table(network = net, label_set = sn, source = sets[[sn]]$source, n_tumors = length(lab), n_labels = length(tab),
                     labels = paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse = ", "),
                     subtypes_ge50 = length(big), n_edges = nrow(ed),
                     frac_identity = round(mean(ed$class == "identity"), 3), frac_signaling = round(mean(ed$class == "signaling"), 3),
                     frac_weak = round(mean(ed$class == "weak"), 3), frac_subtype_specific = round(mean(ed$subtype_specific), 3),
                     top_identity_regulators = paste(head(rs[n_targets >= 50][order(-frac_identity), sprintf("%s %.2f", regulator, frac_identity)], 8), collapse = ", "))
    res[[sn]] <- list(rs = rs, ns = ns)
  }
  message(net, ": ", paste(names(res), collapse = ", "))
  res
}

nets <- names(reg); if (nzchar(Sys.getenv("ONLY_NETS"))) nets <- intersect(nets, strsplit(Sys.getenv("ONLY_NETS"), ",")[[1]])
all <- mclapply(nets, function(n) tryCatch(run_net(n), error = function(e) { message(n, " ERROR ", conditionMessage(e)); NULL }),
                mc.cores = N_CORES, mc.preschedule = FALSE)
flat <- unlist(all, recursive = FALSE)
rs <- rbindlist(lapply(flat, `[[`, "rs")); ns <- rbindlist(lapply(flat, `[[`, "ns"))
setcolorder(rs, c("network", "label_set")); fwrite(rs, file.path(out, "regulon_summary.tsv"), sep = "\t")
fwrite(ns[order(label_set, -frac_identity)], file.path(out, "network_summary.tsv"), sep = "\t")
fwrite(ns[, .(network, label_set, source, n_tumors, labels)], file.path(out, "label_sets.tsv"), sep = "\t")
print(ns[order(label_set, -frac_identity), .(network, label_set, source, n_tumors, n_labels, frac_identity, frac_signaling, frac_weak, frac_subtype_specific)])
message("IDENTITY EDGES DONE: ", nrow(ns), " network x label sets")
