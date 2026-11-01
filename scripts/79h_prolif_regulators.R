#!/usr/bin/env Rscript
# Proliferation-linked regulators of the project ARACNe networks, for the split reversal score (scripts/79i).
# Why (user 2026-10-12, caveat C3 / C4 of docs/notes/drug_side_decisions_and_caveats.md): every tumor has proliferation
# regulators up vs normal and every cytostatic drug turns them down, so proliferation reversal ranks the same drugs
# high in every disease; in tumor-derived networks lineage regulators (FOXA1 / GATA3 in breast) also carry
# proliferation genes as (negative) targets, so growth arrest is read as their activation. Both kinds are marked here.
# Rule (fixed before looking at drug results): proliferation genes = HALLMARK_E2F_TARGETS + HALLMARK_G2M_CHECKPOINT
# (MSigDB 2025.1, symbols -> GENCODE v23 IDs via results/step5/id_map.tsv). Per network x regulator: one-sided
# hypergeometric test of proliferation genes among its targets (either sign; universe = all targets of the network),
# BH over the network's regulators. prolif_linked = FDR < FDR_MAX and fold enrichment >= MIN_FOLD.
# Output: results/lincs/prolif_regulators.tsv (net, regulator (symbol), n_targets, n_prolif, n_prolif_pos,
#         n_prolif_neg, fold, p, fdr, prolif_linked)
# VERSION=v2 [REGSET=tf]: v2 regulons (cross-mappable edges removed) -> prolif_regulators_v2[_tf].tsv
# REGULON=collectri (D75 consensus, 2026-10-28): the CollecTRI regulon as one network "ALL" (universe = all its targets)
#   -> prolif_regulators_collectri.tsv
suppressPackageStartupMessages(library(data.table))
FDR_MAX <- 0.05; MIN_FOLD <- 2
# REGSET=tf_base | tf_sig1 | tf_sig2 (D83): the comparison target lists -> prolif_regulators_v2_<REGSET>.tsv
RS <- Sys.getenv("REGSET", "full"); stopifnot(RS %in% c("full", "tf", "tf_base", "tf_sig1", "tf_sig2"))
SFX <- paste0(if (identical(Sys.getenv("VERSION", "v1"), "v2")) "_v2" else "", if (RS != "full") paste0("_", RS) else "")
stopifnot(SFX != "_tf")
REGULON <- Sys.getenv("REGULON", "aracne_project"); stopifnot(REGULON %in% c("aracne_project", "collectri"))
if (REGULON == "collectri") SFX <- "_collectri"
arac <- if (REGULON == "collectri") list(ALL = readRDS("results/step5/regulons_core.rds")$collectri) else
  readRDS(paste0("data/processed/step5/regulons_aracne_project", SFX, ".rds"))
idm <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
gmt <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
names(gmt) <- vapply(gmt, `[`, "", 1)
psym <- unique(unlist(lapply(gmt[c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")], `[`, -(1:2))))
pgen <- unique(idm[symbol %in% psym, gene_id])
message(sprintf("proliferation genes: %d symbols, %d mapped", length(psym), length(pgen)))

out <- rbindlist(lapply(names(arac), function(net) {
  r <- arac[[net]]
  univ <- unique(unlist(lapply(r, function(x) names(x$tfmode))))
  K <- sum(univ %in% pgen); N <- length(univ)
  d <- rbindlist(lapply(names(r), function(g) {
    tm <- r[[g]]$tfmode; ip <- names(tm) %in% pgen
    data.table(regulator = g, n_targets = length(tm), n_prolif = sum(ip),
               n_prolif_pos = sum(ip & tm > 0.5), n_prolif_neg = sum(ip & tm <= 0.5))
  }))
  d[, `:=`(net = net, fold = round((n_prolif / n_targets) / (K / N), 3),
           p = phyper(n_prolif - 1, K, N - K, n_targets, lower.tail = FALSE))]
  d[, fdr := p.adjust(p, "BH")]
  d
}))
out[, `:=`(prolif_linked = fdr < FDR_MAX & fold >= MIN_FOLD)]
setcolorder(out, c("net", "regulator"))
fwrite(out, paste0("results/lincs/prolif_regulators", SFX, ".tsv"), sep = "\t")
print(out[, .(regulators = .N, prolif_linked = sum(prolif_linked)), by = net])
