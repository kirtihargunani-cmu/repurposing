#!/usr/bin/env Rscript
# Step 5e: project ARACNe networks (scripts/43_run_aracne.sh) -> VIPER regulon objects.
#   viper::aracne2regulon on the consolidated network.txt (3-col: regulator, target, MI) and the same
#   expression matrix the network was built from (cohort tumors, log2(TPM+1)); the mode of regulation (sign)
#   and edge likelihood come from that matrix. Full regulator set, all targets (no pruning) - PI-approved
#   2026-10-04. Regulators renamed Ensembl -> HGNC symbol (as in regulons_core.rds; the
#   few without a current symbol keep their Ensembl ID); targets stay Ensembl.
# Usage: scripts/44_aracne_regulons.R [NET ...]   (no args = every network with a finished network.txt)
# Output: data/processed/step5/regulons_aracne_project.rds (list keyed by network, e.g. BRCA, COADREAD;
#         rebuilt from all finished networks on every run), checksum in results/step5/regulon_checksums.sha256
suppressPackageStartupMessages({library(data.table); library(viper)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
IN <- "data/processed/step5/aracne_inputs"; NETDIR <- "data/processed/step5/aracne"
OUT <- "data/processed/step5/regulons_aracne_project.rds"

nets <- commandArgs(trailingOnly = TRUE)
if (!length(nets)) nets <- basename(dirname(Sys.glob(file.path(NETDIR, "*", "network.txt"))))
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
ens2sym <- setNames(idmap$symbol, idmap$gene_id)

old <- if (file.exists(OUT)) readRDS(OUT) else list()
for (net in nets) {
  f <- file.path(NETDIR, net, "network.txt"); stopifnot(file.exists(f))
  x <- fread(file.path(IN, paste0(net, ".exp.tsv")))
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene
  tmp <- tempfile(fileext = ".txt")  # aracne2regulon's 3col reader has no header handling (header -> MI all NA)
  fwrite(fread(f)[, 1:3], tmp, sep = "\t", col.names = FALSE)
  reg <- aracne2regulon(tmp, m, format = "3col", verbose = FALSE); unlink(tmp)
  sym <- ens2sym[names(reg)]; sym[is.na(sym)] <- names(reg)[is.na(sym)]  # no current HGNC symbol: keep Ensembl ID
  stopifnot(!anyDuplicated(sym))
  names(reg) <- sym
  old[[net]] <- reg
  sz <- vapply(reg, function(r) length(r$tfmode), integer(1))
  message(sprintf("%s: %d regulators, %d edges, targets per regulator median %d (range %d-%d)",
                  net, length(reg), sum(sz), as.integer(median(sz)), min(sz), max(sz)))
}
saveRDS(old, OUT)
ck <- "results/step5/regulon_checksums.sha256"
lines <- grep(OUT, readLines(ck), fixed = TRUE, invert = TRUE, value = TRUE)
writeLines(c(lines, system2("sha256sum", OUT, stdout = TRUE)), ck)
message("ARACNE REGULONS DONE: ", paste(names(old), collapse = ", "))
