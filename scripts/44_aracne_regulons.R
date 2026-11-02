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
# VERSION=v2: inputs/networks/output with suffix _v2; REGSET=tf (v2 only) uses the TF-only networks (aracne_v2_tf ->
#   regulons_aracne_project_v2_tf.rds). v2 drops edges between cross-mappable gene pairs (symmetric cross-mappability
#   > 0, Saha & Battle 2018; data/raw/crossmap/) BEFORE aracne2regulon (D46/D51, check C10), as the edge likelihood is
#   scaled by the network's largest MI; counts per network in results/step5/crossmap_edges_removed<sfx>.tsv.
# OUT_TAG=<tag>: write regulons_aracne_project<sfx><tag>.rds instead, e.g. OUT_TAG=_brca_controls for the BRCA_LUM
#   size-matched controls (D81), which stay out of the production regulon file.
suppressPackageStartupMessages({library(data.table); library(viper)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
SFX <- if (identical(Sys.getenv("VERSION", "v1"), "v2")) "_v2" else ""
REGSET <- Sys.getenv("REGSET", "full"); stopifnot(REGSET %in% c("full", "tf"), REGSET == "full" || SFX == "_v2")
TFX <- if (REGSET == "tf") "_tf" else ""
IN <- paste0("data/processed/step5/aracne_inputs", SFX); NETDIR <- paste0("data/processed/step5/aracne", SFX, TFX)
OUT <- paste0("data/processed/step5/regulons_aracne_project", SFX, TFX, Sys.getenv("OUT_TAG"), ".rds")
CROSSMAP <- SFX == "_v2"
if (CROSSMAP) {
  cm <- fread(cmd = "zcat data/raw/crossmap/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz", header = FALSE,
              col.names = c("a", "b", "s"))[s > 0]
  cm[, `:=`(a = sub("\\..*", "", a), b = sub("\\..*", "", b), s = NULL)]
  cm <- unique(rbind(cm, cm[, .(a = b, b = a)])); setkey(cm, a, b)
  rm_log <- list()
}

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
  ed <- fread(f)[, 1:3]; setnames(ed, c("a", "b", "mi"))
  if (CROSSMAP) {
    cross <- !is.na(cm[ed, on = .(a, b), mult = "first", which = TRUE])
    rm_log[[net]] <- data.table(network = net, edges = nrow(ed), cross_mappable = sum(cross), frac = round(mean(cross), 4))
    message(sprintf("%s: %d of %d edges cross-mappable (%.1f%%), removed", net, sum(cross), nrow(ed), 100 * mean(cross)))
    ed <- ed[!cross]
  }
  fwrite(ed, tmp, sep = "\t", col.names = FALSE)
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
if (CROSSMAP) {
  lf <- paste0("results/step5/crossmap_edges_removed", SFX, TFX, Sys.getenv("OUT_TAG"), ".tsv")
  prev <- if (file.exists(lf)) fread(lf)[!network %in% names(rm_log)] else NULL
  fwrite(rbind(prev, rbindlist(rm_log))[order(network)], lf, sep = "\t")
}
ck <- "results/step5/regulon_checksums.sha256"
lines <- grep(OUT, readLines(ck), fixed = TRUE, invert = TRUE, value = TRUE)
writeLines(c(lines, system2("sha256sum", OUT, stdout = TRUE)), ck)
message("ARACNE REGULONS DONE: ", paste(names(old), collapse = ", "))
