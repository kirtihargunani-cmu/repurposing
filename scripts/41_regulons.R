#!/usr/bin/env Rscript
# Step 5b: regulons mapped into the Step 4 (GENCODE v23 Ensembl) ID space.
#
# Resources
#   collectri   OmniPath REST, dataset=collectri; NF-kB / AP-1 complexes collapsed to NFKB / AP1  open
#   dorothea    OmniPath REST, dataset=dorothea, levels A-C                 academic/non-commercial
#   (decoupleR::get_collectri / get_dorothea fail with OmnipathR 4.0 against the current server, so the
#    REST tables are pulled directly and decoupleR's processing is reproduced here.)
#   Sign: mor = -1 if consensus_inhibition and not consensus_stimulation, else +1 (CollecTRI default).
#   aracne_<x>  aracne.networks precomputed TCGA regulons (Entrez)                          open
# OmniPath content changes between releases, so the first run saves raw pulls to
# data/raw/regulons/<resource>_<date>.tsv and later runs reuse the newest saved pull.
#
# Outputs
#   results/step5/regulons_<collectri|dorothea_ABC>.tsv   source, source_gene_id, target, target_gene_id, mor, ...
#   results/step5/regulons_core.rds          viper-style regulon objects (CollecTRI, DoRothEA) on Ensembl IDs
#   data/processed/step5/regulons_aracne_<cohort>.tsv, regulons_aracne.rds   (large; not in git)
#   results/step5/regulon_checksums.sha256   sha256 of both .rds files; Steps 6-16 must load these exact bytes
#   regulon_coverage.tsv         per resource x TF: n_targets, n_targets_mapped, frac_mapped, source mapped,
#                                in Lambert TF list -> the regulon-coverage statistic
suppressPackageStartupMessages({library(data.table); library(aracne.networks)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
raw_dir <- "data/raw/regulons"; dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/step5", showWarnings = FALSE)
dir.create("data/processed/step5", recursive = TRUE, showWarnings = FALSE)
out_dir <- function(nm) if (startsWith(nm, "aracne_")) "data/processed/step5" else "results/step5"

idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
tfs <- fread("results/step5/tf_universe.tsv")
sym2ens <- setNames(idmap$gene_id, idmap$symbol)
ent2ens <- setNames(idmap$gene_id, as.character(idmap$entrez_id))
COMPLEXES <- list(AP1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2"),
                  NFKB = c("NFKB1", "NFKB2", "RELA", "RELB", "REL"))

pull <- function(name, fetch) {
  saved <- sort(Sys.glob(file.path(raw_dir, paste0(name, "_*.tsv"))), decreasing = TRUE)
  if (length(saved)) { message(name, ": using ", saved[1]); return(fread(saved[1])) }
  x <- as.data.table(fetch())
  f <- file.path(raw_dir, sprintf("%s_%s.tsv", name, format(Sys.Date(), "%Y%m%d")))
  fwrite(x, f, sep = "\t"); message(name, ": fetched ", nrow(x), " rows -> ", f)
  x
}

OP <- "https://omnipathdb.org/interactions?genesymbols=yes&organisms=9606&license=academic&format=tsv"
regs <- list()
ct <- pull("collectri", function() fread(paste0(OP, "&datasets=collectri&fields=curation_effort,references")))
complex_name <- function(sym) {
  m <- strsplit(sym, "_", fixed = TRUE)
  vapply(m, function(x) if (all(x %in% COMPLEXES$NFKB)) "NFKB" else if (all(x %in% COMPLEXES$AP1)) "AP1"
         else paste(x, collapse = "_"), "")
}
ct[, source := fifelse(startsWith(source, "COMPLEX:"), complex_name(source_genesymbol), source_genesymbol)]
ct[, mor := fifelse(consensus_inhibition & !consensus_stimulation, -1, 1)]
regs$collectri <- ct[, .(mor = fifelse(sum(mor) < 0, -1, 1), likelihood = 1,
                         n_references = max(curation_effort)), by = .(source, target = target_genesymbol)]
do <- pull("dorothea_ABC", function() fread(paste0(OP, "&datasets=dorothea&dorothea_levels=A,B,C&fields=dorothea_level")))
do[, mor := fifelse(consensus_inhibition & !consensus_stimulation, -1, 1)]
regs$dorothea_ABC <- do[, .(mor = fifelse(sum(mor) < 0, -1, 1), likelihood = 1,
                            confidence = min(dorothea_level)), by = .(source = source_genesymbol, target = target_genesymbol)]

# aracne.networks: one regulon object per TCGA cohort, Entrez-keyed
aracne_objs <- grep("^regulon", data(package = "aracne.networks")$results[, "Item"], value = TRUE)
for (obj in aracne_objs) {
  e <- new.env(); data(list = obj, package = "aracne.networks", envir = e); r <- e[[obj]]
  regs[[paste0("aracne_", sub("^regulon", "", obj))]] <- rbindlist(lapply(names(r), function(tf)
    data.table(source = tf, target = names(r[[tf]]$tfmode), mor = r[[tf]]$tfmode,
               likelihood = r[[tf]]$likelihood)))
}

source_gene <- function(s, entrez) {
  if (entrez) return(unname(ent2ens[s]))
  out <- unname(sym2ens[s])
  out[s %in% names(COMPLEXES)] <- "complex"
  out
}

regulons_rds <- list(); cov <- list()
for (nm in names(regs)) {
  r <- copy(regs[[nm]])
  entrez <- startsWith(nm, "aracne_")
  r[, source := as.character(source)][, target := as.character(target)]
  r[, source_gene_id := source_gene(source, entrez)]
  r[, target_gene_id := if (entrez) unname(ent2ens[target]) else unname(sym2ens[target])]
  if (entrez) {  # readable names for Entrez-keyed networks
    ens2sym <- setNames(idmap$symbol, idmap$gene_id)
    r[, source := fifelse(is.na(source_gene_id), source, ens2sym[source_gene_id])]
  }
  in_lambert <- function(g, s) g %in% tfs$gene_id | (s %in% names(COMPLEXES))
  cov[[nm]] <- r[, .(n_targets = .N, n_targets_mapped = sum(!is.na(target_gene_id)),
                     source_mapped = !is.na(source_gene_id[1]),
                     in_lambert = in_lambert(source_gene_id[1], source[1])),
                 by = .(resource = rep(nm, nrow(r)), source, source_gene_id)]
  fwrite(r, sprintf("%s/regulons_%s.tsv", out_dir(nm), nm), sep = "\t", na = "NA")
  # viper regulon object on Ensembl IDs; sources keep their name (complexes have no single gene)
  m <- r[!is.na(target_gene_id)]
  m <- m[, .(mor = mean(mor), likelihood = max(likelihood)), by = .(source, target_gene_id)]
  regulons_rds[[nm]] <- lapply(split(m, m$source), function(x)
    list(tfmode = setNames(x$mor, x$target_gene_id), likelihood = x$likelihood))
}
cov <- rbindlist(cov)[, frac_mapped := round(n_targets_mapped / n_targets, 3)]
fwrite(cov, "results/step5/regulon_coverage.tsv", sep = "\t", na = "NA")
aracne <- startsWith(names(regulons_rds), "aracne_")
rds <- c(core = "results/step5/regulons_core.rds", aracne = "data/processed/step5/regulons_aracne.rds")
saveRDS(regulons_rds[!aracne], rds["core"])
saveRDS(regulons_rds[aracne], rds["aracne"])
writeLines(system2("sha256sum", rds, stdout = TRUE), "results/step5/regulon_checksums.sha256")

print(cov[, .(TFs = .N, sources_mapped = sum(source_mapped), in_lambert = sum(in_lambert),
              edges = sum(n_targets), edges_mapped = sum(n_targets_mapped),
              median_frac_mapped = median(frac_mapped),
              TFs_lt_10_mapped = sum(n_targets_mapped < 10)), by = resource])
message("REGULONS DONE")
