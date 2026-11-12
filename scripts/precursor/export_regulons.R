# Precursor (reported only): export the regulon target weights for the breast networks, so the gene-level
# restricted-reversal check can run in python. Writes regulator x target rows with tfmode and likelihood.
# Nothing production is read beyond the regulon objects and nothing production is written.
suppressPackageStartupMessages(library(data.table))
NETS <- c("BRCA_BAS", "BRCA", "BRCA_LUMA", "BRCA_LUMB")       # BRCA-Basal / Her2 / LumA / LumB, per 79e's net map
OUT  <- "data/processed/lincs/precursor_regulons"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

dump <- function(reg, tag) {
  rows <- rbindlist(lapply(names(reg), function(n) {
    r <- reg[[n]]
    data.table(regulator = n, target = names(r$tfmode), tfmode = as.numeric(r$tfmode),
               likelihood = as.numeric(r$likelihood))
  }))
  fwrite(rows, file.path(OUT, paste0(tag, ".tsv.gz")), sep = "\t")
  cat(sprintf("%-12s %5d regulators %8d edges\n", tag, length(reg), nrow(rows)))
}

a <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
for (n in NETS) dump(a[[n]], n)
core <- readRDS("results/step5/regulons_core.rds")
cat("core sets:", paste(names(core), collapse = ", "), "\n")
if ("collectri" %in% names(core)) dump(core[["collectri"]], "collectri")
