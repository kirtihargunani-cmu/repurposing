#!/usr/bin/env Rscript
# D125 A1 step 2: convert the 79n NES dump (D125_NES_DIR, one .rds per set x network) into the project's f32 +
# names convention so the sweep can read it from python. Reads and writes only inside the dump directory.
# Per tag: <tag>.nes.f32 (regulators x signatures, column-major float32), <tag>.regulators.txt, <tag>.sigs.txt,
#          <tag>.targets.tsv.gz (regulator, target: the pruned target sets as used for the NES).
suppressPackageStartupMessages(library(data.table))
d <- Sys.getenv("D125_NES_DIR", "data/processed/lincs/d125_a1")
f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
stopifnot(length(f) > 0)
for (p in f) {
  tag <- sub("\\.rds$", "", basename(p))
  if (file.exists(file.path(d, paste0(tag, ".nes.f32")))) next
  x <- readRDS(p)
  v <- x$nes
  writeBin(as.vector(as.numeric(v)), file.path(d, paste0(tag, ".nes.f32")), size = 4)
  writeLines(rownames(v), file.path(d, paste0(tag, ".regulators.txt")))
  writeLines(colnames(v), file.path(d, paste0(tag, ".sigs.txt")))
  tg <- x$targets[rownames(v)]                    # regulons in NES row order; prep() already pruned them
  fwrite(data.table(regulator = rep(names(tg), lengths(tg)), target = unlist(tg, use.names = FALSE)),
         file.path(d, paste0(tag, ".targets.tsv.gz")), sep = "\t")
  message(sprintf("%s: %d regulators x %d signatures, %d regulon edges", tag, nrow(v), ncol(v), length(unlist(tg))))
}
message("A1 EXPORT DONE")
