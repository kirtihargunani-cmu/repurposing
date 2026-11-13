#!/usr/bin/env Rscript
# Export an existing VIPER/ARACNe regulon object's signed TF->target edges.
# Usage: Rscript scripts/77b_export_regulon_edges.R REGULONS.rds NETWORK OUT.tsv
# NETWORK may be '-' if the RDS directly contains a single named TF list.
# Output columns: regulator, target, sign, confidence.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3L) {
  stop('Usage: Rscript 77b_export_regulon_edges.R <regulons.rds> <network-name-or-> <output.tsv>')
}
obj <- readRDS(args[[1]])
net <- args[[2]]
if (net != '-') {
  if (!is.list(obj) || is.null(names(obj)) || !net %in% names(obj)) {
    stop('Network not present in RDS: ', net, '. Top-level names: ',
         paste(head(names(obj), 30L), collapse=', '))
  }
  obj <- obj[[net]]
}
# For a viper regulon, each named TF has a list containing named vectors
# tfmode (signed mode of regulation) and likelihood (edge confidence).
valid_entry <- function(x) is.list(x) && !is.null(x$tfmode) &&
  !is.null(x$likelihood) && !is.null(names(x$tfmode))
if (!is.list(obj) || is.null(names(obj)) || !all(vapply(obj, valid_entry, logical(1)))) {
  stop('Selected object is not a named VIPER regulon list with tfmode and likelihood. ',
       'Inspect str(readRDS(file)) and select the correct network object.')
}
parts <- lapply(names(obj), function(tf) {
  x <- obj[[tf]]
  genes <- intersect(names(x$tfmode), names(x$likelihood))
  mor <- as.numeric(x$tfmode[genes])
  conf <- as.numeric(x$likelihood[genes])
  ok <- is.finite(mor) & is.finite(conf) & mor != 0 & conf > 0
  data.frame(regulator=tf, target=genes[ok], sign=sign(mor[ok]),
             confidence=abs(mor[ok]) * conf[ok], stringsAsFactors=FALSE)
})
out <- do.call(rbind, parts)
if (!nrow(out)) stop('No nonzero signed edges to export')
dir.create(dirname(args[[3]]), recursive=TRUE, showWarnings=FALSE)
write.table(out, args[[3]], sep='\t', quote=FALSE, row.names=FALSE)
cat('Exported', nrow(out), 'edges from', length(obj), 'regulators to', args[[3]], '\n')
