#!/usr/bin/env Rscript
# Step 5 (D83): signaling-only target lists from the TF-only regulons, using the D82 edge classes (48_identity_edges.R).
#   sig1 (main)  drop edges classed identity under the histology labels (networks without histology labels unchanged)
#   sig2         drop edges classed identity under the histology OR the molecular labels
#   base         the current lists unchanged (same file contents; a separate name so the comparison runs side by side)
# Kept edges keep their mode (tfmode) and likelihood. No size filter here: VIPER's own minimum (10 targets present)
# applies to every arm alike; 'lost' = regulators with >= 10 targets before and < 10 after.
# Outputs: data/processed/step5/regulons_aracne_project_v2_tf_{base,sig1,sig2}.rds;
#          results/step5/identity_edges_v2_tf/signaling_lists.tsv (per network x arm: edges before / after, regulators lost)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_TARGETS <- 10
IE <- "results/step5/identity_edges_v2_tf"
reg <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
ident <- function(net, sets) {
  ff <- file.path(IE, sprintf("edges_%s_%s.tsv.gz", net, sets)); ff <- ff[file.exists(ff)]
  if (!length(ff)) return(NULL)
  unique(rbindlist(lapply(ff, function(f) fread(f, select = c("regulator", "target", "class"))[class == "identity", .(regulator, target)])))
}
log <- list()
for (arm in c("base", "sig1", "sig2")) {
  sets <- switch(arm, base = character(), sig1 = "histology", sig2 = c("histology", "molecular"))
  out <- lapply(setNames(names(reg), names(reg)), function(net) {
    r <- reg[[net]]; id <- if (length(sets)) ident(net, sets) else NULL
    if (is.null(id) || !nrow(id)) {
      log[[paste(arm, net)]] <<- data.table(arm = arm, network = net, regulators = length(r), edges = sum(lengths(lapply(r, `[[`, "tfmode"))),
                                            regulators_after = length(r), edges_after = sum(lengths(lapply(r, `[[`, "tfmode"))), lost = "", n_lost = 0L)
      return(r)
    }
    drop <- split(id$target, id$regulator)
    r2 <- lapply(setNames(names(r), names(r)), function(g) {
      x <- r[[g]]; k <- !names(x$tfmode) %in% drop[[g]]
      list(tfmode = x$tfmode[k], likelihood = x$likelihood[k])
    })
    n2 <- vapply(r2, function(x) length(x$tfmode), integer(1))
    n1 <- vapply(r, function(x) length(x$tfmode), integer(1))[names(n2)]
    lost <- names(n2)[n1 >= MIN_TARGETS & n2 < MIN_TARGETS]
    r2 <- r2[n2 > 0]; class(r2) <- class(r)
    log[[paste(arm, net)]] <<- data.table(arm = arm, network = net, regulators = length(r), edges = sum(lengths(lapply(r, `[[`, "tfmode"))),
                                          regulators_after = length(r2), edges_after = sum(n2), lost = paste(lost, collapse = ","), n_lost = length(lost))
    r2
  })
  saveRDS(out, sprintf("data/processed/step5/regulons_aracne_project_v2_tf_%s.rds", arm))
}
lg <- rbindlist(log); fwrite(lg, file.path(IE, "signaling_lists.tsv"), sep = "\t")
print(lg[arm != "base" & edges_after < edges, .(arm, network, edges, edges_after, pct_kept = round(100 * edges_after / edges, 1), n_lost,
                                               lost = substr(lost, 1, 70))])
