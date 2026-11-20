#!/usr/bin/env Rscript
# Flat regulon membership for the FULL ARACNe set (TFs + kinases + signalling), LINCS gene space.
# INTERPRETER: /home/kirti/miniforge3/envs/repurposing/bin/Rscript
suppressPackageStartupMessages(library(data.table))
setwd("/home/kirti/projects/repurposing")
R <- readRDS("data/processed/step5/regulons_aracne_project_v2.rds")
g <- fread("data/processed/lincs/drug_line_profiles_genes.tsv")$gene_id
out <- rbindlist(lapply(names(R), function(n) {
  r <- R[[n]]
  rbindlist(lapply(names(r), function(x) {
    tg <- names(r[[x]]$tfmode)
    k <- tg %in% g
    if (!any(k)) return(NULL)
    data.table(network = n, regulator = x, gene_id = tg[k],
               tfmode = as.numeric(r[[x]]$tfmode)[k], likelihood = as.numeric(r[[x]]$likelihood)[k])
  }))
}))
f <- "data/processed/lincs/trace_aracne/aracne_full_regulon_edges.tsv.gz"
fwrite(out, f, sep = "\t")
cat(sprintf("%s: %d edges, %d networks, %d regulators (LINCS gene space)\n",
            f, nrow(out), uniqueN(out$network), uniqueN(out$regulator)))
