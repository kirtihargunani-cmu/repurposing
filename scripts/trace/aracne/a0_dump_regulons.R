#!/usr/bin/env Rscript
# TRACE-ARACNe stage 0: signed ARACNe regulons -> one flat table in the LINCS gene space.
#
# Source: data/processed/step5/regulons_aracne_project_v2_tf.rds (scripts/44_aracne_regulons.R,
#   viper::aracne2regulon on the per-cohort TCGA tumour matrix; tfmode = signed mode of regulation,
#   likelihood = edge weight). Regulators are HGNC symbols, targets Ensembl.
# This is the ARACNe analogue of results/lincs/d128/a1/regulon_edges.parquet, which is CollecTRI-only.
#
# Output: data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz
#   network regulator gene_id tfmode likelihood sign in_lincs_gene_space feature_space
suppressPackageStartupMessages({library(data.table)})
setwd("/home/kirti/projects/repurposing")
IN  <- "data/processed/step5/regulons_aracne_project_v2_tf.rds"
OUT <- "data/processed/lincs/trace_aracne/aracne_regulon_edges.tsv.gz"

g <- fread("data/processed/lincs/drug_line_profiles_genes.tsv")        # gene_id symbol feature_space
setkey(g, gene_id)
R <- readRDS(IN)
cat(sprintf("networks: %d\n", length(R)))

out <- rbindlist(lapply(names(R), function(net) {
  reg <- R[[net]]
  rbindlist(lapply(names(reg), function(tf) {
    tm <- reg[[tf]]$tfmode; lk <- reg[[tf]]$likelihood
    data.table(network = net, regulator = tf, gene_id = names(tm),
               tfmode = as.numeric(tm), likelihood = as.numeric(lk))
  }))
}))
out[, sign := ifelse(tfmode >= 0, 1L, -1L)]
out[, in_lincs_gene_space := gene_id %in% g$gene_id]
out[g, feature_space := i.feature_space, on = "gene_id"]
out[is.na(feature_space), feature_space := "not measured"]

fwrite(out, OUT, sep = "\t")
cat(sprintf("wrote %s: %d edges, %d networks, %d regulators, %.1f%% in LINCS gene space\n",
            OUT, nrow(out), uniqueN(out$network), uniqueN(out$regulator),
            100 * mean(out$in_lincs_gene_space)))
print(out[, .(regulators = uniqueN(regulator), edges = .N,
              in_lincs = sum(in_lincs_gene_space)), by = network][order(network)][1:10])
