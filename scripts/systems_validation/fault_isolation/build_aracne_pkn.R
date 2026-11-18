#!/usr/bin/env Rscript
# F-3b: export the FULL ARACNe regulator set (TFs + kinases + signalling) as a CARNIVAL PKN.
# Frozen by docs/FREEZE_F3b_aracne_arm.md: UNPRUNED, sign = sign(tfmode).
# INTERPRETER: /home/kirti/miniforge3/envs/repurposing/bin/Rscript
suppressPackageStartupMessages(library(data.table))
setwd("/home/kirti/projects/repurposing")
OUT <- "results/systems_validation/pkn"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
R <- readRDS("data/processed/step5/regulons_aracne_project_v2.rds")
idm <- fread("results/step5/id_map.tsv")[primary_for_hgnc==TRUE]; e2s <- setNames(idm$symbol, idm$gene_id)
for (net in commandArgs(trailingOnly=TRUE)) {
  r <- R[[net]]; stopifnot(!is.null(r))
  ed <- rbindlist(lapply(names(r), function(g) data.table(
        source=g, interaction=ifelse(r[[g]]$tfmode>=0,1L,-1L),
        target=unname(e2s[names(r[[g]]$tfmode)]),
        tfmode=as.numeric(r[[g]]$tfmode), likelihood=as.numeric(r[[g]]$likelihood))))
  ed <- ed[!is.na(target) & source!=target]
  ed <- ed[order(-likelihood)][!duplicated(paste(source,target))]
  f <- file.path(OUT, sprintf("aracne_pkn_%s.tsv", net))
  fwrite(ed[, .(source, interaction, target)], f, sep="\t")
  fwrite(ed, file.path(OUT, sprintf("aracne_pkn_%s_full.tsv", net)), sep="\t")
  cat(sprintf("%-10s %7d edges, %5d regulators, %6d nodes, %.0f%% activating -> %s\n",
      net, nrow(ed), uniqueN(ed$source), uniqueN(c(ed$source,ed$target)),
      100*mean(ed$interaction==1), f))
}
