#!/usr/bin/env Rscript
# Step 7, descriptive companion to the BRCA_LUM test (66_brca_lum_test.R; NOT part of the D81 rule): what do the luminal
# regulators' regulons contain on each network? D81 tests whether activity scores change; a regulon can shift from
# luminal-identity genes to estrogen-responsive genes while the regulator's score keeps its sign. This reports that shift.
# Networks: BRCA (pooled), BRCA_LUM, BRCA_SUB1-3 (size-matched controls). Regulators: ESR1 (headline), FOXA1, GATA3, XBP1.
# Per regulator x network, on the full regulon (before the per-group filters of 50_viper.R; cross-mappable edges removed):
#   n_targets, n_pos (targets with positive mode of regulation)
#   overlap with pooled BRCA: shared targets, Jaccard, fraction of this network's targets also in the pooled regulon
#   estrogen response: Hallmark ESTROGEN_RESPONSE_EARLY + _LATE (MSigDB h.all v2025.1, union; symbols -> Ensembl),
#     restricted to the network's gene universe (genes in its ARACNe input). n_er, frac_er, fold enrichment over the
#     universe, Fisher p (greater), and frac_er_pos (fraction of those ER genes with positive mode; ER mostly activates)
# Reading: the controls show how much the regulon moves with size alone; BRCA_LUM's change beyond that is the luminal effect.
# Output: results/step7/brca_lum_test/regulon_content.tsv (+ regulon_content_targets.tsv: ESR1 targets per network)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
REGS <- c("ESR1", "FOXA1", "GATA3", "XBP1")
out <- "results/step7/brca_lum_test"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

prod <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
cf <- "data/processed/step5/regulons_aracne_project_v2_tf_brca_controls.rds"
ctrl <- if (file.exists(cf)) readRDS(cf) else list()
nets <- c(list(BRCA = prod$BRCA, BRCA_LUM = prod$BRCA_LUM), ctrl[intersect(paste0("BRCA_SUB", 1:3), names(ctrl))])
if (length(nets) < 5) warning("controls missing: ", paste(setdiff(paste0("BRCA_SUB", 1:3), names(nets)), collapse = ", "))

idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
gmt <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
er_sym <- unique(unlist(lapply(gmt[vapply(gmt, `[`, "", 1) %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE")], `[`, -(1:2))))
er <- unique(idmap[symbol %in% er_sym, gene_id])
message("Hallmark estrogen response: ", length(er_sym), " symbols, ", length(er), " mapped to Ensembl")

universe <- function(net) fread(file.path("data/processed/step5/aracne_inputs_v2", paste0(net, ".exp.tsv")), select = 1)[[1]]
tg <- function(net, r) { x <- nets[[net]][[r]]; if (is.null(x)) NULL else x$tfmode }

res <- rbindlist(lapply(REGS, function(r) rbindlist(lapply(names(nets), function(net) {
  m <- tg(net, r); if (is.null(m)) return(data.table(regulator = r, network = net, n_targets = 0L))
  u <- universe(net); t <- intersect(names(m), u); erU <- intersect(er, u); pooled <- names(tg("BRCA", r))
  hit <- intersect(t, erU)
  ft <- fisher.test(matrix(c(length(hit), length(t) - length(hit), length(erU) - length(hit),
                             length(u) - length(t) - length(erU) + length(hit)), 2), alternative = "greater")
  data.table(regulator = r, network = net, n_targets = length(m), n_pos = sum(m > 0),
             shared_with_pooled = length(intersect(names(m), pooled)),
             jaccard_pooled = round(length(intersect(names(m), pooled)) / length(union(names(m), pooled)), 3),
             frac_in_pooled = round(mean(names(m) %in% pooled), 3),
             n_er = length(hit), frac_er = round(length(hit) / length(t), 4),
             fold_er = round((length(hit) / length(t)) / (length(erU) / length(u)), 2), fisher_p = signif(ft$p.value, 3),
             frac_er_pos = round(mean(m[hit] > 0), 3))
}))), fill = TRUE)
fwrite(res, file.path(out, "regulon_content.tsv"), sep = "\t")

sym <- setNames(idmap$symbol, idmap$gene_id)
tl <- rbindlist(lapply(names(nets), function(net) { m <- tg(net, "ESR1"); if (is.null(m)) return(NULL)
  data.table(network = net, target = names(m), symbol = sym[names(m)], mode = round(m, 3), hallmark_er = names(m) %in% er) }))
fwrite(tl, file.path(out, "regulon_content_targets.tsv"), sep = "\t")
print(res)
