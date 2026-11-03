suppressPackageStartupMessages(library(data.table)); SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"; setwd("~/projects/repurposing")
H <- fread(file.path(SCR, "hidden_factor_edges.tsv.gz"))
rg <- H[, .(n = .N, removed5 = mean(cls_k5 == "removed"), removed10 = mean(cls_k10 == "removed")), by = .(net, regulator)][n >= 10]
cat("per-regulator share removed (k=5): quantiles", round(quantile(rg$removed5, c(.1, .5, .9)), 2), "\n")
print(rg[regulator %in% c("ESR1", "GATA3", "FOXA1", "XBP1", "MYC", "E2F1", "SOX10", "MITF", "CDX2", "NKX2-1", "SP1", "TP53")][order(net, -removed5)], nrows = 50)
# vs 79t drug-side disagreement
ag <- fread("results/lincs/regulon_agreement_v2_tf.tsv")[network %in% unique(rg$net), .(net = network, regulator, rho_all)]
a <- merge(rg, ag, by = c("net", "regulator"))
cat("\n79t: Spearman(removed5, ARACNe-CollecTRI rho) =", round(cor(a$removed5, a$rho_all, method = "s"), 3), " n", nrow(a), "\n")
print(a[, .(n = .N, median_rho = round(median(rho_all), 3), share_disagree = round(mean(rho_all < 0), 3)), by = .(b = cut(removed5, c(0, .4, .6, .8, 1)))][order(b)])
# vs D84 benchmark (V0 on these nets, per TF)
o <- fread(file.path(SCR, "variant_eval_hf_sigs.tsv.gz"))[variant == "V0", .(auc = 1 - mean(fpct), correct = mean(correct), n = .N), by = cmap_name]
rr <- rg[, .(removed5 = mean(removed5)), by = regulator]
b <- merge(o, rr, by.x = "cmap_name", by.y = "regulator")
cat("\nD84 benchmark: Spearman(removed5, per-TF AUC) =", round(cor(b$removed5, b$auc, method = "s"), 3), " n", nrow(b), "\n")
print(b[, .(n = .N, auc = round(mean(auc), 3), correct = round(mean(correct), 3)), by = .(bin = cut(removed5, c(0, .4, .6, .8, 1)))][order(bin)])
