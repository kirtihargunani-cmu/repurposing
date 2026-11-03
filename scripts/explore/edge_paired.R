suppressPackageStartupMessages({library(data.table)})
setwd("~/projects/repurposing")
SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
E <- readRDS(file.path(SCR, "edges_features.rds")); setkey(E, net, regulator, target)
R <- fread(file.path(SCR, "edge_rows.tsv.gz"))
R <- merge(R, E[, .(net, regulator, target, mode, recur)], by = c("net", "regulator", "target"))
sg <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[set == "aracne_v2_tf" & own_ok == TRUE, .(sig_id, aracne_matched, cell_iname)]
R <- merge(R, sg, by = "sig_id")
R[, rec := fifelse(recur >= 8, "recur>=8", fifelse(recur <= 2, "recur<=2", "mid"))]
R[, agree := y > 0]
# within-signature contrast: agreement of recurrent minus non-recurrent edges, per signature with both
w <- dcast(R[rec != "mid", .(a = mean(agree), n = .N), by = .(sig_id, regulator, technology, aracne_matched, rec)],
           sig_id + regulator + technology + aracne_matched ~ rec, value.var = c("a", "n"))
w <- w[`n_recur>=8` >= 5 & `n_recur<=2` >= 5]
w[, d := `a_recur>=8` - `a_recur<=2`]
cat("signatures with >= 5 edges in both bins:", nrow(w), " TFs:", uniqueN(w$regulator), "\n")
cat("within-signature diff (recur>=8 minus <=2): mean", round(mean(w$d), 3), " median", round(median(w$d), 3),
    " share > 0:", round(mean(w$d > 0), 3), " sign test p:", signif(binom.test(sum(w$d > 0), sum(w$d != 0))$p.value, 3), "\n")
print(w[, .(sigs = .N, mean_d = round(mean(d), 3), share_pos = round(mean(d > 0), 3)), by = .(technology, aracne_matched)])
# per TF (average over its signatures)
t <- w[, .(d = mean(d)), by = regulator]
cat("per TF: n", nrow(t), " mean", round(mean(t$d), 3), " share > 0", round(mean(t$d > 0), 3), "\n")
# does the effect hold when the signature's TF is barely moved (weak) vs strongly read? use per-sig mean y as strength
st <- R[, .(s = mean(y)), by = sig_id]; w <- merge(w, st, by = "sig_id")
w[, strength := cut(s, quantile(s, 0:3 / 3), include.lowest = TRUE, labels = c("weak", "mid", "strong"))]
print(w[, .(sigs = .N, mean_d = round(mean(d), 3), share_pos = round(mean(d > 0), 3)), by = strength][order(strength)])
# mode sign within signature
w2 <- dcast(R[, .(a = mean(agree), n = .N), by = .(sig_id, ms = fifelse(mode > 0, "pos", "neg"))], sig_id ~ ms, value.var = c("a", "n"))[n_neg >= 5 & n_pos >= 5]
cat("\nwithin-signature pos minus neg mode: sigs", nrow(w2), " mean", round(mean(w2$a_pos - w2$a_neg), 3), " share > 0", round(mean(w2$a_pos > w2$a_neg), 3), "\n")
# ESR1 BRCA: recurrence of its targets and anti-ER response (from earlier scratch file, symbols)
d <- fread(file.path(SCR, "esr1_pos_targets.tsv"))
gn <- fread("data/processed/lincs/level5_bing_genes.tsv"); s2e <- setNames(gn$ensembl_id, gn$gene_symbol)
d[, target := s2e[gene]]; d <- merge(d, E[net == "BRCA" & regulator == "ESR1", .(target, recur)], by = "target")
d[, rb := cut(recur, c(0, 2, 7, 28))]
cat("\nESR1 (BRCA) positive targets under anti-estrogens in MCF7 by recurrence:\n")
print(d[, .(n = .N, mean_change = round(mean(anti), 2), share_up = round(mean(anti > 0), 2)), by = rb][order(rb)])
print(d[order(-anti)][1:10, .(gene, anti = round(anti, 2), recur)])
cat("\nedges per network by recurrence (all regulators):\n"); print(E[, .(share = round(.N / nrow(E), 3)), by = .(rb = cut(recur, c(0, 1, 2, 3, 7, 15, 28)))][order(rb)])
