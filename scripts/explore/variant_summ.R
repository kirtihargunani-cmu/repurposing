suppressPackageStartupMessages(library(data.table)); SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
o <- fread(file.path(SCR, "variant_eval_sigs.tsv.gz"))
cat("overall (each variant on its own scored signatures):\n")
print(o[, .(sigs = .N, tfs = uniqueN(cmap_name), correct = round(mean(correct), 3), auc = round(1 - mean(fpct), 3),
            correct_kd = round(mean(correct[direction < 0]), 3), correct_oe = round(mean(correct[direction > 0]), 3), n_reg = as.numeric(median(n_reg))), by = variant])
# common signatures (scored by all variants)
com <- o[, .N, by = sig_id][N == 4, sig_id]; oc <- o[sig_id %in% com]
cat("\ncommon signatures:", length(com), "TFs", uniqueN(oc$cmap_name), "\n")
print(oc[, .(correct = round(mean(correct), 3), auc = round(1 - mean(fpct), 3)), by = variant])
t <- oc[, .(c = mean(correct), a = 1 - mean(fpct), n = .N), by = .(cmap_name, variant)]
w <- dcast(t, cmap_name + n ~ variant, value.var = "a")
for (v in c("V1", "V2", "V3")) { d <- w[[v]] - w$V0
  cat(sprintf("%s vs V0 per-TF AUC: mean diff %+.3f, better %d, worse %d, Wilcoxon p %.2g\n", v, mean(d), sum(d > 0), sum(d < 0), wilcox.test(w[[v]], w$V0, paired = TRUE)$p.value)) }
# pass counts (D84 rule) per variant on common signatures
p <- oc[, .(n = .N, mp = median(fpct), fc = mean(correct)), by = .(cmap_name, variant)][n >= 3, .(tested = .N, pass = sum(mp <= .25 & fc >= .6)), by = variant]; print(p)
# split-half CV by TF: choose best variant on half A by AUC, report gain over V0 on half B
set.seed(1); tfs <- unique(t$cmap_name)
cv <- rbindlist(lapply(1:200, function(i) { A <- sample(tfs, length(tfs) %/% 2)
  sa <- t[cmap_name %in% A, .(a = mean(a)), by = variant]; best <- sa[which.max(a), variant]
  sb <- t[!cmap_name %in% A, .(a = mean(a)), by = variant]
  data.table(best = best, gain = sb[variant == best, a] - sb[variant == "V0", a]) }))
cat("\nsplit-half CV (200 splits): best variant chosen", paste(names(table(cv$best)), table(cv$best), collapse = ", "),
    "; held-out AUC gain vs V0 mean", round(mean(cv$gain), 4), " 5-95%:", paste(round(quantile(cv$gain, c(.05, .95)), 4), collapse = " to "), "\n")
print(oc[, .(auc = round(1 - mean(fpct), 3)), by = .(variant, technology)][order(technology, variant)])
