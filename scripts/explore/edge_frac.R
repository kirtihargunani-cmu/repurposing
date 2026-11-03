suppressPackageStartupMessages({library(data.table)})
setwd("~/projects/repurposing"); SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
E <- readRDS(file.path(SCR, "edges_features.rds"))
E[, n_nets_reg := uniqueN(net), by = regulator]; E[, rfrac := recur / n_nets_reg]
R <- merge(fread(file.path(SCR, "edge_rows.tsv.gz")), E[, .(net, regulator, target, mode, recur, rfrac, n_nets_reg)], by = c("net", "regulator", "target"))
R[, agree := y > 0]
cat("regulators by number of networks with a regulon:\n"); print(table(cut(unique(E[, .(regulator, n_nets_reg)])$n_nets_reg, c(0, 5, 10, 20, 28))))
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = .(fb = cut(rfrac, c(0, .1, .2, .35, .5, .75, 1)))][order(fb)])
m1 <- glm(agree ~ log(recur), data = R, family = binomial()); m2 <- glm(agree ~ rfrac, data = R, family = binomial())
m3 <- glm(agree ~ log(recur) + rfrac + log(n_nets_reg), data = R, family = binomial())
cat("AIC recur", round(AIC(m1)), " rfrac", round(AIC(m2)), " both", round(AIC(m3)), "\n"); print(round(summary(m3)$coefficients, 4))
cat("\nregulators with regulons in <= 10 networks: agreement by recurrence fraction\n")
print(R[n_nets_reg <= 10, .(n = .N, agree = round(mean(agree), 3)), by = .(fb = cut(rfrac, c(0, .2, .5, 1)))][order(fb)])
