#!/usr/bin/env Rscript
# Step 5f: bootstrap support for every edge of the consolidated project ARACNe networks.
#   n_boot   number of the 100 bootstrap networks (scripts/43_run_aracne.sh) that contain the edge
#   mi_boot_mean / mi_boot_sd   MI across the bootstraps that contain it
# Support measures how stable an edge is when patients are resampled; it does not say which patients drive it.
# Usage: scripts/45_edge_support.R [NET ...]   (no args = every network with a finished network.txt)
# Outputs: data/processed/step5/aracne/<NET>/edge_support.tsv.gz (network.txt + support columns)
#          results/step5/aracne_edge_support.tsv (per network: support distribution)
# VERSION=v2 [REGSET=tf]: networks in aracne_v2[_tf]/, summary results/step5/aracne_edge_support_v2[_tf].tsv
#   (support over all consolidated edges, before the cross-mappable filter of scripts/44)
suppressPackageStartupMessages({library(data.table); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
SFX <- if (identical(Sys.getenv("VERSION", "v1"), "v2")) "_v2" else ""
TFX <- if (identical(Sys.getenv("REGSET", "full"), "tf")) "_tf" else ""; stopifnot(TFX == "" || SFX == "_v2")
NETDIR <- paste0("data/processed/step5/aracne", SFX, TFX)
nets <- commandArgs(trailingOnly = TRUE)
if (!length(nets)) nets <- basename(dirname(Sys.glob(file.path(NETDIR, "*", "network.txt"))))

run <- function(net) {
  d <- file.path(NETDIR, net); bf <- Sys.glob(file.path(d, "bootstrapNetwork_*.txt"))
  boots <- rbindlist(lapply(bf, fread, select = 1:3))
  sup <- boots[, .(n_boot = .N, mi_boot_mean = mean(MI), mi_boot_sd = sd(MI)), by = .(Regulator, Target)]
  nw <- fread(file.path(d, "network.txt"))
  nw <- merge(nw, sup, by = c("Regulator", "Target"), all.x = TRUE, sort = FALSE)
  nw[is.na(n_boot), n_boot := 0L]
  fwrite(nw[order(-n_boot, -MI)], file.path(d, "edge_support.tsv.gz"), sep = "\t")
  q <- quantile(nw$n_boot, c(0.1, 0.25, 0.5, 0.75, 0.9))
  data.table(net = net, n_bootstraps = length(bf), edges = nrow(nw), edges_seen_any_bootstrap = nrow(sup),
             support_min = min(nw$n_boot), support_p10 = q[[1]], support_p25 = q[[2]], support_median = q[[3]],
             support_p75 = q[[4]], support_p90 = q[[5]], support_max = max(nw$n_boot),
             frac_ge_50 = round(mean(nw$n_boot >= 50), 3), frac_ge_80 = round(mean(nw$n_boot >= 80), 3),
             rho_support_mi = round(cor(nw$n_boot, nw$MI, method = "spearman"), 3))
}
res <- rbindlist(mclapply(nets, function(n) tryCatch(run(n), error = function(e) { message(n, ": ", conditionMessage(e)); NULL }),
                          mc.cores = 8, mc.preschedule = FALSE))
fwrite(res[order(net)], paste0("results/step5/aracne_edge_support", SFX, TFX, ".tsv"), sep = "\t")
print(as.data.frame(res[order(net)]), row.names = FALSE)
message("EDGE SUPPORT DONE: ", nrow(res), " networks")
