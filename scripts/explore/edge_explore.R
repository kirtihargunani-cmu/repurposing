## Exploratory: which ARACNe edges are confirmed by TF perturbations (D84 knockdown + overexpression)?
suppressPackageStartupMessages({library(data.table)})
setwd("~/projects/repurposing")
SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
proc <- "data/processed/lincs"
genes <- fread(file.path(proc, "level5_bing_genes.tsv"))
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
ld <- function(f, ids) { X <- matrix(readBin(file.path(proc, f), "numeric", size = 4, n = nrow(genes) * length(ids)), nrow = nrow(genes))
  X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- ids; X }
mk <- fread(file.path(proc, "knockdown_bing_sigs.tsv.gz")); mo <- fread(file.path(proc, "overexpression_bing_sigs.tsv.gz"))
X <- cbind(ld("knockdown_bing.f32", mk$sig_id), ld("overexpression_bing.f32", mo$sig_id))
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])
landmark <- genes$ensembl_id[ok_g][genes$feature_space[ok_g] == "landmark"]

sg <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[set == "aracne_v2_tf" & own_ok == TRUE]
cat("signatures:", nrow(sg), "TFs:", uniqueN(sg$cmap_name), "\n")
generic <- rowMeans(abs(pmin(pmax(X, -5), 5)))      # how much a gene moves under any perturbation

arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds"); arac$BRCA_LUM <- NULL
nets <- names(arac)
E <- rbindlist(lapply(nets, function(n) rbindlist(lapply(names(arac[[n]]), function(r)
  data.table(net = n, regulator = r, target = names(arac[[n]][[r]]$tfmode), mode = arac[[n]][[r]]$tfmode, lik = arac[[n]][[r]]$likelihood)))))
E[, hub := .N, by = .(net, target)]                                        # regulators listing this target in the net
E[, hub_pct := frank(hub) / .N, by = net]
E[, recur := .N, by = .(regulator, target, sgn = sign(mode))]              # networks with the same signed edge
E[, target_is_reg := target %in% sym2ens[unique(regulator)]]
col <- readRDS("results/step5/regulons_core.rds")$collectri
C <- rbindlist(lapply(names(col), function(r) data.table(regulator = r, target = names(col[[r]]$tfmode), cmode = col[[r]]$tfmode)))
E <- merge(E, C, by = c("regulator", "target"), all.x = TRUE)
E[, curated := fifelse(is.na(cmode), "not_in_collectri", fifelse(sign(cmode) == sign(mode), "collectri_same_sign", "collectri_opposite"))]
idf <- list.files("results/step5/identity_edges_v2_tf", "^edges_.*_histology.tsv.gz$", full.names = TRUE)
ID <- rbindlist(lapply(idf, function(f) { d <- fread(f, select = c("regulator", "target", "class", "retention"))
  d[, net := sub("edges_(.*)_histology.tsv.gz", "\\1", basename(f))] }))
E <- merge(E, ID, by = c("net", "regulator", "target"), all.x = TRUE)
gmt <- strsplit(readLines("data/raw/msigdb/h.all.v2025.1.Hs.symbols.gmt"), "\t")
prolif <- sym2ens[unique(unlist(lapply(gmt[sapply(gmt, `[`, 1) %in% c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")], `[`, -(1:2))))]
E[, prolif := target %in% prolif]
E[, generic := generic[target]]
setkey(E, net, regulator)

# edge x signature outcomes
rows <- rbindlist(lapply(seq_len(nrow(sg)), function(i) { s <- sg[i]; ns <- strsplit(s$nets, ";")[[1]]
  e <- E[.(ns, s$cmap_name), nomatch = 0L][target %in% rownames(X) & target != sym2ens[s$cmap_name]]
  if (!nrow(e)) return(NULL)
  z <- X[e$target, s$sig_id]
  e[, `:=`(sig_id = s$sig_id, technology = s$technology, z = z, y = pmin(pmax(z, -5), 5) * sign(mode) * s$direction,
           tf_nes = s$nes * s$direction)] }))
rows[, agree := y > 0]
cat("edge x signature rows:", nrow(rows), "; overall agreement", round(mean(rows$agree), 3), "mean y", round(mean(rows$y), 3), "\n")
fwrite(rows[, .(net, regulator, target, sig_id, technology, y)], file.path(SCR, "edge_rows.tsv.gz"), sep = "\t")

summ <- function(by) rows[, .(n = .N, agree = round(mean(agree), 3), mean_y = round(mean(y), 3)), by = by][order(get(by[1]))]
rows[, mode_sign := fifelse(mode > 0, "pos", "neg")]
rows[, absmode_bin := cut(abs(mode), c(0, .3, .5, .7, .9, 1), include.lowest = TRUE)]
rows[, lik_bin := cut(lik, quantile(lik, 0:4 / 4), include.lowest = TRUE)]
rows[, hub_bin := cut(hub_pct, c(0, .5, .75, .9, .97, 1), include.lowest = TRUE)]
rows[, recur_bin := cut(recur, c(0, 1, 3, 7, 15, 28))]
rows[, generic_bin := cut(generic, quantile(generic, 0:4 / 4), include.lowest = TRUE)]
rows[, landmark := target %in% landmark]
rows[, strong_tf := tf_nes > 0 & FALSE]  # placeholder
for (b in list("mode_sign", "absmode_bin", "lik_bin", "hub_bin", "recur_bin", "target_is_reg", "curated", "class",
               "prolif", "generic_bin", "landmark", "technology")) { cat("\n==", b, "\n"); print(summ(b)) }
cat("\n== mode_sign x technology\n"); print(rows[, .(n = .N, agree = round(mean(agree), 3)), by = .(mode_sign, technology)][order(mode_sign)])
cat("\n== mode_sign x hub_bin\n"); print(dcast(rows[, .(a = round(mean(agree), 3)), by = .(mode_sign, hub_bin)], hub_bin ~ mode_sign))
cat("\n== mode_sign x recur_bin\n"); print(dcast(rows[, .(a = round(mean(agree), 3)), by = .(mode_sign, recur_bin)], recur_bin ~ mode_sign))
# joint model: which features carry independent signal
m <- glm(agree ~ mode_sign + abs(mode) + log(lik) + hub_pct + log(recur) + target_is_reg + curated + prolif + generic + landmark + technology,
         data = rows, family = binomial())
print(round(summary(m)$coefficients, 4))
saveRDS(E, file.path(SCR, "edges_features.rds"))
