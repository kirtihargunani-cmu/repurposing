## Exploratory: do CollecTRI edges move as predicted under their TF's knockdown / overexpression (D84 signatures)?
suppressPackageStartupMessages(library(data.table)); setwd("~/projects/repurposing")
proc <- "data/processed/lincs"; genes <- fread(file.path(proc, "level5_bing_genes.tsv"))
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
ld <- function(f, ids) { X <- matrix(readBin(file.path(proc, f), "numeric", size = 4, n = nrow(genes) * length(ids)), nrow = nrow(genes))
  X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- ids; X }
mk <- fread(file.path(proc, "knockdown_bing_sigs.tsv.gz")); mo <- fread(file.path(proc, "overexpression_bing_sigs.tsv.gz"))
X <- cbind(ld("knockdown_bing.f32", mk$sig_id), ld("overexpression_bing.f32", mo$sig_id))
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g]); landmark <- genes$ensembl_id[ok_g][genes$feature_space[ok_g] == "landmark"]
sg <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[own_ok == TRUE]
col <- readRDS("results/step5/regulons_core.rds")$collectri
C <- rbindlist(lapply(names(col), function(r) data.table(regulator = r, target = names(col[[r]]$tfmode), mode = col[[r]]$tfmode)))
C[, size := .N, by = regulator]
ar <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds"); ar$BRCA_LUM <- NULL
A <- unique(rbindlist(lapply(ar, function(n) rbindlist(lapply(names(n), function(r) data.table(regulator = r, target = names(n[[r]]$tfmode)))))))
A[, in_aracne := TRUE]; C <- merge(C, A, by = c("regulator", "target"), all.x = TRUE); C[is.na(in_aracne), in_aracne := FALSE]
setkey(C, regulator)
one <- function(st) { s <- unique(sg[set == st, .(sig_id, cmap_name, technology, direction)])
  rbindlist(lapply(seq_len(nrow(s)), function(i) { e <- C[.(s$cmap_name[i]), nomatch = 0L][target %in% rownames(X) & target != sym2ens[s$cmap_name[i]]]
    if (!nrow(e)) return(NULL)
    e[, `:=`(sig_id = s$sig_id[i], technology = s$technology[i], y = pmin(pmax(X[target, s$sig_id[i]], -5), 5) * sign(mode) * s$direction[i])] })) }
R <- one("collectri"); R[, agree := y > 0]
cat("CollecTRI edge x signature rows:", nrow(R), " signatures", uniqueN(R$sig_id), " TFs", uniqueN(R$regulator), "\n")
cat("overall agreement", round(mean(R$agree), 3), " mean y", round(mean(R$y), 3), "\n")
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = .(sign = fifelse(mode > 0, "activation", "repression"))])
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = technology])
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = .(also_in_aracne = in_aracne)])
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = .(regulon_size = cut(size, c(0, 25, 100, 300, 5000)))][order(regulon_size)])
print(R[, .(n = .N, agree = round(mean(agree), 3)), by = .(landmark = target %in% landmark)])
t <- R[, .(agree = mean(agree), n = .N), by = regulator][n >= 50]
cat("per TF (>= 50 edge x signature rows):", nrow(t), "TFs; agreement quantiles", round(quantile(t$agree, c(.1, .25, .5, .75, .9)), 2),
    "; TFs below 0.5:", sum(t$agree < .5), "\n")
print(t[regulator %in% c("ESR1", "MYC", "TP53", "E2F1", "SP1", "SP3", "AR", "NR3C1", "STAT3", "RELA", "JUN", "HIF1A")][order(-agree)])
