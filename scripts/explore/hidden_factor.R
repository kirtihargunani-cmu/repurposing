## Exploratory: edges that survive removal of the top k hidden factors (PCs) of the tumor expression, and their
## agreement with TF perturbations. Networks BRCA, LUAD, COADREAD, SKCM; k = 5, 10, 20; 50% retention (D82 cutoff).
suppressPackageStartupMessages({library(data.table); library(parallel)})
setwd("~/projects/repurposing"); SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
NETS <- c("BRCA", "LUAD", "COADREAD", "SKCM"); KS <- c(5, 10, 20); RETAIN <- 0.5; WEAK <- 0.1
reg <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
cor_rows <- function(a, B) as.vector(B %*% a) / (sqrt(sum(a^2)) * sqrt(rowSums(B^2)))

run_net <- function(net) {
  t0 <- Sys.time()
  x <- fread(file.path("data/processed/step5/aracne_inputs_v2", paste0(net, ".exp.tsv")))
  m <- as.matrix(x[, -1]); rownames(m) <- x$gene
  R <- t(apply(m, 1, rank)); R <- R - rowMeans(R)              # genes x tumors, centred ranks
  R <- R[rowSums(R^2) > 0, ]
  S <- crossprod(R / sqrt(rowSums(R^2)))                          # tumors x tumors, genes standardised
  ev <- eigen(S, symmetric = TRUE)
  varexp <- ev$values / sum(ev$values)
  rg <- reg[[net]]
  pairs <- rbindlist(lapply(names(rg), function(r) { ge <- if (startsWith(r, "ENSG")) r else sym2ens[[r]]
    tg <- intersect(names(rg[[r]]$tfmode), rownames(R))
    if (is.null(ge) || is.na(ge) || !ge %in% rownames(R) || !length(tg)) return(NULL)
    data.table(regulator = r, reg_gene = ge, target = tg, mode = rg[[r]]$tfmode[tg]) }))
  res <- pairs[, .(net = net, regulator, target, mode)]
  idx <- split(seq_len(nrow(pairs)), pairs$regulator)
  rc <- function(M) { v <- numeric(nrow(pairs)); for (i in idx) v[i] <- cor_rows(M[pairs$reg_gene[i[1]], ], M[pairs$target[i], , drop = FALSE]); v }
  res[, rho_p := rc(R)]
  for (k in KS) {
    U <- ev$vectors[, 1:k, drop = FALSE]
    Rr <- R - (R %*% U) %*% t(U)                                  # remove k factors from every gene
    v <- rc(Rr)
    res[[paste0("rho_k", k)]] <- v
  }
  message(sprintf("%s: %d tumors, %d edges; variance explained by PC1-5 %.3f, 1-10 %.3f, 1-20 %.3f (%.1f min)", net, ncol(R), nrow(res),
                  sum(varexp[1:5]), sum(varexp[1:10]), sum(varexp[1:20]), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  res
}
H <- rbindlist(mclapply(NETS, run_net, mc.cores = 4))
for (k in KS) { r <- H[[paste0("rho_k", k)]]
  H[[paste0("ret_k", k)]] <- r / H$rho_p
  H[[paste0("cls_k", k)]] <- fifelse(abs(H$rho_p) < WEAK, "weak", fifelse(sign(r) == sign(H$rho_p) & abs(r) >= RETAIN * abs(H$rho_p), "survives", "removed")) }
fwrite(H, file.path(SCR, "hidden_factor_edges.tsv.gz"), sep = "\t")
cat("\nshare of edges per class:\n")
print(rbindlist(lapply(KS, function(k) H[, .(k = k, survives = round(mean(get(paste0("cls_k", k)) == "survives"), 3),
  removed = round(mean(get(paste0("cls_k", k)) == "removed"), 3), weak = round(mean(get(paste0("cls_k", k)) == "weak"), 3)), by = net])))

## perturbation agreement
R <- fread(file.path(SCR, "edge_rows.tsv.gz"))[net %in% NETS]
R <- merge(R, H, by = c("net", "regulator", "target"))
E <- readRDS(file.path(SCR, "edges_features.rds"))[net %in% NETS]; E[, n_nets_reg := uniqueN(net), by = regulator]
E2 <- readRDS(file.path(SCR, "edges_features.rds")); E2[, n_nets_reg := uniqueN(net), by = regulator]; E2[, rfrac := recur / n_nets_reg]
R <- merge(R, E2[net %in% NETS, .(net, regulator, target, rfrac, class_d82 = class)], by = c("net", "regulator", "target"))
R[, agree := y > 0]
cat("\nedge x signature rows on these networks:", nrow(R), " signatures", uniqueN(R$sig_id), " TFs", uniqueN(R$regulator), "; overall agreement", round(mean(R$agree), 3), "\n")
for (k in KS) { cl <- paste0("cls_k", k); cat("\n== k =", k, "\n")
  print(R[, .(n = .N, agree = round(mean(agree), 3), mean_y = round(mean(y), 3)), by = c(cl)])
  rb <- cut(R[[paste0("ret_k", k)]], c(-Inf, 0, .25, .5, .75, 1, Inf))
  print(R[abs(rho_p) >= WEAK, .(n = .N, agree = round(mean(agree), 3)), by = .(retention = rb[abs(R$rho_p) >= WEAK])][order(retention)])
  w <- dcast(R[get(cl) != "weak", .(a = mean(agree), n = .N), by = c("sig_id", cl)], as.formula(paste("sig_id ~", cl)), value.var = c("a", "n"))
  w <- w[n_survives >= 5 & n_removed >= 5]
  cat(sprintf("within-signature survives minus removed: sigs %d, mean %+.3f, share > 0 %.2f, sign test p %.2g\n", nrow(w),
      mean(w$a_survives - w$a_removed), mean(w$a_survives > w$a_removed), binom.test(sum(w$a_survives > w$a_removed), sum(w$a_survives != w$a_removed))$p.value))
  f <- as.formula(paste("agree ~ rfrac +", paste0("I(", cl, " == 'survives')")))
  print(round(summary(glm(f, data = R[get(cl) != "weak"], family = binomial()))$coefficients, 4))
  cat("recurrence by class: median rfrac", paste(R[, .(m = round(median(rfrac), 3)), by = c(cl)][, paste(get(cl), m)], collapse = "; "), "\n")
}
cat("\n== by network (k = 10)\n"); print(dcast(R[, .(a = round(mean(agree), 3)), by = .(net, cls_k10)], net ~ cls_k10))
cat("\n== D82 histology class vs k10 class (BRCA, COADREAD)\n"); print(R[!is.na(class_d82), .N, by = .(class_d82, cls_k10)][order(class_d82)])
## ESR1 in BRCA: which of its edges survive, and do the anti-estrogen risers survive?
d <- fread(file.path(SCR, "esr1_pos_targets.tsv")); gn <- fread("data/processed/lincs/level5_bing_genes.tsv"); s2e <- setNames(gn$ensembl_id, gn$gene_symbol)
d[, target := s2e[gene]]; d <- merge(d, H[net == "BRCA" & regulator == "ESR1"], by = "target")
cat("\nESR1 (BRCA) positive targets: anti-estrogen change in MCF7 by class\n")
for (k in KS) print(d[, .(k = k, n = .N, mean_change = round(mean(anti), 2), share_up = round(mean(anti > 0), 2)), by = c(paste0("cls_k", k))])
print(d[order(-anti)][1:10, .(gene, anti = round(anti, 2), rho_p = round(rho_p, 2), rho_k5 = round(rho_k5, 2), rho_k10 = round(rho_k10, 2), rho_k20 = round(rho_k20, 2))])
