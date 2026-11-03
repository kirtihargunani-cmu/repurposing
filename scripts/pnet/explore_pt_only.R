## Exploratory (not D86): perturbation network (P) vs tumor ARACNe (T) on regulators both score (CollecTRI not required),
## D84 knockdown / overexpression signatures in MCF7, PC3, A375, landmark genes. Also P vs CollecTRI on regulators both score.
suppressPackageStartupMessages({library(data.table); library(viper)}); setwd("~/projects/repurposing")
MIN_T <- 10; LINES <- c("MCF7", "PC3", "A375"); TNET <- c(MCF7 = "BRCA", PC3 = "PRAD", A375 = "SKCM")
R <- readRDS("data/processed/lincs/pnet/regulons.rds")
genes <- fread("data/processed/lincs/level5_bing_genes.tsv"); ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
ld <- function(f, ids) { X <- matrix(readBin(file.path("data/processed/lincs", f), "numeric", size = 4, n = nrow(genes) * length(ids)), nrow = nrow(genes))
  X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- ids; X }
mk <- fread("data/processed/lincs/knockdown_bing_sigs.tsv.gz"); mo <- fread("data/processed/lincs/overexpression_bing_sigs.tsv.gz")
X <- cbind(ld("knockdown_bing.f32", mk$sig_id), ld("overexpression_bing.f32", mo$sig_id))
lm <- genes$ensembl_id[ok_g & genes$feature_space == "landmark"]; X <- X[lm, ]; sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])
sg <- unique(fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[own_ok == TRUE & cell_iname %in% LINES, .(sig_id, cmap_name, cell_iname, technology, direction)])
prep <- function(rg) { r <- lapply(names(rg), function(n) { tm <- rg[[n]]$tfmode; lk <- rg[[n]]$likelihood; e <- sym2ens[n]
  k <- names(tm) %in% lm & (is.na(e) | names(tm) != e); if (sum(k) < MIN_T) return(NULL); list(tfmode = tm[k], likelihood = lk[k]) })
  names(r) <- names(rg); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r }
nes_of <- function(reg, M, tf) { n <- ncol(M); if (n == 1) M <- cbind(M, M); v <- viper(M, reg, method = "none", minsize = MIN_T, verbose = FALSE); v[cbind(match(tf, rownames(v)), seq_len(n))] }
Tn <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds"); C_ <- prep(readRDS("results/step5/regulons_core.rds")$collectri)
out <- list()
for (l in LINES) { T_ <- prep(Tn[[TNET[l]]]); for (a in c("raw", "specific")) { P_ <- prep(R[[sprintf("%s_%s/all", l, a)]])
  for (cmp in c("T", "C")) { O <- if (cmp == "T") T_ else C_; cm <- intersect(names(P_), names(O)); s <- sg[cell_iname == l & cmap_name %in% cm]
    if (!nrow(s)) next; M <- X[, s$sig_id, drop = FALSE]
    out[[paste(l, a, cmp)]] <- data.table(line = l, arm = a, vs = cmp, s, ok_P = sign(nes_of(P_, M, s$cmap_name)) == s$direction, ok_O = sign(nes_of(O, M, s$cmap_name)) == s$direction) } } }
S <- rbindlist(out)
fwrite(S, "results/lincs/pnet/explore_pt_only_signatures.tsv", sep = "\t")
r <- S[, { b <- sum(ok_P & !ok_O); w <- sum(!ok_P & ok_O); .(signatures = .N, tfs = uniqueN(cmap_name), P = round(mean(ok_P), 3), other = round(mean(ok_O), 3),
          P_minus_other = round(mean(ok_P) - mean(ok_O), 3), P_only = b, other_only = w, sign_p = signif(if (b + w) binom.test(b, b + w)$p.value else 1, 2)) }, by = .(vs, arm)]
print(r); print(S[, .(sigs = .N, P = round(mean(ok_P), 3), other = round(mean(ok_O), 3)), by = .(vs, arm, line)][order(vs, arm, line)])
