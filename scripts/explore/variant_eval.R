## Exploratory: regulon-level effect of recurrence corrections on the D84 perturbation benchmark (ARACNe v2 TF)
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
setwd("~/projects/repurposing"); SCR <- "/tmp/claude-1014/-home-kirti/814325cb-1e21-49ee-b5b6-d8f7ba4623d2/scratchpad"
MIN_TARGETS <- 10; MIN_FRAC <- 0.5; NC <- 28
proc <- "data/processed/lincs"
genes <- fread(file.path(proc, "level5_bing_genes.tsv"))
ok_g <- !is.na(genes$ensembl_id) & genes$ensembl_id != "" & !duplicated(genes$ensembl_id)
ld <- function(f, ids) { X <- matrix(readBin(file.path(proc, f), "numeric", size = 4, n = nrow(genes) * length(ids)), nrow = nrow(genes))
  X <- X[ok_g, ]; rownames(X) <- genes$ensembl_id[ok_g]; colnames(X) <- ids; X }
mk <- fread(file.path(proc, "knockdown_bing_sigs.tsv.gz")); mo <- fread(file.path(proc, "overexpression_bing_sigs.tsv.gz"))
X <- cbind(ld("knockdown_bing.f32", mk$sig_id), ld("overexpression_bing.f32", mo$sig_id))
sym2ens <- setNames(genes$ensembl_id[ok_g], genes$gene_symbol[ok_g])
sg <- fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[set == "aracne_v2_tf" & own_ok == TRUE]
X <- X[, sg$sig_id]
E <- readRDS(file.path(SCR, "edges_features.rds")); E[, n_nets_reg := uniqueN(net), by = regulator]; E[, rfrac := recur / n_nets_reg]
setkey(E, net, regulator)

build <- function(variant) {
  e <- switch(variant, V0 = E, V1 = copy(E)[, lik := lik * rfrac], V2 = E[recur > 1], V3 = E[rfrac >= 0.2])
  lapply(split(e, by = "net"), function(d) { r <- lapply(split(d, by = "regulator"), function(x) list(tfmode = setNames(x$mode, x$target), likelihood = x$lik)); r })
}
prep <- function(rg, orig) {
  r <- lapply(names(rg), function(n) { tm <- rg[[n]]$tfmode; lk <- rg[[n]]$likelihood; e <- sym2ens[n]
    k <- names(tm) %in% rownames(X) & (is.na(e) | names(tm) != e)
    if (sum(k) < MIN_TARGETS || sum(k) / length(tm) < MIN_FRAC) return(NULL)
    list(tfmode = tm[k], likelihood = lk[k]) })
  names(r) <- names(rg); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r
}
score_net <- function(reg, net) {
  cols <- sg[grepl(paste0("(^|;)", net, "(;|$)"), nets), sig_id]
  if (!length(cols)) return(NULL)
  v <- viper(X[, cols, drop = FALSE], reg, method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  m <- sg[match(cols, sig_id)]; i <- match(m$cmap_name, rownames(v)); ok <- !is.na(i)
  pct <- apply(v, 2, function(x) (rank(x) - 1) / (length(x) - 1))
  data.table(net = net, sig_id = cols[ok], nes = v[cbind(i[ok], which(ok))], pct = pct[cbind(i[ok], which(ok))], n_reg = nrow(v))
}
out <- rbindlist(lapply(c("V0", "V1", "V2", "V3"), function(vn) { t0 <- Sys.time(); A <- build(vn)
  r <- rbindlist(mclapply(names(A), function(n) score_net(prep(A[[n]]), n), mc.cores = NC))
  r <- r[, .(nes = mean(nes), pct = mean(pct), n_reg = round(mean(n_reg))), by = sig_id]
  message(sprintf("%s: %d signatures scored, %.1f min", vn, nrow(r), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  cbind(variant = vn, r) }))
out <- merge(out, sg[, .(sig_id, cmap_name, technology, direction)], by = "sig_id")
out[, `:=`(fpct = fifelse(direction < 0, pct, 1 - pct), correct = sign(nes) == direction)]
fwrite(out, file.path(SCR, "variant_eval_sigs.tsv.gz"), sep = "\t")
