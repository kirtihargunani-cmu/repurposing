#!/usr/bin/env Rscript
# Exploration (2026-11-03): identity vs drug-response links. Networks on the 978 landmark genes and the 81 measured
# regulators (same ARACNe settings): within-line drug networks (D86, line-specific arm, all drugs), DepMap baseline
# networks across cell lines (ALL 1,588 lines; BREAST 69 lines; scripts/pnet/depmap_baseline_inputs.py) and the tumor
# ARACNe v2 TF networks of the matching lineage (MCF7 -> BRCA, PC3 -> PRAD, A375 -> SKCM), restricted to the same genes.
# Tumor links are classed by where else they appear: drug network of the line ("drug response"), DepMap baseline
# ("identity": genes co-vary across lines), both, neither. Each class is checked against the D84 knockdown /
# overexpression signatures of the line (link agrees if the target moves as the link's sign x perturbation direction).
# Outputs results/lincs/pnet/baseline_compare_{overlap,classes,agreement,regulators}.tsv
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
res <- "results/lincs/pnet"; pn <- "data/processed/lincs/pnet"; pd <- "data/processed/lincs/pnet_depmap"
LINES <- data.table(line = c("MCF7", "PC3", "A375"), ModelID = c("ACH-000019", "ACH-000090", "ACH-000219"), tumor = c("BRCA", "PRAD", "SKCM"))
idm <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]; e2s <- setNames(idm$symbol, idm$gene_id); s2e <- setNames(idm$gene_id, idm$symbol)
regs <- readLines(file.path(pd, "regulators.txt")); genes <- fread(file.path(pn, "MCF7_specific/all.exp.tsv"), select = 1)$gene
N <- length(regs) * (length(genes) - 1)
rd <- function(f) { x <- fread(f); setnames(x, c("regulator", "target", "mi", "p")); x[regulator %in% regs & target %in% genes & regulator != target, .(regulator, target, mi)] }
nets <- c(list(DepMap_ALL = rd(file.path(pd, "ALL/net_all/network.txt")), DepMap_BREAST = rd(file.path(pd, "BREAST/net_all/network.txt"))),
          setNames(lapply(LINES$line, function(l) rd(file.path(pn, paste0(l, "_specific/net_all/network.txt")))), paste0("drug_", LINES$line)))
tum <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
for (i in seq_len(nrow(LINES))) {
  rg <- tum[[LINES$tumor[i]]]
  nets[[paste0("tumor_", LINES$tumor[i])]] <- rbindlist(lapply(names(rg), function(r) { e <- if (startsWith(r, "ENSG")) r else s2e[r]
    if (is.na(e) || !e %in% regs) return(NULL); tm <- rg[[r]]$tfmode; k <- names(tm) %in% genes & names(tm) != e
    if (!any(k)) return(NULL); data.table(regulator = e, target = names(tm)[k], mode = sign(tm[k])) }))
}
key <- function(x) paste(x$regulator, x$target)
ov <- rbindlist(lapply(names(nets), function(a) rbindlist(lapply(names(nets), function(b) {
  A <- key(nets[[a]]); B <- key(nets[[b]]); k <- length(intersect(A, B))
  data.table(a = a, b = b, edges_a = length(A), edges_b = length(B), shared = k, share_of_a = round(k / length(A), 3),
             fold_over_chance = round(k / (length(A) * length(B) / N), 2)) }))))
fwrite(ov, file.path(res, "baseline_compare_overlap.tsv"), sep = "\t")

## tumor links by class, and their agreement with the line's knockdown / overexpression signatures
gl <- fread("data/processed/lincs/level5_bing_genes.tsv"); ok_g <- !is.na(gl$ensembl_id) & gl$ensembl_id != "" & !duplicated(gl$ensembl_id)
rdX <- function(f, n) { X <- matrix(readBin(file.path("data/processed/lincs", f), "numeric", size = 4, n = nrow(gl) * n), nrow = nrow(gl)); X <- X[ok_g, , drop = FALSE]; rownames(X) <- gl$ensembl_id[ok_g]; X }
mk <- fread("data/processed/lincs/knockdown_bing_sigs.tsv.gz"); mo <- fread("data/processed/lincs/overexpression_bing_sigs.tsv.gz")
X <- cbind(rdX("knockdown_bing.f32", nrow(mk)), rdX("overexpression_bing.f32", nrow(mo))); colnames(X) <- c(mk$sig_id, mo$sig_id)
sg <- unique(fread("results/lincs/perturbation_benchmark_d84_signatures.tsv.gz")[own_ok == TRUE, .(sig_id, cmap_name, direction, cell_iname)])
cls <- rbindlist(lapply(seq_len(nrow(LINES)), function(i) {
  t <- copy(nets[[paste0("tumor_", LINES$tumor[i])]]); dk <- key(nets[[paste0("drug_", LINES$line[i])]])
  t[, `:=`(line = LINES$line[i], in_drug = key(t) %in% dk, in_all = key(t) %in% key(nets$DepMap_ALL),
           in_breast = if (LINES$line[i] == "MCF7") key(t) %in% key(nets$DepMap_BREAST) else NA)]
  t[, class := fifelse(in_drug & in_all, "both", fifelse(in_drug, "drug response only", fifelse(in_all, "identity only (DepMap)", "neither")))]
  t }))
fwrite(cls[, .(line, regulator = e2s[regulator], target = e2s[target], mode, in_drug, in_all, in_breast, class)], file.path(res, "baseline_compare_classes.tsv"), sep = "\t")
ag <- rbindlist(lapply(seq_len(nrow(LINES)), function(i) {
  s <- sg[cell_iname == LINES$line[i] & sig_id %in% colnames(X)]; s[, regulator := s2e[cmap_name]]
  m <- merge(cls[line == LINES$line[i]], s, by = "regulator", allow.cartesian = TRUE)[target %in% rownames(X)]
  if (!nrow(m)) return(NULL)
  m[, z := X[cbind(match(target, rownames(X)), match(sig_id, colnames(X)))]]
  m[, agree := sign(z) == mode * direction]; m }))
agt <- rbind(ag[, .(rows = .N, signatures = uniqueN(sig_id), regulators = uniqueN(regulator), agree = round(mean(agree), 3)), by = .(line, class)],
             ag[, .(line = "pooled", rows = .N, signatures = uniqueN(sig_id), regulators = uniqueN(regulator), agree = round(mean(agree), 3)), by = class])
# within-signature contrast: drug-response links vs identity-only links of the same perturbation
w <- dcast(ag[class %in% c("drug response only", "identity only (DepMap)", "neither")], sig_id ~ class, value.var = "agree", fun.aggregate = mean)
setnames(w, make.names(names(w)))
fwrite(agt[order(line, class)], file.path(res, "baseline_compare_agreement.tsv"), sep = "\t")

## per regulator in MCF7: where its tumor links sit
rg <- cls[line == "MCF7", .(tumor_links = .N, drug = round(mean(in_drug), 2), depmap_all = round(mean(in_all), 2), depmap_breast = round(mean(in_breast), 2)), by = regulator][order(-tumor_links)]
rg[, regulator := e2s[regulator]]
fwrite(rg, file.path(res, "baseline_compare_regulators.tsv"), sep = "\t")

options(width = 200)
cat("== links per network and overlap (share of row network's links found in column network; x chance)\n")
print(dcast(ov[a != b], a ~ b, value.var = "share_of_a")); print(dcast(ov[a != b], a ~ b, value.var = "fold_over_chance"))
cat("\n== tumor links by class\n"); print(cls[, .N, by = .(line, class)][, share := round(N / sum(N), 3), by = line][order(line, class)])
cat("\n== agreement with the line's knockdown / overexpression signatures\n"); print(agt[order(line, class)])
for (cc in c("drug.response.only", "identity.only..DepMap.")) if (all(c(cc, "neither") %in% names(w))) {
  x <- w[!is.na(get(cc)) & !is.na(neither)]; d <- x[[cc]] - x$neither
  cat(sprintf("within signature, %s minus neither: %d signatures, mean %+.3f, share > 0 %.2f, sign test p %.2g\n", cc, nrow(x), mean(d), mean(d > 0),
              binom.test(sum(d > 0), sum(d != 0))$p.value)) }
cat("\n== MCF7: regulators with the most tumor (BRCA) links on landmark genes\n"); print(head(rg, 25))
