## Exploratory (2026-11-03): does the luminal-only breast network (BRCA_LUM, D65 / D81) read anti-estrogens in MCF7
## differently from the pooled BRCA network? Same profiles (79c line-specific and the 6 h early signature), same pruning
## as 79d; ESR1 / FOXA1 / GATA3 / XBP1 NES and regulon sizes.
suppressPackageStartupMessages({library(data.table); library(viper)})
setwd("~/projects/repurposing")
arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds"); cat("networks:", paste(grep("BRCA", names(arac), value = TRUE), collapse = ", "), "\n")
pg <- fread("data/processed/lincs/drug_line_profiles_genes.tsv"); idx <- fread("data/processed/lincs/drug_line_profiles_index.tsv")
P <- matrix(readBin("data/processed/lincs/drug_line_profiles_bing_specific.f32", "numeric", size = 4, n = nrow(pg) * nrow(idx)), nrow = nrow(pg)); rownames(P) <- pg$gene_id
cmp <- fread("results/lincs/compounds.tsv")[, .(pert_id, cmap_name)]
ae <- c("tamoxifen", "fulvestrant", "raloxifene", "toremifene", "bazedoxifene", "clomifene")
I <- merge(idx[, .(pert_id, ModelID, col = .I)], cmp, by = "pert_id")[ModelID == "ACH-000019" & cmap_name %in% ae]
prune <- function(rg) { r <- lapply(rg, function(x) { k <- names(x$tfmode) %in% rownames(P); if (sum(k) < 10 || mean(k) < 0.5) return(NULL)
  list(tfmode = x$tfmode[k], likelihood = x$likelihood[k]) }); r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r }
R <- c("ESR1", "FOXA1", "GATA3", "XBP1", "SPDEF", "PGR", "MYB")
for (n in c("BRCA", "BRCA_LUM")) { rg <- prune(arac[[n]])
  v <- viper(P[, I$col], rg, method = "none", minsize = 10, verbose = FALSE)
  cat(sprintf("\n%s (%d regulators): targets %s\n", n, length(rg), paste(sprintf("%s %d", intersect(R, names(rg)), sapply(rg[intersect(R, names(rg))], function(x) length(x$tfmode))), collapse = ", ")))
  print(round(apply(v[intersect(R, rownames(v)), , drop = FALSE], 1, function(x) c(median = median(x), neg = mean(x < 0))), 2)) }
