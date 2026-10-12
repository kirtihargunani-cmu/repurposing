#!/usr/bin/env Rscript
# Step 9a: process the colorectal replication cohorts into patient-level regulator activity, matching the TCGA COAD
# Step 8 input as closely as each platform allows. Rectal clusters cannot be replicated here (GSE39582 is colon only;
# CPTAC-COAD has 1 rectal tumor, dropped), so this targets the TCGA-COAD 2-cluster result.
#   GSE39582  (Affymetrix U133 Plus 2, author-normalized log2 values from the series matrix; 566 tumors + 19 normals)
#             probes -> genes: probes annotated to exactly one Entrez gene (GPL570 annot), Entrez -> GENCODE v23 Ensembl
#             (results/step5/id_map.tsv); one probe per gene = highest mean expression ("MaxMean", Miller et al. 2011)
#   CPTAC-COAD (GDC STAR counts, GENCODE v36; 105 colon tumors): unversioned Ensembl matched to v23 protein-coding genes,
#             edgeR filterByExpr + TMM logCPM
#   both      ESTIMATE scores (21_estimate.R approach; platform "affymetrix" / "illumina"); gene-level correction as in
#             61_quality_correct.R: tumors regressed on ESTIMATE score (+ the long-minus-short RNA-quality index for
#             CPTAC; not for the array, where degradation acts through probe position, not gene length); residual + mean
#   regulons  TCGA COADREAD ARACNe regulons with the mappability filter (genes with mappability < 0.8 removed; edges
#             between cross-mappable genes removed) and the regulators kept for TCGA-COAD (50_viper.R filter)
#   reference PI decision 2026-10-09: each patient vs the cohort's average tumor (viper method = "scale"); clustering
#             compares patients with each other, so this matches the TCGA module scores after standardization
# Outputs: data/processed/step9/<cohort>/{expression,patient_nes}.tsv.gz (+ expression_raw / counts_raw for CMS); results/step9/<cohort>_samples.tsv;
#          results/step9/processing_summary.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(viper); library(estimate)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(90)
RAW <- "data/raw/replication"; MIN_TARGETS <- 10
dir.create("results/step9", recursive = TRUE, showWarnings = FALSE)
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol); ens2sym <- setNames(idmap$symbol, idmap$gene_id)
pc <- fread("results/step3/genes_protein_coding_v23.tsv")$gene_id
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))

# ---- regulons: TCGA COADREAD + mappability filter, TCGA-COAD regulator set
gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
lowmap <- sub("\\..*", "", gm[!is.na(map) & map < 0.8, id])
cm <- fread(cmd = "zcat data/raw/crossmap/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz", header = FALSE,
            col.names = c("a", "b", "s"))[s > 0]
cm[, `:=`(a = sub("\\..*", "", a), b = sub("\\..*", "", b), s = NULL)]; cm <- unique(rbind(cm, cm[, .(a = b, b = a)])); setkey(cm, a, b)
reg0 <- readRDS("data/processed/step5/regulons_aracne_project.rds")$COADREAD
keep_regs <- fread("results/step6/aracne_project/regulator_filter.tsv")[group_id == "TCGA-COAD" & keep == TRUE, regulator]
keep_regs <- keep_regs[!(keep_regs %in% lowmap) & !(sym2ens[keep_regs] %in% lowmap)]
ens <- vapply(keep_regs, function(r) if (startsWith(r, "ENSG")) r else unname(sym2ens[r]), "")
n <- vapply(reg0[keep_regs], function(r) length(r$tfmode), 1L)
ed <- data.table(rn = rep(keep_regs, n), a = rep(ens, n), b = unlist(lapply(reg0[keep_regs], function(r) names(r$tfmode)), use.names = FALSE))
ed[, cross := !is.na(cm[ed, on = .(a, b), mult = "first", which = TRUE])]; rm(cm)
drop <- split(ed$b[ed$cross], ed$rn[ed$cross])
regs <- setNames(lapply(keep_regs, function(rn) { r <- reg0[[rn]]; k <- !(names(r$tfmode) %in% c(drop[[rn]], lowmap))
  list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) }), keep_regs)
message(sprintf("regulons: %d TCGA-COAD regulators after mappability filter; %d cross-mappable edges removed", length(regs), sum(ed$cross)))

estimate_scores <- function(lx, platform) {  # lx: genes (v23 Ensembl) x samples, log scale
  sym <- ens2sym[rownames(lx)]; data(common_genes, package = "estimate", envir = environment())
  ok <- !is.na(sym) & sym %in% common_genes$GeneSymbol; m <- lx[ok, , drop = FALSE]; rownames(m) <- sym[ok]
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  gct <- tempfile(fileext = ".gct"); out <- tempfile(fileext = ".gct"); ids <- sprintf("S%04d", seq_len(ncol(m)))
  writeLines(c("#1.2", paste(nrow(m), ncol(m), sep = "\t")), gct)
  fwrite(data.table(NAME = rownames(m), Description = rownames(m), setNames(as.data.table(m), ids)), gct, sep = "\t", append = TRUE, col.names = TRUE)
  capture.output(estimateScore(gct, out, platform = platform))
  s <- fread(out, skip = 2, header = TRUE); v <- as.matrix(s[, -(1:2)]); rownames(v) <- s$NAME
  message(sprintf("ESTIMATE (%s): %d of %d genes", platform, nrow(m), nrow(common_genes)))
  data.table(sample = colnames(lx), estimate_stromal = v["StromalScore", ], estimate_immune = v["ImmuneScore", ],
             estimate_score = v["ESTIMATEScore", ])
}
qindex <- function(m) { m <- m[apply(m, 1, median) >= 1 & rownames(m) %in% names(span), , drop = FALSE]
  z <- (m - rowMeans(m)) / pmax(apply(m, 1, sd), 1e-6); q <- quantile(span[rownames(z)], c(0.2, 0.8))
  colMeans(z[span[rownames(z)] >= q[2], , drop = FALSE]) - colMeans(z[span[rownames(z)] <= q[1], , drop = FALSE]) }
correct <- function(m, X) t(qr.resid(qr(X), t(m))) + rowMeans(m)
score <- function(cohort, ex, samples) {  # ex: corrected tumor expression (genes x patients)
  ex <- ex[!(rownames(ex) %in% lowmap), ]
  reg <- lapply(regs, function(r) { k <- names(r$tfmode) %in% rownames(ex); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  reg <- reg[vapply(reg, function(r) length(r$tfmode) >= MIN_TARGETS, logical(1))]; class(reg) <- "regulon"
  nes <- viper(ex, reg, method = "scale", minsize = MIN_TARGETS, verbose = FALSE)
  d <- file.path("data/processed/step9", cohort); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(gene_id = rownames(ex), round(ex, 4)), file.path(d, "expression.tsv.gz"), sep = "\t")
  fwrite(data.table(regulator = rownames(nes), round(nes, 3)), file.path(d, "patient_nes.tsv.gz"), sep = "\t")
  fwrite(samples, file.path("results/step9", paste0(cohort, "_samples.tsv")), sep = "\t", na = "NA")
  data.table(cohort = cohort, n_tumors = ncol(ex), n_genes = nrow(ex), n_regulators_scored = nrow(nes),
             frac_tcga_regulators = round(nrow(nes) / length(regs), 3))
}

# ---- GSE39582
sm <- readLines(gzfile(file.path(RAW, "gse39582/GSE39582_series_matrix.txt.gz")))
b <- grep("^!series_matrix_table_begin", sm); e <- grep("^!series_matrix_table_end", sm)
ex <- fread(text = sm[(b + 1):(e - 1)]); probes <- ex$ID_REF; X <- as.matrix(ex[, -1]); rownames(X) <- probes
geo <- gsub('"', "", strsplit(sm[grep("^!Sample_geo_accession", sm)], "\t")[[1]][-1])
chars <- lapply(grep("^!Sample_characteristics_ch1", sm, value = TRUE), function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
cl <- data.table(sample = geo)
for (v in chars) { k <- sub(":.*", "", v[which(grepl(":", v))[1]]); if (is.na(k)) next
  val <- ifelse(grepl(":", v), trimws(sub("^[^:]*:", "", v)), NA); val[val %in% c("N/A", "")] <- NA
  set(cl, j = make.names(k), value = val) }
ann <- fread(cmd = paste("zcat", file.path(RAW, "gse39582/GPL570.annot.gz"), "| awk 'f; /^ID\\t/{f=1; print}'"), sep = "\t", quote = "")
ann <- ann[!grepl("///", `Gene ID`) & nzchar(`Gene ID`), .(probe = ID, entrez = as.integer(`Gene ID`))]
ann[, gene_id := idmap$gene_id[match(entrez, idmap$entrez_id)]]; ann <- ann[!is.na(gene_id) & gene_id %in% pc & probe %in% probes]
ann[, mean_expr := rowMeans(X[probe, , drop = FALSE])]
best <- ann[order(-mean_expr)][!duplicated(gene_id)]
G <- X[best$probe, ]; rownames(G) <- best$gene_id
message(sprintf("GSE39582: %d single-gene probes -> %d genes", nrow(ann), nrow(G)))
tum <- cl[dataset != "Non Tumoral", sample]
es <- estimate_scores(G, "affymetrix"); cl <- merge(cl, es, by = "sample")
Gt <- G[, tum]; pur <- cl[match(tum, sample), estimate_score]
Gc <- correct(Gt, cbind(1, pur - mean(pur)))
dir.create("data/processed/step9/gse39582", recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(gene_id = rownames(G), round(G, 4)), "data/processed/step9/gse39582/expression_raw.tsv.gz", sep = "\t")  # uncorrected, all samples (CMS)
s1 <- score("gse39582", Gc, cl[, role := fifelse(dataset == "Non Tumoral", "normal", "tumor")])

# ---- CPTAC-COAD
man <- fread(file.path(RAW, "cptac_coad/gdc_manifest.tsv"))[primary_site == "colon"]
cnt <- sapply(man$file_name, function(f) { x <- fread(file.path(RAW, "cptac_coad/star_counts", f), skip = 1)[grepl("^ENSG", gene_id)]
  setNames(x$unstranded, sub("\\..*", "", x$gene_id)) })
colnames(cnt) <- man$case; cnt <- cnt[rownames(cnt) %in% pc, ]
cnt <- cnt[!duplicated(rownames(cnt)), ]
dir.create("data/processed/step9/cptac_coad", recursive = TRUE, showWarnings = FALSE)
fwrite(data.table(gene_id = rownames(cnt), cnt), "data/processed/step9/cptac_coad/counts_raw.tsv.gz", sep = "\t")  # unfiltered counts (CMS)
y <- DGEList(cnt); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
L <- cpm(y, log = TRUE, prior.count = 1)
cs <- fread(cmd = paste("grep -v '^#'", file.path(RAW, "cptac_coad/cbioportal/data_clinical_sample.txt")))
cp <- fread(cmd = paste("grep -v '^#'", file.path(RAW, "cptac_coad/cbioportal/data_clinical_patient.txt")))
cl2 <- merge(data.table(sample = colnames(L)), merge(cs, cp, by = "PATIENT_ID", all.x = TRUE), by.x = "sample", by.y = "PATIENT_ID", all.x = TRUE)
es2 <- estimate_scores(L, "illumina"); cl2 <- merge(cl2, es2, by = "sample"); qi <- qindex(L); cl2[, quality_index := qi[sample]]
pur2 <- cl2[match(colnames(L), sample), estimate_score]; q2 <- qi[colnames(L)]
Lc <- correct(L, cbind(1, pur2 - mean(pur2), q2 - mean(q2)))
message(sprintf("CPTAC-COAD: %d colon tumors x %d genes (filterByExpr)", ncol(L), nrow(L)))
s2 <- score("cptac_coad", Lc, cl2[, role := "tumor"])

summ <- rbind(s1, s2)
summ[, msi_like := c(sprintf("dMMR %d / pMMR %d", cl[role == "tumor" & mmr.status == "dMMR", .N], cl[role == "tumor" & mmr.status == "pMMR", .N]),
                     sprintf("MSI-H %d / MSS %d", cl2[MSI_STATUS == "MSI-H", .N], cl2[MSI_STATUS == "MSS", .N]))]
fwrite(summ, "results/step9/processing_summary.tsv", sep = "\t"); print(summ)
message("STEP 9 CRC PROCESSING DONE")
