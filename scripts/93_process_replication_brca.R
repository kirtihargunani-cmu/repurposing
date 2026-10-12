#!/usr/bin/env Rscript
# Step 9c: process the breast replication cohorts into patient-level regulator activity per PAM50 subtype, matching the
# TCGA BRCA Step 8 input as closely as each platform allows (same approach as 90_process_replication_crc.R).
#   METABRIC   (Illumina HT-12 array, cBioPortal brca_metabric log2 intensities; datahub commit in DATAHUB_COMMIT)
#              rows -> genes: Entrez -> GENCODE v23 Ensembl (results/step5/id_map.tsv), protein-coding; one row per gene =
#              highest mean expression (MaxMean). Subtype = CLAUDIN_SUBTYPE (PAM50 + claudin-low); Basal / Her2 / LumA /
#              LumB used (claudin-low, Normal, NC not in the TCGA groups).
#   CPTAC-BRCA (GDC STAR counts, GENCODE v36; primary tumors): unversioned Ensembl matched to v23 protein-coding genes,
#              edgeR filterByExpr + TMM logCPM (per subtype). Subtype = cBioPortal PAM50 (brca_cptac_2020); tumors without
#              a PAM50 call or Normal-like dropped.
#   both       ESTIMATE ("affymetrix" setting is for Affymetrix; METABRIC is Illumina, so platform "illumina" for both);
#              gene-level correction within each subtype: tumors regressed on ESTIMATE score (+ the long-minus-short
#              RNA-quality index for CPTAC, not for the array); residual + mean
#   regulons   TCGA BRCA ARACNe network with the mappability filter (genes with mappability < 0.8 removed; edges between
#              cross-mappable genes removed), restricted to the regulators kept for the matching TCGA subtype group
#   reference  each patient vs the cohort-subtype's average tumor (viper method "scale"), as for colorectal
#              (PI 2026-10-09); module scores are standardized within cohort-subtype downstream
# Outputs: data/processed/step9/<cohort>/<group>/{expression,patient_nes}.tsv.gz; results/step9/<cohort>_samples.tsv;
#          results/step9/processing_summary_brca.tsv
suppressPackageStartupMessages({library(data.table); library(edgeR); library(viper); library(estimate)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root); set.seed(93)
RAW <- "data/raw/replication"; MIN_TARGETS <- 10
SUB <- c(Basal = "BRCA-Basal", Her2 = "BRCA-Her2", LumA = "BRCA-LumA", LumB = "BRCA-LumB")
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
sym2ens <- setNames(idmap$gene_id, idmap$symbol); ens2sym <- setNames(idmap$symbol, idmap$gene_id)
pc <- fread("results/step3/genes_protein_coding_v23.tsv")$gene_id
pm <- fread("data/raw/solid_toil/gencode.v23.annotation.gene.probemap")
span <- setNames(pm$chromEnd - pm$chromStart + 1, sub("\\..*", "", pm$id))

# ---- regulons: TCGA BRCA + mappability filter (subtype regulator sets applied in score())
gm <- fread("data/raw/crossmap/hg38_gene_mappability.txt.gz", header = FALSE, col.names = c("id", "map"))
lowmap <- sub("\\..*", "", gm[!is.na(map) & map < 0.8, id])
cm <- fread(cmd = "zcat data/raw/crossmap/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz", header = FALSE,
            col.names = c("a", "b", "s"))[s > 0]
cm[, `:=`(a = sub("\\..*", "", a), b = sub("\\..*", "", b), s = NULL)]; cm <- unique(rbind(cm, cm[, .(a = b, b = a)])); setkey(cm, a, b)
reg0 <- readRDS("data/processed/step5/regulons_aracne_project.rds")$BRCA
rf <- fread("results/step6/aracne_project/regulator_filter.tsv")[group_id %in% SUB & keep == TRUE]
all_regs <- unique(rf$regulator); all_regs <- all_regs[!(all_regs %in% lowmap) & !(sym2ens[all_regs] %in% lowmap) & all_regs %in% names(reg0)]
ens <- vapply(all_regs, function(r) if (startsWith(r, "ENSG")) r else unname(sym2ens[r]), "")
n <- vapply(reg0[all_regs], function(r) length(r$tfmode), 1L)
ed <- data.table(rn = rep(all_regs, n), a = rep(ens, n), b = unlist(lapply(reg0[all_regs], function(r) names(r$tfmode)), use.names = FALSE))
ed[, cross := !is.na(cm[ed, on = .(a, b), mult = "first", which = TRUE])]; rm(cm)
drop <- split(ed$b[ed$cross], ed$rn[ed$cross])
regs <- setNames(lapply(all_regs, function(rn) { r <- reg0[[rn]]; k <- !(names(r$tfmode) %in% c(drop[[rn]], lowmap))
  list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) }), all_regs)
message(sprintf("regulons: %d BRCA regulators (union of subtypes) after mappability filter; %d cross-mappable edges removed", length(regs), sum(ed$cross)))

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
score <- function(cohort, gid, ex) {  # ex: corrected tumor expression (genes x patients) of one subtype
  ex <- ex[!(rownames(ex) %in% lowmap), ]
  keep <- intersect(rf[group_id == gid, regulator], names(regs))
  reg <- lapply(regs[keep], function(r) { k <- names(r$tfmode) %in% rownames(ex); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  reg <- reg[vapply(reg, function(r) length(r$tfmode) >= MIN_TARGETS, logical(1))]; class(reg) <- "regulon"
  nes <- viper(ex, reg, method = "scale", minsize = MIN_TARGETS, verbose = FALSE)
  d <- file.path("data/processed/step9", cohort, gid); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(gene_id = rownames(ex), round(ex, 4)), file.path(d, "expression.tsv.gz"), sep = "\t")
  fwrite(data.table(regulator = rownames(nes), round(nes, 3)), file.path(d, "patient_nes.tsv.gz"), sep = "\t")
  data.table(cohort = cohort, group_id = gid, n_tumors = ncol(ex), n_genes = nrow(ex), n_regulators_scored = nrow(nes),
             frac_tcga_regulators = round(nrow(nes) / length(keep), 3))
}

# ---- METABRIC
mx <- fread(file.path(RAW, "metabric/data_mrna_illumina_microarray.txt"))
X <- as.matrix(mx[, -(1:2)]); gid_row <- idmap$gene_id[match(mx$Entrez_Gene_Id, idmap$entrez_id)]
ok <- !is.na(gid_row) & gid_row %in% pc & rowSums(is.na(X)) == 0
X <- X[ok, ]; gid_row <- gid_row[ok]; o <- order(-rowMeans(X)); keep <- o[!duplicated(gid_row[o])]
G <- X[keep, ]; rownames(G) <- gid_row[keep]; rm(X, mx)
mp <- fread(cmd = paste("grep -v '^#'", file.path(RAW, "metabric/data_clinical_patient.txt")))
cl <- data.table(sample = colnames(G))[mp, on = .(sample = PATIENT_ID), nomatch = 0]
cl[, subtype := CLAUDIN_SUBTYPE]
message(sprintf("METABRIC: %d genes x %d tumors; subtypes: %s", nrow(G), ncol(G), paste(names(table(cl$subtype)), table(cl$subtype), collapse = ", ")))
es <- estimate_scores(G, "illumina"); cl <- merge(cl, es, by = "sample")
s1 <- rbindlist(lapply(names(SUB), function(st) { s <- cl[subtype == st, sample]; pur <- cl[match(s, sample), estimate_score]
  score("metabric", SUB[[st]], correct(G[, s], cbind(1, pur - mean(pur)))) }))
fwrite(cl[, .(sample, subtype, ER_IHC, INTCLUST, AGE_AT_DIAGNOSIS, CELLULARITY, OS_MONTHS, OS_STATUS, RFS_MONTHS, RFS_STATUS,
              estimate_stromal, estimate_immune, estimate_score)], "results/step9/metabric_samples.tsv", sep = "\t", na = "NA")

# ---- CPTAC-BRCA
man <- fread(file.path(RAW, "cptac_brca/gdc_manifest.tsv"))[sample_types == "Primary Tumor"]
cnt <- sapply(man$file_name, function(f) { x <- fread(file.path(RAW, "cptac_brca/star_counts", f), skip = 1)[grepl("^ENSG", gene_id)]
  setNames(x$unstranded, sub("\\..*", "", x$gene_id)) })
colnames(cnt) <- man$case; cnt <- cnt[rownames(cnt) %in% pc, ]; cnt <- cnt[!duplicated(rownames(cnt)), ]
cs <- fread(cmd = paste("grep -v '^#'", file.path(RAW, "cptac_brca/cbioportal/data_clinical_sample.txt")))
cp <- fread(cmd = paste("grep -v '^#'", file.path(RAW, "cptac_brca/cbioportal/data_clinical_patient.txt")))
cc <- merge(cs, cp, by = "PATIENT_ID", all.x = TRUE)[, case := sub("^X", "", PATIENT_ID)]
cl2 <- data.table(sample = colnames(cnt))[cc, on = .(sample = case), nomatch = 0][, subtype := PAM50]
message(sprintf("CPTAC-BRCA: %d of %d tumors with clinical data; subtypes: %s", nrow(cl2), ncol(cnt),
                paste(names(table(cl2$subtype)), table(cl2$subtype), collapse = ", ")))
s2 <- rbindlist(lapply(names(SUB), function(st) { s <- cl2[subtype == st, sample]
  y <- DGEList(cnt[, s]); y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]; y <- calcNormFactors(y)
  L <- cpm(y, log = TRUE, prior.count = 1); e <- estimate_scores(L, "illumina"); qi <- qindex(L)
  cl2[match(e$sample, sample), `:=`(estimate_stromal = e$estimate_stromal, estimate_immune = e$estimate_immune,
                                     estimate_score = e$estimate_score, quality_index = qi[e$sample])]
  pur <- e$estimate_score
  score("cptac_brca", SUB[[st]], correct(L, cbind(1, pur - mean(pur), qi - mean(qi)))) }))
fwrite(cl2[, .(sample, subtype, NMF_CLUSTER, TUMOR_STAGE, ER_UPDATED_CLINICAL_STATUS, ERBB2_UPDATED_CLINICAL_STATUS,
               ESTIMATE_TUMORPURITY, estimate_stromal, estimate_immune, estimate_score, quality_index)],
       "results/step9/cptac_brca_samples.tsv", sep = "\t", na = "NA")
summ <- rbind(s1, s2); fwrite(summ, "results/step9/processing_summary_brca.tsv", sep = "\t"); print(summ)
message("STEP 9 BRCA PROCESSING DONE")
