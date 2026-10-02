#!/usr/bin/env Rscript
# Step 6b: compare group-level regulator activity (msviper NES) across groups at a chosen FDR.
# REGULON=collectri (default); FDR_CUT=0.10 (default).
# Outputs in results/step6/<REGULON>/:
#   fdr<cut>_group_summary.tsv   per group: regulators tested / passing (up, down), |NES| needed to pass,
#                                max |NES|, analytic-null comparison, group sizes
#   fdr<cut>_nes_matrix.tsv      regulators passing in >= 1 group x groups (NES; FDR in a parallel table)
#   fdr<cut>_fdr_matrix.tsv
#   controls.tsv                 known-biology checks (expected direction vs observed)
# Significance uses the permutation (normal-split) null only. ANALYTIC=1 adds msviper's analytic-null NES
# (assumes independent genes; decoupleR's default) as an internal sensitivity check - off by default and not
# used for any call, because regulon targets are co-expressed (PI decision 2026-10-02).
suppressPackageStartupMessages({library(data.table); library(viper)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
REGULON <- Sys.getenv("REGULON", "collectri"); CUT <- as.numeric(Sys.getenv("FDR_CUT", "0.10"))
d <- file.path("results/step6", REGULON); tag <- sprintf("fdr%02d", round(CUT * 100))

mr <- fread(file.path(d, "group_mr.tsv"))
filt <- fread(file.path(d, "regulator_filter.tsv"))[keep == TRUE]
gs <- fread(file.path(d, "group_summary.tsv"))
regs_all <- readRDS("results/step5/regulons_core.rds")[[REGULON]]
zmat <- list(gtex = fread("results/step4/z_tumor_vs_gtex.tsv"), adjacent = fread("results/step4/z_tumor_vs_adjacent.tsv"))

ANALYTIC <- Sys.getenv("ANALYTIC", "0") == "1"
analytic <- if (!ANALYTIC) data.table(group_id = character(), regulator = character(), nes_analytic = numeric(),
                                      p_analytic = numeric()) else rbindlist(lapply(seq_len(nrow(gs)), function(i) {
  gid <- gs$group_id[i]; z <- zmat[[gs$contrast[i]]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]
  sig <- setNames(z[[gid]], z$gene_id)
  reg <- lapply(regs_all[filt[group_id == gid, regulator]], function(r) {
    k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  a <- msviper(sig, reg, minsize = 10, verbose = FALSE)
  data.table(group_id = gid, regulator = names(a$es$nes), nes_analytic = a$es$nes, p_analytic = a$es$p.value)
}))
analytic[, fdr_analytic := p.adjust(p_analytic, "BH"), by = group_id]
mr <- merge(mr, analytic, by = c("group_id", "regulator"), all.x = TRUE)

summ <- mr[, .(tested = .N, pass = sum(fdr < CUT), up = sum(fdr < CUT & nes > 0), down = sum(fdr < CUT & nes < 0),
               nes_needed = if (any(fdr < CUT)) round(min(abs(nes[fdr < CUT])), 2) else NA_real_,
               max_abs_nes = round(max(abs(nes)), 2),
               pass_analytic = if (ANALYTIC) sum(fdr_analytic < CUT) else NA_integer_,
               rho_perm_vs_analytic = if (ANALYTIC) round(cor(nes, nes_analytic, method = "spearman"), 3) else NA_real_),
            by = group_id]
summ <- merge(summ, gs[, .(group_id, contrast, n_tumor, n_normal)], by = "group_id")
setorder(summ, -pass)
fwrite(summ, file.path(d, paste0(tag, "_group_summary.tsv")), sep = "\t")
cat(sprintf("Groups with >= 1 regulator at FDR < %.2f: %d of %d\n", CUT, summ[pass > 0, .N], nrow(summ)))
cat(sprintf("Spearman(passing count, n_tumor) = %.2f; (passing count, n_normal) = %.2f\n",
            cor(summ$pass, summ$n_tumor, method = "spearman"), cor(summ$pass, summ$n_normal, method = "spearman")))

hits <- mr[fdr < CUT, unique(regulator)]
wide <- function(v) dcast(mr[regulator %in% hits], regulator ~ group_id, value.var = v)
nm <- wide("nes"); fm <- wide("fdr")
nm[, n_groups_pass := rowSums(as.matrix(fm[, -1]) < CUT, na.rm = TRUE)]
nm[, n_up := rowSums(as.matrix(fm[, -1]) < CUT & as.matrix(nm[, 2:(ncol(fm))]) > 0, na.rm = TRUE)]
setorder(nm, -n_groups_pass)
fwrite(nm[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], file.path(d, paste0(tag, "_nes_matrix.tsv")), sep = "\t", na = "NA")
fwrite(fm[, lapply(.SD, function(x) if (is.numeric(x)) signif(x, 3) else x)], file.path(d, paste0(tag, "_fdr_matrix.tsv")), sep = "\t", na = "NA")
cat(sprintf("%d regulators pass FDR < %.2f in >= 1 group\n", length(hits), CUT))

# known-biology checks: expected sign of activity in tumor vs normal
ctl <- rbindlist(list(
  data.table(regulator = c("E2F1", "FOXM1", "MYBL2"), expected = 1, why = "proliferation",
             group_id = list(unique(mr$group_id))),
  data.table(regulator = c("ESR1", "FOXA1", "GATA3"), expected = 1, why = "luminal breast identity",
             group_id = list(c("BRCA-LumA", "BRCA-LumB"))),
  data.table(regulator = c("MYC", "FOXM1"), expected = 1, why = "basal breast", group_id = list("BRCA-Basal")),
  data.table(regulator = c("EPAS1", "HIF1A"), expected = 1, why = "VHL loss / hypoxia", group_id = list("TCGA-KIRC")),
  data.table(regulator = c("HNF4A", "HNF1A"), expected = -1, why = "loss of hepatocyte identity", group_id = list("TCGA-LIHC")),
  data.table(regulator = c("AR", "ERG"), expected = 1, why = "androgen signaling / TMPRSS2-ERG", group_id = list("TCGA-PRAD")),
  data.table(regulator = "MYC", expected = 1, why = "WNT-MYC", group_id = list(c("TCGA-COAD", "TCGA-READ")))
))
ctl <- ctl[, .(group_id = unlist(group_id)), by = .(regulator, expected, why)]
ctl <- merge(ctl, mr[, .(group_id, regulator, nes = round(nes, 2), fdr = signif(fdr, 2))],
             by = c("group_id", "regulator"), all.x = TRUE)
ctl[, verdict := fifelse(is.na(nes), "not tested (filtered)", fifelse(sign(nes) != expected, "wrong direction",
                 fifelse(fdr < CUT, "correct, passes", "correct direction, below cutoff")))]
fwrite(ctl, file.path(d, "controls.tsv"), sep = "\t", na = "NA")
print(ctl[, .N, by = verdict])
