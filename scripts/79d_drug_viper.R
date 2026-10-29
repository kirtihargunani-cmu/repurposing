#!/usr/bin/env Rscript
# Regulator activity of every drug x cell-line profile, with the project ARACNe regulons (user, 2026-10-12).
#
# Same regulons and rules as the disease side (scripts/50_viper.R, REGULON=aracne_project): one network per cohort
# (COAD / READ -> COADREAD), full regulator set (TFs, kinases, signaling proteins). For each network used by a group
# with a disease signature, the regulators scored are those kept for at least one of its groups
# (results/step6/aracne_project/regulator_filter.tsv), re-checked on the LINCS gene universe with the same size rules
# (>= MIN_TARGETS targets present and >= MIN_FRAC of the mapped targets present); targets pruned to LINCS genes.
# Input: drug x line consensus profiles on landmark + best-inferred genes (scripts/79c_drug_profiles.py; MODZ z-scores,
# so viper is run with method = "none", no further scaling). PROFILE=specific (default: line generic response removed
# at the gene level, see 79c) | raw. ONLY_NETS (comma-separated) limits the run; N_CORES networks in parallel.
# Outputs: data/processed/lincs/drug_viper_<PROFILE>/<net>.f32 (float32, regulators x profiles, column-major = profile-major)
#          + <net>.regulators.txt; results/lincs/drug_viper_regulators.tsv (per network x regulator: targets in
#          ARACNe, in LINCS, kept)
# REGULON=collectri | dorothea_ABC (user 2026-10-28, ESR1 finding): curated, context-free regulons (results/step5/
#   regulons_core.rds, as the disease side) scored once for all diseases: network "ALL", every regulator with >= MIN_TARGETS
#   and >= MIN_FRAC targets on LINCS genes (disease-side expression filters act in 79e through the tested MRs).
#   Outputs data/processed/lincs/drug_viper_<PROFILE>_<REGULON>/ALL.f32, results/lincs/drug_viper_regulators_<REGULON>.tsv
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)
MIN_TARGETS <- 10; MIN_FRAC <- 0.5
PROFILE <- Sys.getenv("PROFILE", "specific"); N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
SFX <- if (identical(Sys.getenv("VERSION", "v1"), "v2")) "_v2" else ""
REGSET <- Sys.getenv("REGSET", "full"); stopifnot(REGSET %in% c("full", "tf"), REGSET == "full" || SFX == "_v2")
TFX <- if (REGSET == "tf") "_tf" else ""
REGULON <- Sys.getenv("REGULON", "aracne_project"); stopifnot(REGULON %in% c("aracne_project", "collectri", "dorothea_ABC"))
TAG <- if (REGULON == "aracne_project") paste0(SFX, TFX) else paste0("_", REGULON)  # curated: not versioned
ONLY_NETS <- Filter(nzchar, strsplit(Sys.getenv("ONLY_NETS", ""), ",")[[1]])
proc <- "data/processed/lincs"; out_dir <- file.path(proc, paste0("drug_viper_", PROFILE, TAG))
dir.create(out_dir, showWarnings = FALSE)

genes <- fread(file.path(proc, "drug_line_profiles_genes.tsv"))
idx <- fread(file.path(proc, "drug_line_profiles_index.tsv"))
pfile <- if (PROFILE == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32"
P <- matrix(readBin(file.path(proc, pfile), "numeric", size = 4,
                    n = nrow(genes) * nrow(idx)), nrow = nrow(genes))
dimnames(P) <- list(genes$gene_id, paste(idx$pert_id, idx$ModelID, sep = "|"))
cat(sprintf("profiles (%s): %d genes x %d drug-line profiles\n", PROFILE, nrow(P), ncol(P)))

if (REGULON == "aracne_project") {
  arac <- readRDS(paste0("data/processed/step5/regulons_aracne_project", SFX, TFX, ".rds"))
  grp <- fread("results/step1/disease_groups_toil.tsv")
  net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
  groups <- fread(paste0("results/step4", SFX, "/signature_summary.tsv"))[status == "ok", unique(group_id)]
  filt <- fread(paste0("results/step6/aracne_project", SFX, TFX, "/regulator_filter.tsv"))[keep == TRUE & group_id %in% groups]
  filt[, net := net_of[group_id]]
} else {
  arac <- list(ALL = readRDS("results/step5/regulons_core.rds")[[REGULON]])
  filt <- data.table(net = "ALL", regulator = names(arac$ALL))
}

nets <- sort(unique(filt$net)); if (length(ONLY_NETS)) nets <- intersect(nets, ONLY_NETS)
reg_tab <- mclapply(nets, function(nt) {
  net <- nt  # data.table: "net" inside [ ] is the column
  t0 <- Sys.time()
  regs_all <- arac[[net]]
  cand <- intersect(unique(filt[net == nt, regulator]), names(regs_all))
  st <- rbindlist(lapply(cand, function(r) {
    tg <- names(regs_all[[r]]$tfmode); n_in <- sum(tg %in% rownames(P))
    data.table(net = net, regulator = r, n_targets = length(tg), n_lincs = n_in, frac_lincs = round(n_in / length(tg), 3))
  }))
  st[, kept := n_lincs >= MIN_TARGETS & frac_lincs >= MIN_FRAC]
  reg <- lapply(regs_all[st[kept == TRUE, regulator]], function(r) {
    k <- names(r$tfmode) %in% rownames(P); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  class(reg) <- "regulon"
  nes <- viper(P, reg, method = "none", minsize = MIN_TARGETS, verbose = FALSE)  # regulators x profiles
  stopifnot(identical(colnames(nes), colnames(P)))
  writeBin(as.numeric(nes), file.path(out_dir, paste0(net, ".f32")), size = 4)
  writeLines(rownames(nes), file.path(out_dir, paste0(net, ".regulators.txt")))
  st[, scored := regulator %in% rownames(nes)]
  message(sprintf("%-9s candidates %4d, kept on LINCS genes %4d, scored %4d (%.1f min)", net, nrow(st), sum(st$kept),
                  nrow(nes), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  st
}, mc.cores = N_CORES, mc.preschedule = FALSE)
bad <- vapply(reg_tab, inherits, logical(1), "try-error"); if (any(bad)) stop("failed: ", paste(nets[bad], collapse = ", "))
if (!length(ONLY_NETS)) fwrite(rbindlist(reg_tab), paste0("results/lincs/drug_viper_regulators", TAG, ".tsv"), sep = "\t")
