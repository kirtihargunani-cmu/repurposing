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
# REGSET=tf_base | tf_sig1 | tf_sig2 (D83): the comparison target lists (main 49_signaling_regulons.R); outputs drug_viper_specific_v2_<REGSET>
REGSET <- Sys.getenv("REGSET", "full"); stopifnot(REGSET %in% c("full", "tf", "tf_base", "tf_sig1", "tf_sig2"), REGSET == "full" || SFX == "_v2")
TFX <- if (REGSET == "full") "" else paste0("_", REGSET)
REGULON <- Sys.getenv("REGULON", "aracne_project"); stopifnot(REGULON %in% c("aracne_project", "collectri", "dorothea_ABC"))
TAG <- if (REGULON == "aracne_project") paste0(SFX, TFX) else paste0("_", REGULON)  # curated: not versioned
ONLY_NETS <- Filter(nzchar, strsplit(Sys.getenv("ONLY_NETS", ""), ",")[[1]])
# RUN_TAG (D85, 2026-11-02): appended to the output folder (e.g. drug_viper_raw_v2_tf_d85); the regulator table then goes
#   into that folder instead of results/lincs, so production outputs are not touched
RUN_TAG <- Sys.getenv("RUN_TAG", "")
proc <- "data/processed/lincs"; out_dir <- file.path(proc, paste0("drug_viper_", PROFILE, TAG, RUN_TAG))
dir.create(out_dir, showWarnings = FALSE)

genes <- fread(file.path(proc, "drug_line_profiles_genes.tsv"))
idx <- fread(file.path(proc, "drug_line_profiles_index.tsv"))
pfile <- if (PROFILE == "raw") "drug_line_profiles_bing.f32" else "drug_line_profiles_bing_specific.f32"
# PROFILE_FILE (D90, 2026-11-03): another profile file with the 79c layout (same index and genes), path from the repo root
if (nzchar(Sys.getenv("PROFILE_FILE"))) { pfile <- normalizePath(Sys.getenv("PROFILE_FILE")); proc_p <- "" } else proc_p <- proc
P <- matrix(readBin(if (nzchar(proc_p)) file.path(proc_p, pfile) else pfile, "numeric", size = 4,
                    n = nrow(genes) * nrow(idx)), nrow = nrow(genes))
dimnames(P) <- list(genes$gene_id, paste(idx$pert_id, idx$ModelID, sep = "|"))
cat(sprintf("profiles (%s): %d genes x %d drug-line profiles\n", PROFILE, nrow(P), ncol(P)))

if (REGULON == "aracne_project") {
  arac <- readRDS(paste0("data/processed/step5/regulons_aracne_project", SFX, TFX, ".rds"))
  grp <- fread("results/step1/disease_groups_toil.tsv")
  net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
  groups <- fread(paste0("results/step4", SFX, "/signature_summary.tsv"))[status == "ok", unique(group_id)]
  filt <- fread(paste0("results/step6/aracne_project", SFX, TFX, "/regulator_filter.tsv"))[keep == TRUE & group_id %in% groups]
  # EXTRA_RDS / NET_MAP / MR_OVERRIDE (D92, 2026-11-03): extra networks, group -> network map, and the mapped groups'
  #   regulator filter from another Step 6 folder. Defaults: off.
  if (nzchar(Sys.getenv("EXTRA_RDS"))) { ex <- readRDS(Sys.getenv("EXTRA_RDS")); arac <- c(arac, ex[setdiff(names(ex), names(arac))]) }
  # D97 (2026-11-03): production network map config/group_networks.tsv is the default for the v2 TF-only set; NET_MAP overrides
  nmf <- Sys.getenv("NET_MAP", if (SFX == "_v2" && REGSET == "tf" && file.exists("config/group_networks.tsv")) "config/group_networks.tsv" else "")
  if (nzchar(nmf)) {
    nm <- fread(nmf); net_of[nm$group_id] <- nm$net
    if (nzchar(Sys.getenv("MR_OVERRIDE"))) filt <- rbind(filt[!group_id %in% nm$group_id],
      fread(file.path(Sys.getenv("MR_OVERRIDE"), "regulator_filter.tsv"))[keep == TRUE & group_id %in% nm$group_id], fill = TRUE)
  }
  filt[, net := net_of[group_id]]
  # D119: subgroup regulators (Step 6 strata filter) are candidates in their parent's network too
  sff <- paste0("results/step6/aracne_project", SFX, TFX, "_strata/regulator_filter.tsv")
  if (file.exists(sff) && file.exists(paste0("results/step10", SFX, "/strata.tsv"))) {
    st <- unique(fread(paste0("results/step10", SFX, "/strata.tsv"))[, .(stratum, group_id)])
    sf <- fread(sff)[keep == TRUE][st, on = .(group_id = stratum), nomatch = 0]
    sf[, net := net_of[i.group_id]]
    n0 <- nrow(unique(filt[, .(net, regulator)]))
    filt <- rbind(filt, sf[, .(regulator, group_id, net)], fill = TRUE)
    cat(sprintf("D119: subgroup filter adds %d network x regulator candidates\n", nrow(unique(filt[, .(net, regulator)])) - n0))
  }
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
if (!length(ONLY_NETS)) fwrite(rbindlist(reg_tab), if (nzchar(RUN_TAG)) file.path(out_dir, "regulators.tsv") else
                                 paste0("results/lincs/drug_viper_regulators", TAG, ".tsv"), sep = "\t")
