#!/usr/bin/env Rscript
# Per-target breakdown of every disease master-regulator NES (user 2026-11-22, D179): the nested array
# "regulator NES -> the regulon's targets, with their expression change and how much each adds to the NES".
#
# 50_viper.R keeps only regulator, nes, p, size. This step does not re-score anything. It rebuilds the exact
# msviper / aREA input (same signature, same regulator filter, same target pruning, same gene universe) and splits
# the enrichment score (ES) into one additive term per target. Read-only on every production file.
#
# aREA (viper 1.44, method "loop"), for one regulator with targets j (mode m_j in [-1, 1], likelihood l_j):
#   genes ranked over the universe U (msviper ges.filter: signature genes that are a target of ANY regulator passed in)
#   t2_j = qnorm(rank_j / (|U| + 1))                                  two-tailed (signed) score
#   t1_j = qnorm(|2 * rank_j/(|U|+1) - 1| + (1 - max)/2)             one-tailed (|change|) score
#   a_j = m_j l_j t2_j ;  b_j = (1 - |m_j|) l_j t1_j ;  S1 = sum a ;  S2 = sum b ;  s = sign(S1) (0 -> +1)
#   ES = s (|S1| + S2 [S2 > 0]) / sum l
# so ES = sum_j c_j exactly, with  c_j = (a_j + s b_j [S2 > 0]) / sum l.
# The production NES is aecdf(ES) against the permutation null (msviper nullmodel): a monotone, regulator-specific
# rescaling that the null permutations do not let us replay exactly. contrib_nes = c_j * NES / ES allocates the NES
# in proportion (it sums to the NES exactly; it is an allocation, not a separate measurement). share = c_j / ES.
#
# Checks (written to <out>/regulon_contributions_check.tsv, the run stops if (1) or (2) fails):
#   (1) sum_j c_j == ES from viper:::aREA on the same input (|diff| < 1e-9)
#   (2) the number of targets per regulator == `size` in group_mr.tsv
#   (3) Spearman(ES * w, stored NES) per unit, w = sqrt(sum (l/max l)^2) (aREA's analytic NES): reported, not gated,
#       because the stored NES used the permutation null.
#
# Env: VERSION=v2 REGSET=tf (production, D160) | STRATA=1 (Step 10 subgroups) | N_CORES (default 6) | ONLY_GROUPS
#      EDGES=1 also writes every regulon edge (net, regulator, target, tfmode, likelihood) for the drug side (79e2).
# Outputs: results/step6/aracne_project<SFX><TFX>[_strata]/regulon_contributions.tsv.gz
#            group_id regulator nes es target target_symbol tfmode likelihood z t2 contrib_es contrib_nes share
#          .../regulon_contributions_check.tsv
#          data/processed/step5/regulons_aracne_project<SFX><TFX>_edges.tsv.gz (EDGES=1)
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
root <- normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), ".."))
setwd(root)

SFX <- if (identical(Sys.getenv("VERSION", "v2"), "v2")) "_v2" else ""
REGSET <- Sys.getenv("REGSET", "tf"); TFX <- if (REGSET == "full") "" else paste0("_", REGSET)
ST <- identical(Sys.getenv("STRATA"), "1")
N_CORES <- as.integer(Sys.getenv("N_CORES", "6"))
tag <- paste0("aracne_project", SFX, TFX)
mr_dir <- file.path("results/step6", paste0(tag, if (ST) "_strata"))

arac <- readRDS(paste0("data/processed/step5/regulons_", tag, ".rds"))
if (identical(Sys.getenv("EDGES"), "1")) {
  ed <- rbindlist(lapply(names(arac), function(nt) rbindlist(lapply(names(arac[[nt]]), function(r) {
    x <- arac[[nt]][[r]]; data.table(net = nt, regulator = r, target = names(x$tfmode), tfmode = x$tfmode, likelihood = x$likelihood) }))))
  f <- paste0("data/processed/step5/regulons_", tag, "_edges.tsv.gz"); fwrite(ed, f, sep = "\t")
  cat(sprintf("edges: %d rows over %d networks -> %s\n", nrow(ed), length(arac), f))
}

grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
nmf <- if (SFX == "_v2" && REGSET == "tf" && file.exists("config/group_networks.tsv")) "config/group_networks.tsv" else ""
if (nzchar(nmf)) { nm <- fread(nmf); net_of[nm$group_id] <- nm$net }   # same map as 50_viper.R (D97)
sig_dir <- if (ST) "results/step10_v2" else paste0("results/step4", SFX)
summ <- fread(file.path(sig_dir, "signature_summary.tsv"))[status == "ok"]
contrast_for <- summ[, .(contrast = if ("gtex" %in% contrast) "gtex" else "adjacent"), by = group_id]
zmat <- list(gtex = fread(file.path(sig_dir, "z_tumor_vs_gtex.tsv")), adjacent = fread(file.path(sig_dir, "z_tumor_vs_adjacent.tsv")))
parent_of <- if (ST) { sp <- unique(fread(file.path(sig_dir, "strata.tsv"))[, .(stratum, group_id)]); setNames(sp$group_id, sp$stratum) } else NULL
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc == TRUE]
ens2sym <- setNames(idmap$symbol, idmap$gene_id)
filt <- fread(file.path(mr_dir, "regulator_filter.tsv"))[keep == TRUE]
mr <- fread(file.path(mr_dir, "group_mr.tsv"))

units <- intersect(contrast_for$group_id, unique(mr$group_id))
only <- Filter(nzchar, strsplit(Sys.getenv("ONLY_GROUPS"), ",")[[1]]); if (length(only)) units <- intersect(units, only)

run_unit <- function(gid) {
  par <- if (ST) parent_of[[gid]] else gid
  regs_all <- arac[[net_of[[par]]]]
  contrast <- contrast_for[group_id == gid, contrast]
  z <- zmat[[contrast]][, c("gene_id", gid), with = FALSE][!is.na(get(gid))]
  sig <- setNames(z[[gid]], z$gene_id)
  reg <- lapply(regs_all[filt[group_id == gid, regulator]], function(r) {          # as 50_viper.R run_group()
    k <- names(r$tfmode) %in% names(sig); list(tfmode = r$tfmode[k], likelihood = r$likelihood[k]) })
  U <- names(sig)[names(sig) %in% unique(c(names(reg), unlist(lapply(reg, function(x) names(x$tfmode)), use.names = FALSE)))]
  g <- sig[U]                                                                        # msviper ges.filter
  t2 <- rank(g) / (length(g) + 1); t1 <- abs(t2 - 0.5) * 2; t1 <- t1 + (1 - max(t1)) / 2
  t1 <- qnorm(t1); t2 <- qnorm(t2)
  m <- mr[group_id == gid]
  out <- rbindlist(lapply(m$regulator, function(r) {
    x <- reg[[r]]; tg <- names(x$tfmode); mo <- unname(x$tfmode); l <- x$likelihood
    a <- mo * l * t2[tg]; b <- (1 - abs(mo)) * l * t1[tg]
    S1 <- sum(a); S2 <- sum(b); s <- if (S1 == 0) 1 else sign(S1)
    c_ <- (a + s * b * (S2 > 0)) / sum(l)
    data.table(regulator = r, target = tg, tfmode = mo, likelihood = l, z = unname(g[tg]), t2 = unname(t2[tg]),
               contrib_es = unname(c_), w_nes = sqrt(sum((l / max(l))^2)))
  }))
  es_ref <- viper:::aREA(matrix(g, ncol = 1, dimnames = list(U, gid)), reg[m$regulator], method = "loop", minsize = 0)$es[, 1]
  out[, es := sum(contrib_es), by = regulator]
  chk <- out[, .(es = es[1], n = .N, w = w_nes[1]), by = regulator][m, on = "regulator"]
  chk[, es_ref := es_ref[regulator]]
  out <- m[, .(regulator, nes)][out, on = "regulator"]
  out[, `:=`(group_id = gid, contrib_nes = contrib_es * nes / es, share = contrib_es / es,
             target_symbol = unname(ens2sym[target]), w_nes = NULL)]
  setcolorder(out, c("group_id", "regulator", "nes", "es", "target", "target_symbol", "tfmode", "likelihood", "z", "t2",
                     "contrib_es", "contrib_nes", "share"))
  list(out = out, check = data.table(group_id = gid, contrast = contrast, net = net_of[[par]], n_universe = length(U),
       n_regulators = nrow(m), max_abs_es_diff = max(abs(chk$es - chk$es_ref)), n_size_mismatch = sum(chk$n != chk$size),
       rho_es_nes = cor(chk$es * chk$w, chk$nes, method = "spearman"), sign_agree = mean(sign(chk$es) == sign(chk$nes))))
}

res <- mclapply(units, function(u) tryCatch(run_unit(u), error = function(e) list(err = data.table(group_id = u, error = conditionMessage(e)))),
                mc.cores = N_CORES, mc.preschedule = FALSE)
err <- rbindlist(lapply(res, `[[`, "err")); if (nrow(err)) { print(err); stop("units failed") }
chk <- rbindlist(lapply(res, `[[`, "check")); out <- rbindlist(lapply(res, `[[`, "out"))
sfx <- Sys.getenv("OUT_TAG")
fwrite(chk, file.path(mr_dir, paste0("regulon_contributions_check", sfx, ".tsv")), sep = "\t")
print(chk[, .(units = .N, regulators = sum(n_regulators), max_es_diff = max(max_abs_es_diff), size_mismatch = sum(n_size_mismatch),
              min_rho = round(min(rho_es_nes), 3), median_rho = round(median(rho_es_nes), 3), min_sign_agree = round(min(sign_agree), 3))])
stopifnot(max(chk$max_abs_es_diff) < 1e-9, sum(chk$n_size_mismatch) == 0)
num <- c("nes", "es", "tfmode", "likelihood", "z", "t2", "contrib_es", "contrib_nes", "share")
out[, (num) := lapply(.SD, signif, 6), .SDcols = num]
f <- file.path(mr_dir, paste0("regulon_contributions", sfx, ".tsv.gz")); fwrite(out, f, sep = "\t")
cat(sprintf("%d rows (%d units x regulator x target) -> %s\n", nrow(out), length(units), f))
