#!/usr/bin/env Rscript
# D88 (rule fixed 2026-11-03 before running, docs/DECISIONS.md): score the time / dose arms (scripts/timedose/td_profiles.py)
# with tumor ARACNe v2 TF (regulators and pruning as 79d) and CollecTRI (reported), then apply the rule.
#   (1) ligand profiles (D84 list): lineage-matched networks else all holding the TF (as 79n / 79v), share in the
#       expected direction, paired by drug x TF x line; sign test over discordant rows
#   (2) growth-arrest coupling: per network, median over regulators of |r(activity, proliferation score)|; median over networks
# Outputs: results/lincs/timedose/{ligands.tsv, coupling.tsv, summary.tsv, decision.tsv, esr1_mcf7.tsv}
suppressPackageStartupMessages({library(data.table); library(viper); library(parallel)})
args <- commandArgs(trailingOnly = FALSE)
setwd(normalizePath(file.path(dirname(sub("--file=", "", args[grep("--file=", args)])), "../..")))
MIN_TARGETS <- 10; MIN_FRAC <- 0.5; GAIN <- 0.05; P_MAX <- 0.05
N_CORES <- min(as.integer(Sys.getenv("N_CORES", "20")), 24L)
proc <- "data/processed/lincs"; td <- file.path(proc, "timedose"); res <- "results/lincs/timedose"; dir.create(res, FALSE, TRUE)
genes <- fread(file.path(proc, "drug_line_profiles_genes.tsv"))$gene_id
COMP <- list(time = c("A0", "T6", "T24"), dose = c("A0", "DL", "DH")); VERS <- c("specific", "raw")

## regulons: ARACNe v2 TF with the 79d regulator filter, CollecTRI
arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
grp <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD", "TCGA-READ"), "COADREAD", sub("^[^-]+-", "", grp$cohort)), grp$group_id)
groups <- fread("results/step4_v2/signature_summary.tsv")[status == "ok", unique(group_id)]
filt <- fread("results/step6/aracne_project_v2_tf/regulator_filter.tsv")[keep == TRUE & group_id %in% groups][, net := net_of[group_id]]
prune <- function(rg, cand) { r <- lapply(rg[intersect(cand, names(rg))], function(x) { k <- names(x$tfmode) %in% genes
  if (sum(k) < MIN_TARGETS || mean(k) < MIN_FRAC) return(NULL); list(tfmode = x$tfmode[k], likelihood = x$likelihood[k]) })
  r <- Filter(Negate(is.null), r); class(r) <- "regulon"; r }
REG <- lapply(setNames(sort(unique(filt$net)), sort(unique(filt$net))), function(n) prune(arac[[n]], unique(filt[net == n, regulator])))
cl <- readRDS("results/step5/regulons_core.rds")$collectri; REG$CollecTRI_ALL <- prune(cl, names(cl))

## lineage matching (as 79n)
keys <- fread(cmd = "grep -v '^#' config/celligner_tissue_keys.tsv")
ln <- unique(keys[group_id %in% names(net_of), .(lineage = unlist(strsplit(lineages, ";")), net = net_of[group_id]), by = group_id][, .(lineage, net)])
lin <- setNames(fread("data/raw/depmap_24q4/Model.csv")[, OncotreeLineage], fread("data/raw/depmap_24q4/Model.csv")[, ModelID])

## ligands
lig <- unique(fread("results/lincs/d84_ligands.tsv")[profile == "specific", .(drug, tf, direction, ModelID)])
cmp <- fread("results/lincs/compounds.tsv")[, .(drug = cmap_name, pert_id)]
lig <- merge(lig, cmp, by = "drug", allow.cartesian = TRUE)

jobs <- CJ(comp = names(COMP), arm = c("A0", "T6", "T24", "DL", "DH"), ver = VERS, net = names(REG))[mapply(function(c, a) a %in% COMP[[c]], comp, arm)]
cache <- file.path(td, "scored.rds")
if (file.exists(cache)) out <- readRDS(cache) else {
message(sprintf("%d scoring jobs on %d cores", nrow(jobs), N_CORES))
out <- mclapply(seq_len(nrow(jobs)), function(j) {
  J <- jobs[j]; I <- fread(file.path(td, paste0(J$comp, "_index.tsv")))
  M <- matrix(readBin(file.path(td, sprintf("%s_%s_%s.f32", J$comp, J$arm, J$ver)), "numeric", size = 4, n = length(genes) * nrow(I)), nrow = length(genes))
  rownames(M) <- genes; colnames(M) <- paste(I$pert_id, I$ModelID, sep = "|")
  pr <- fread(file.path(td, paste0(J$comp, "_prolif.tsv")))[[paste0(J$arm, "_", J$ver)]]
  v <- viper(M, REG[[J$net]], method = "none", minsize = MIN_TARGETS, verbose = FALSE)
  cp <- data.table(J, regulator = rownames(v), abs_r = abs(as.vector(cor(t(v), pr))))
  L <- merge(lig, I[, .(pert_id, ModelID, col = .I)], by = c("pert_id", "ModelID"))[tf %in% rownames(v)]
  L <- if (nrow(L)) cbind(J, L[, .(drug, tf, direction, ModelID, nes = v[cbind(match(tf, rownames(v)), col)])]) else NULL
  list(cp = cp, lig = L)
}, mc.cores = N_CORES, mc.preschedule = FALSE)
bad <- vapply(out, inherits, logical(1), "try-error"); if (any(bad)) stop(out[[which(bad)[1]]])
saveRDS(out, cache) }
CP <- rbindlist(lapply(out, `[[`, "cp")); LG <- rbindlist(lapply(out, `[[`, "lig"))
LG[, set := fifelse(net == "CollecTRI_ALL", "collectri", "aracne")]
LG[set == "aracne", matched := paste(lin[ModelID], net) %in% ln[, paste(lineage, net)]]
LG <- rbind(LG[set == "aracne", if (any(matched)) .(nes = mean(nes[matched]), matched = TRUE) else .(nes = mean(nes), matched = FALSE),
               by = .(comp, arm, ver, set, drug, tf, direction, ModelID)],
            LG[set == "collectri", .(comp, arm, ver, set, drug, tf, direction, ModelID, nes, matched = NA)])
LG <- LG[, .(nes = mean(nes), matched = matched[1]), by = .(comp, arm, ver, set, drug, tf, direction, ModelID)]   # several IDs per drug
w <- dcast(LG, comp + ver + drug + tf + direction + ModelID ~ set + arm, value.var = "nes")
w[, consensus_A0 := rowMeans(cbind(aracne_A0, collectri_A0), na.rm = TRUE)]
for (a in c("T6", "T24", "DL", "DH")) if (paste0("aracne_", a) %in% names(w)) w[[paste0("consensus_", a)]] <- rowMeans(cbind(w[[paste0("aracne_", a)]], w[[paste0("collectri_", a)]]), na.rm = TRUE)
fwrite(w, file.path(res, "ligands.tsv"), sep = "\t")
CPs <- CP[, .(med = median(abs_r, na.rm = TRUE)), by = .(comp, arm, ver, set = fifelse(net == "CollecTRI_ALL", "collectri", "aracne"), net)][
  , .(coupling = round(median(med), 3), nets = .N), by = .(comp, arm, ver, set)]
fwrite(CPs, file.path(res, "coupling.tsv"), sep = "\t")

cmpare <- function(x, a, b, cmp_, ver_, set_) {
  x <- x[comp == cmp_ & ver == ver_]; ca <- paste0(set_, "_", a); cb <- paste0(set_, "_", b)
  x <- x[!is.na(get(ca)) & !is.na(get(cb))]
  ha <- sign(x[[ca]]) == x$direction; hb <- sign(x[[cb]]) == x$direction; wn <- sum(ha & !hb); ls <- sum(!ha & hb)
  ca2 <- CPs[comp == sub(" .*", "", cmp_) & ver == ver_ & set == fifelse(set_ == "consensus", "aracne", set_)]
  data.table(comp = cmp_, ver = ver_, set = set_, arm = a, vs = b, rows = nrow(x), tfs = uniqueN(x$tf), correct_arm = round(mean(ha), 4), correct_vs = round(mean(hb), 4),
             diff = round(mean(ha) - mean(hb), 4), wins = wn, losses = ls, p_sign = if (wn + ls) signif(binom.test(wn, wn + ls)$p.value, 3) else NA_real_,
             coupling_arm = ca2[arm == a, coupling], coupling_vs = ca2[arm == b, coupling])
}
SM <- rbindlist(lapply(VERS, function(v) rbindlist(lapply(c("aracne", "collectri", "consensus"), function(s) rbind(
  cmpare(w, "T6", "T24", "time", v, s), cmpare(w, "T6", "A0", "time", v, s), cmpare(w, "T24", "A0", "time", v, s),
  cmpare(w, "DL", "DH", "dose", v, s), cmpare(w, "DL", "A0", "dose", v, s), cmpare(w, "DH", "A0", "dose", v, s), fill = TRUE)))), fill = TRUE)
## reported: T6 on pairs active at 6 h
It <- fread(file.path(td, "time_index.tsv"))[tas_6h >= 0.4, .(pert_id, ModelID)]
act <- merge(w[comp == "time"], unique(merge(lig[, .(drug, ModelID, pert_id)], It, by = c("pert_id", "ModelID"))[, .(drug, ModelID)]), by = c("drug", "ModelID"))
SM <- rbind(SM, rbindlist(lapply(VERS, function(v) rbind(cmpare(act, "T6", "T24", "time", v, "aracne"), cmpare(act, "T6", "A0", "time", v, "aracne"))))[, comp := "time (active at 6 h)"], fill = TRUE)
fwrite(SM, file.path(res, "summary.tsv"), sep = "\t")
fwrite(w[tf == "ESR1" & ModelID == "ACH-000019"], file.path(res, "esr1_mcf7.tsv"), sep = "\t")

g <- function(a, b, cmp_) SM[comp == cmp_ & ver == "specific" & set == "aracne" & arm == a & vs == b]
t1 <- g("T6", "T24", "time"); t2 <- g("T6", "A0", "time"); d1 <- g("DL", "DH", "dose")
pass <- function(x) x$diff >= GAIN && isTRUE(x$p_sign < P_MAX) && x$wins > x$losses && x$coupling_arm <= x$coupling_vs
txt <- function(x) sprintf("%s vs %s: %.1f%% vs %.1f%% (%+.1f points, %d rows, %d wins / %d losses, p = %.3g); coupling %.3f vs %.3f",
                           x$arm, x$vs, 100 * x$correct_arm, 100 * x$correct_vs, 100 * x$diff, x$rows, x$wins, x$losses, x$p_sign, x$coupling_arm, x$coupling_vs)
dec <- data.table(question = c("time: T6 vs T24", "time: T6 vs A0", "dose: DL vs DH"), result = c(txt(t1), txt(t2), txt(d1)), met = c(pass(t1), pass(t2), pass(d1)))
dec <- rbind(dec, data.table(question = c("DECISION time", "DECISION dose"),
  result = c(if (pass(t1) && pass(t2)) "T6 goes forward to a ranking-level rule" else "time stays descriptive",
             if (pass(d1)) "DL goes forward to a ranking-level rule" else "dose stays descriptive"), met = c(pass(t1) && pass(t2), pass(d1))))
fwrite(dec, file.path(res, "decision.tsv"), sep = "\t")
options(width = 220); print(SM); print(CPs[order(comp, ver, set, arm)])
print(w[tf == "ESR1" & ModelID == "ACH-000019" & ver == "specific" & direction < 0, .(comp, drug, aracne_A0, aracne_T6, aracne_T24, aracne_DL, aracne_DH)])
print(dec)
