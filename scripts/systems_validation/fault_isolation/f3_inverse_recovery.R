#!/usr/bin/env Rscript
# F-3: inverse CARNIVAL driver recovery + alternative-optima enumeration, per the FROZEN
# registration in docs/FREEZE_F2_driver_recovery.md / f2_registration.json (commit 9cbc91e).
#
# INTERPRETER: /home/kirti/miniforge3/envs/sysval/bin/Rscript
# Usage: f3_inverse_recovery.R <cohort> [betaWeight]
#
# Reads NOTHING from F-4's held-out perturbations, no outcome label, no d128/a4.
suppressPackageStartupMessages({library(CARNIVAL); library(data.table); library(jsonlite)})
setwd("/home/kirti/projects/repurposing")
args <- commandArgs(trailingOnly = TRUE)
COHORT <- args[1]; BETA <- if (length(args) > 1) as.numeric(args[2]) else 0.2
PKNF <- if (length(args) > 2) args[3] else "results/systems_validation/pkn/signed_pkn.tsv"
TAG  <- if (length(args) > 3) args[4] else ""            # "" = curated arm (F-3), else F-3b
TLIM <- if (nzchar(TAG)) 1800 else 600                    # frozen in FREEZE_F3b section 6
REG <- fromJSON("results/systems_validation/fault_isolation/f2_registration.json")
stopifnot(COHORT %in% names(REG$pilots))
DRIVER <- REG$pilots[[COHORT]]$driver
OUT <- "results/systems_validation/fault_isolation"
D <- file.path("results/systems_validation/fault_isolation/work",
               paste0("f3", TAG, "_", COHORT, "_b", BETA))   # persistent: enumeration is slow, do not redo it
dir.create(D, recursive = TRUE, showWarnings = FALSE)
CBC <- "/home/kirti/miniforge3/envs/sysval/bin/cbc"

pkn <- fread(PKNF)[, .(source, interaction, target)]
tf <- fread(file.path(OUT, paste0("tf_activity_", COHORT, ".tsv")))
meas_df <- tf[in_pkn == TRUE & fdr < 0.05]                      # FROZEN selection rule
meas <- setNames(as.numeric(meas_df$nes), meas_df$regulator)    # FROZEN: continuous NES, unaltered
cat(sprintf("%s | driver %s | betaWeight %.2f | %d measured TFs | PKN %d edges\n",
            COHORT, DRIVER, BETA, length(meas), nrow(pkn)))

o <- defaultCbcSolveCarnivalOptions()
o$solverPath <- CBC; o$outputFolder <- D; o$workdir <- D
o$keepLPFiles <- TRUE; o$cleanTmpFiles <- FALSE
o$betaWeight <- BETA; o$timelimit <- TLIM; o$threads <- 1

t0 <- Sys.time()
gen <- tryCatch(generateLpFileCarnival(measurements = meas, priorKnowledgeNetwork = as.data.frame(pkn),
                                       carnivalOptions = o),
                error = function(e) structure(list(msg = conditionMessage(e)), class = "err"))
if (inherits(gen, "err")) {
  cat("LP GENERATION FAILED:", gen$msg, "\n")
  write(toJSON(list(cohort = COHORT, driver = DRIVER, betaWeight = BETA, status = "lp_generation_failed",
                    message = gen$msg), auto_unbox = TRUE, pretty = TRUE),
        file.path(OUT, sprintf("f3%s_%s_b%s.json", TAG, COHORT, BETA)))
  quit(status = 0)
}
gen_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
e <- new.env(); load(gen[["parsedDataFile"]], envir = e)
nodesDf <- get("variables", envir = e)$nodesDf
cat(sprintf("  LP generated in %.0fs: %d nodes, %d LP variables\n", gen_s, nrow(nodesDf), nrow(nodesDf) * 4))

# map LP binaries -> gene names, so enumeration output is interpretable
map <- rbind(data.table(var = nodesDf$nodesUpVars, node = nodesDf$nodes, dir = "up"),
             data.table(var = nodesDf$nodesDownVars, node = nodesDf$nodes, dir = "down"))
fwrite(map, file.path(D, "varmap.tsv"), sep = "\t")

enum_json <- file.path(D, "enum.json")
if (file.exists(enum_json)) {
  cat("  reusing cached enumeration:", enum_json, "\n")
  enum_s <- NA_real_
} else {
cmd <- c("scripts/systems_validation/fault_isolation/enumerate_optima.py",
         "--lp", gen[["lpFile"]], "--cbc", CBC, "--max-solutions", "25",
         "--seconds", as.character(TLIM), "--out", enum_json)
t1 <- Sys.time()
system2("/home/kirti/miniforge3/envs/sysval/bin/python", cmd, stdout = file.path(D, "enum.log"),
        stderr = file.path(D, "enum.err"))
enum_s <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
}

if (!file.exists(enum_json)) {
  cat("ENUMERATION FAILED; see", file.path(D, "enum.err"), "\n")
  cat(paste(utils::tail(readLines(file.path(D, "enum.err")), 5), collapse = "\n"), "\n")
  quit(status = 0)
}
en <- fromJSON(enum_json, simplifyVector = FALSE)
nopt <- en$n_alternative_optima
cat(sprintf("  enumeration: %d alternative optima in %.0fs (status %s)\n", nopt, enum_s, en$status))

# ---- recovery call, per the FROZEN definitions -------------------------------------------------
m <- setNames(map$node, map$var); mdir <- setNames(map$dir, map$var)
per_sol <- lapply(en$solutions, function(s) {
  a <- unlist(s$active)
  data.table(idx = s$index, node = unname(m[a]), dir = unname(mdir[a]))
})
driver_in <- vapply(per_sol, function(d) DRIVER %in% d$node, logical(1))
driver_dir <- vapply(per_sol, function(d) { x <- d[node == DRIVER]; if (nrow(x)) x$dir[1] else NA_character_ },
                     character(1))
expected <- REG$recovery$expected_sign[[DRIVER]]
sign_ok <- driver_dir == ifelse(expected == "+", "up", "down")
call <- if (en$status == "truncated_at_max_solutions") {
  "NOT EVALUABLE (cap hit, not exhausted)"
} else if (nopt > 0 && all(driver_in) && all(sign_ok, na.rm = TRUE)) {
  "STRONG"
} else if (any(driver_in) && any(sign_ok, na.rm = TRUE)) {
  "WEAK"
} else {
  "NOT RECOVERED"
}
core <- Reduce(intersect, lapply(per_sol, function(d) d$node))
uni <- Reduce(union, lapply(per_sol, function(d) d$node))
cat(sprintf("  driver %s in %d of %d optima; sign expected %s, observed %s\n", DRIVER,
            sum(driver_in), nopt, expected, paste(unique(na.omit(driver_dir)), collapse = "/")))
cat(sprintf("  RECOVERY CALL: %s\n", call))
cat(sprintf("  core nodes in every optimum: %d; union: %d\n", length(core), length(uni)))

res <- list(cohort = COHORT, driver = DRIVER, betaWeight = BETA, arm = ifelse(nzchar(TAG),"aracne_full","curated"), pkn_file = PKNF,
            status = REG$pilots[[COHORT]]$status,
            n_measured_tfs = length(meas), pkn_edges = nrow(pkn),
            lp_seconds = round(gen_s, 1), enum_seconds = round(enum_s, 1),
            n_alternative_optima = nopt, enum_status = en$status, optimum = en$optimum,
            driver_in_n_optima = sum(driver_in), driver_direction = unique(na.omit(driver_dir)),
            expected_sign = expected, recovery_call = call,
            n_core_nodes = length(core), n_union_nodes = length(uni),
            median_pairwise_jaccard = en$median_pairwise_jaccard,
            core_nodes = head(sort(core), 200),
            frozen_by = "docs/FREEZE_F2_driver_recovery.md",
            reads_no_outcome_labels = TRUE)
write(toJSON(res, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(OUT, sprintf("f3%s_%s_b%s.json", TAG, COHORT, BETA)))
cat("wrote", file.path(OUT, sprintf("f3%s_%s_b%s.json", TAG, COHORT, BETA)), "\n")
