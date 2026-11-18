#!/usr/bin/env Rscript
# SPEC_systems_theoretic_mechanism_validation.md, Gate F0 (F-MATH):
#   "the installed package reproduces its standard and inverse toy examples with auditable solutions."
#
# This runs the PUBLISHED CARNIVAL package's own bundled fixtures - it does not reimplement or
# relabel anything. Both modes are exercised, on both available open solvers, and the vanilla
# result is compared against the package's own expected-result fixture.
#
# Env: /home/kirti/miniforge3/envs/sysval (isolated; production envs untouched).
# Output: results/systems_validation/synthetic/gate_f0_carnival.json
suppressPackageStartupMessages({library(CARNIVAL); library(jsonlite)})
setwd("/home/kirti/projects/repurposing")
LIB <- system.file(package = "CARNIVAL")
OUT <- "results/systems_validation/synthetic"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
TMP <- file.path(tempdir(), "carnival_f0"); dir.create(TMP, showWarnings = FALSE)

ld <- function(f) { e <- new.env(); load(file.path(LIB, f), envir = e); get(ls(e)[1], envir = e) }
perturb  <- ld("toy_perturbations_ex1.RData")
meas1    <- ld("toy_measurements_ex1.RData")
net1     <- ld("toy_network_ex1.RData")
meas2    <- ld("toy_measurements_ex2.RData")
net2     <- ld("toy_network_ex2.RData")
net2_fix <- net2
expected <- ld("toy_result_expected_ex1.RData")

cat(sprintf("CARNIVAL %s | solvers: %s\n", packageVersion("CARNIVAL"),
            paste(unlist(getSupportedSolvers()), collapse = ", ")))
cat(sprintf("ex1 vanilla : %d perturbations, %d measurements, %d network edges\n",
            length(perturb), length(meas1), nrow(net1)))
cat(sprintf("ex2 inverse : %d measurements, %d network edges\n", length(meas2), nrow(net2)))

CBC <- "/home/kirti/miniforge3/envs/sysval/bin/cbc"   # Sys.which() fails: the env bin is not on
                                                      # Rscript's PATH when invoked by absolute path
opts <- function(solver) {
  o <- if (solver == "cbc") defaultCbcSolveCarnivalOptions() else defaultLpSolveCarnivalOptions()
  o$outputFolder <- TMP; o$workdir <- TMP
  if (solver == "cbc") o$solverPath <- CBC
  o$timelimit <- 180
  o
}

# FINDING (recorded, not worked around silently): the package's own ex2 (inverse) fixture ships with
# columns Source/Effect/Target, while the installed 2.20.0 validator requires source/interaction/target
# - the ex1 fixture already uses the new names. The fixture predates the rename. Renaming the columns
# is an input-format adaptation, NOT a change to CARNIVAL's maths; the edges and signs are untouched.
net2_raw_cols <- NULL
run <- function(label, f, ...) {
  t0 <- Sys.time()
  r <- tryCatch(f(...), error = function(e) structure(list(error = conditionMessage(e)), class = "err"))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "err")) { cat(sprintf("  %-28s FAILED: %s\n", label, r$error))
                            return(list(status = "error", message = r$error, seconds = dt)) }
  nn <- if (!is.null(r$nodesAttributes)) nrow(r$nodesAttributes) else NA
  ne <- if (!is.null(r$weightedSIF)) nrow(r$weightedSIF) else NA
  cat(sprintf("  %-28s OK  nodes=%s edges=%s  %.1fs\n", label, nn, ne, dt))
  list(status = "ok", n_nodes = nn, n_edges = ne, seconds = dt,
       nodes = if (!is.na(nn)) as.data.frame(r$nodesAttributes) else NULL,
       sif = if (!is.na(ne)) as.data.frame(r$weightedSIF) else NULL)
}

net2_raw_cols <- colnames(net2_fix)
if (!all(c("source","interaction","target") %in% colnames(net2_fix))) {
  colnames(net2_fix)[match(c("Source","Effect","Target"), colnames(net2_fix))] <- c("source","interaction","target")
  cat(sprintf("\nex2 fixture columns renamed %s -> %s (format only)\n",
              paste(net2_raw_cols, collapse="/"), paste(colnames(net2_fix), collapse="/")))
}
if (is.data.frame(meas2) && nrow(meas2) == 1) meas2 <- unlist(meas2[1, , drop = TRUE])

cat("\nGate F0 runs:\n")
res <- list()
for (s in c("lpSolve", "cbc")) {
  res[[paste0("vanilla_", s)]] <- run(paste("vanilla (standard) /", s), runVanillaCarnival,
      perturbations = perturb, measurements = meas1, priorKnowledgeNetwork = net1,
      carnivalOptions = opts(s))
  res[[paste0("inverse_", s)]] <- run(paste("inverse /", s), runInverseCarnival,
      measurements = meas2, priorKnowledgeNetwork = net2_fix, carnivalOptions = opts(s))
}

# --- comparison against the package's OWN expected result (the actual gate) -------------------
cmp <- list()
exp_nodes <- as.data.frame(expected$nodesAttributes)
for (k in c("vanilla_lpSolve", "vanilla_cbc")) {
  r <- res[[k]]
  if (r$status != "ok") { cmp[[k]] <- list(comparable = FALSE, reason = "run failed"); next }
  got <- r$nodes
  key <- intersect(names(got), names(exp_nodes))
  m <- merge(exp_nodes, got, by = "Node", suffixes = c("_exp", "_got"))
  same_sign <- if (all(c("AvgAct_exp", "AvgAct_got") %in% names(m)))
                 mean(sign(as.numeric(m$AvgAct_exp)) == sign(as.numeric(m$AvgAct_got))) else NA
  cmp[[k]] <- list(comparable = TRUE, expected_nodes = nrow(exp_nodes), got_nodes = nrow(got),
                   matched_nodes = nrow(m), shared_fields = key,
                   node_activity_sign_agreement = same_sign)
  cat(sprintf("\n%s vs package expected fixture: %d expected nodes, %d returned, %d matched, sign agreement %s\n",
              k, nrow(exp_nodes), nrow(got), nrow(m),
              if (is.na(same_sign)) "n/a" else sprintf("%.3f", same_sign)))
}

ok <- all(vapply(res, function(x) identical(x$status, "ok"), logical(1)))
manifest <- list(gate = "F0 (CARNIVAL published-toy replication)",
                 spec = "docs/SPEC_systems_theoretic_mechanism_validation.md",
                 carnival_version = as.character(packageVersion("CARNIVAL")),
                 r_version = R.version.string, env = "/home/kirti/miniforge3/envs/sysval",
                 cbc = CBC, ex2_fixture_columns_renamed = net2_raw_cols, lpSolve = as.character(packageVersion("lpSolve")),
                 api_note = paste("installed 2.20.0 exposes runVanillaCarnival/runInverseCarnival;",
                                  "the legacy runCARNIVAL(inputObj=,measObj=,netObj=) also exists.",
                                  "The spec's F-MATH section names only the legacy form."),
                 runs = lapply(res, function(x) x[c("status","n_nodes","n_edges","seconds","message")]),
                 expected_fixture_comparison = cmp,
                 all_runs_ok = ok)
write(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null"),
      file.path(OUT, "gate_f0_carnival.json"))
cat(sprintf("\nGATE F0: %s\n", if (ok) "all four runs completed" else "AT LEAST ONE RUN FAILED"))
cat("wrote", file.path(OUT, "gate_f0_carnival.json"), "\n")
