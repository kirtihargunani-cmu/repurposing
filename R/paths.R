# R/paths.R - single source of truth for build paths. Sourced after setwd(root).
# Override the active build with PIPELINE_BUILD (NOT BUILD: conda sets that).
`%||%` <- function(a, b) if (is.null(a) || !nzchar(a)) b else a

.build_cfg <- function() {
  f <- "config/builds.yaml"
  if (!file.exists(f)) stop("missing ", f, " (run from the project root)")
  yaml::read_yaml(f)
}

build <- function(name = NULL) {
  cfg <- .build_cfg()
  name <- name %||% Sys.getenv("PIPELINE_BUILD") %||% cfg$active_build
  b <- cfg$builds[[name]]
  if (is.null(b)) stop("unknown build '", name, "'; known: ", paste(names(cfg$builds), collapse = ", "))
  b$name <- name
  b
}

p_signatures       <- function(b) b$signatures
p_regulons         <- function(b) b$regulons
p_step6            <- function(b) b$step6
p_patient_nes      <- function(b) b$patient_nes
p_group_mr         <- function(b) file.path(b$step6, "group_mr.tsv")
p_regulator_filter <- function(b) file.path(b$step6, "regulator_filter.tsv")
p_step8            <- function(b) b$step8
p_step8_features   <- function(b) b$step8_features

# groups with an ok Step 4 signature for this build
p_groups <- function(b) {
  s <- data.table::fread(file.path(b$signatures, "signature_summary.tsv"))
  sort(unique(s[status == "ok", group_id]))
}

# fail loudly rather than silently resolving to stale data
require_build <- function(b, need = character()) {
  miss <- character()
  for (k in need) {
    p <- b[[k]]
    if (is.null(p) || !file.exists(p)) miss <- c(miss, sprintf("%s -> %s", k, if (is.null(p)) "<undeclared>" else p))
  }
  if (length(miss)) stop("build '", b$name, "' missing inputs:\n  ", paste(miss, collapse = "\n  "))
  message(sprintf("build: %s", b$name))
  invisible(TRUE)
}
