#!/usr/bin/env Rscript
# Gate F1 — the spec's adversarial F-MATH fixtures, run against the INSTALLED CARNIVAL 2.20.0.
# These test BEHAVIOURAL claims ("reports non-identifiability", "never invent a path"), which can
# only be confirmed empirically. A claim the package does not honour is a spec defect.
# Output: results/systems_validation/synthetic/gate_f1_adversarial.json
suppressPackageStartupMessages({library(CARNIVAL); library(jsonlite)})
setwd("/home/kirti/projects/repurposing")
OUT <- "results/systems_validation/synthetic"; dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
TMP <- file.path(tempdir(), "f1"); dir.create(TMP, showWarnings=FALSE)
o <- defaultLpSolveCarnivalOptions(); o$outputFolder <- TMP; o$workdir <- TMP; o$timelimit <- 120
net <- function(...) { d <- do.call(rbind, list(...)); data.frame(source=d[,1],
          interaction=as.numeric(d[,2]), target=d[,3], stringsAsFactors=FALSE) }
res <- list()
run <- function(lab, f, ...) {
  r <- tryCatch(f(...), error=function(e) structure(list(msg=conditionMessage(e)), class="err"),
                warning=function(w) structure(list(msg=conditionMessage(w)), class="warn"))
  if (inherits(r,"err")) { cat(sprintf("  %-46s ERROR: %s\n", lab, substr(r$msg,1,70)))
                           return(list(status="error", message=r$msg)) }
  na <- r$nodesAttributes; sif <- r$weightedSIF
  cat(sprintf("  %-46s nodes=%-3s edges=%-3s\n", lab,
              if(is.null(na)) "0" else nrow(na), if(is.null(sif)) "0" else nrow(sif)))
  list(status="ok", n_nodes=if(is.null(na)) 0 else nrow(na), n_edges=if(is.null(sif)) 0 else nrow(sif),
       nodes=if(is.null(na)) NULL else as.data.frame(na), sif=if(is.null(sif)) NULL else as.data.frame(sif))
}
sgn <- function(r, node) { if (is.null(r$nodes)) return(NA)
  i <- match(node, r$nodes$Node); if (is.na(i)) return(NA); sign(as.numeric(r$nodes$AvgAct[i])) }

cat("F1-a: four-node signed chain X->A(+), A->B(-), B->C(+); activate X\n")
chain <- net(c("X",1,"A"), c("A",-1,"B"), c("B",1,"C"))
res$chain_pos <- run("perturb X=+1, measure A=+1,B=-1,C=-1", runVanillaCarnival,
  perturbations=c(X=1), measurements=c(A=1,B=-1,C=-1), priorKnowledgeNetwork=chain, carnivalOptions=o)
res$chain_neg <- run("perturb X=-1, measure A=-1,B=+1,C=+1 (reversed)", runVanillaCarnival,
  perturbations=c(X=-1), measurements=c(A=-1,B=1,C=1), priorKnowledgeNetwork=chain, carnivalOptions=o)
cat(sprintf("    signs under X=+1: A=%s B=%s C=%s  (spec: +, -, -)\n",
            sgn(res$chain_pos,"A"), sgn(res$chain_pos,"B"), sgn(res$chain_pos,"C")))
cat(sprintf("    signs under X=-1: A=%s B=%s C=%s  (spec: reversed)\n",
            sgn(res$chain_neg,"A"), sgn(res$chain_neg,"B"), sgn(res$chain_neg,"C")))

cat("\nF1-b: two sources with indistinguishable measured effects -> must NOT invent a unique driver\n")
two <- net(c("S1",1,"H"), c("S2",1,"H"), c("H",1,"T1"), c("H",1,"T2"))
res$two_sources <- run("inverse, measure T1=+1,T2=+1", runInverseCarnival,
  measurements=c(T1=1,T2=1), priorKnowledgeNetwork=two, carnivalOptions=o)
if (!is.null(res$two_sources$nodes)) {
  nd <- res$two_sources$nodes
  s <- nd[nd$Node %in% c("S1","S2"), c("Node","AvgAct")]
  cat("    source nodes returned:\n"); print(s, row.names=FALSE)
  res$two_sources$both_sources_present <- all(c("S1","S2") %in% nd$Node)
  res$two_sources$exactly_one_active <- sum(abs(as.numeric(s$AvgAct)) > 0) == 1
}

cat("\nF1-c: contradictory observations -> report disagreement/infeasibility, never invent a path\n")
res$contradiction <- run("chain, measure A=+1 AND B=+1 (sign-impossible)", runVanillaCarnival,
  perturbations=c(X=1), measurements=c(A=1,B=1), priorKnowledgeNetwork=chain, carnivalOptions=o)
if (!is.null(res$contradiction$nodes)) {
  cat(sprintf("    returned A=%s B=%s - an inconsistent pair cannot BOTH be satisfied through A-|B\n",
              sgn(res$contradiction,"A"), sgn(res$contradiction,"B")))
}

cat("\nF1-d: measured TF disconnected from the network -> must not fabricate reachability\n")
res$disconnected <- run("measure Z, absent from the PKN", runVanillaCarnival,
  perturbations=c(X=1), measurements=c(A=1,Z=1), priorKnowledgeNetwork=chain, carnivalOptions=o)
if (!is.null(res$disconnected$nodes)) {
  cat(sprintf("    Z present in solution: %s\n", "Z" %in% res$disconnected$nodes$Node))
  res$disconnected$Z_in_solution <- "Z" %in% res$disconnected$nodes$Node
}

write(toJSON(list(gate="F1 (adversarial F-MATH fixtures)",
                  carnival=as.character(packageVersion("CARNIVAL")), results=res),
             auto_unbox=TRUE, pretty=TRUE, null="null"),
      file.path(OUT,"gate_f1_adversarial.json"))
cat("\nwrote", file.path(OUT,"gate_f1_adversarial.json"), "\n")
