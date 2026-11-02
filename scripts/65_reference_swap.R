
#!/usr/bin/env Rscript
# Reference-swap test (written and run 2026-10-27 on the TF arm; header added 2026-10-29): how much do the group
# master regulators depend on which normal reference is used? For every group with both a GTEx and an adjacent-normal
# Step 4 v2 signature, msviper (TF-only v2 regulons, 50_viper.R filters) three ways:
#   A  GTEx signature,     null from GTEx normals      (the v2 setup)
#   B  adjacent signature, null from adjacent normals
#   C  adjacent signature, null from GTEx normals      (separates the signature's effect from the null's)
# Reported per group and variant: NES spread, FDR < 0.05 count and direction, Spearman of NES between A and B / A and C,
# and of the two gene signatures. Caveat: the GTEx null here uses all of the group's GTEx normals, as 50_viper.R did
# before the 2026-10-28 fix (D73); a rerun should use the Step 4a reference set.
# Outputs: results/step6/reference_swap/reference_swap_stats.tsv, reference_swap_nes.tsv.gz (not in git)
suppressPackageStartupMessages({library(data.table); library(edgeR); library(limma); library(viper); library(parallel)})
setwd("/home/kirti/projects/repurposing")
MIN_LOGCPM <- 1; MIN_TARGETS <- 10; MIN_FRAC <- 0.5; N_PERM <- 1000
N_CORES <- as.integer(Sys.getenv("N_CORES","12")); set.seed(1)
COMPLEXES <- list(AP1=c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2"), NFKB=c("NFKB1","NFKB2","RELA","RELB","REL"))

arac <- readRDS("data/processed/step5/regulons_aracne_project_v2_tf.rds")
grp  <- fread("results/step1/disease_groups_toil.tsv")
net_of <- setNames(ifelse(grp$cohort %in% c("TCGA-COAD","TCGA-READ"),"COADREAD", sub("^[^-]+-","",grp$cohort)), grp$group_id)
idmap <- fread("results/step5/id_map.tsv")[primary_for_hgnc==TRUE]; sym2ens <- setNames(idmap$gene_id, idmap$symbol)
summ <- fread("results/step4_v2/signature_summary.tsv")[status=="ok"]
w <- dcast(summ, group_id~contrast, value.var="n_normal"); GROUPS <- w[!is.na(gtex)&!is.na(adjacent), group_id]
zg <- fread("results/step4_v2/z_tumor_vs_gtex.tsv"); za <- fread("results/step4_v2/z_tumor_vs_adjacent.tsv")
cov <- fread("results/step3/covariates_toil.tsv")[primary_per_patient==TRUE]
counts <- readRDS("data/processed/step4/counts_protein_coding.rds")
cat("groups:", length(GROUPS), "\n")

perm_null <- function(nor, genes){
  p <- sapply(seq_len(N_PERM), function(i){
    lab <- sample(rep(c(TRUE,FALSE), length.out=ncol(nor)))
    a <- nor[,lab,drop=FALSE]; b <- nor[,!lab,drop=FALSE]
    t <- (rowMeans(a)-rowMeans(b))/sqrt(apply(a,1,var)/ncol(a)+apply(b,1,var)/ncol(b)+1e-6)
    zscoreT(pmin(pmax(t,-50),50), df=ncol(a)+ncol(b)-2) })
  rownames(p) <- genes; p }

run_one <- function(gid){
  tryCatch({
  g1 <- zg[!is.na(get(gid)), .(gene_id, z=get(gid))]; g2 <- za[!is.na(get(gid)), .(gene_id, z=get(gid))]
  genes <- intersect(intersect(g1$gene_id, g2$gene_id), rownames(counts))
  sig_g <- setNames(g1[match(genes,gene_id), z], genes); sig_a <- setNames(g2[match(genes,gene_id), z], genes)
  d <- cov[group_id==gid & role %in% c("tumor","gtex_normal","adjacent_normal") & sample %in% colnames(counts)]
  y <- calcNormFactors(DGEList(counts[genes, d$sample])); lc <- cpm(y, log=TRUE, prior.count=1)
  tum <- lc[, d[role=="tumor", sample], drop=FALSE]
  nor_g <- lc[, d[role=="gtex_normal", sample], drop=FALSE]; nor_a <- lc[, d[role=="adjacent_normal", sample], drop=FALSE]
  med <- apply(tum,1,median)
  expressed <- function(src){ if (src %in% names(COMPLEXES)){ gg<-sym2ens[COMPLEXES[[src]]]; return(any(med[gg[!is.na(gg)]]>=MIN_LOGCPM,na.rm=TRUE)) }
    gg <- if (startsWith(src,"ENSG")) src else sym2ens[src]; !is.na(gg) && !is.na(med[gg]) && med[gg]>=MIN_LOGCPM }
  regs_all <- arac[[net_of[[gid]]]]
  keep <- vapply(names(regs_all), function(src){ tg<-names(regs_all[[src]]$tfmode); n<-sum(tg %in% genes)
    expressed(src) && n>=MIN_TARGETS && n/length(tg)>=MIN_FRAC }, logical(1))
  reg <- lapply(regs_all[names(keep)[keep]], function(r){ k<-names(r$tfmode) %in% genes; list(tfmode=r$tfmode[k], likelihood=r$likelihood[k]) })
  class(reg) <- "regulon"
  ng <- perm_null(nor_g, genes); na_ <- perm_null(nor_a, genes)
  fit <- function(s,n){ m <- msviper(s, reg, nullmodel=n, minsize=MIN_TARGETS, verbose=FALSE)
    data.table(regulator=names(m$es$nes), nes=as.numeric(m$es$nes), p=as.numeric(m$es$p.value))[, fdr:=p.adjust(p,"BH")] }
  A <- fit(sig_g, ng); B <- fit(sig_a, na_); C <- fit(sig_a, ng)
  long <- rbindlist(list(cbind(variant="A_gtex_sig_gtex_null",A), cbind(variant="B_adj_sig_adj_null",B), cbind(variant="C_adj_sig_gtex_null",C)))
  long[, group_id:=gid]
  st <- long[, .(n_reg=.N, nes_median=round(median(nes),3), nes_mad=round(mad(nes),3), nes_mean=round(mean(nes),3),
                 frac_pos=round(mean(nes>0),3), n_fdr05=sum(fdr<0.05),
                 frac_up_fdr05=round(if (sum(fdr<0.05)>0) mean(nes[fdr<0.05]>0) else NA_real_,3)), by=.(group_id,variant)]
  m <- dcast(long, regulator~variant, value.var="nes")
  st[, rho_A_B := round(cor(m$A_gtex_sig_gtex_null, m$B_adj_sig_adj_null, method="spearman"),3)]
  st[, rho_A_C := round(cor(m$A_gtex_sig_gtex_null, m$C_adj_sig_gtex_null, method="spearman"),3)]
  st[, rho_sig  := round(cor(sig_g, sig_a, method="spearman"),3)]
  st[, n_normal_gtex := ncol(nor_g)][, n_normal_adj := ncol(nor_a)][, n_tumor := ncol(tum)][, n_genes := length(genes)]
  list(stats=st, long=long)
  }, error=function(e) list(stats=data.table(group_id=gid, variant="ERROR", err=conditionMessage(e)), long=NULL))
}
res <- mclapply(GROUPS, run_one, mc.cores=N_CORES)
stats <- rbindlist(lapply(res,`[[`,"stats"), fill=TRUE)
long  <- rbindlist(Filter(Negate(is.null), lapply(res,`[[`,"long")), fill=TRUE)
dir.create("results/step6/reference_swap", recursive=TRUE, showWarnings=FALSE)
fwrite(stats, "results/step6/reference_swap/reference_swap_stats.tsv", sep="\t")
fwrite(long,  "results/step6/reference_swap/reference_swap_nes.tsv.gz", sep="\t")
cat("ERRORS:", sum(stats$variant=="ERROR"), "\n"); print(stats[variant=="ERROR"])
print(stats[variant!="ERROR", .(group_id,variant,nes_median,nes_mad,n_fdr05,frac_up_fdr05,rho_A_B)][order(group_id,variant)])
cat("DONE\n")
