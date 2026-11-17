#!/usr/bin/env Rscript
# ============================================================================
# compareFONSE_PANSE.R -- compare a FONSE fit's estimates against PANSE (and ROC)
# on the S. cerevisiae S288c genome.
#
# Compares:
#   1. phi (synthesis rate):  FONSE vs PANSE,  FONSE vs ROC   (join by GeneID)
#   2. selection vs nonsense error rate: FONSE dOmega vs PANSE NSERate (by Codon)
#   3. mutation bias:         FONSE dM vs ROC dM             (join by Codon)
#
# Usage: compareFONSE_PANSE.R <fonse_fit_dir> [out_pdf]
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
fonse_dir <- if (length(args) >= 1) args[1] else stop("need fonse fit dir")
out_pdf   <- if (length(args) >= 2) args[2] else file.path(fonse_dir, "fonse_vs_panse.pdf")

# Repository root: this script lives in <root>/03_r_scripts/. Set YEAST_ROOT to override;
# when sourced interactively (no --file argument) the working directory is used.
script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
YR <- Sys.getenv("YEAST_ROOT", unset = if (length(script_path)) normalizePath(file.path(dirname(script_path[1]), "..")) else getwd())
panse_dir <- file.path(YR, "01_results/00_panse_fits",
  "2026-01-07_Weinberg_etal_2016_all_frames_filter_genes_shorter_than_225_codons_200_ramp_updated_posterior_variable_nse",
  "final_restart/Parameter_est")
roc_dir <- file.path(YR, "00_data/01_correlates/ROC_Estimates")

g1 <- function(dir, pat) { f <- list.files(dir, pattern = pat, full.names = TRUE); if(!length(f)) NA_character_ else f[1] }

# ---- phi ----
fonse_phi <- read.csv(file.path(fonse_dir, "Parameter_est/gene_expression.txt"), stringsAsFactors = FALSE)
panse_phi <- read.csv(file.path(panse_dir, "gene_expression.txt"), stringsAsFactors = FALSE)
roc_phi   <- read.csv(file.path(roc_dir, "saccharomyces_cerevisiae.max.cds_Phi.csv"), stringsAsFactors = FALSE)
idcol <- function(d) names(d)[1]
mphi  <- function(d) d[["Mean"]]

ph <- function(a, b, la, lb) {
  m <- merge(data.frame(GeneID=a[[idcol(a)]], A=mphi(a)),
             data.frame(GeneID=b[[idcol(b)]], B=mphi(b)), by="GeneID")
  m <- m[is.finite(m$A) & is.finite(m$B) & m$A>0 & m$B>0, ]
  list(m=m, n=nrow(m),
       pearson=cor(log10(m$A), log10(m$B)),
       spearman=cor(m$A, m$B, method="spearman"), la=la, lb=lb)
}
phi_fp <- ph(fonse_phi, panse_phi, "FONSE phi", "PANSE phi")
phi_fr <- ph(fonse_phi, roc_phi,   "FONSE phi", "ROC phi")

# ---- codon-specific ----
fonse_sel <- read.csv(g1(file.path(fonse_dir,"Parameter_est"), "Selection\\.csv$"), stringsAsFactors=FALSE)
fonse_mut <- read.csv(g1(file.path(fonse_dir,"Parameter_est"), "Mutation\\.csv$"),  stringsAsFactors=FALSE)
panse_nse <- read.csv(g1(panse_dir, "NSERate\\.csv$"), stringsAsFactors=FALSE)
roc_mut   <- read.csv(file.path(roc_dir, "saccharomyces_cerevisiae.max.cds_Mutation.csv"), stringsAsFactors=FALSE)

cc <- function(a, b, la, lb, log_b=FALSE) {
  m <- merge(data.frame(Codon=a$Codon, A=a$Mean), data.frame(Codon=b$Codon, B=b$Mean), by="Codon")
  m <- m[is.finite(m$A) & is.finite(m$B), ]
  list(m=m, n=nrow(m),
       pearson=cor(m$A, m$B), spearman=cor(m$A, m$B, method="spearman"),
       la=la, lb=lb, log_b=log_b)
}
sel_nse <- cc(fonse_sel, panse_nse, "FONSE dOmega", "PANSE NSERate")
mut_mut <- cc(fonse_mut, roc_mut,   "FONSE dM",     "ROC dM")

# ---- report ----
cat(sprintf("\n=== FONSE vs PANSE/ROC on S288c (FONSE fit: %s) ===\n", fonse_dir))
rep1 <- function(x) cat(sprintf("  %-14s vs %-14s  n=%-5d  Pearson(log)=%.3f  Spearman=%.3f\n",
                                 x$la, x$lb, x$n, x$pearson, x$spearman))
rep2 <- function(x) cat(sprintf("  %-14s vs %-14s  n=%-3d  Pearson=%.3f  Spearman=%.3f\n",
                                 x$la, x$lb, x$n, x$pearson, x$spearman))
cat("phi (per gene):\n");   rep1(phi_fp); rep1(phi_fr)
cat("codon-specific:\n");   rep2(sel_nse); rep2(mut_mut)

pdf(out_pdf, width=10, height=8); par(mfrow=c(2,2))
sp1 <- function(x){ plot(log10(x$m$A), log10(x$m$B), pch=20, col="#00000055",
   xlab=paste0("log10 ", x$la), ylab=paste0("log10 ", x$lb),
   main=sprintf("%s vs %s\nr=%.3f (Sp %.3f), n=%d", x$la, x$lb, x$pearson, x$spearman, x$n)); abline(0,1,col="red",lty=2)}
sp2 <- function(x){ plot(x$m$A, x$m$B, pch=20, col="#0000aaaa",
   xlab=x$la, ylab=x$lb, main=sprintf("%s vs %s\nPearson=%.3f, Spearman=%.3f, n=%d", x$la, x$lb, x$pearson, x$spearman, x$n)) }
sp1(phi_fp); sp1(phi_fr); sp2(sel_nse); sp2(mut_mut)
dev.off()
cat("plots ->", out_pdf, "\n")
