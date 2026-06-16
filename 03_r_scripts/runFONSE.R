#!/usr/bin/env Rscript
# ============================================================================
# Script:       runFONSE.R
# Version:      1.0.0
# Created:      2026-06-16
# Revised:      2026-06-16
# Purpose:      Fit the AnaCoDa FONSE model to a coding-sequence FASTA and
#               write phi (synthesis) + codon-specific (dM / dOmega) estimates,
#               for comparison against PANSE estimates on the same genome.
#
# Mirrors the group's rocAnalysis.R conventions (sphi init, stop-adapting at
# 50% of samples, getCSPEstimates / getExpressionEstimates outputs), but for
# model = "FONSE". a1 (initiation) and a2 (elongation) cost coefficients are
# held FIXED at 4 by default (the new AnaCoDa default); pass --est_a1 / --est_a2
# to estimate them.
#
# Usage:        runFONSE.R -i CDS.fasta -o OUTDIR [-s SAMPLES] [-t THIN]
#                          [-a ADAPT] [-n THREADS] [--est_a1] [--est_a2]
#                          [--lib /path/to/AnaCoDa/lib]
# ============================================================================

suppressMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("-i","--input", type="character", help="CDS FASTA")
parser$add_argument("-o","--output", type="character", help="Output directory")
parser$add_argument("-s","--samp", type="integer", default=2000, help="Number of samples")
parser$add_argument("-t","--thin", type="integer", default=10, help="Thinning")
parser$add_argument("-a","--adapt", type="integer", default=50, help="Adaptive width (in samples)")
parser$add_argument("-n","--threads", type="integer", default=1, help="MCMC threads")
parser$add_argument("--est_a1", action="store_true", help="Estimate initiation cost a1 (default: fixed at 4)")
parser$add_argument("--est_a2", action="store_true", help="Estimate elongation cost a2 (default: fixed at 4)")
parser$add_argument("--lib", type="character", default=NULL, help="AnaCoDa lib.loc (.libPaths prepend)")
parser$add_argument("--init_phi", type="character", default=NULL, help="Warm-start phi CSV (GeneID,phi)")
parser$add_argument("--init_dm", type="character", default=NULL, help="Warm-start dM CSV (AA,Codon,Mean)")
parser$add_argument("--init_domega", type="character", default=NULL, help="Warm-start dOmega CSV (AA,Codon,Mean)")
args <- parser$parse_args()
print(args)

if (!is.null(args$lib)) .libPaths(c(args$lib, .libPaths()))
suppressMessages(library(AnaCoDa))
cat("AnaCoDa:", as.character(packageVersion("AnaCoDa")), "from", find.package("AnaCoDa"), "\n")

input        <- args$input
directory    <- args$output
samples      <- args$samp
thinning     <- args$thin
adaptiveWidth<- args$adapt
num_threads  <- args$threads

genome <- initializeGenomeObject(file = input)
size   <- length(genome)
cat("genome size:", size, "genes\n")

numMixtures     <- 1
geneAssignment  <- rep(1, size)
sphi_init       <- rep(1.5, numMixtures)   # yeast-scale default; refined by data via adaptation
percent.to.keep <- 0.5                     # stop adapting at 50%, estimate from last 50%

# Warm-start phi (aligned to genome gene order) if supplied.
init_phi <- NULL
if (!is.null(args$init_phi)) {
  ptab <- read.csv(args$init_phi, stringsAsFactors = FALSE)
  init_phi <- setNames(ptab[[2]], ptab[[1]])[getNames(genome)]
  sphi_init <- rep(sd(log(init_phi)), numMixtures)
  cat("warm-start phi: loaded", sum(!is.na(init_phi)), "genes; sphi_init =", round(sphi_init,3), "\n")
}

parameter <- initializeParameterObject(genome, sphi_init, numMixtures, geneAssignment,
                                        split.serine = TRUE,
                                        mixture.definition = "selectionShared",
                                        model = "FONSE",
                                        initial.expression.values = init_phi,
                                        init.initiation.cost = 4,
                                        init.elongation.cost = 4)
# Warm-start codon-specific parameters (dM, dOmega) if supplied; fix=FALSE so they are estimated.
if (!is.null(args$init_dm))     { parameter$initMutationCategories(args$init_dm, 1, FALSE);  cat("warm-start dM loaded\n") }
if (!is.null(args$init_domega)) { parameter$initSelectionCategories(args$init_domega, 1, FALSE); cat("warm-start dOmega loaded\n") }
if (args$est_a1) { parameter$estimateInitiationCost(); cat("a1: ESTIMATED\n") } else cat("a1: fixed at 4\n")
if (args$est_a2) { parameter$estimateElongationCost(); cat("a2: ESTIMATED\n") } else cat("a2: fixed at 4\n")

steps.to.adapt <- (samples * thinning) * (1 - percent.to.keep)

dir.create(directory, recursive = TRUE, showWarnings = FALSE)
for (sub in c("Graphs","Restart_files","Parameter_est","R_objects"))
  dir.create(file.path(directory, sub), showWarnings = FALSE)

mcmc <- initializeMCMCObject(samples = samples, thinning = thinning,
                             adaptive.width = adaptiveWidth,
                             est.expression = TRUE, est.csp = TRUE, est.hyper = TRUE)
mcmc$setStepsToAdapt(steps.to.adapt)

model <- initializeModelObject(parameter, "FONSE")
setRestartSettings(mcmc, file.path(directory, "Restart_files/rstartFile.rst"), adaptiveWidth, FALSE)

cat("=== starting FONSE MCMC:", samples, "samples x thin", thinning, "=",
    samples*thinning, "iterations on", num_threads, "threads ===\n")
sys.runtime <- system.time(runMCMC(mcmc, genome, model, num_threads, 0))
cat("=== MCMC done in", round(sys.runtime["elapsed"]/60, 1), "min ===\n")
write.table(data.frame(Value = names(sys.runtime), Time = as.vector(sys.runtime)),
            file = file.path(directory, "mcmc_runtime.csv"), sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)

# --- Write the estimates FIRST (these are the comparison outputs) ---
# Codon-specific estimates: Mutation (dM) + Selection (dOmega) from last 50%.
getCSPEstimates(parameter, file.path(directory, "Parameter_est", basename(input)),
                1, samples * percent.to.keep, report.original.ref = TRUE)

# Per-gene phi.
expressionValues <- getExpressionEstimates(parameter, c(1:size), samples * percent.to.keep, genome = genome)
write.table(expressionValues, file = file.path(directory, "Parameter_est/gene_expression.txt"),
            sep = ",", col.names = TRUE, quote = FALSE, row.names = FALSE)

# --- Then persist the R objects (best-effort; do not lose estimates on failure) ---
tryCatch(writeParameterObject(parameter, file.path(directory, "R_objects/parameter.Rda")),
         error = function(e) cat("WARN: writeParameterObject failed:", conditionMessage(e), "\n"))
tryCatch(writeMCMCObject(mcmc, file.path(directory, "R_objects/mcmc.Rda")),
         error = function(e) cat("WARN: writeMCMCObject failed:", conditionMessage(e), "\n"))

# a1 / a2 traces (constant at 4 unless estimated).
trace <- parameter$getTraceObject()
a1t <- unlist(trace$getInitiationCostTrace())
a2t <- unlist(trace$getElongationCostTrace())
cat("a1 last:", a1t[length(a1t)], " a2 last:", a2t[length(a2t)], "\n")

cat("FONSE FIT COMPLETE ->", directory, "\n")
