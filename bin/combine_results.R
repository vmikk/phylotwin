#!/usr/bin/env Rscript

## Combine `estimate_diversity` results with Biodiverse-based results,
## and export results in various formats (GeoJSON, GeoPackage, tab-delimited)

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("h3")
load_pckg("sf")


cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-e", "--estdiv"), type = "character", default = "diversity_estimates.qs",    help = "Results from `estimate_diversity` process (QS format)"),
    make_option(c("-b", "--biodiv"), type = "character", default = "Biodiverse_results.txt.gz", help = "Results from Biodiverse subworkflow (tab-delimited format)"),
    make_option(c("-r", "--resolution"), action="store", default = 4, type='integer', help="H3 resolution (e.g., 4)"),
    make_option(c("-o", "--output"), type = "character", default = "diversity_estimates", help = "Output prefix")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))
# print(opt)

## Validation of arguments

## At least one of input files must be provided
if(is.na(opt$estdiv) && is.na(opt$biodiv)){
  cat("At least one of the input files with diversity estimates must be provided.\n", file=stderr()); stop()
}

if(is.na(opt$resolution)){ cat("H3 resolution is not specified.\n", file=stderr()); stop() }
if(opt$resolution < 1 || opt$resolution > 15){ cat("H3 resolution must be between 1 and 15.\n", file=stderr()); stop() }


## Input parameters
ESTDIV     <- opt$estdiv
BIODIV     <- opt$biodiv
RESOLUTION <- opt$resolution
OUTPUT     <- opt$output

cat("\nInput-output parameters:\n")
cat("  Diversity estimates:",   ESTDIV, "\n")
cat("  Biodiverse results:",    BIODIV, "\n")
cat("  H3 resolution:",         RESOLUTION, "\n")
cat("  Output prefix:",         OUTPUT, "\n")


