#!/usr/bin/env Rscript

## Estimate and visualize species diversity

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("PhyloMeasures")
load_pckg("ape")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")

cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-i", "--input"),  type = "character", default = NA, help = "File with aggregated species occurrences (Parquet, long format)"),
    make_option(c("-t", "--tree"),   type = "character", default = NA, help = "Phylogenetic tree (Newick format)"),
    make_option(c("-d", "--div"),    type = "character", default = NA, help = "Diversity metrics (comma-separated list)"),
    make_option(c("-o", "--output"), type = "character", default = NA, help = "Output prefix"),
    make_option("--threads",         type = "integer",   default = 4,  help = "Number of CPUs to use")

)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)



