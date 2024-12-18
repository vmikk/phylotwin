#!/usr/bin/env Rscript

## Estimate and visualize species diversity

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("ape")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")

