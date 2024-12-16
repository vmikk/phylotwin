#!/usr/bin/env Rscript

## Prepare a data subset for diversity estimation

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("sf")
load_pckg("geos")
load_pckg("wk")
load_pckg("ape")
load_pckg("duckdb")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")
load_pckg("h3")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-i", "--inpdir"), type = "character", default = NULL, help = "Input - Directory with pre-processed species occurrence counts (Parquet format)"),
    make_option(c("-o", "--output"), type = "character", default = NULL, help = "Output prefix"),

    ## Taxonomy filters
    make_option("--tree",        action="store", default=NA, type='character', help="Phylogenetic tree (Newick format) defining the 'major' group of interest"),
    make_option("--phylum",      action="store", default=NA, type='character', help="Comma-separated list of phyla to select"),
    make_option("--class",       action="store", default=NA, type='character', help="Comma-separated list of classes to select"),
    make_option("--order",       action="store", default=NA, type='character', help="Comma-separated list of orders to select"),
    make_option("--family",      action="store", default=NA, type='character', help="Comma-separated list of families to select"),
    make_option("--genus",       action="store", default=NA, type='character', help="Comma-separated list of genera to select"),
    make_option("--specieskeys", action="store", default=NA, type='character', help="File with user-supplied GBIF specieskeys"),

    ## Spatial filters
    make_option("--resolution", action="store", default=4, type='integer', help="H3 resolution (e.g., 4)"),
    make_option("--country", action="store", default=NA, type='character', help="Comma-separated list of country codes (e.g., AU,CA), ISO 3166 format"),
    make_option("--latmin",  action="store", default=NA, type='double',    help="Minimum latitude"),
    make_option("--latmax",  action="store", default=NA, type='double',    help="Maximum latitude"),
    make_option("--lonmin",  action="store", default=NA, type='double',    help="Minimum longitude"),
    make_option("--lonmax",  action="store", default=NA, type='double',    help="Maximum longitude"),
    make_option("--polygon", action="store", default=NA, type='character', help="Custom area of interest (a file with polygons in WKT, GeoPackage, or GeoJSON format)"),

    ## Additional filters
    make_option("--minyear", action="store", default=1945, type='integer', help="Minimum year of occurrence (default, 1945)"),
    make_option("--maxyear", action="store", default=NA,   type='integer', help="Maximum year of occurrence"),
    make_option("--basisofrecord", action="store", default=NA, type='character', help="Basis of record to include from the data"),

    make_option("--data", action="store", default="./data", type='character', help="Path to the internal data of the pipeline"),
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}


## Input parameters
INPDIR <- opt$inpdir
OUTPUT <- opt$output

TREE   <- to_na( opt$tree )
PHYLUM <- to_na( opt$phylum )
CLASS  <- to_na( opt$class )
ORDER  <- to_na( opt$order )
FAMILY <- to_na( opt$family )
GENUS  <- to_na( opt$genus )
SPECIESKEYS  <- to_na( opt$specieskeys )

RESOLUTION <- opt$resolution
COUNTRY <- to_na( opt$country )
LATMIN  <- as.numeric( to_na(opt$latmin) )
LATMAX  <- as.numeric( to_na(opt$latmax) )
LONMIN  <- as.numeric( to_na(opt$lonmin) )
LONMAX  <- as.numeric( to_na(opt$lonmax) )
POLYGON <- to_na( opt$polygon )

MINYEAR <- as.numeric(to_na( opt$minyear) )
MAXYEAR <- as.numeric(to_na( opt$maxyear) )
BASISOFRECORD <- to_na( opt$basisofrecord )

DATA <- opt$data


