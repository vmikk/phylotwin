
cat("Hypothesis testing\n\n")

## TODO:
# - treat polygons separately (splitref, splittest)
# - duckdb H3 extensions path


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("sf")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("h3")
load_pckg("PhyloMeasures")
load_pckg("phyloregion")
load_pckg("ape")

cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-1", "--polygons_reference"), type = "character", default = NA, help = "Reference polygons (GeoJSON or GeoPackage)"),
    make_option(c("-2", "--polygons_test"),      type = "character", default = NA, help = "Test polygons (GeoJSON or GeoPackage)"),
    # make_option("--splitref",  action="store_true", default=TRUE, help="Treat reference polygons separately"),
    # make_option("--splittest", action="store_true", default=TRUE, help="Treat test polygons separately"),
    make_option(c("-o", "--occurrences"),        type = "character", default = NA, help = "Species occurrences (Parquet, long format)"),
    make_option(c("-t", "--tree"),               type = "character", default = NA, help = "Phylogenetic tree (Nexus format)"),
    make_option("--resolution", action="store", default=4, type='integer', help="H3 resolution (e.g., 4)"),
    make_option(c("-r", "--results"),            type = "character", default = NA, help = "Results prefix")
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


## Parameters validation
if(is.na(opt$polygons_reference)) { cat("Reference polygons are not specified.\n",  file=stderr()); stop() }
if(is.na(opt$polygons_test))      { cat("Test polygons are not specified.\n",       file=stderr()); stop() }
if(is.na(opt$occurrences))        { cat("Species occurrences are not specified.\n", file=stderr()); stop() }
if(is.na(opt$tree))               { cat("Phylogenetic tree is not specified.\n",    file=stderr()); stop() }
if(is.na(opt$resolution))         { cat("H3 resolution is not specified.\n",        file=stderr()); stop() }
if(opt$resolution < 1 || opt$resolution > 15){ cat("H3 resolution must be between 1 and 15.\n", file=stderr()); stop() }


## Input parameters
POLYGONS_REFERENCE <- opt$polygons_reference
POLYGONS_TEST      <- opt$polygons_test
# SPLIT_REFERENCE    <- as.logical(opt$splitref)
# SPLIT_TEST         <- as.logical(opt$splittest)
OCCURRENCES        <- opt$occurrences
TREE               <- opt$tree
RESULTS            <- opt$results
RESOLUTION         <- opt$resolution

cat("\nParameters parsed:\n")
cat("  Reference polygons:",  POLYGONS_REFERENCE, "\n")
cat("  Test polygons:",       POLYGONS_TEST, "\n")
# cat("  Split reference:",     SPLIT_REFERENCE, "\n")
# cat("  Split test:",          SPLIT_TEST, "\n")
cat("  Species occurrences:", OCCURRENCES, "\n")
cat("  Phylogenetic tree:",   TREE, "\n")
cat("  H3 resolution:",       RESOLUTION, "\n")
cat("  Results prefix:",      RESULTS, "\n")

