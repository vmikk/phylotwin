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




## Input parameters
OCC     <- opt$input
TREE    <- opt$tree
DIV     <- opt$div
THREADS <- opt$threads

cat("\nParameters parsed:\n")
cat("  Input:", OCC, "\n")
cat("  Tree:", TREE, "\n")
cat("  Diversity metrics:", DIV, "\n")
cat("  Number of threads:", THREADS, "\n")

## Specify the number of threads for data.table
setDTthreads(threads = THREADS)


##########################################################
########################################################## Load and reshape data
##########################################################

cat("\n\n-------- Loading data --------\n\n")

## Load phylogenetic tree
cat("Loading phylogenetic tree: ", TREE, "\n")
tree <- read.tree(TREE)

## Load aggregated species occurrences
cat("Loading aggregated species occurrences\n")
occ <- read_parquet(OCC)
setDT(occ)

if(any(!occ$total_records >= 1)){
  cat("..some records have missing values\n")
  occ <- occ[ total_records >= 1 ]
}

## Scale species occurrences to presence/absence
occ[ , PA := 1 ]

## Reshape species occurrences to wide format
cat("Reshaping species occurrences to wide format\n")
datt <- dcast(
  data = occ,
  formula = H3 ~ specieskey,
  value.var = "PA",            # value.var = "total_records",
  fun.aggregate = sum
)

# any(datt[,-1] > 1)

## Subset phylogenetic tree to species present in the data
cat("Subsetting phylogenetic tree\n")
tree <- keep.tip(tree, intersect(tree$tip.label, colnames(datt)[-1]))

## Data summary
cat("Data summary:\n")
cat("  Number of species:", ncol(datt) - 1, "\n")
cat("  Number of H3 cells:", nrow(datt), "\n")


## Convert data to matrix
cat("Converting data to matrix\n")
datm <- as.matrix(datt[,-1])
rownames(datm) <- datt$H3


##########################################################
########################################################## Estimate species diversity
##########################################################

cat("\n\n-------- Estimating species diversity --------\n\n")

## Parse diversity metrics to estimate
MEASURES <- sort(unique( strsplit(DIV, ",")[[1]] ))
cat("Number of diversity metrics to estimate:", length(MEASURES), "\n")
cat("Metrics: ", paste(MEASURES, collapse = ", "), "\n")

## Initialize results
RES <- vector("list")


######## Number of records per H3 cell

cat("..Estimating number of records per H3 cell\n")
num_records <- occ[ , .(NumRecords = sum(total_records, na.rm = TRUE)), by = "H3" ]

######## Traditional indices

cat("..Estimating species richness\n")
RES <- c(RES, list(Richness = rowSums(datt[,-1]) ))



######## PhyloMeasure-based indices

if("PD" %in% MEASURES){
  cat("..Estimating PD\n")
  RES <- c(RES, list(PD = PhyloMeasures::pd.query(tree = tree, matrix = datm) ))
}
if("MPD" %in% MEASURES){
  cat("..Estimating MPD\n")
  RES <- c(RES, list(MPD = PhyloMeasures::mpd.query(tree = tree, matrix = datm) ))
}
if("MNTD" %in% MEASURES){
  cat("..Estimating MNTD\n")
  RES <- c(RES, list(MNTD = PhyloMeasures::mntd.query(tree = tree, matrix = datm) ))
}

if("SES.PD" %in% MEASURES){
  cat("..Estimating SES.PD\n")
  RES <- c(RES, list(SES.PD = PhyloMeasures::pd.query(tree = tree, matrix = datm, standardize = T) ))
}
if("SES.MPD" %in% MEASURES){
  cat("..Estimating SES.MPD\n")
  RES <- c(RES, list(SES.MPD = PhyloMeasures::mpd.query(tree = tree, matrix = datm, standardize = T) ))
}
if("SES.MNTD" %in% MEASURES){
  cat("..Estimating SES.MNTD\n")
  RES <- c(RES, list(SES.MNTD = PhyloMeasures::mntd.query(tree = tree, matrix = datm, standardize = T) ))
}



