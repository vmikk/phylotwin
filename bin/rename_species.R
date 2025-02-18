#!/usr/bin/env Rscript

## In Biodiverse, there is an issue with handling species names with whitespaces
## Putting them in quotes does not work well, so we need to rename them

## Usage:
# Rscript bin/rename_species.R \
#   --tree /path/to/tree.nex.gz \
#   --occurrences /path/to/occurrences.csv \
#   --output_tree tree_renamed.nex \
#   --output_occurrences occurrences_renamed.csv


cat("Renaming species in tree and occurrences\n\n")

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("ape")


cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-t", "--tree"),         type = "character", default = "phylogenetic_tree.nex",   help = "Phylogenetic tree in Nexus format"),
    make_option(c("-o", "--occurrences"),  type = "character", default = "aggregated_counts.csv",   help = "Occurrences data (CSV format)"),
    make_option(c("--output_tree"),        type = "character", default = "tree_renamed.nwk",        help = "Output tree (Newick format)"),
    make_option(c("--output_occurrences"), type = "character", default = "occurrences_renamed.csv", help = "Output occurrences data (CSV format)")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))
# print(opt)

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


## Input parameters
TREE         <- opt$tree
OCCURRENCES  <- opt$occurrences
OUTPUT_TREE  <- opt$output_tree
OUTPUT_OCC   <- opt$output_occurrences

cat("\nInput-output parameters:\n")
cat("  Tree:",               TREE, "\n")
cat("  Occurrences:",        OCCURRENCES, "\n")
cat("  Output tree:",        OUTPUT_TREE, "\n")
cat("  Output occurrences:", OUTPUT_OCC, "\n")



## For debugging
# TREE         <- "phylogenetic_tree.nex"
# OCCURRENCES  <- "aggregated_counts.csv"
# OUTPUT_TREE  <- "tree_renamed.nwk"
# OUTPUT_OCC   <- "occurrences_renamed.csv"

##########################################################
########################################################## Load and merge tables
##########################################################

cat("\n\n-------- Processing data --------\n\n")

cat(".. Loading tree\n")
tree <- ape::read.nexus(TREE)

cat(".. Loading occurrences\n")
occ  <- fread(OCCURRENCES)

## Create a table with original species names
cat(".. Creating a table with species names\n")
SPP <- data.table(
  original = sort(unique(c(
    tree$tip.label,
    unique(occ$species)
  ))))

SPP[ , renamed := ape::makeLabel(original) ]

## Renaming
cat(".. Renaming species in tree\n")
tree$tip.label <- SPP$renamed[ match(x = tree$tip.label, table = SPP$original) ]

cat(".. Renaming species in occurrences\n")
SPP <- merge(x = occ, y = SPP, by.x = "species", by.y = "original", all.x = TRUE)
SPP[ , species := NULL ]
setnames(x = SPP, old = "renamed", new = "species")
setorder(x = SPP, H3, species)
setcolorder(x = SPP, neworder = c("H3", "Latitude", "Longitude", "species", "total_records"))

## Export data
cat("Exporting data\n")

cat(".. Exporting tree\n")
ape::write.tree(phy = tree, file = OUTPUT_TREE)

cat(".. Exporting occurrences\n")
fwrite(x = SPP, file = OUTPUT_OCC, sep = ",", quote = FALSE)

cat("\nDone!\n")

