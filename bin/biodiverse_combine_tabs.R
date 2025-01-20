
## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")

cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-i", "--inpdir"), type = "character", default = ".",   help = "Input directory"),
    make_option(c("-p", "--prefix"), type = "character", default = "RND", help = "Biodiverse prefix of intput files"),
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))
# print(opt)

## Input parameters
INPDIR     <- opt$inpdir
PREFIX     <- opt$prefix

##########################################################
########################################################## Taxonomy-based filters
##########################################################



cat("\n\n-------- Loading Biodiverse results --------\n\n")



