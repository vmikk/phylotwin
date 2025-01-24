#!/usr/bin/env Rscript


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("h3")

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

## Function to read Biodiverse results
read_bd <- function(file){
  
  res <- fread(file, sep = ",")
  
  ## ELEMENT = concatenated coords 20.4069805389076:68.992705038507
  res[, ELEMENT := NULL]
  
  setnames(res,
    old = c("Axis_1", "Axis_0"),
    new = c("Latitude", "Longitude"))
  
  res[ , H3 := h3::geo_to_h3(res[, .(Latitude, Longitude)], res = RESOLUTION) ]
  res[ , c("Latitude", "Longitude") := NULL ]
  setcolorder(res, "H3")

  return(res)
}




