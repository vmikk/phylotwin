#!/usr/bin/env Rscript

## Combine `estimate_diversity` results with Biodiverse-based results,
## and export results in various formats (GeoJSON, GeoPackage, tab-delimited)

cat("Combining results and exporting spatial data\n\n")

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
    make_option(c("-e", "--estdiv"), type = "character", default = "estimate_diversity_results.qs", help = "Results from `estimate_diversity` process (QS format)"),
    make_option(c("-b", "--biodiv"), type = "character", default = "Biodiverse_results.txt.gz",     help = "Results from Biodiverse subworkflow (tab-delimited format)"),
    make_option(c("-r", "--resolution"), action="store", default = 4, type='integer', help="H3 resolution (e.g., 4)"),
    make_option(c("-o", "--output"), type = "character", default = "diversity_estimates", help = "Output prefix")
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


#### Validation of arguments

## At least one of input files must be provided
if(is.na(opt$estdiv) && is.na(opt$biodiv)){
  cat("At least one of the input files with diversity estimates must be provided.\n", file=stderr()); stop()
}
if(! is.na(opt$estdiv) && ! file.exists(opt$estdiv)){
  cat("Input file with diversity estimates does not exist.\n", file=stderr()); stop()
}
if(! is.na(opt$biodiv) && ! file.exists(opt$biodiv)){
  cat("Input file with Biodiverse results does not exist.\n", file=stderr()); stop()
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




#### For debugging
# ESTDIV <- "estimate_diversity_results.qs"
# BIODIV <- "Biodiverse_results.txt.gz"
# RESOLUTION <- 4
# OUTPUT <- "diversity_estimates"


##########################################################
########################################################## Load and merge tables
##########################################################

cat("\n\n-------- Loading and merging tables --------\n\n")

if(!is.na(ESTDIV)){
  cat("..Loading `estimate_diversity` results\n")
  ESTDIVt <- qs::qread(ESTDIV)
}

if(!is.na(BIODIV)){
  cat("..Loading Biodiverse results\n")
  BIODIVt <- fread(BIODIV)
}

## Merge tables
if(!is.na(ESTDIV) && !is.na(BIODIV)){
  cat("...Merging tables\n")
  RES <- merge(x = ESTDIVt, y = BIODIVt, by = "H3", all.x = TRUE)
} else if(!is.na(ESTDIV)){
  cat("...Using `estimate_diversity` results\n")
  RES <- ESTDIVt
} else if(!is.na(BIODIV)){
  cat("...Using Biodiverse results\n")
  RES <- BIODIVt
} else {
  cat("...No input tables found.\n", file=stderr()); stop()
}


##########################################################
########################################################## Prepare spatial data
##########################################################

cat("\n\n-------- Preparing spatial data --------\n\n")

cat("..Preparing H3 polygons\n")
H3_poly <- h3_to_geo_boundary_sf(RES$H3)

## Break polygons at the anti-meridian/date line (180/-180 degrees)
cat("..Break antimeridian\n")
H3_poly <- st_wrap_dateline(
  x = H3_poly,
  options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
  quiet = TRUE)

cat("..Adding diversity estimates to H3 polygons\n")
vars_to_exclude <- c("H3", "MostRangeRestricted", "MostPhylogeneticallyOriginal")
vars <- colnames(RES)[! colnames(RES) %in% vars_to_exclude]
H3_poly <- cbind(H3_poly, RES[, ..vars])

cat("..Exporting polygons with divsity estimates in GeoPackage format\n")
st_write(
  obj   = H3_poly,
  dsn   = paste0(OUTPUT, ".gpkg"),
  layer = "diversity_estimates")

cat("..Exporting polygons with divsity estimates in GeoJSON format\n")
st_write(
  obj   = H3_poly,
  dsn   = paste0(OUTPUT, ".geojson"))


##########################################################
########################################################## Export results
##########################################################


cat("\n\n-------- Exporting tab-delimited file --------\n\n")

## Add grid cell coordinates
cat("..Adding geo-coordinates for grid cell centers\n")
RES[, c("Latitude", "Longitude") := as.data.table(h3::h3_to_geo(H3)) ]
 
## Reorder columns
cat("..Reordering columns\n")
setcolorder(RES, c("H3", "Latitude", "Longitude"))

cat("..Exporting tab-delimited file\n")
fwrite(x = RES, file = paste0(OUTPUT, ".txt"), sep = "\t")



cat("\n\n-------- All done --------\n\n")
