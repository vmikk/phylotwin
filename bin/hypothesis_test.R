
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

##########################################################

# ## Debugging
# POLYGONS_REFERENCE <- "poly1.geojson"
# POLYGONS_TEST      <- "poly2.geojson"
# # SPLIT_REFERENCE  <- TRUE
# # SPLIT_TEST       <- TRUE
# OCCURRENCES        <- "aggregated_counts.parquet"
# TREE               <- "phylogenetic_tree.nex"
# RESOLUTION         <- 4
# RESULTS            <- "tests"
# Sys.setenv(PATH=paste("~/.nextflow/assets/vmikk/phylotwin/bin", Sys.getenv("PATH"), sep=":"))



##########################################################
########################################################## Load data
##########################################################

## Function to load polygons
load_polygons <- function(poly, separate = FALSE){
  # separate - if TRUE, return a list of polygons
  
  if(!file.exists(poly)){
    cat("ERROR: File with polygons not found: ", poly, "\n", file=stderr())
    stop()
  }
  
  res <- try( st_read(poly, quiet = TRUE) )
  
  if("try-error" %in% class(res)){
    cat("ERROR: Failed to read polygons from the file: ", poly, "\n", file=stderr())
    stop()
  }
  
  ## Check geometry types
  geom_types <- unique(st_geometry_type(res))
  if(! all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))){
    unsupported_geoms <- geom_types[ ! geom_types %in% c("POLYGON", "MULTIPOLYGON") ]
    unsupported_geoms <- paste(unsupported_geoms, collapse = ", ")
    cat("....provided polygon has unsupported geometry type: ", unsupported_geoms, "\n")
    cat("....only POLYGON or MULTIPOLYGON are supported\n")
    cat("....unsupported geometries will be ignored\n")
    res <- res[ st_geometry_type(res) %in% c("POLYGON", "MULTIPOLYGON"), ]
    cat("....number of geometries after filtering: ", nrow(res), "\n")
  }

  if(nrow(res) == 0){
    cat("ERROR: No polygons found in the file: ", poly, "\n", file=stderr())
    stop()
  } else {
    cat("...Loaded ", nrow(res), " geometries\n")
  }

  ## Check CRS
  if(st_crs(res)$epsg != 4326){
    cat("...provided polygon(s) has CRS: EPSG ", st_crs(res)$epsg, "\n")
    cat("...reprojecting to EPSG:4326\n")
    res <- st_transform(res, crs = 4326)
  }

  ## TODO: --- split geometries if requested --> if multipolygons, split them independently within each geometry

  ## Convert MULTIPOLYGONs to POLYGONs
  if("MULTIPOLYGON" %in% st_geometry_type(res)){
    cat("...converting MULTIPOLYGONs to POLYGONs\n")
    res <- st_cast(x = res, to = "POLYGON", warn = FALSE, do_split = TRUE)
  }

  return(res)
} # end of `load_polygons`


cat("\n\n-------- Loading data --------\n\n")


## Load species occurrences
cat("..Loading species occurrences\n")
occ <- open_dataset(OCCURRENCES)

## Validate H3 resolution
tmp <- occ %>% head(1) %>% collect()
if(h3_get_resolution(tmp$H3) != RESOLUTION){
  cat("ERROR: H3 resolution in the species occurrences does not match the specified resolution.\n", file=stderr())
  stop()
}
rm(tmp)

## Load polygons
cat("..Loading reference polygons\n")
polygons_reference <- load_polygons(POLYGONS_REFERENCE)

cat("..Loading test polygons\n")
polygons_test <- load_polygons(POLYGONS_TEST)


## Load phylogenetic tree
cat("..Loading phylogenetic tree\n")
# tree <- ape::read.tree(TREE)               # Newick format
tree <- ape::read.nexus(file = TREE)         # Nexus format



