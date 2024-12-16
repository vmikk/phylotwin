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


## Initialize the grid-cells of interest
HEXES <- vector(mode = "character")

## Get grid-cells from a polygon
if(! is.na(POLYGON)) {

  cat("Loading user-supplied polygon\n")
  polygon_ext <- tolower( tools::file_ext(POLYGON) )
  
  ## Supported extensions
  ext_wkt     <- c("wkt", "wkb", "ewkt", "ewkb", "agf")
  ext_gpkg    <- c("gpkg", "geopackage")
  ext_geojson <- c("geojson", "json", "geojsonp", "geojsons", "geojsonsl", "geojsonss", "geojsonsz")

  if(polygon_ext %in% ext_wkt){
    cat("..based on file extension, WKT format detected\n")
    
    # POLY <- geos::geos_read_wkt(wkt = POLYGON, fix_structure = TRUE)
    # POLY <- wk::read_wkt(POLY)
    # cat("..WKT summary:\n")
    # print(wk::wk_count(POLY))

    # WKT <- geos::geos_write_wkt(POLY)   # fixed in case of errors

    # https://github.com/cran/rgeos/blob/master/R/rgeos_wkt.R

    cat("Currently, support for WKT is not implemented\n")
    stop()
  }

  if(!polygon_ext %in% c(ext_gpkg, ext_geojson)){
    cat("The specified polygon file has an unsupported extension.\n", file=stderr())
    cat("Supported extensions (case-insensitive): ", paste(c(ext_gpkg, ext_geojson), collapse = ", "), "\n")
    stop()
  }

  if(polygon_ext %in% c("gpkg", "geopackage")){
    cat("..based on file extension, GeoPackage format detected\n")
  }
  if(polygon_ext %in% ext_geojson){
    cat("..based on file extension, GeoJSON format detected\n")
  }

  ## Read the user-supplied polygon file
  cat("...reading the polygon file\n")
  POLY <- st_read(POLYGON, promote_to_multi = TRUE)

  ## Check geometry types
  geom_types <- unique(st_geometry_type(POLY))
  if(! all(geom_types %in% c("POLYGON", "MULTIPOLYGON"))){
    unsupported_geoms <- geom_types[ ! geom_types %in% c("POLYGON", "MULTIPOLYGON") ]
    unsupported_geoms <- paste(unsupported_geoms, collapse = ", ")
    cat("...provided polygon has unsupported geometry type: ", unsupported_geoms, "\n")
    cat("...only POLYGON or MULTIPOLYGON are supported\n")
    cat("...unsupported geometries will be ignored\n")
    POLY <- POLY[ st_geometry_type(POLY) %in% c("POLYGON", "MULTIPOLYGON"), ]
    cat("...number of geometries after filtering: ", nrow(POLY), "\n")
    if(nrow(POLY) == 0){
      cat("...no valid geometries left\n")
      stop()
    }
  }

  ## Check CRS
  if(st_crs(POLY)$epsg != 4326){
    cat("...provided polygon has CRS: EPSG ", st_crs(POLY)$epsg, "\n")
    cat("...reprojecting to EPSG:4326\n")
    POLY <- st_transform(POLY, crs = 4326)
  }

  ## Convert MULTIPOLYGONs to POLYGONs
  cat("...converting MULTIPOLYGONs to POLYGONs\n")
  POLY <- st_cast(x = POLY, to = "POLYGON", warn = FALSE, do_split = TRUE)

  ## Convert to WKT format
  cat("...converting to WKT format\n")
  WKT <- st_as_text(POLY$geom)
  
  ## Get all hexagons with CENTERS contained in a given polygon
  # POLYGON_HEXES <- h3::polyfill(POLY, res = RESOLUTION)

  ## Export to file
  cat("...exporting WKT to file\n")
  tmp_wkt <- tempfile(pattern = "polygon_wkt", fileext = ".wkt")
  fwrite(
    x = data.table(Polygon = WKT),
    file = tmp_wkt,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## Get H3 cells that intersect with the polygon
  cat("...getting H3 cells\n")
  tmp_h3 <- tempfile(pattern = "polygon_h3", fileext = ".txt")
  
  system2(
    command = system(command = "which wkt_polygon_to_h3_cells.sh", intern = TRUE),
    args = c(
      "-i", tmp_wkt, 
      "-o", tmp_h3, 
      "-r", RESOLUTION, 
      "-a", "experimental", 
      "-c", "CONTAINMENT_OVERLAPPING"),
    stdout = "", stderr = "",             # output stdout and stderr to R console
    wait = TRUE)

  ## Load the H3 cells
  POLYGON_HEXES <- fread(tmp_h3, header = TRUE, sep = "\t")
  POLYGON_HEXES <- POLYGON_HEXES$h3_cell

  ## Add to the main list of grid-cells
  HEXES <- c(HEXES, POLYGON_HEXES)
}



