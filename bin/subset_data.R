#!/usr/bin/env Rscript

## Prepare a data subset for diversity estimation

## Notes:
# - Taxonomy filters are applied with AND condtion for the taxonomic ranks,
#   then, user-supplied species keys are added to the list
#   if no constraints are specified, all species keys from a phylogenetic tree are selected



## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("optparse")
load_pckg("data.table")
load_pckg("sf")
load_pckg("ape")
load_pckg("duckdb")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")
load_pckg("glue")
load_pckg("openxlsx")
# load_pckg("geos")
# load_pckg("wk")
# load_pckg("h3")

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

## Validation of the required arguments
if(is.na(opt$inpdir)){ cat("Input directory with pre-processed species occurrence counts in Parquet format is not specified.\n", file=stderr()); stop() }
if(is.na(opt$output)){ cat("Output prefix is not specified.\n", file=stderr()); stop() }
if(!is.na(opt$specieskeys) && !file.exists(opt$specieskeys)){ cat("The specified file with specieskeys does not exist.\n", file=stderr()); stop() }

if(!is.na(opt$latmin) && (opt$latmin < -90 || opt$latmin > 90)){ cat("The minimum latitude is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$latmax) && (opt$latmax < -90 || opt$latmax > 90)){ cat("The maximum latitude is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$latmin) && !is.na(opt$latmax) && opt$latmin > opt$latmax){ cat("Spatial filter: The minimum latitude is greater than the maximum latitude.\n", file=stderr()); stop() }

if(!is.na(opt$lonmin) && (opt$lonmin < -180 || opt$lonmin > 180)){ cat("The minimum longitude is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$lonmax) && (opt$lonmax < -180 || opt$lonmax > 180)){ cat("The maximum longitude is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$lonmin) && !is.na(opt$lonmax) && opt$lonmin > opt$lonmax){ cat("Spatial filter: The minimum longitude is greater than the maximum longitude.\n", file=stderr()); stop() }

## Temporary (may be changed later) - check if all coordinates are specified for a full bounding box (or none are specified)
coords <- c(opt$lonmin, opt$latmin, opt$lonmax, opt$latmax)
if( ! (all(sapply(coords, is.na)) || all(!sapply(coords, is.na))) ) {
  stop("You must specify the bounding box coordinates together (lonmin, latmin, lonmax, latmax).")
}

if(!is.na(opt$minyear) && (opt$minyear < 0 || opt$minyear > 3000)){ cat("The minimum year is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$maxyear) && (opt$maxyear < 0 || opt$maxyear > 3000)){ cat("The maximum year is out of range.\n", file=stderr()); stop() }
if(!is.na(opt$minyear) && !is.na(opt$maxyear) && opt$minyear > opt$maxyear){ cat("Additional filter: The minimum year is greater than the maximum year.\n", file=stderr()); stop() }

if(!is.na(opt$country)){
  ## December 2024: 249 current officially assigned ISO 3166-1 alpha-2 codes
  ## In Natural Earth v.5.1.2, there are 237 ISO-A2 codes
  ## ++ see https://www.statoids.com/wab.html

  valid_countries <- c(
    "AD", "AE", "AF", "AG", "AI", "AL", "AM", "AO", "AQ", "AR", "AS", "AT", "AU", "AW", "AX", "AZ", 
    "BA", "BB", "BD", "BE", "BF", "BG", "BH", "BI", "BJ", "BL", "BM", "BN", "BO", "BQ", "BR", "BS", "BT", "BV", "BW", "BY", "BZ", 
    "CA", "CC", "CD", "CF", "CG", "CH", "CI", "CK", "CL", "CM", "CN", "CO", "CR", "CU", "CV", "CW", "CX", "CY", "CZ",
    "DE", "DJ", "DK", "DM", "DO", "DZ", 
    "EC", "EE", "EG", "EH", "ER", "ES", "ET", 
    "FI", "FJ", "FK", "FM", "FO", "FR", 
    "GA", "GB", "GD", "GE", "GF", "GG", "GH", "GI", "GL", "GM", "GN", "GP", "GQ", "GR", "GS", "GT", "GU", "GW", "GY",
    "HK", "HM", "HN", "HR", "HT", "HU", 
    "ID", "IE", "IL", "IM", "IN", "IO", "IQ", "IR", "IS", "IT", 
    "JE", "JM", "JO", "JP", 
    "KE", "KG", "KH", "KI", "KM", "KN", "KP", "KR", "KW", "KY", "KZ", 
    "LA", "LB", "LC", "LI", "LK", "LR", "LS", "LT", "LU", "LV", "LY", 
    "MA", "MC", "MD", "ME", "MF", "MG", "MH", "MK", "ML", "MM", "MN", "MO", "MP", "MQ", "MR", "MS", "MT", "MU", "MV", "MW", "MX", "MY", "MZ", 
    "NA", "NC", "NE", "NF", "NG", "NI", "NL", "NO", "NP", "NR", "NU", "NZ", 
    "OM", "PA", "PE", "PF", "PG", "PH", "PK", "PL", "PM", "PN", "PR", "PS", "PT", "PW", "PY", 
    "QA", "RE", "RO", "RS", "RU", "RW", 
    "SA", "SB", "SC", "SD", "SE", "SG", "SH", "SI", "SJ", "SK", "SL", "SM", "SN", "SO", "SR", "SS", "ST", "SV", "SX", "SY", "SZ", 
    "TC", "TD", "TF", "TG", "TH", "TJ", "TK", "TL", "TM", "TN", "TO", "TR", "TT", "TV", "TW", "TZ", 
    "UA", "UG", "UM", "US", "UY", "UZ", "VA", "VC", "VE", "VG", "VI", "VN", "VU", 
    "WF", "WS", "XK", "YE", "YT", "ZA", "ZM", "ZW")

  # valid_countries[ ! valid_countries %in% NE$ISO_A2_EH ]    # not_supported_codes
  # NE[ ! NE$ISO_A2_EH %in% valid_countries, ]

  ## Codes missing in the Natural Earth data
  not_supported_codes <- c("BQ", "BV", "CC", "CX", "GF", "GI", "GP", "MQ", "RE", "SJ", "TK", "UM", "YT")

  valid_countries <- setdiff(valid_countries, not_supported_codes)

  ## Check if the country codes are valid
  countries <- toupper(unique(strsplit(x = opt$country, split = ",")[[1]]))
  if(any(! countries %in% valid_countries)){
    cat("The following country codes are not valid: ", 
        paste(countries[ !countries %in% valid_countries], collapse = ", "), "\n",
        file = stderr())
    stop()
  }
}



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
cat("\nInput-output parameters:\n")
cat("  Input directory:", INPDIR, "\n")
cat("  Output prefix:", OUTPUT, "\n")

cat("\nTaxonomy filters:\n")
cat("  Tree:",   TREE, "\n")
cat("  Phylum:", PHYLUM, "\n")
cat("  Class:",  CLASS, "\n")
cat("  Order:",  ORDER, "\n")
cat("  Family:", FAMILY, "\n")
cat("  Genus:",  GENUS, "\n")
cat("  Specieskeys:", SPECIESKEYS, "\n")

cat("\nSpatial filters:\n")
cat("  Country:", COUNTRY, "\n")
cat("  Latitude:", LATMIN, "-", LATMAX, "\n")
cat("  Longitude:", LONMIN, "-", LONMAX, "\n")
cat("  Polygon:", POLYGON, "\n")

cat("\nAdditional filters:\n")
cat("  Minimum year:", MINYEAR, "\n")
cat("  Maximum year:", MAXYEAR, "\n")
cat("  Basis of record:", BASISOFRECORD, "\n")

cat("\nInternal data of the pipeline: ", DATA, "\n")


##########################################################
########################################################## Taxonomy-based filters
##########################################################

## Initialize species keys
SPECIESKEYS_SELECTED <- vector(mode = "character")

## Based on the specified phylogenetic tree, load corresponding taxonomy table
cat("\nLoading taxonomic data for the specified phylogenetic tree\n")
  
## Get the name of taxonomy table (the same prefix as for the tree file)
curated_trees <- c(
  "Ferns_FTOL_1-7-0.nwk.gz", "FishTree.nwk.gz", "Henriquez-Piskulich_2024_BeetreeOfLife.nwk.gz", "Kawahara_2023_Butterflies.nwk.gz", 
  "Tietje_2023_Seed_plants_TACT.nwk.gz", "Varga_2019_Mushrooms.nwk.gz", "VertLife_Amphibians.nwk.gz", 
  "VertLife_Birds.nwk.gz", "VertLife_Mammals.nwk.gz", "VertLife_Squamates.nwk.gz")

if(basename(TREE) %in% curated_trees){
  cat("...using one of the curated phylogenies\n")

  tax_table <- file.path(DATA, "TaxonomyTables", "CuratedTrees_TaxonomyTable.parquet")

  if(!file.exists(tax_table)){
    cat("File with taxonomy table does not exist (", tax_table, ").\n", file=stderr())
    stop()
  }

  tree_id <- sub(pattern = "\\.nwk.gz$", replacement = "", x = basename(TREE))
}

## Load taxonomy table
# tax_table <- qs::qread(tax_table)
# tax_table <- arrow::read_parquet(tax_table)

tax_table <- open_dataset(tax_table) |> 
  filter(TreeID == tree_id) |> 
  collect() |> setDT()

cat("...number of records in taxonomy table:", nrow(tax_table), "\n")

## Apply taxonomic filters if specified
if(!is.na(PHYLUM)) {
  cat("...filtering by Phylum\n")
  phyla <- unique(strsplit(PHYLUM, ",")[[1]])
  tax_table <- tax_table[tax_table$phylum %in% phyla, ]
}
if(!is.na(CLASS)) {
  cat("...filtering by Class\n")
  classes <- unique(strsplit(CLASS, ",")[[1]])
  tax_table <- tax_table[tax_table$class %in% classes, ]
}
if(!is.na(ORDER)) {
  cat("...filtering by Order\n")
  orders <- unique(strsplit(ORDER, ",")[[1]])
  tax_table <- tax_table[tax_table$order %in% orders, ]
}
if(!is.na(FAMILY)) {
  cat("...filtering by Family\n")
  families <- unique(strsplit(FAMILY, ",")[[1]])
  tax_table <- tax_table[tax_table$family %in% families, ]
}
if(!is.na(GENUS)) {
  cat("...filtering by Genus\n")
  genera <- unique(strsplit(GENUS, ",")[[1]])
  tax_table <- tax_table[tax_table$genus %in% genera, ]
}

if(is.na(PHYLUM) && is.na(CLASS) && is.na(ORDER) && is.na(FAMILY) && is.na(GENUS)){
  cat("...no taxonomic filters specified, selecting all species keys from the tree\n")
} else {
  cat("...number of records in taxonomy table after filtering:", nrow(tax_table), "\n")
}

## Get species keys from filtered taxonomy
SPECIESKEYS_SELECTED <- c(SPECIESKEYS_SELECTED, unique(tax_table$specieskey))

## User-supplied species keys
if(! is.na(SPECIESKEYS)) {
  cat("\nLoading user-supplied species keys\n")
  USER_SPECIESKEYS <- fread(SPECIESKEYS, header = FALSE)$V1
  cat("...number of records in user-supplied species keys:", length(USER_SPECIESKEYS), "\n")
  
  SPECIESKEYS_SELECTED <- unique(c(SPECIESKEYS_SELECTED, USER_SPECIESKEYS))
}


cat("The total number of species keys selected:", length(SPECIESKEYS_SELECTED), "\n")

## Validate
if(length(SPECIESKEYS_SELECTED) == 0) {
  cat("No species selected for the analysis. Check your taxonomic filters or species keys file.\n", file=stderr())
  stop()
}


##########################################################
########################################################## Grid-cell-based filters
##########################################################
# - Country
# - Bounding box
# - Polygon

## Function to fetch H3 cells from WKT
get_h3_cells_from_wkt <- function(wkt, h3res) {
  # wkt - vectort of polygons in WKT format
  # h3res - H3 resolution

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
  hexes <- fread(tmp_h3, header = TRUE, sep = "\t")

  return(hexes)
}


## Initialize the grid-cells of interest
HEXES <- vector(mode = "character")

## Load grid-cells per country 
if(! is.na(COUNTRY)) {

  cat("Loading grid-cells at country level and the specified H3 resolution\n")
  COUNTRY <- toupper(unique(strsplit(x = COUNTRY, split = ",")[[1]]))
  country_path <- file.path(DATA, "Countries_H3", paste0("h3_", RESOLUTION), paste0(COUNTRY, ".txt.gz"))
  country_h3 <- alply(.data = country_path, .margins = 1, .fun = fread)
  country_h3 <- rbindlist(country_h3)
  
  cat("..number of grid-cell IDs loaded: ", nrow(country_h3))

  COUNTRY_HEXES <- unique(country_h3$h3_cell)
  rm(country_h3)

  ## Add to the main list of grid-cells
  HEXES <- c(HEXES, COUNTRY_HEXES)
}


## Get grid-cells from a bounding box
if(! is.na(LONMIN) && ! is.na(LONMAX) && ! is.na(LATMIN) && ! is.na(LATMAX)) {
  cat("Finding grid-cells from a bounding box\n")
  boundingbox <- c(LONMIN, LATMIN, LONMAX, LATMAX)
  
  ## Convert to WKT format
  BBOX <- sprintf('POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))',
    boundingbox[1], boundingbox[2],  # bottom-left
    boundingbox[3], boundingbox[2],  # bottom-right
    boundingbox[3], boundingbox[4],  # top-right
    boundingbox[1], boundingbox[4],  # top-left
    boundingbox[1], boundingbox[2]   # close polygon by returning to start
  )

  ## Get H3 cells that intersect with the bounding box
  BBOX_HEXES <- get_h3_cells_from_wkt(wkt = BBOX, h3res = RESOLUTION)
  cat("..number of grid-cells from the bounding box: ", nrow(BBOX_HEXES), "\n")

  ## Add to the main list of grid-cells
  HEXES <- c(HEXES, BBOX_HEXES$h3_cell)
}


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

  ## Get H3 cells that intersect with the polygon
  POLYGON_HEXES <- get_h3_cells_from_wkt(wkt = WKT, h3res = RESOLUTION)
  cat("..number of grid-cells from the polygon: ", nrow(POLYGON_HEXES), "\n")

  ## Add to the main list of grid-cells
  HEXES <- c(HEXES, POLYGON_HEXES$h3_cell)
}


## Overall list of grid-cells
HEXES <- unique(HEXES)
if(length(HEXES) > 0){
  cat("Total number of grid-cells to extract from the data: ", length(HEXES), "\n")
}




##########################################################
########################################################## Prepare a request to the database
##########################################################

cat("Preparing the database query\n")

## Export files with H3 cells and species keys
if(!length(HEXES) > 0){
  stop("No H3 cells selected for the analysis. Check your spatial filters.\n")
}
if(!length(SPECIESKEYS_SELECTED) > 0){
  stop("No species keys selected for the analysis. Check your taxonomic filters or species keys file.\n")
}

cat("Exporting H3 cells and species keys to files\n")

fwrite( 
  x = data.table(h3_cell = sort(HEXES)),
  file = "h3_cells.txt.gz",
  sep = "\t", quote = FALSE, col.names = TRUE, compress = "gzip")

fwrite(
  x = data.table(specieskey = sort(SPECIESKEYS_SELECTED)),
  file = "species_keys.txt.gz",
  sep = "\t", quote = FALSE, col.names = TRUE, compress = "gzip")


cat("\nInitializing DuckDB\n")

## Initialize DuckDB in-memory database
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")


cat("Preparing the query for DuckDB\n")

## Prepare the query for DuckDB
main_query <- glue("
-- Create temporary tables for filter values
CREATE TEMP TABLE h3_cells AS 
SELECT h3_cell AS h3_cell 
FROM read_csv('{h3_cells}', header=true, delim='\\t');

CREATE TEMP TABLE species_keys AS 
SELECT CAST(specieskey AS BIGINT) AS specieskey 
FROM read_csv('{species_keys}', header=true, delim='\\t');

-- Main filtered aggregation
COPY (
    WITH filtered_data AS (
        SELECT *
        FROM read_parquet('{INPDIR}/*.parquet')
        WHERE H3 IN (SELECT h3_cell FROM h3_cells)
        AND specieskey IN (SELECT specieskey FROM species_keys)",
  h3_cells     = "h3_cells.txt.gz",
  species_keys = "species_keys.txt.gz")

if(! is.na(MINYEAR)){ main_query <- glue(main_query, " AND year >= {MINYEAR}") }
if(! is.na(MAXYEAR)){ main_query <- glue(main_query, " AND year <= {MAXYEAR}") }

main_query <- glue(
  main_query, "
    )
    SELECT 
        H3,
        specieskey,
        SUM(record_count) as total_records
    FROM filtered_data
    GROUP BY H3, specieskey
) TO '{output_parquet_long}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 5);",
  output_parquet_long = "aggregated_counts.parquet")

main_query <- glue(
  main_query, "
-- Dataset keys with record counts
COPY (
    WITH filtered_data AS (
        SELECT *
        FROM read_parquet('{INPDIR}/*.parquet')
        WHERE H3 IN (SELECT h3_cell FROM h3_cells)
        AND specieskey IN (SELECT specieskey FROM species_keys)")

if(! is.na(MINYEAR)){ main_query <- glue(main_query, " AND year >= {MINYEAR}") }
if(! is.na(MAXYEAR)){ main_query <- glue(main_query, " AND year <= {MAXYEAR}") }

main_query <- glue(
  main_query, "
    ),
    unnested_keys AS (
        SELECT 
            UNNEST(dataset_keys) as datasetKey,
            1 as num_records,  -- count each dataset key once
            record_count       -- keep the pre-aggregated count (if there are multiple datasets already aggregated, we do not know the actual number of records per dataset key)
        FROM filtered_data
    )
    SELECT 
        datasetKey,
        SUM(num_records) as num_records,                    -- actual count of occurrences
        SUM(record_count) as dataset_records_approximate    -- sum of pre-aggregated counts
    FROM unnested_keys
    GROUP BY datasetKey
    ORDER BY dataset_records_approximate DESC
) TO '{output_datasetkeys}' (HEADER, DELIMITER '\\t');",
  output_datasetkeys = "dataset_keys.tsv")




