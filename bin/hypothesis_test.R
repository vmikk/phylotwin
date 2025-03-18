#!/usr/bin/env Rscript

## TODO:
# - treat polygons separately (splitref, splittest)
# - duckdb H3 extensions path

## Output files:
# - diversity table       (prefix + "_diversity.txt", e.g. `HypTest_diversity.txt`)
# - species originalities (prefix + "_species_originalities.txt", e.g. `HypTest_species_originalities.txt`)


cat("Hypothesis testing\n\n")

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
load_pckg("jsonlite")

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
    make_option("--resolution", action="store",  type='integer',     default=4,    help="H3 resolution (e.g., 4)"),
    make_option(c("-r", "--results"),            type = "character", default = NA, help = "Results prefix"),
    make_option("--duckdb_extdir", action="store", default=NA, type='character', help="Path to the DuckDB extensions directory")
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
# SPLIT_REFERENCE  <- as.logical(opt$splitref)
# SPLIT_TEST       <- as.logical(opt$splittest)
OCCURRENCES        <- opt$occurrences
TREE               <- opt$tree
RESULTS            <- opt$results
RESOLUTION         <- opt$resolution
DUCKDB_EXT_DIR     <- ifelse( is.na(opt$duckdb_extdir), yes = "", no = opt$duckdb_extdir)

cat("\nParameters parsed:\n")
cat("  Reference polygons:",  POLYGONS_REFERENCE, "\n")
cat("  Test polygons:",       POLYGONS_TEST, "\n")
# cat("  Split reference:",     SPLIT_REFERENCE, "\n")
# cat("  Split test:",          SPLIT_TEST, "\n")
cat("  Species occurrences:", OCCURRENCES, "\n")
cat("  Phylogenetic tree:",   TREE, "\n")
cat("  H3 resolution:",       RESOLUTION, "\n")
cat("  Results prefix:",      RESULTS, "\n")
cat("  DuckDB extensions directory:", DUCKDB_EXT_DIR, "\n")
cat("  Working directory:",   getwd(), "\n")

##########################################################

# ## Debugging
# POLYGONS_REFERENCE <- "poly1.geojson"
# POLYGONS_TEST      <- "poly2.geojson"
# # SPLIT_REFERENCE  <- TRUE
# # SPLIT_TEST       <- TRUE
# OCCURRENCES        <- "aggregated_counts.parquet"
# TREE               <- "phylogenetic_tree.nex"
# RESOLUTION         <- 4
# RESULTS            <- "HypTest"

## Local
# Sys.setenv(PATH=paste(normalizePath(file.path(Sys.getenv("HOME"), ".nextflow/assets/vmikk/phylotwin/bin")), Sys.getenv("PATH"), sep=":"))
## Sys.getenv("PATH")
# DUCKDB_EXT_DIR <- ""

## Docker
# Sys.setenv(PATH=paste(normalizePath(file.path(Sys.getenv("HOME"), ".nextflow/assets/vmikk/phylotwin/bin")), Sys.getenv("PATH"), sep=":"))
# DUCKDB_EXT_DIR <- "/usr/local/bin/duckdb_ext"


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



##########################################################
########################################################## Get H3 cells of polygons
##########################################################

cat("\n\n-------- Converting polygons to H3 cells --------\n\n")

## Function to fetch H3 cells from WKT
get_h3_cells_from_polygons <- function(poly, h3res = RESOLUTION) {
  # poly  - simple feature collection with polygons
  # h3res - H3 resolution

  ## Get data in WKT format
  cat("...converting to WKT format\n")
  WKT <- st_as_text(poly$geom)

  ## Export to file
  cat("...exporting WKT to file\n")
  tmp_wkt <- tempfile(pattern = "polygon_wkt", fileext = ".wkt", tmpdir = getwd())
  fwrite(
    x = data.table(Polygon = WKT),
    file = tmp_wkt,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  ## Get H3 cells that intersect with the polygon
  cat("...getting H3 cells\n")
  tmp_h3 <- tempfile(pattern = "polygon_h3", fileext = ".txt", tmpdir = getwd())
  
  ## Arguments for the script
  ARGS <- c(
      "-i", tmp_wkt, 
      "-o", tmp_h3, 
      "-r", h3res, 
      "-a", "experimental", 
      "-c", "CONTAINMENT_OVERLAPPING")
  
  if(DUCKDB_EXT_DIR != ""){
    ARGS <- c(ARGS, "-e", DUCKDB_EXT_DIR)
  }

  system2(
    command = system(command = "which wkt_polygon_to_h3_cells.sh", intern = TRUE),
    args = ARGS,
    stdout = "", stderr = "",             # output stdout and stderr to R console
    wait = TRUE)

  ## Load the H3 cells
  hexes <- fread(tmp_h3, header = TRUE, sep = "\t")

  ## Clean up
  unlink(tmp_wkt)
  unlink(tmp_h3)

  return(hexes)
}



## Get all H3 cells that intersect with the polygon
cat("..Getting H3 cells from the reference polygons\n")

h3_reference <- get_h3_cells_from_polygons(
    poly = polygons_reference,
    h3res = RESOLUTION)

cat("...number of grid-cells from reference polygons: ", nrow(h3_reference), "\n")


## Get H3 cells from the test polygons
cat("..Getting H3 cells from the test polygons\n")

h3_test <- get_h3_cells_from_polygons(
    poly = polygons_test,
    h3res = RESOLUTION)

cat("...number of grid-cells from test polygons: ", nrow(h3_test), "\n")


cat("\nCounting the number of grid cells\n")

# num_grid_cells <- rowwiseDT(   # supported only in newer data.table versions
#   Geometry=, GridCells=,
#   "EntireArea", occ %>% select(H3) %>% unique() %>% collect() %>% nrow(),
#   "Reference",  nrow(h3_reference),
#   "Test",       nrow(h3_test))

num_grid_cells <- data.table(
  Geometry     = c("EntireArea", "Reference", "Test"),
  GridCells    = c(
    occ %>% select(H3) %>% unique() %>% collect() %>% nrow(),
    nrow(h3_reference),
    nrow(h3_test)
  ))



##########################################################
########################################################## Get occurrences per geometry
##########################################################

cat("\n\n-------- Extracting species occurrences per geometry --------\n\n")

## Load occurrences
cat("..Reading species occurrences per geometry and summarizing species abundances\n")

occ_glob <- occ %>%
  group_by(species) %>%
  summarise(Abundance = sum(total_records)) %>%
  collect() %>%
  setDT()

occ_reference <- occ %>% 
  filter(H3 %in% h3_reference$h3_cell) %>% 
  group_by(species) %>%
  summarise(Abundance = sum(total_records)) %>%
  collect() %>%
  setDT()

occ_test <- occ %>% 
  filter(H3 %in% h3_test$h3_cell) %>% 
  group_by(species) %>%
  summarise(Abundance = sum(total_records)) %>%
  collect() %>%
  setDT()


## Reshape to wide format
cat("..Reshaping to wide format\n")

tab <- rbind(
    data.table(Geometry = "EntireArea", occ_glob),
    data.table(Geometry = "Reference",  occ_reference),
    data.table(Geometry = "Test",       occ_test))

tabw <- dcast(
    data = tab,
    formula = Geometry ~ species,
    value.var = "Abundance",
    fun.aggregate = sum)

## Convert to matrix
datm <- as.matrix(tabw[ , -1])
datm[ datm > 0 ] <- 1
rownames(datm) <- tabw$Geometry


##########################################################
########################################################## Species originality estimation
##########################################################

cat("\n\n-------- Estimating species originality --------\n\n")

## Estimate fair proportion index (== FP == ED == Evolutionary distinctiveness)
cat("..Estimating fair proportion index\n")
FP <- phyloregion::evol_distinct(tree, type = "fair.proportion")
FP <- data.table(species = names(FP), FP = FP)

## Estimate the inverse of the range size of each taxon. Sum(RR) = Weighted endemism
cat("..Estimating the inverse of the range size of each taxon\n")
RR <- occ %>%
  group_by(species) %>%
  summarise(RangeSize = n()) %>%
  collect() %>%
  setDT()

## Combine data
SPEC <- merge(x = FP, y = RR, by = "species", all = TRUE)

## Top-N species
fpp <- quantile(x = SPEC$FP,        probs = 0.05)
rrp <- quantile(x = SPEC$RangeSize, probs = 0.05)

SPEC[ , TopPhylogeneticallyDistinct := fifelse(FP <= fpp,        "Yes", "No") ]
SPEC[ , TopRangeRestricted          := fifelse(RangeSize <= rrp, "Yes", "No") ]

SPEC[ , In_EntireArea := TRUE ]
SPEC[ , In_Reference  := species %in% occ_reference$species ]
SPEC[ , In_Test       := species %in% occ_test$species ]


## Count number of "top" species per area
TOPSP <- rbind(
    data.table(
      Geometry = "EntireArea",
      PhylogeneticallyDistinctSpecies = sum(SPEC[ TopPhylogeneticallyDistinct %in% "Yes" ]$In_EntireArea),
      RangeRestrictedSpecies          = sum(SPEC[ TopRangeRestricted          %in% "Yes" ]$In_EntireArea)),
    data.table(
      Geometry = "Reference",
      PhylogeneticallyDistinctSpecies = sum(SPEC[ TopPhylogeneticallyDistinct %in% "Yes" ]$In_Reference),
      RangeRestrictedSpecies          = sum(SPEC[ TopRangeRestricted          %in% "Yes" ]$In_Reference)),
    data.table(
      Geometry = "Test",
      PhylogeneticallyDistinctSpecies = sum(SPEC[ TopPhylogeneticallyDistinct %in% "Yes" ]$In_Test),
      RangeRestrictedSpecies          = sum(SPEC[ TopRangeRestricted          %in% "Yes" ]$In_Test))
)



##########################################################
########################################################## Diversity estimation
##########################################################

cat("\n\n-------- Estimating diversity --------\n\n")

## Subset tree to species present in the matrix
cat("..Subsetting tree to species present in the matrix\n")
tree <- ape::drop.tip(tree, tip = setdiff(tree$tip.label, colnames(datm)))

cat("..Estimating diversity\n")
DIV <- data.table(
    Geometry              = rownames(datm),
    TotalSpecies          = rowSums(datm),
    PD                    = PhyloMeasures::pd.query(tree = tree, matrix = datm),
    SES.PD                = PhyloMeasures::pd.query(tree = tree, matrix = datm, standardize = T),
    PhyloEndemismWeighted = phyloregion::phylo_endemism(x = datm, phy = tree, weighted = TRUE),
    PhyloEndemismStrict   = phyloregion::phylo_endemism(x = datm, phy = tree, weighted = FALSE)
    )

DIV[ , PD_proportion := round(PD / DIV[ Geometry == "EntireArea" ]$PD * 100, 2) ]

## Add number of grid cells
DIV <- merge(DIV, num_grid_cells, by = "Geometry", all.x = TRUE)

## Conservation value (as per Cadotte et al. 2010, doi:10.1111/j.1472-4642.2010.00650.x)
DIV[ , ConservationValue := PD / GridCells ]

## Add number of "top" species
DIV <- merge(x = DIV, y = TOPSP, by = "Geometry", all.x = TRUE)

## Reorder columns
setcolorder(DIV, c(
  "Geometry",
  "GridCells",
  "TotalSpecies", "PhylogeneticallyDistinctSpecies", "RangeRestrictedSpecies",
  "PD", "PD_proportion", "SES.PD",
  "PhyloEndemismWeighted", "PhyloEndemismStrict",
  "ConservationValue"
  ))

## Export results
cat("..Exporting results\n")

cat("...writing diveristy table\n")
fwrite(x = DIV, file = paste0(RESULTS, "_diversity.txt"), sep = "\t")

cat("...writing species originalities\n")
fwrite(x = SPEC, file = paste0(RESULTS, "_species_originalities.txt"), sep = "\t")

##########################################################
########################################################## Beutify the results for web-GUI
##########################################################

cat("\n\n-------- Beutifying results for web-GUI --------\n\n")

cat("..Formatting results\n")

DIVN <- copy(DIV)
# DIVN[ , PD          := round(PD, 2) ]
DIVN[ , PD_proportion := round(PD_proportion, 2) ]
DIVN[ , SES.PD        := round(SES.PD, 2) ]
DIVN[ , PhyloEndemismWeighted := round(PhyloEndemismWeighted, 2) ]
DIVN[ , PhyloEndemismStrict   := round(PhyloEndemismStrict, 2) ]
DIVN[ , ConservationValue     := round(ConservationValue, 2) ]

DIVN[ Geometry %in% "EntireArea", Geometry := "Entire area" ]
DIVN[ Geometry %in% "EntireArea", SES.PD := NA ]
DIVN[ Geometry %in% "EntireArea", ConservationValue := NA ]

## Rename columns
cat("..Renaming columns\n")

cat("..Adding variable metadata\n")

## Create metadata for each column with colors and descriptions
metadata <- list(
  Geometry = list(
    color = "#1f77b4",
    description = "Area."
  ),
  GridCells = list(
    color = "#ff7f0e",
    description = "Number of grid cells in the area."
  ),
  Species = list(
    color = "#2ca02c",
    description = "Total number of species in the area."
  ),
  PhylogeneticallyDistinctSpecies = list(
    color = "#d62728",
    description = "Count of the most phylogenetically distinct species (phylogenetic distinctiveness <= 5th percentile)."
  ),
  RangeRestrictedSpecies = list(
    color = "#9467bd",
    description = "Count of the most range-restricted species (range size <= 5th percentile)."
  ),
  PD = list(
    color = "#8c564b",
    description = "Phylogenetic Diversity (PD) metric."
  ),
  PD_proportion = list(
    color = "#e377c2",
    description = "Proportion of the total phylogenetic diversity."
  ),
  SES.PD = list(
    color = "#7f7f7f",
    description = "Standardized Effect Size (SES)for PD. Negative values suggest lower diversity than expected (clustering), while positive values indicate higher diversity (overdispersion)."
  ),
  PhyloEndemismWeighted = list(
    color = "#bcbd22",
    description = "Weighted phylogenetic endemism measures spatial uniqueness by summing, for each branch at a site, its length multiplied by the inverse of its range."
  ),
  PhyloEndemismStrict = list(
    color = "#17becf",
    description = "Strict measure of phylogenetic endemism (the total amount of branch length found only in this area)."
  ),
  ConservationValue = list(
    color = "#000000",
    description = "Indicator of the conservation value represented as the average phylogenetic diversity per grid cell, reflecting the spatial concentration of evolutionary history (as per Cadotte et al. 2010)."
  )
)

cat("..Beutifying results\n")

## Combine the raw data and metadata into one list
DIVB <- list(
  data = DIVN,
  metadata = metadata
)

## Convert the output to JSON
DIVB_json <- toJSON(DIVB, pretty = TRUE, auto_unbox = TRUE)
write(DIVB_json, file = paste0(RESULTS, "_diversity.json"))

