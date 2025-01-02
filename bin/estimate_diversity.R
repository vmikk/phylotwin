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
load_pckg("phyloregion")
load_pckg("canaper")
load_pckg("sf")
load_pckg("h3")
load_pckg("ape")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")
load_pckg("future")

cat("\n Parsing command line arguments\n")

## Define the option parser
option_list <- list(
    ## Input-output parameters
    make_option(c("-i", "--input"),  type = "character", default = NA, help = "File with aggregated species occurrences (Parquet, long format)"),
    make_option(c("-t", "--tree"),   type = "character", default = NA, help = "Phylogenetic tree (Newick format)"),
    make_option(c("-d", "--div"),    type = "character", default = NA, help = "Diversity metrics (comma-separated list)"),
    make_option(c("-o", "--output"), type = "character", default = NA, help = "Output prefix"),
    make_option(c("-r", "--randomizations"), type = "integer", default = 100, help = "Number of randomizations"),
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


## Parameters validation
if(is.na(opt$input)) { cat("Input is not specified.\n", file=stderr()); stop() }
if(is.na(opt$output)){ cat("Output is not specified.\n", file=stderr()); stop() }
if(is.na(opt$tree))  { cat("No phylogenetic tree specified.\n", file=stderr()); stop() }

## Default diversity metrics
if(is.null(opt$div)){ opt$div <- "PD,SES.PD" }


## Input parameters
OCC     <- opt$input
TREE    <- opt$tree
DIV     <- opt$div
RANDOMIZATIONS <- as.integer(opt$randomizations)
THREADS <- opt$threads
OUTPUT  <- opt$output

cat("\nParameters parsed:\n")
cat("  Input:", OCC, "\n")
cat("  Output prefix:", OUTPUT, "\n")
cat("  Tree:", TREE, "\n")
cat("  Diversity metrics:", DIV, "\n")
cat("  Number of randomizations:", RANDOMIZATIONS, "\n")
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
  cat("..Estimating PD (phylogenetic diversity)\n")
  RES <- c(RES, list(PD = PhyloMeasures::pd.query(tree = tree, matrix = datm) ))
}
if("MPD" %in% MEASURES){
  cat("..Estimating MPD (mean pairwise distance)\n")
  RES <- c(RES, list(MPD = PhyloMeasures::mpd.query(tree = tree, matrix = datm) ))
}
if("MNTD" %in% MEASURES){
  cat("..Estimating MNTD (mean nearest taxon distance)\n")
  RES <- c(RES, list(MNTD = PhyloMeasures::mntd.query(tree = tree, matrix = datm) ))
}

if("SES.PD" %in% MEASURES){
  cat("..Estimating SES.PD (standardized effect size of PD)\n")
  RES <- c(RES, list(SES.PD = PhyloMeasures::pd.query(tree = tree, matrix = datm, standardize = T) ))
}
if("SES.MPD" %in% MEASURES){
  cat("..Estimating SES.MPD (standardized effect size of MPD)\n")
  RES <- c(RES, list(SES.MPD = PhyloMeasures::mpd.query(tree = tree, matrix = datm, standardize = T) ))
}
if("SES.MNTD" %in% MEASURES){
  cat("..Estimating SES.MNTD (standardized effect size of MNTD)\n")
  RES <- c(RES, list(SES.MNTD = PhyloMeasures::mntd.query(tree = tree, matrix = datm, standardize = T) ))
}


######## phyloregion-based metrics


if("WeightedEndemism" %in% MEASURES){
  cat("..Estimating WeightedEndemism (weighted endemism)\n")
  RES <- c(RES, list(WeightedEndemism = phyloregion::weighted_endemism(x = datm) ))
}

if("PhyloEndemismWeighted" %in% MEASURES){
  cat("..Estimating PhyloEndemismWeighted (phylogenetic endemism weighted)\n")
  RES <- c(RES, list(PhyloEndemismWeighted = phyloregion::phylo_endemism(x = datm, phy = tree, weighted = TRUE) ))
}
if("PhyloEndemismStrict" %in% MEASURES){
  cat("..Estimating PhyloEndemismStrict (phylogenetic endemism strict)\n")
  RES <- c(RES, list(PhyloEndemismStrict = phyloregion::phylo_endemism(x = datm, phy = tree, weighted = FALSE) ))
}


########## CANAPE

if("CANAPE" %in% MEASURES){
  cat("..Estimating CANAPE (Categorical Analysis of Neo- And Paleo-Endemism)\n")

  cat("... Preparing data\n")
  dattt <- copy(datt)
  setnames(dattt, 
    old = colnames(dattt)[ ! colnames(dattt) %in% "H3" ],
    new = paste0("sp_", colnames(dattt)[ ! colnames(dattt) %in% "H3" ]) )
  setDF(dattt)
  rownames(dattt) <- dattt$H3
  dattt$H3 <- NULL

  treee <- tree
  treee$tip.label <- paste0("sp_", treee$tip.label)

  cat("... Generating a set of random communities\n")
  rnd <- cpr_rand_test(
    comm = dattt,
    phy = treee,
    metrics = c("pe", "rpe"),    # "pd", "rpd", 
    n_reps = RANDOMIZATIONS,
    null_model = "curveball",    # sequential algo for binary data (Strona et al., 2014), preserves row frequencies
    n_iterations = 100000,       # for sequential null models
    thin = 1,
    tbl_out = TRUE)


  cat("... Classifying endemism\n")
  canape_res <- canaper::cpr_classify_endem(df = rnd)
  if(any(canape_res$site != datt$H3)){
    stop("Site names are not matching --- a fix should be implemented in the future\n")
  }
  RES <- c(RES, list(CANAPE = canape_res$endem_type ))
}







######## Combine results

cat("Combining results\n")
RES <- do.call("cbind", RES)
RES <- data.table(H3 = rownames(datm), RES)

cat("..Adding number of records per H3 cell\n")
RES <- merge(x = RES, y = num_records, by = "H3", all.x = TRUE)

### Estimate sampling redundancy for each cell (how well the are is sampled) 
###  = 1 - (Richness / Number of specimens)
### (see Mishler et al., 2020; DOI: 10.1111/jse.12590) 
cat("..Estimating redundancy index\n")
RES[, Redundancy := ( 1 - (Richness / NumRecords) )]





##########################################################
########################################################## Prepare spatial data
##########################################################

cat("\n\n-------- Preparing spatial data --------\n\n")

cat("..Preparing H3 polygons\n")
H3_poly <- h3_to_geo_boundary_sf(RES$H3)

cat("..Adding diversity estimates to H3 polygons\n")
vars <- colnames(RES)[! colnames(RES) %in% "H3" ]
H3_poly <- cbind(H3_poly, RES[, ..vars])     

cat("..Exporting polygons with divsity estimates in GeoPackage format\n")
st_write(
  obj   = H3_poly,
  dsn   = paste0(OUTPUT, ".gpkg"),
  layer = "diversity_estimates")      


##########################################################
########################################################## Export results
##########################################################

cat("\n\n-------- Exporting results --------\n\n")

cat("Exporting QS file\n")
qs::qsave(RES, paste0(OUTPUT, ".qs"), preset = "custom", algorithm = "zstd", compress_level = 14L, nthreads = THREADS)

cat("Exporting tab-delimited file\n")

## Add grid cell coordinates
cat("..Adding geo-coordinates for grid cell centers\n")
RES[, c("Latitude", "Longitude") := as.data.table(h3::h3_to_geo(H3)) ]
 
## Reorder columns
cat("..Reordering columns\n")
setcolorder(RES, c("H3", "Latitude", "Longitude"))

cat("..Exporting tab-delimited file\n")
fwrite(x = RES, file = paste0(OUTPUT, ".txt"), sep = "\t")

cat("..Exporting a list of estimated diversity indices\n")
inds <- colnames(RES)[! colnames(RES) %in% c("H3", "Latitude", "Longitude") ]
fwrite(
  x = data.table(DivIndices = inds),
  file = "Div_indices.txt",
  sep = "\t", col.names = FALSE)

