#!/usr/bin/env Rscript

## Estimate diversity


## Notes:
# - Richness, NumRecords, and Redundancy index are always estimated

cat("Estimating diversity\n\n")

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
load_pckg("ape")
load_pckg("arrow")
load_pckg("dplyr")
load_pckg("qs")
# load_pckg("canaper")
# load_pckg("future")  # currently, used only for CANAPE

# load_pckg("SpadeR")
# load_pckg("BiodiverseR")     # https://github.com/biogeospatial/BiodiverseR

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
########################################################## For debugging
##########################################################

# OCC  <- "aggregated_counts.parquet"
# TREE <- "Ferns_FTOL_1-7-0.nwk.gz"
# DIV  <- "PD,MNTD,SES.PD,RPD,PhyloEndemismWeighted,CANAPE"
# RANDOMIZATIONS <- 20
# THREADS <- 4
# OUTPUT <- "diversity_estimates"


##########################################################
########################################################## Load and reshape data
##########################################################

cat("\n\n-------- Loading data --------\n\n")

## Load phylogenetic tree
cat("Loading phylogenetic tree: ", TREE, "\n")

## Function to load phylogenetic tree
load_tree <- function(tree_file){
  ## Check if file is gzip-compressed and get the base file name
  is_gz <- grepl("\\.gz$", tree_file, ignore.case = TRUE)
  base_file <- if (is_gz) sub("\\.gz$", "", tree_file, ignore.case = TRUE) else tree_file
  
  ## Extract the file extension (in lower case)
  ext <- tolower(tools::file_ext(base_file))

  ## Load tree in Nexus format
  if (ext %in% c("nex", "nexus")) {
    TRE <- ape::read.nexus(file = tree_file)
  } else if (ext %in% c("tre", "nwk", "newick")) {
    TRE <- ape::read.tree(file = tree_file)
  } else {
    stop("Unsupported file extension: ", ext)
  }
  return(TRE)
}

tree <- load_tree(TREE)

## Load aggregated species occurrences
cat("Loading aggregated species occurrences\n")
occ <- try( read_parquet(OCC) )

if("try-error" %in% class(occ) || nrow(occ) == 0){
  cat("..NO SPECIES OCCURRENCES FOUND; exiting with code 99\n")
  quit(status = 99)
}

## Convert to data.table
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
  formula = H3 ~ species,
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

if("RPD" %in% MEASURES){
  cat("..Estimating RPD (relative phylogenetic diversity)\n")
  RES <- c(RES, list(RPD = phyloregion::RPD(x = datm, phy = tree) ))
}

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


########## CANAPE (canaper-based) -- in the pipeline, we use "native" Biodiverse implementation

# if("CANAPE" %in% MEASURES){
#   cat("..Estimating CANAPE (Categorical Analysis of Neo- And Paleo-Endemism)\n")
# 
#   if(nrow(datt) <= 5){
#     cat("... WARNING: randomizations might be unreliable for small datasets\n")
#   }
# 
#   ## Enable parallel backend for some metrics
#   if(THREADS > 1){
#     cat("... Enabling parallel backend\n")
#     plan(multicore, workers = THREADS)
#   }
# 
#   cat("... Preparing data\n")
#   dattt <- copy(datt)
#   setnames(dattt, 
#     old = colnames(dattt)[ ! colnames(dattt) %in% "H3" ],
#     new = paste0("sp_", colnames(dattt)[ ! colnames(dattt) %in% "H3" ]) )
#   setDF(dattt)
#   rownames(dattt) <- dattt$H3
#   dattt$H3 <- NULL
# 
#   treee <- tree
#   treee$tip.label <- paste0("sp_", treee$tip.label)
# 
#   cat("... Generating a set of random communities\n")
#   rnd <- cpr_rand_test(
#     comm = dattt,
#     phy = treee,
#     metrics = c("pe", "rpe"),    # "pd", "rpd", 
#     n_reps = RANDOMIZATIONS,
#     null_model = "curveball",    # sequential algo for binary data (Strona et al., 2014), preserves row frequencies
#     n_iterations = 100000,       # for sequential null models
#     thin = 1,
#     tbl_out = TRUE)
# 
#   ## Switch back to non-parallel mode
#   if(THREADS > 1){
#     cat("... Switching back to non-parallel mode\n")
#     plan(sequential)
#   }
# 
#   cat("... Classifying endemism\n")
#   canape_res <- canaper::cpr_classify_endem(df = rnd)
#   if(any(canape_res$site != datt$H3)){
#     stop("Site names are not matching --- a fix should be implemented in the future\n")
#   }
#   RES <- c(RES, list(CANAPE = canape_res$endem_type ))
# }





######## Combine results

cat("Combining results\n")


RES <- do.call(what = data.table:::cbind.data.table, RES)
RES[ , H3 := rownames(datm) ]
setcolorder(RES, "H3")

cat("..Adding number of records per H3 cell\n")
RES <- merge(x = RES, y = num_records, by = "H3", all.x = TRUE)

### Estimate sampling redundancy for each cell (how well the are is sampled) 
###  = 1 - (Richness / Number of specimens)
### (see Mishler et al., 2020; DOI: 10.1111/jse.12590) 
cat("..Estimating redundancy index\n")
RES[, Redundancy := ( 1 - (Richness / NumRecords) )]




##########################################################
########################################################## Species originality and range
##########################################################

cat("\n\n-------- Estimating species originalities and ranges --------\n\n")

## Estimate FP index (fair proportion). Sum(FP) = PD
cat("..Estimating FP index (fair proportion)\n")
FP <- phyloregion::evol_distinct(tree, type = "fair.proportion")
FP <- data.table(species = names(FP), FP = FP)

## Estimate the inverse of the range size of each taxon. Sum(RR) = Weighted endemism
cat("..Estimating the inverse of the range size of each taxon\n")
RR <- occ[ , .(RangeSize = .N), by = "species" ]
RR[ , RangeSizeInv := 1 / RangeSize ]

## Combine data
SPEC <- merge(x = FP, y = RR, by = "species", all = TRUE)

## Export data
cat("..Exporting species originalities and ranges\n")
fwrite(x = SPEC, file = "Species_originalities_and_ranges.txt", sep = "\t")

## Report top-N species
if(!is.na(TOPNSP) && TOPNSP > 0){
  cat("..Reporting top-", TOPNSP, "species per grid cell\n")
  
  ## Add species originalities and ranges to occurrences
  occ_spec <- merge(x = occ, y = SPEC, by = "species", all.x = TRUE)

  ## FP
  setorder(occ_spec, H3, -FP)
  top_fp  <- occ_spec[, head(.SD, TOPNSP), by = "H3", .SDcols = "species"]
  top_fpm <- top_fp[ , .(MostPhylogeneticallyOriginal = paste0(species, collapse = ", ")), by = "H3" ]

  ## Ranges
  setorder(occ_spec, H3, -RangeSizeInv)
  top_rr  <- occ_spec[, head(.SD, TOPNSP), by = "H3", .SDcols = "species"]
  top_rrm <- top_rr[ , .(MostRangeRestricted = paste0(species, collapse = ", ")), by = "H3" ]

  ## Combine
  top_sp <- merge(x = top_fpm, y = top_rrm, by = "H3", all = TRUE)

  ## Export originalities and ranges into a separate file
  # fwrite(x = top_sp, file = "Top_species_per_grid_cell.txt", sep = "\t")

  ## Add top species to diversity estimates
  cat("..Adding top species to diversity estimates\n")
  RES <- merge(x = RES, y = top_sp, by = "H3", all.x = TRUE)

}

# Top driving species in phyloregions
# phyloregion::indicators


##########################################################
########################################################## Export results
##########################################################

cat("\n\n-------- Exporting results --------\n\n")

cat("Exporting QS file\n")
qs::qsave(RES, paste0(OUTPUT, ".qs"),
  preset = "custom", algorithm = "zstd", compress_level = 8L, nthreads = THREADS)


# cat("..Exporting a list of estimated diversity indices\n")
# inds <- colnames(RES)[! colnames(RES) %in% c("H3", "Latitude", "Longitude") ]
# fwrite(
#   x = data.table(DivIndices = inds),
#   file = "Div_indices.txt",
#   sep = "\t", col.names = FALSE)

