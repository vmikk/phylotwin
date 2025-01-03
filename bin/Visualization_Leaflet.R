#!/usr/bin/env Rscript

## Visualization of diversity estimates using intaractive maps (Leaflet-based)

cat("Leaflet-based visualization of biodiversity estimates\n")


## Usage example:
# Visualization_Leaflet.R \
#   --divestimates ${divestimates} \
#   --variables    ${index} \
#   --shortid      TRUE \
#   --antimeridianfix TRUE

#   --observed  "RND_SPATIAL_RESULTS.csv" \
#   --sesscores "RND_rand--z_scores--SPATIAL_RESULTS.csv" \
#   --sigscores "RND_rand--SPATIAL_RESULTS.csv" \
#   --canape    "RND_rand--CANAPE--.csv" \
#   --hurlbert  "RND_HURLBERT_ES.csv" \
#   --reccounts "Record_counts_H3.RData" \
#   --variables "RICHNESS_ALL,PD,SES_PD,PD_P,ENDW_WE,SES_ENDW_WE,PE_WE,SES_PE_WE,CANAPE,Redundancy" \
#   --isolate_maps FALSE \
#   --palette "quantile" \
#   --color "RdYlBu" \
#   --bins 5 \
#   --output "03.Plots/Choropleth.html"

## For a list of indices available in Biodiverse, see:
# https://github.com/shawnlaffan/biodiverse/wiki/IndicesDevVersion

## For Standardized-Effect-Size-based variables, add `SES_` prefix (e.g., `SES_PD`)
## For Hurlbert ES indices, add `ES_` prefix (e.g., `ES_50`)

## For CANAPE (categorical analysis of neo- and paleoendemism; Mishler et al., 2014), add `CANAPE` to the variables list
# CANAPE is able to distinguish different types of centres of endemism, and can thus give insights 
# into different evolutionary and ecological processes that may be responsible for these patterns. 
# - The centres of paleo-endemism indicate places where there are over-representation of long branches that are rare across the landscape.
# - The centres of neo-endemism indicate an area where there is an over-representation of short branches that are rare on the landscape.
# - Mixture of both paleo-endemism and neo-endemism
# - Super-endemic sites

## H3 grid cell index is represented as a 15-character hexadecimal string (e.g., `830021fffffffff`)

## Output:
# - Choropleth.html = Leaflet-based interactive map
# - Biodiverse_results_merged.txt = Diversity metrics (per grid cell) in a tabular format
# - Diversity_estimates.gpkg = Polygons with diversity metrics in GeoPackage format

## TO DO:
# - SES for Hurlbert ES indices?
# - handle NA values (e.g., if Richness == 1)


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(".. ", paste(pkg, packageVersion(pkg), "\n"))
}

load_pckg("optparse")

cat("\n Parsing command line arguments\n")

## Parse arguments
option_list <- list(
  make_option(c("-d", "--divestimates"), action="store", default=NA,   type='character', help="Input file with diversity estimates (tab-delimited or QS)"),
  make_option(c("-v", "--variable"),     action="store", default="PD", type='character', help="Name of the diversity variable to plot"),
  make_option("--variable2",             action="store", default=NULL, type='character', help="Name of the second diversity variable to plot"),
  make_option(c("-j", "--redundancy"),   action="store", default=0,    type='double',    help="Redundancy threshold for hiding the grid cells with low number of records (disabled by default)"),
  make_option(c("--shortid"), action="store", default=TRUE, type='logical', help="Shorten H3 index name of grid cell labels on the map"),
  make_option(c("--antimeridianfix"), action="store", default=TRUE, type='logical', help="Fix H3 polygons that cross the antimeridian"),
  make_option(c("--saveqs"), action="store", default=FALSE, type='logical', help="Save the Leaflet data in QS format")
)
opt <- parse_args(OptionParser(option_list=option_list))

## Function to convert text "NA"s to NA
to_na <- function(x){
  if(x %in% c("NA", "null", "Null")){ x <- NA }
  return(x)
}

## Replaces "null"s from Nextflow with NA
opt <- lapply(X = opt, FUN = to_na)


cat("\n-------- Loading packages --------\n")
load_pckg("data.table")
load_pckg("sf")
load_pckg("h3")
load_pckg("plyr")
load_pckg("leaflet")
if(!is.null(opt$variable2)){ load_pckg("leaflet.extras2") }
load_pckg("mapview")
load_pckg("qs")


## Validation of the required argiments
cat("\n-------- Parameters validation --------\n")
if(is.na(opt$divestimates)){
  stop("Input file with diversity estimates is not specified.\n")
}

## Assign variables
INPUT          <- opt$divestimates
VARIABLE       <- opt$variable
VARIABLE2      <- opt$variable2
REDUNDANCYTRSH <- as.numeric(to_na( opt$redundancy ))
SHORTID        <- as.logical( opt$shortid )
ANTIFIX        <- as.logical( opt$antimeridianfix )
SAVEQS         <- as.logical( opt$saveqs )

## Hardcoded parameters
PALETTE     <- "quantile"
COLOR       <- "RdYlBu"
BINS        <- 5
COLORSES    <- "threat"

## Check the redundancy range
if(!is.na(REDUNDANCYTRSH)){
  if(! (0 <= REDUNDANCYTRSH & REDUNDANCYTRSH <= 1) ){
    stop("Redundacy threshold should be in the [0,1] range.\n")
  }
} else {
  REDUNDANCYTRSH <- 0
}

## Log assigned variables
cat(paste("Diversity estimates: ",          INPUT,          "\n", sep=""))
cat(paste("Index to plot: ",                VARIABLE,       "\n", sep=""))
if(!is.null(VARIABLE2)){
  cat(paste("Second index to plot: ",       VARIABLE2,      "\n", sep=""))
}
cat(paste("Redundancy threshold: ",         REDUNDANCYTRSH, "\n", sep=""))
cat(paste("Display short H3 index names: ", SHORTID,        "\n", sep=""))
cat(paste("Antimeridian fix: ",             ANTIFIX,        "\n", sep=""))
cat(paste("Save QS data: ",                 SAVEQS,         "\n", sep=""))


##########################################################
########################################################## Debug
##########################################################

# INPUT          <- "diversity_estimates.qs"
# VARIABLE       <- "PD"
# VARIABLE2      <- "Redundancy"  # NULL
# REDUNDANCYTRSH <- 0
# SHORTID        <- TRUE
# ANTIFIX        <- TRUE
# SAVEQS         <- TRUE
# PALETTE        <- "quantile"
# COLOR          <- "RdYlBu"
# BINS           <- 5
# COLORSES       <- "threat"

##########################################################
########################################################## Load and prepare data
##########################################################

cat("\n\n-------- Preparing data --------\n\n")

## If there are multiple variables selected - split them
if(any(grepl(pattern = ",", x = VARIABLES))){
  VARIABLES <- unique( strsplit(x = VARIABLES, split = ",")[[1]] )
  if(length(VARIABLES) > 2){
    stop("More than two variables selected (currently not supported)\n")
  }
}

## Load input data
cat("..Loading diversity estimates\n")

tab_ext <- tolower( tools::file_ext(INPUT) )
if(tab_ext %in% "qs"){
  res <- qs::qread(INPUT)
} else if(tab_ext %in% c("txt", "tsv", "txt.gz", "tsv.gz", "gz")){
  res <- fread(INPUT, sep = "\t")
} else {
  cat("The specified diversity estimates file has an unsupported extension.\n", file=stderr())
  cat("Supported extensions (case-insensitive): ", paste(c("qs", "txt"), collapse = ", "), "\n")
  stop()
}

## Check if the selected index is in the tables
if(any(!VARIABLES %in% colnames(res))){
  cat("The selected indices are not present in the diversity estimates table!\n")
  cat("Please check the spelling of index names\n")
  stop()
}

## Subset data
cat("..Subsetting data\n")
clz <- unique(c("H3", "Redundancy", VARIABLES))
clz <- clz[ clz %in% colnames(res) ]
res <- res[, ..clz]

## Remove cells with the redundancy index below the specified threshold
if(REDUNDANCYTRSH > 0 & "Redundancy" %in% colnames(res)){
  cat("Removing grid cells with low redundancy index\n")
  cat("..The specified threshold is ", REDUNDANCYTRSH, "\n")
  cat("..There are ", sum(res$Redundancy <= REDUNDANCYTRSH), " grid cells to be removed\n")
  cat("..And ", sum(res$Redundancy > REDUNDANCYTRSH), " grid cells will be preserved\n")

  if(any(res$Redundancy <= REDUNDANCYTRSH)){
    res <- res[ Redundancy > REDUNDANCYTRSH ]
  }
}

## Check the data
if(!nrow(res) > 0){
  stop("There are no grid cells to display!\n")
}


##########################################################
########################################################## Visualization
##########################################################

cat("\n\n-------- Visualization --------\n\n")

## Prepare spatial polygons
cat("..Preparing spatial data\n")
H3_poly <- h3_to_geo_boundary_sf(res$H3)

cat("..Adding diversity estimates to polygons\n")
vars <- colnames(res)[! colnames(res) %in% "H3" ]
H3_poly <- cbind(H3_poly, res[, ..vars])     

## Assign rownames as H3 grid cell IDs
cat("..Adding H3 cell names\n")
if(SHORTID == TRUE){
  ## Truncate fff-tail of index name (e.g., `830021fffffffff` -> `830021f`)
  rownames(H3_poly) <- gsub(pattern = "f+$", replacement = "f", x = res$H3)
} else {
  rownames(H3_poly) <- res$H3
}

## Fix H3 polygons that cross the antimeridian by cutting them in two
if(ANTIFIX == TRUE){
  cat("..Fixing antimeridian issue\n")
  H3_poly <- st_wrap_dateline(H3_poly, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
}

# plot(H3_poly)


############################################## Leaflet choropleths

cat("..Creating leaflet map\n")

## Set value for the minZoom and maxZoom settings
leaflet(options = leafletOptions(minZoom = 0, maxZoom = 18))

## Function to create a label with single variable
single_label <- function(num, name){

  if( is.numeric(num) ){
    res <- sprintf(
      "<strong>%s:</strong> %.3g<br/>",
      rep(name, times = length(num)),
      num
      )
  } else {
    res <- sprintf(
      "<strong>%s:</strong> %s<br/>",
      rep(name, times = length(num)),
      num
      )
  }
  return(res)
}
# single_label(H3_poly$PD, name = "PD")           # numeric labels
# single_label(H3_poly$CANAPE, name = "CANAPE")   # factor labels

## Create labels for all variables
cat("..Generating polygon labels\n")
labels <- alply(.data = VARIABLES, .margins = 1, .fun = function(v){
  single_label(num = H3_poly[[ v ]], name = v)
  })

## Add H3 index to the labels
labels <- c(labels, 
  list( single_label(num = rownames(H3_poly), name = "H3 index") )
  )


## Concatenate labels from all variables
labels <- do.call(paste0, labels)

## Mark labels as HTML
labels <- labels %>% lapply(htmltools::HTML)


## Create a map, add OpenStreetMap map tiles as background
cat("..Building basemap\n")
m <- leaflet() %>% addTiles()


## Color palette
# pal <- colorQuantile(palette = "RdYlBu", domain = H3_poly$PD, n = 5, reverse = TRUE)
#
# bins <- c(0, 10, 20, 50, 100, 200, 500, 1000, Inf)
# pal <- colorBin("YlOrRd", domain = H3_poly$PD, bins = 7, reverse = TRUE) # , bins = bins)
#
# pal <- colorNumeric(palette = "YlGnBu", domain = H3_poly$PD, reverse = TRUE)



cat("..Preparing color palettes\n")

## Function to create color palette
gen_color_palette <- function(x, type = "quantile", col = "RdYlBu", colses = "threat", nbins = 5, rev = TRUE){

  ## Bin numeric data via the quantile function
  if(type %in% "quantile"){

    ## To fix `'breaks' are not unique` error, check the number of potential bins in advance
    testbins <- quantile(x,
      probs = seq(0, 1, length.out = nbins + 1),
      na.rm = TRUE, names = FALSE)
    
    nbins <- length(unique(testbins)) - 1

    pal <- colorQuantile(palette = col, domain = x, n = nbins, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- nbins
  }

  ## Equal-interval binning (based on `cut` function)
  if(type %in% "equal"){
    pal <- colorBin(palette = col, domain = x, bins = nbins, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- nbins
  }

  ## Simple linear mapping from continuous numeric data to an interpolated palette
  if(type %in% "continuous"){
    pal <- colorNumeric(palette = col, domain = x, reverse = rev, na.color = "#808080")
    attr(pal, which = "newbins") <- 999
  }

  ## SES-score mapping (symmetric around zero)
  if(type %in% "ses"){

    ## Hotspot-type colors (high SES values = red, low = blue)
    if(colses %in% "hotspots"){
    ses_colors <- c(
      HighlyNegative = "#27408B",
      Negative       = "#4876FF",
      NotSignificant = "#FAFAD2",
      Positive       = "#FF0000",
      HighlyPositive = "#8B0000"
      )
    }

    ## Threat-type colors (high SES values = blue, low = red), as in Mishler et al., 2014
    if(colses %in% "threat"){
    ses_colors <- c(
      HighlyNegative = "#8B0000",
      Negative       = "#FF0000",
      NotSignificant = "#FAFAD2",
      Positive       = "#4876FF",
      HighlyPositive = "#27408B"
      )
    }

    ses_bins <-c(-1000, -2.58, -1.96, 1.96, 2.58, 1000)

    pal <- colorBin(palette = ses_colors, bins = ses_bins,
      na.color = "#808080", domain = c(-1000, 1000))
    
    attr(pal, which = "newbins") <- 5
    # plot(-4:4, col = pal(-4:4), cex = 2, pch = 16) # test
  }

  ## CANAPE-style mapping (neo / paleo / mixed / super)
  if(type %in% "canape"){

    canape_colors <- c(
      Neo_endemism   = "#FF0000",
      Paleo_endemism = "#4876FF",
      NotSignificant = "#FAFAD2",
      Mixed_endemism = "#CB7FFF",
      Super_endemism = "#9D00FF"
      )

    pal <- colorFactor(
      palette = canape_colors,
      levels = names(canape_colors),
      na.color = "#808080")

    attr(pal, which = "newbins") <- 5
    # plot(1:length(canape_colors), col = pal(names(canape_colors)), cex = 2, pch = 16) # test
  }

  ## Redundancy-style mapping (yellow-brown, [0,1])
  if(type %in% "redundancy"){
    pal <- colorNumeric(palette = "YlOrBr", domain = seq(0, 1, 0.01), reverse = FALSE, na.color = "#808080")
    attr(pal, which = "newbins") <- 999
    # plot(seq(0,1,0.05), col = pal(seq(0,1,0.05)), cex = 2, pch = 16) # test
  }

  return(pal)
}
# e.g., gen_color_palette(1:10, type = "quantile", nbins = 5)
#       gen_color_palette(1:10, type = "equal",    nbins = 5)
#       gen_color_palette(1:10, type = "continuous")

## Make color palette for all variables
pals <- list()

VARIABLES_ses <- grep(pattern = "^SES_", x = VARIABLES, value = TRUE)
VARIABLES_raw <- grep(pattern = "^SES_", x = VARIABLES, value = TRUE, invert = TRUE)
VARIABLES_raw <- VARIABLES_raw[ ! VARIABLES_raw %in% c("CANAPE", "Redundancy") ]

## Colors for "raw" variables
if(length(VARIABLES_raw) > 0){

  pals_raw <- alply(.data = VARIABLES_raw, .margins = 1,
    .fun = function(v, ...){ gen_color_palette(x = H3_poly[[ v ]], ...) }, 
    type = PALETTE, col = COLOR, nbins = BINS, rev = TRUE)

  names(pals_raw) <- VARIABLES_raw
  pals <- c(pals, pals_raw)
  rm(pals_raw)
}


## Colors of SES-scores
if(length(VARIABLES_ses) > 0){

  pals_ses <- alply(.data = VARIABLES_ses, .margins = 1,
    .fun = function(v, ...){ gen_color_palette(x = H3_poly[[ v ]], ...) }, 
    type = "ses", colses = COLORSES)

  names(pals_ses) <- VARIABLES_ses
  pals <- c(pals, pals_ses)
  rm(pals_ses)
}


cat("..Adding polygons\n")

# ## Add PD polygons to the map
# m <- m %>% 
#   addPolygons(data = H3_poly,
#     fillColor = ~pal(PD),
#     group = "PD",
#     opacity = 0.8,
#     fillOpacity = 0.8,
#     weight = 0.3, color = "white", dashArray = "1",
#     highlightOptions = highlightOptions(
#       weight = 2, color = "#777777", dashArray = "1",
#       opacity = 0.8, bringToFront = TRUE),
#     label = labels,
#     labelOptions = labelOptions(
#       style = list("font-weight" = "normal", padding = "3px 8px"),
#       textsize = "10px", direction = "auto")
#     ) %>%
#   addLegend("bottomright", pal = pal, values = H3_poly$PD,
#     title = "PD", group = "PD",  opacity = 1
#     )



## Shortcut to add polygons with legend to the map
add_polygons_with_legend <- function(m, v, pal, ses_labels = FALSE){
  res <- m %>% 
      addPolygons(data = H3_poly,
        fillColor = ~ pal( H3_poly[[v]] ),
        group = v,
        opacity = 0.8,
        fillOpacity = 0.8,
        weight = 0.3, color = "white", dashArray = "1",
        highlightOptions = highlightOptions(
          weight = 2, color = "#777777", dashArray = "1",
          opacity = 0.8, bringToFront = TRUE),
        label = labels,
        labelOptions = labelOptions(
          style = list("font-weight" = "normal", padding = "3px 8px"),
          textsize = "10px", direction = "auto")
      ) %>%
      addLegend("bottomright", pal = pal, values = H3_poly[[v]],
        title = v, group = v,  opacity = 1)

  ## Relable upper bounds of effect sizes
  ## (unfortunately, addLegend(labels =...) does not work with `pal` argument)
  if(ses_labels == TRUE){
    
    newlabs <- c("< -2.58", "< -1.96", "-1.96 - 1.96", "> 1.96", "> 2.58")
    newlabs <- gsub(pattern = " - ", replacement = " &ndash; ", x = newlabs)

    slotid <- length(res$x$calls)
    if(
      "labels" %in% names(res$x$calls[[ slotid ]]$args[[1]]) &                # `labels` slot is present
      length(res$x$calls[[ slotid ]]$args[[1]]$labels) == length(newlabs)     # and has the same number of categories
      ){
      res$x$calls[[ slotid ]]$args[[1]]$labels <- newlabs
    }
  }

  return(res)
}
## Example:
# add_polygons_with_legend(m, "PD", pal = pals[[ "PD" ]])


## Loop throug all variables
## Except "CANAPE" (it has a categorical color scheme) and "Redundancy"
for(v in VARIABLES[ ! VARIABLES %in% c("CANAPE", "Redundancy") ]){

  cat("... ", v, "\n")

  if(v %in% VARIABLES_raw & PALETTE %in% "quantile"){
    if(attr(pals[[v]], "newbins") != BINS){
    cat(".... number of bins was adjusted to ", attr(pals[[v]], "newbins"), "\n")
    }
  }

  ## SES lables?
  if(v %in% VARIABLES_ses){
    ses_lab <- TRUE
  } else {
    ses_lab <- FALSE
  }

  tmp <- try( add_polygons_with_legend(m = m, v = v, pal = pals[[v]], ses_labels = ses_lab) )

  ## Quantile palette may fail
  if(("try-error" %in% class(tmp) | attr(pals[[v]], "newbins") == 1) & PALETTE %in% "quantile"){
   cat(".... Warning: quantile palette failed, trying continuous palette.\n") 

   ## Generate new color palette (continuous)
   tmppal <- gen_color_palette(x = H3_poly[[ v ]],
    type = "continuous", col = COLOR, nbins = BINS, rev = TRUE)

   tmp <- add_polygons_with_legend(m = m, v = v, pal = tmppal)

   rm(tmppal)
  } # end of `try-error`

  m <- tmp
  rm(tmp, ses_lab)

}   # end of loop
rm(v)


## Add CANAPE to the plot
if("CANAPE" %in% VARIABLES){

  cat("... CANAPE\n")

  canape_pal <- gen_color_palette(
    x = H3_poly[[ "CANAPE" ]],
    type = "canape")

  m <- m %>% 
    addPolygons(data = H3_poly,
      fillColor = ~ canape_pal( H3_poly[[ "CANAPE" ]] ),
      group = "CANAPE",
      opacity = 0.8,
      fillOpacity = 0.8,
      weight = 0.3, color = "white", dashArray = "1",
      highlightOptions = highlightOptions(
        weight = 2, color = "#777777", dashArray = "1",
        opacity = 0.8, bringToFront = TRUE),
      label = labels,
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "10px", direction = "auto")
    ) %>%
    addLegend("bottomright", pal = canape_pal, values = H3_poly[[ "CANAPE" ]],
      title = "CANAPE", group = "CANAPE",  opacity = 1)

} # end of CANAPE


## Add Redundancy index to the plot
if("Redundancy" %in% VARIABLES){

  cat("... Redundancy\n")

  redundancy_pal <- gen_color_palette(
    x = H3_poly[[ "Redundancy" ]],
    type = "redundancy")

  m <- m %>% 
    addPolygons(data = H3_poly,
      fillColor = ~ redundancy_pal( H3_poly[[ "Redundancy" ]] ),
      group = "Redundancy",
      opacity = 0.8,
      fillOpacity = 0.8,
      weight = 0.3, color = "white", dashArray = "1",
      highlightOptions = highlightOptions(
        weight = 2, color = "#777777", dashArray = "1",
        opacity = 0.8, bringToFront = TRUE),
      label = labels,
      labelOptions = labelOptions(
        style = list("font-weight" = "normal", padding = "3px 8px"),
        textsize = "10px", direction = "auto")
    ) %>%
    addLegend("bottomright", pal = redundancy_pal, values = H3_poly[[ "Redundancy" ]],
      title = "Sampling redundancy", group = "Redundancy",  opacity = 1)

}


## Add variable selector
cat("..Adding variable selector\n")
m <- m %>%
  addLayersControl(
    overlayGroups = c(VARIABLES),
    options = layersControlOptions(collapsed = FALSE)
    )

## Hide all vars except the first one
cat("..Hiding variables\n")
m <- m %>% 
  hideGroup(VARIABLES[-1])


cat("..Exporting the results\n")

## Use mapdeck as the rendering platform instead of leaflet
# mapviewOptions(platform = "mapdeck")

cat("...Exporting Choropleth map in HTML\n")

## Save a leaflet map as .html
mapshot(
  m,
  url  = paste0("Choropleth_", opt$variables, ".html"),   # NB. use non-splitted variable name as provided by the user
  file = NULL,
  selfcontained = TRUE)

# remove_controls = c("zoomControl", "layersControl",
#    "homeButton", "scaleBar", "drawToolbar", "easyButton")

# library(htmlwidgets)
# saveWidget(m, file="m.html")


if(SAVEQS == TRUE){
  cat("...Exporting Leaflet object\n")

  attr(m, which = "VARIABLES") <- VARIABLES

  # saveRDS(object = m, file = "Leaflet_object.RData", compress = "xz")

  qs::qsave(m, paste0("Leaflet_", opt$variables, ".qs"),
    preset = "custom", algorithm = "zstd", compress_level = 14L)
}

cat("Plotting finished\n")
