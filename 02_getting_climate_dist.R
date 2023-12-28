# rm(list=ls())

setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(data.table)
library(ape)
library(phytools)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

source("00_utility_functions.R")

#-------------------------------
# Getting climate data
#-------------------------------
# Load cleaned points:
# This file is very large, so it's on gitignore (not in the repo)
# Data is available at:
# Dorey, et al. (2023). A globally synthesised and flagged bee occurrence dataset and cleaning workflow. Scientific Data, 10(1), 747.
all_cleaned_points <- fread("OutputData/05_cleaned_database.csv")

# Filter according to species that are in the phylogeny (i.e. keeping only those that are sampled)
bee_tree <- read.tree("original_data/ML_beetree.tre")
all_cleaned_points$species <- gsub(" ","_",all_cleaned_points$species)
all_cleaned_points <- subset(all_cleaned_points, all_cleaned_points$species %in% bee_tree$tip.label)

# 1. Thinning occurrence data (this will keep only a few records per grid cell; it's a workaround collecting bias)
thinned_points <- Thinning(all_cleaned_points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1)

# Added step! Remove suspicious and exotic species:
colnames(thinned_points) <- c("species","lat","lon")
exotic <- read.csv("original_data/exotic_species_to_remove.csv")[,1] # table provided by Aline after visual inspection of points
suspicious <- read.csv("original_data/species_with_suspicious_distribution.csv")[,1] # table provided by Aline after visual inspection of points
remove <- c(exotic, suspicious)
thinned_points <- subset(thinned_points, !thinned_points$species %in% remove)

# Saving filtered points
save(thinned_points, file="curated_data/thinned_points_res1.Rsave")

# load("curated_data/thinned_points_res1.Rsave")

#--------------------------------------------

# 2. Now getting summary statistics of climatic variables for each species
all_layers <- list.files("env_layers", paste(c(".tif$",".bil$"),collapse = "|")) 
# note: env_layers is also on gitignore for now due to the size of the layer files
labels <- gsub(paste(c(".tif$",".bil$"),collapse = "|"),"", all_layers)
all_layers <- lapply(paste0("env_layers/",all_layers), raster)
names(all_layers) <- labels

# Now extracting data and saving curated datasets:
for(i in 1:length(all_layers)){
  one_layer <- all_layers[[i]]
  one_label <- names(all_layers)[i]
  allpoints <- DataFromPoints(thinned_points, one_layer)
  #write.csv(allpoints, file=paste0("climate_data/",one_label,"_allpoints.csv"), row.names=F)
  summstats <- GetClimateSummStats_custom(allpoints, type="raw")
  write.csv(summstats, file=paste0("curated_data/",one_label,"_climate_summstats.csv"), row.names=F)
  cat(one_label, "done.", "\n")
}

