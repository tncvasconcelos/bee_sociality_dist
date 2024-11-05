
#-------------------------------------------------------------------------------
#------------------------------------SETUP--------------------------------------
rm(list=ls())

# Libraries
library(sp)
library(sf)
library(raster)
library(ggplot2)
library(gridExtra)
library(maptools)
library(dplyr)
library(rnaturalearth)

# Directories
getwd()                                       
wd <- "/Users/lenarh/Desktop/bee_sociality_dist/proportion_heatmaps"
setwd(wd)
data_wd <- paste0(wd, "/data")
data_wd

# Source code
functions_path <- file.path(wd, "example code", "00_utility_functions_synthesis.R") 
source(functions_path)

# Load and check data
tree_data <- read.csv(file.path(data_wd, "bees_traits.csv")) 
head(tree_data)
colnames(tree_data)[3] <- "species" # to be consistent with all_bees_data
colnames(tree_data) 
nrow(tree_data)

# tree_data contains trait information (sociality, nesting) for phylogeny tip species

library(data.table)
all_bees_data <- fread(file.path(data_wd, "05_cleaned_database.csv")) # fast read; relies on data.table package
head(all_bees_data)
colnames(all_bees_data)
nrow(all_bees_data)

# all_bees_data contains distribution data for all bees with data in GBIF

# Make species columns and names consistent between tree_data and all_bees_data
tree_data$tips <- gsub("_", " ", tree_data$tips) # remove underscores so species names match
#colnames(tree_data) <- c("family","tribe","species","sociality","nesting") # rename "tips" column to "species" (consistent with all_bees_data)
head(tree_data)
colnames(tree_data)
nrow(tree_data)


#-------------------------------------------------------------------------------
#-----------------MAKE DATAFRAME OF SPECIES STILL TO SCORE----------------------

# Will use dplyr to find species in all_bees_data that aren't in tree_data
library(dplyr)

# Perform anti-join
bee_species_not_in_tree <- anti_join(all_bees_data, tree_data, by = "species")
head(bee_species_not_in_tree) # contains rows of all_bees_data where "species" column doesn't have a match in tree_data

# Investigate species not in tree
nrow(bee_species_not_in_tree) # 761200 rows
head(bee_species_not_in_tree$species) # There are multiple observations per species
length(unique(bee_species_not_in_tree$species)) # 7655 unique species (in all_bees_data, but not tree_data) - this is what we want

# Select unique species, keeping family and subfamily columns
unique_species_not_in_tree <- unique(bee_species_not_in_tree[c("family", "subfamily", "species")])

# Write unique species names to a CSV file
write.csv(unique_species_not_in_tree, "unique_species_not_in_tree.csv", row.names = FALSE)

# Load in that data
unique_df <- read.csv("unique_species_not_in_tree.csv")
head(unique_df)
unique(unique_df$family) # Unique species from all 7 families

# Let's see how many species there are per family that I need to score
unique_families <- unique(unique_species_not_in_tree$family)

# Loop through each unique family
for (family in unique_families) {
  # Filter the dataframe to include only rows where the "family" column is equal to the current family
  family_rows <- unique_species_not_in_tree[unique_species_not_in_tree$family == family, ]
  
  # Count the number of rows in the filtered subset
  num_family_rows <- nrow(family_rows)
  
  # Print
  cat("Number of species to score for family", family, ":", num_family_rows, "\n")
}

# Weird things: no Apis, some subfamilies missing from this dataset. All spp. are in both datasets?


#-------------------------------------------------------------------------------
#---------------------------------SUBSET DATAFRAME------------------------------

# Try using a subset of all_bees_data for heatmaps; finish scoring sociality later
# Note: make sure to transform "tips" column in tree_data to "species" before creating this subset
subset_all_bees <- all_bees_data[all_bees_data$species %in% tree_data$species, 
                                 c("species", "decimalLatitude", "decimalLongitude")]
colnames(subset_all_bees) <- c("species", "lat", "lon")
colnames(subset_all_bees)
nrow(subset_all_bees)
head(subset_all_bees)
# this contains occurrence (lat, lon) information for phylogeny tip species

# Save CSV to working directory 
write.csv(subset_all_bees, file.path(data_wd, "subset_all_bees.csv"), row.names = FALSE)

length(unique(subset_all_bees$species)) # 3952 species, so some species in tree_data aren't in all_bees_data
# Or their species names are formatted inconsistently. 
# re-ran this on 10/30/24: now 3748 species

length(unique(tree_data$species)) # 4538 species... so not all of our tip species have GBIF observations in all_bees data.
# re-ran on 10/30/24: now 4293 species
length(unique(all_bees_data$species)) # was originally 11,607

# Now to use this subset to make plots.
# Once we finish scoring, we can use all_bees_data instead.


#-------------------------------------------------------------------------------
#-------------------------SPECIES RICHNESS PER AREA-----------------------------

# First, we have to generate the total species richness for each area (i.e., count # of species within each area from global map)
# DON'T RUN THIS FUNCTION unless you need to re-make species lists - results are loaded below
organize.bubble.plot2 <- function(points, twgd_data) {
  focal_areas <- as.character(twgd_data$LEVEL3_COD)  # Names of areas in global map shapefile
  species <- points[,1] # Vector of species names
  points <- points[,c(3,2)] # Vector of lat and long
  sp::coordinates(points) <- ~ lon + lat # Transforming into sp format coordinates
  results <- matrix(nrow=0, ncol=4) # Create empty results matrix
  list_result1 <- list()
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    if(nrow(area_plus_buffer)>0) {
      res <- sp::over(points, area_plus_buffer)
      if(any(!is.na(res$LEVEL1_COD))) { # Should we use LEVEL3?
        sp_rich <- unique(species[which(!is.na(res$LEVEL1_COD))])
        list_result1[[i]] <- sp_rich
        names(list_result1)[i] <- one_area
        n_points <- length(sp_rich)
        coords <- sp::coordinates(area_plus_buffer)
        centroid_x <- mean(coords[, 1])
        centroid_y <- mean(coords[, 2])
        centroids <- sp::SpatialPoints(matrix(c(centroid_x, centroid_y), ncol = 2)) # Making centroid
        lon <- raster::extent(centroids)[1]
        lat <- raster::extent(centroids)[3]
        results <- rbind(results, cbind(n_points, one_area, lon, lat)) # Adds row for each iteration that is species richness in one_area for that iteration 
      }
      cat(i, "\r")
    }
  }
  results <- as.data.frame(results)
  results$n_points <- as.numeric(results$n_points)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  save(list_result1, file="list_result1.Rsave")
  return(results)
}

load("list_results_bees.Rsave") # Result of running the above function
list_result1 # Loading the data - list of species present in each country shapefile 
# Includes species not in tree_data that we haven't scored yet
# And also the names are formatted differently than tree_data
# Will have to harmonize the names 


#-------------------------------------------------------------------------------
#------------------------PLOTTING: SPECIES RICHNESS-----------------------------
# Here, we're using ALL bees (not just the subset we scored for sociality and nesting)

twgd_path = "TWDG/wgsrpd-master/level3/level3.shp"

# Load in map data
twgd_data <- st_read(twgd_path)
class(twgd_data) # sf dataframe

# Convert map data to sp object to prepare for organize.bubble.plot2 function
# and set CRS to NA so it is consistent with points
twgd_data <- as(twgd_data, "Spatial")
class(twgd_data) # SpatialPolygonsDataFrame
proj4string(twgd_data) <- CRS("") # set CRS to NA
proj4string(twgd_data) # NA

# Load in bee occurrence points
load("thinned_points_res1.Rsave") # Use thinned points instead of all from GBIF
colnames(thinned_points) <- c("species", "lat","lon")
head(thinned_points). # ASK THAIS ABOUT THIS

# Generate species richness per area 
richness_per_area <- organize.bubble.plot2(points = thinned_points, twgd_data) # 369
write.csv(richness_per_area, file.path(data_wd, "richness_per_area.csv"), row.names = FALSE) # Write to CSV so you don't have to run every time

# Prepare data for mapping
twgd_data_sf <- sf::st_as_sf(twgd_data) # convert twgd_data into sf object
colnames(twgd_data_sf) 
#colnames(results)
twgd_data_bees <- merge(twgd_data_sf, richness_per_area, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_bees <- subset(twgd_data_bees , twgd_data_bees$LEVEL1_COD%in%c(7,8)) # Just the Americas (7,8)
colnames(twgd_data_bees)
head(twgd_data_bees) # contains polygons with # of bee species per polygon

# Mapping
spp_rich_heatmap <- ggplot(data = twgd_data_bees) +
  geom_sf(aes(fill = n_points)) +
  scale_fill_viridis_c(option = "viridis") +
  theme_classic()

spp_rich_heatmap

#-------------------------------------------------------------------------------
#-----------------------PLOTTING: PROPORTION SOCIAL-----------------------------

# For this, we'll use subset_all_bees and tree_data
#   subset_all_bees is latitude and longitude observations for tip species we have scored
#   tree_data is trait scoring for tip species

# First, we need to calculate proportion social species for each region in TWGD level 3
# And ultimately get a dataset that has one row per TWGD level 3 region, 
#   the geometry that describes that region in MULTIPOLYGON format,
#   and proportion of species in that region which are social.

# Loading in bee occurrence data
subset_all_bees <- read.csv(file.path(wd, "subset_all_bees.csv"))
head(subset_all_bees)

# Load in map data
twgd_path <- "TWDG/wgsrpd-master/level3/level3.shp"
twgd_data <- st_read(twgd_path)
class(twgd_data) # sf dataframe
st_geometry_type(twgd_data) # MULTIPOLYGON geometry
st_is_valid(twgd_data) # contains invalid geometries...

# Running subset_all_bees through organize.bubble.plot2() to generate total species richness
twgd_data <- as(twgd_data, "Spatial") # convert map data to sp object to prepare for organize.bubble.plot2 function
class(twgd_data) # SpatialPolygonsDataFrame
proj4string(twgd_data) <- CRS("") # set CRS to NA to be consistent with points in subset_all_bees
proj4string(twgd_data) # NA
results_total_richness <- organize.bubble.plot2(points = subset_all_bees, twgd_data = twgd_data, colnames = "spp_rich") # generate species richness per area 

head(results_total_richness) 
nrow(results_total_richness) # 314
summary(results_total_richness$spp_rich)

# Merge datasets into one that contains both occurrences + sociality scoring
merged_data <- merge(subset_all_bees, tree_data, by = "species") 
sociality_occurrence <- merged_data[ , c("species", "sociality", "lat", "lon")] # select relevant columns from merged dataset
# sociality_occurrence <- st_as_sf(sociality_occurrence, coords = c("lon", "lat")) # convert into sf object
head(sociality_occurrence)

# Subset sociality_occurrence to include only social species
#   then re-run organize.bubble.plot2() to get species richness per polygon
#   but just for social species
just_social_spp <- sociality_occurrence[sociality_occurrence$sociality == "social", ]
head(just_social_spp)
results_social_richness <- organize.bubble.plot2(points = just_social_spp, twgd_data, colname = "social_rich") 
head(results_social_richness)
nrow(results_social_richness) # 309 - why less than results_total_richness?

# Join results_social_richness and results_total_richness
all_rich <- merge(results_social_richness, results_total_richness, by = "one_area")
head(all_rich)
nrow(all_rich)

# Calculate proportion social social_rich/spp_rich
#   If social_rich/spp_rich is Inf or NA, make it 0
all_rich$prop_social <- all_rich$social_rich / all_rich$spp_rich
head(all_rich)
all_rich <- all_rich[ , c("one_area", "prop_social")]
head(all_rich) 

# Join proportion sociality data with polygon coordinates that correspond to one_area codes
bee_twgd_sociality <- merge(twgd_data, all_rich, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_sociality <- subset(bee_twgd_sociality, bee_twgd_sociality$LEVEL1_COD%in%c(7,8)) # subset to Americas
bee_twgd_sociality <- bee_twgd_sociality[ , c("LEVEL3_NAM","LEVEL3_COD","prop_social","geometry")]
write.csv(bee_twgd_sociality, "bee_twgd_sociality.csv", row.names = FALSE) # write to CSV so we don't have to run it all again
head(bee_twgd_sociality)

# Plot
sociality_heatmap <- ggplot(data = bee_twgd_sociality) +
  geom_sf(aes(fill = prop_social)) + 
  scale_fill_viridis_c(option = "C",alpha = 0.9,  name = "Proportion Social Bee Species") +
  theme_classic() +
  coord_sf(ylim = c(-60, 90), xlim = c(-170, 0), expand = FALSE)

sociality_heatmap

#-------------------------------------------------------------------------------
#----------------PLOTTING: PROPORTION ABOVE-GROUND NESTING----------------------

# Merge datasets into one that contains both occurrences + nesting scoring
merged_data <- merge(subset_all_bees, tree_data, by = "species") 
nesting_occurrence <- merged_data[ , c("species","nesting","lat","lon")] # select relevant columns from merged dataset
#sociality_occurrence <- st_as_sf(sociality_occurrence, coords = c("lon", "lat")) # convert into sf object
head(nesting_occurrence)

# Generate subset of nesting_occurrence dataset that only contains above-ground nesting species
just_above_ground <- nesting_occurrence[nesting_occurrence$nesting == "above-ground", ]
head(just_above_ground)

# Run above-ground nesting subset through function to get richness of above-ground nesting species per area
twgd_data <- as(twgd_data, "Spatial") # convert map data to sp object to prepare for organize.bubble.plot2 function
class(twgd_data) # SpatialPolygonsDataFrame
proj4string(twgd_data) <- CRS("") # set CRS to NA to be consistent with points in subset_all_bees
proj4string(twgd_data) # NA
results_nesting_richness <- organize.bubble.plot2(points = just_above_ground, twgd_data, colname = "above_ground_rich") 
head(results_nesting_richness)
nrow(results_nesting_richness) # 310
write.csv(results_nesting_richness, "results_nesting_richness.csv", row.names = FALSE)

# Join results_nesting_richness and results_total_richness
all_rich2 <- merge(results_nesting_richness, results_total_richness, by = "one_area")
head(all_rich2)
nrow(all_rich2)

# Calculate proportion above-ground nesters
all_rich2$prop_nest <- all_rich2$above_ground_rich / all_rich2$spp_rich
head(all_rich2)
all_rich2 <- all_rich2[ , c("one_area", "prop_nest")]
head(all_rich2) 

# Join proportion nesting data with polygon coordinates that correspond to one_area codes
bee_twgd_nest <- merge(twgd_data, all_rich2, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_nest <- subset(bee_twgd_nest, bee_twgd_nest$LEVEL1_COD%in%c(7,8)) # subset to Americas
bee_twgd_nest <- bee_twgd_nest[ , c("LEVEL3_NAM","LEVEL3_COD","prop_nest","geometry")]
write.csv(bee_twgd_nest, "bee_twgd_nest.csv", row.names = FALSE) # write to CSV so we don't have to run it all again
head(bee_twgd_nest)

# Plot
nesting_heatmap <- ggplot(data = bee_twgd_nest) +
  geom_sf(aes(fill = prop_nest)) + 
  scale_fill_viridis_c(option = "C", alpha = 0.9,  name = "Proportion Above-Ground Nesting Bee Species") +
  theme_classic() +
  coord_sf(ylim = c(-60, 90), xlim = c(-170, 0), expand = FALSE)

nesting_heatmap


#-------------------------------------------------------------------------------
#-----------------------------PLOTTING: TOGETHER--------------------------------
# First, we need a merged dataset that includes the columns one_area, geometry, prop_social, and prop_nest
#     Since we saved these datasets separately as CSVs, we will load them, join them, then create the faceted plot
#     With a proportion social heatmap as A, and a proportion above-ground nesting heatmap as B

# Load libraries
library(dplyr)
library(ggplot2)
library(sf)
library(gridExtra)

# Load datasets required for heatmaps
bee_twgd_sociality <- read.csv("bee_twgd_sociality.csv", header = TRUE) 
bee_twgd_nest <- read.csv("bee_twgd_nest.csv", header = TRUE)

# Convert to sf objects
bee_twgd_sociality <- st_as_sf(bee_twgd_sociality) 
bee_twgd_nest <- st_as_sf(bee_twgd_nest)

# Arrange plots side by side using grid.arrange from gridExtra
grid.arrange(sociality_heatmap, nesting_heatmap, ncol = 2)

