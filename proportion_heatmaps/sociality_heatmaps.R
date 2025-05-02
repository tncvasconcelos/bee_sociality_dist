# ==============================================================================
# Heatmaps
# ==============================================================================
# Plots heatmaps of the Americas showing proportion of social and above-ground species
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load libraries
# ------------------------------------------------------------------------------

#rm(list=ls())

library(sp)
library(sf)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(data.table)
library(patchwork)

getwd()      
wd <- "/Users/lenarh/Desktop/bee_sociality_dist/proportion_heatmaps"
setwd(wd)
data_wd <- paste0(wd, "/data")
data_wd

# Source code
# functions_path <- file.path(wd, "example code", "00_utility_functions_synthesis.R") 
# source(functions_path)
# Note 11/6/24 - function that is in this document works, but one in 00_utility_functions_synthesis.R doesn't


# ------------------------------------------------------------------------------
# Load trait and occurrence data
# ------------------------------------------------------------------------------

tree_spp_traits <- read.csv(file.path(data_wd, "bees_traits.csv")) 
colnames(tree_spp_traits)[3] <- "species" # to be consistent with gbif_occurrence_data
colnames(tree_spp_traits) 
nrow(tree_spp_traits) # 4293
# tree_spp_traits contains trait information (sociality, nesting) for phylogeny tip species

gbif_occurrence_data <- fread(file.path(data_wd, "05_cleaned_database.csv")) # fast read; relies on data.table package
head(gbif_occurrence_data)
colnames(gbif_occurrence_data)
nrow(gbif_occurrence_data) # 6890148
# gbif_occurrence_data contains distribution data for all bees with data in GBIF

tree_spp_traits$species

# Make species columns and names consistent between tree_spp_traits and gbif_occurrence_data
tree_spp_traits$species <- gsub("_", " ", tree_spp_traits$species) # remove underscores so species names match
tree_spp_traits$species


# ------------------------------------------------------------------------------
# Make dataframe of species still to score
# ------------------------------------------------------------------------------
# Find species in gbif_occurrence_data that aren't in tree_spp_traits
# ------------------------------------------------------------------------------

# # Perform anti-join
# bee_species_not_in_tree <- anti_join(gbif_occurrence_data, tree_spp_traits, by = "species")
# head(bee_species_not_in_tree) # contains rows of gbif_occurrence_data where "species" column doesn't have a match in tree_spp_traits
# 
# # Investigate species not in tree
# nrow(bee_species_not_in_tree) # 761200 rows
# head(bee_species_not_in_tree$species) # There are multiple observations per species
# length(unique(bee_species_not_in_tree$species)) # 7655 unique species (in gbif_occurrence_data, but not tree_spp_traits) - this is what we want
# 
# # Select unique species, keeping family and subfamily columns
# unique_species_not_in_tree <- unique(bee_species_not_in_tree[c("family", "subfamily", "species")])
# unique(unique_species_not_in_tree$family) # Unique species from all 7 families
# 
# # Write unique species names to a CSV file
# write.csv(unique_species_not_in_tree, "unique_species_not_in_tree.csv", row.names = FALSE)
# 
# # Let's see how many species there are per family that I need to score
# unique_families <- unique(unique_species_not_in_tree$family)
# 
# # Loop through each unique family
# for (family in unique_families) {
#   # Filter the dataframe to include only rows where the "family" column is equal to the current family
#   family_rows <- unique_species_not_in_tree[unique_species_not_in_tree$family == family, ]
#   
#   # Count the number of rows in the filtered subset
#   num_family_rows <- nrow(family_rows)
#   
#   # Print
#   cat("Number of species to score for family", family, ":", num_family_rows, "\n")
# }
# 
# # Weird things: no Apis, some subfamilies missing from this dataset. All spp. are in both datasets?
# 


# ------------------------------------------------------------------------------
# Subset gbif occurrence data to only include species scored for nesting/sociality
# ------------------------------------------------------------------------------

# Note: make sure to transform "tips" column in tree_spp_traits to "species" before creating this subset
gbif_subset <- gbif_occurrence_data[gbif_occurrence_data$species %in% tree_spp_traits$species, 
                                 c("species", "decimalLatitude", "decimalLongitude")]
colnames(gbif_subset) <- c("species", "lat", "lon")
colnames(gbif_subset)
nrow(gbif_subset) # 6063857
head(gbif_subset)

# gbif_subset contains occurrence (lat, lon) information for phylogeny tip species

# Save csv to data wd
write.csv(gbif_subset, file.path(data_wd, "gbif_subset.csv"), row.names = FALSE) # ran 1/10/25

length(unique(gbif_subset$species)) # 3748 species (compared to 4293 in tree_spp_traits), 
# so some phylogeny tip species aren't in gbif_occurrence_data/lack occurrence points

length(unique(gbif_occurrence_data$species)) # was originally 11,607

# Once we finish scoring all bees from the Dorey et al. 2023 dataset, we can use gbif_occurrence_data instead

# Now to use this subset to make plots.


# ------------------------------------------------------------------------------
# Species richness per area function
# ------------------------------------------------------------------------------

# First, we have to generate the total species richness for each area 
# (i.e., count # of species within each area from global map)
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
# Includes species not in tree_spp_traits that we haven't scored yet
# And also the names are formatted differently than tree_spp_traits
# Will have to harmonize the names 


# ------------------------------------------------------------------------------
# Plotting: species richness
# ------------------------------------------------------------------------------
# Here, we're using ALL bees (not just the subset we scored for sociality and nesting)
# To look at overall species richness patterns
# Skip to loading the richness_per_area csv (below) which was already created, unless you need to re-run this
# ------------------------------------------------------------------

twgd_path <- "TWDG/wgsrpd-master/level3/level3.shp"

# Load in map data
twgd_data <- st_read(twgd_path)
class(twgd_data) # sf dataframe

# Convert map data to sp object to prepare for organize.bubble.plot2 function
# and set CRS to NA so it is consistent with points
twgd_data <- as(twgd_data, "Spatial")
class(twgd_data) # SpatialPolygonsDataFrame
proj4string(twgd_data) <- CRS("") 
proj4string(twgd_data) 

# Load in bee occurrence points
load("thinned_points_res1.Rsave") # Use thinned points instead of all from GBIF (i.e., gbif_occurrence_data)
colnames(thinned_points) <- c("species", "lat","lon")
head(thinned_points) # filtered GBIF data for phylogeny tip species 
# (only one occurrence point per cell, exotic species removed)

# Generate species richness per area 
richness_per_area <- organize.bubble.plot2(points = thinned_points, twgd_data) # 369
colnames(richness_per_area)[1] <- "spp_rich"
head(richness_per_area)
write.csv(richness_per_area, file.path(data_wd, "richness_per_area.csv"), row.names = FALSE) # Write to CSV so you don't have to run every time
#richness_per_area <- read.csv(file.path(data_wd, "richness_per_area.csv"))

# Prepare data for mapping
twgd_data_sf <- sf::st_as_sf(twgd_data) # convert twgd_data into sf object
colnames(twgd_data_sf) 
twgd_data_bees <- merge(twgd_data_sf, richness_per_area, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_bees <- subset(twgd_data_bees , twgd_data_bees$LEVEL1_COD%in%c(7,8)) # Just the Americas (7,8)
colnames(twgd_data_bees)
head(twgd_data_bees) # contains polygons with # of bee species per polygon

# Mapping
spp_rich_heatmap <- ggplot(data = twgd_data_bees) +
  geom_sf(aes(fill = spp_rich)) +
  scale_fill_viridis_c() +
  theme_classic()

spp_rich_heatmap


# ------------------------------------------------------------------------------
# Plotting: proportion social
# ------------------------------------------------------------------------------
# For this, we'll use the thinned GBIF points and tree_spp_traits
#   thinned_points is latitude and longitude observations for tip species we have scored, with additional filtering
#   tree_spp_traits is trait scoring for tip species

# First, we need to calculate proportion social species for each region in TWGD level 3
# And ultimately get a dataset that has one row per TWGD level 3 region, 
#   the geometry that describes that region in MULTIPOLYGON format,
#   and proportion of species in that region which are social.
# ------------------------------------------------------------------------------

thinned_points$species <- gsub("_"," ", thinned_points$species) # Make species name formatting consistent

# Create subset of thinned points that contains just species that are 1) in the tree and 2) social
social_subset <- subset(thinned_points, thinned_points$species%in%tree_spp_traits$species[which(tree_spp_traits$sociality_binary=="social")])

# Generate richness of social bees per TWGD area
richness_social_per_area <- organize.bubble.plot2(points = social_subset, twgd_data) # 369
head(richness_social_per_area)
colnames(richness_social_per_area)[1] <- "social_rich"
head(richness_social_per_area)

# Join results_social_richness and richness_per_area
all_rich <- merge(richness_social_per_area, richness_per_area, by = "one_area")
head(all_rich)
all_rich <- all_rich[, -c(6, 7)] # Remove the last two columns
colnames(all_rich)[3:4] <- c("lon", "lat") # Rename columns 3 and 4 to "lon" and "lat"
all_rich <- all_rich[, c("one_area", "lat", "lon", "social_rich", "spp_rich")] # Reorder columns
head(all_rich)
nrow(all_rich) # 306... why?

# Calculate proportion social using social_rich/spp_rich
#   If social_rich/spp_rich is Inf or NA, make it 0
all_rich$prop_social <- all_rich$social_rich / all_rich$spp_rich
head(all_rich)
all_rich <- all_rich[ , c("one_area", "prop_social")]
head(all_rich) # now all_rich just contains the area code and proportion social species

# Add polygon information from TWGD
bee_twgd_sociality <- merge(twgd_data, all_rich, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_sociality <- subset(bee_twgd_sociality, bee_twgd_sociality$LEVEL1_COD%in%c(7,8)) # subset to Americas
bee_twgd_sociality <- st_as_sf(bee_twgd_sociality)
bee_twgd_sociality
head(bee_twgd_sociality)

# Save final dataset which contains area codes, polygons, and proportion social species
write.csv(bee_twgd_sociality, file.path(data_wd, "bee_twgd_sociality.csv"), row.names = FALSE) # ran 1/10/25
bee_twgd_sociality <- read.csv(file.path(data_wd, "bee_twgd_sociality.csv")) # error: "more columns than column names"

# Mapping
prop_social_heatmap <- ggplot(data = bee_twgd_sociality) +
  geom_sf(aes(fill = prop_social)) +
  labs(x = "Latitude", y = "Longitude", fill = "Proportion Social") +  # Relabel axes and legend
  scale_fill_viridis_c(option = "viridis") +
  xlim(-175, 0) +
  theme_classic()

prop_social_heatmap


# ------------------------------------------------------------------------------
# Plotting: proportion above-ground nesting
# ------------------------------------------------------------------------------

thinned_points$species <- gsub("_"," ", thinned_points$species) # Make species name formatting consistent

# Create subset of thinned points that contains just species that are 1) in the tree and 2) above-ground nesters
abvgrnd_subset <- subset(thinned_points, thinned_points$species%in%tree_spp_traits$species[which(tree_spp_traits$nest_binary=="aboveground")])

# Generate richness of social bees per TWGD area
richness_abvgrnd_per_area <- organize.bubble.plot2(points = abvgrnd_subset, twgd_data) # 369
head(richness_abvgrnd_per_area)
colnames(richness_abvgrnd_per_area)[1] <- "aboveground_rich"
head(richness_abvgrnd_per_area)

# Join results_abvgrnd_richness and richness_per_area
all_rich2 <- merge(richness_abvgrnd_per_area, richness_per_area, by = "one_area")
head(all_rich2)
all_rich2 <- all_rich2[, -c(6, 7)] # Remove the last two columns
colnames(all_rich2)[3:4] <- c("lon", "lat") # Rename columns 3 and 4 to "lon" and "lat"
all_rich2 <- all_rich2[, c("one_area", "lat", "lon", "aboveground_rich", "spp_rich")] # Reorder columns
head(all_rich2)
nrow(all_rich2) # 309... why?

# Calculate proportion social using aboveground_rich/spp_rich
#   If aboveground_rich/spp_rich is Inf or NA, make it 0
all_rich2$prop_aboveground <- all_rich2$aboveground_rich / all_rich2$spp_rich
head(all_rich2)
all_rich2 <- all_rich2[ , c("one_area", "prop_aboveground")]
head(all_rich2) # now all_rich2 just contains the area code and proportion above-ground nesting species

# Add polygon information from TWGD
bee_twgd_nest <- merge(twgd_data, all_rich2, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_nest <- subset(bee_twgd_nest, bee_twgd_nest$LEVEL1_COD%in%c(7,8)) # subset to Americas
bee_twgd_nest <- st_as_sf(bee_twgd_nest)
bee_twgd_nest

# Save final dataset which contains area codes, polygons, and proportion above-ground nesting species
write.csv(bee_twgd_nest, file.path(data_wd, "bee_twgd_nest.csv"), row.names = FALSE) # ran 1/10/25

# Mapping
prop_aboveground_heatmap <- ggplot(data = bee_twgd_nest) +
  geom_sf(aes(fill = prop_aboveground)) +
  labs(x = "Latitude", y = "Longitude", fill = "Proportion Above-Ground Nesting") +  # Relabel axes and legend
  scale_fill_viridis_c(option = "viridis") +
  xlim(-175, 0) +
  theme_classic()

prop_aboveground_heatmap


# ------------------------------------------------------------------------------
# Plotting: both heatmaps together
# ------------------------------------------------------------------------------

# Combine the two plots
combined_plot <- prop_aboveground_heatmap + prop_social_heatmap +
  plot_layout(ncol = 1) + # Arrange the plots in one column
  plot_annotation(tag_levels = 'A') # This will add "A" and "B" labels automatically

combined_plot # Display combined plot

# Save plot
ggsave("/Users/lenarh/Desktop/bee_sociality_dist/proportion_heatmaps/plots/combined_plot.pdf", 
       plot = combined_plot, width = 8, height = 6)


# ------------------------------------------------------------------------------
# Plotting: now excluding bumblebees (Bombus)
# ------------------------------------------------------------------------------

# Check how many Bombus are in datasets
sum(grepl("Bombus", tree_spp_traits$species)) # 243
sum(grepl("Bombus", thinned_points$species)) # 33283

# Remove Bombus from trait and GBIF datasets
tree_spp_traits_bombus <- tree_spp_traits[!grepl("^Bombus", tree_spp_traits$species), ]
thinned_points_bombus <- thinned_points[!grepl("^Bombus", thinned_points$species), ]

# Confirm Bombus removed
sum(grepl("Bombus", tree_spp_traits_bombus$species)) # 0 
sum(grepl("Bombus", thinned_points_bombus$species)) # 0

# Species richness
twgd_path <- "TWDG/wgsrpd-master/level3/level3.shp"
twgd_data <- st_read(twgd_path)
twgd_data <- as(twgd_data, "Spatial")
class(twgd_data) 
proj4string(twgd_data) <- CRS("") 
proj4string(twgd_data) 

richness_per_area_bombus <- organize.bubble.plot2(points = thinned_points_bombus, twgd_data) # 369
colnames(richness_per_area_bombus)[1] <- "spp_rich"
head(richness_per_area_bombus)
write.csv(richness_per_area_bombus, file.path(data_wd, "richness_per_area_bombus.csv"), row.names = FALSE)

twgd_data_sf <- sf::st_as_sf(twgd_data)
colnames(twgd_data_sf) 
twgd_data_bees_bombus <- merge(twgd_data_sf, richness_per_area_bombus, by.x="LEVEL3_COD", by.y="one_area")
twgd_data_bees_bombus <- subset(twgd_data_bees_bombus, twgd_data_bees_bombus$LEVEL1_COD%in%c(7,8))
colnames(twgd_data_bees_bombus)
head(twgd_data_bees_bombus) 

spp_rich_heatmap_bombus <- ggplot(data = twgd_data_bees_bombus) +
  geom_sf(aes(fill = spp_rich)) +
  scale_fill_viridis_c() +
  theme_classic()

spp_rich_heatmap_bombus

# Above-ground richness
abvgrnd_subset_bombus <- subset(thinned_points_bombus, thinned_points_bombus$species%in%tree_spp_traits_bombus$species[which(tree_spp_traits_bombus$nest_binary=="aboveground")])
richness_abvgrnd_per_area_bombus <- organize.bubble.plot2(points = abvgrnd_subset_bombus, twgd_data) # 369
head(richness_abvgrnd_per_area_bombus)
colnames(richness_abvgrnd_per_area_bombus)[1] <- "aboveground_rich"
head(richness_abvgrnd_per_area_bombus)

all_rich2_bombus <- merge(richness_abvgrnd_per_area_bombus, richness_per_area_bombus, by = "one_area")
head(all_rich2_bombus)
all_rich2_bombus <- all_rich2_bombus[, -c(6, 7)] 
colnames(all_rich2_bombus)[3:4] <- c("lon", "lat") 
all_rich2_bombus <- all_rich2_bombus[, c("one_area", "lat", "lon", "aboveground_rich", "spp_rich")] 
head(all_rich2_bombus)
nrow(all_rich2_bombus) # 305

all_rich2_bombus$prop_aboveground <- all_rich2_bombus$aboveground_rich / all_rich2_bombus$spp_rich
head(all_rich2_bombus)
all_rich2_bombus <- all_rich2_bombus[ , c("one_area", "prop_aboveground")]
head(all_rich2_bombus)

bee_twgd_nest_bombus <- merge(twgd_data, all_rich2_bombus, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_nest_bombus <- subset(bee_twgd_nest_bombus, bee_twgd_nest_bombus$LEVEL1_COD%in%c(7,8)) 
bee_twgd_nest_bombus <- st_as_sf(bee_twgd_nest_bombus)
bee_twgd_nest_bombus

write.csv(bee_twgd_nest_bombus, file.path(data_wd, "bee_twgd_nest_bombus.csv"), row.names = FALSE)

prop_aboveground_heatmap_bombus <- ggplot(data = bee_twgd_nest_bombus) +
  geom_sf(aes(fill = prop_aboveground)) +
  labs(x = "Latitude", y = "Longitude", fill = "Proportion Above-Ground Nesting") +  
  scale_fill_viridis_c(option = "viridis") +
  xlim(-175, 0) +
  theme_classic()

prop_aboveground_heatmap_bombus # Greenland is still showing up as 1.0 proportion above-ground nesting.

# Social richness
social_subset_bombus <- subset(thinned_points_bombus, thinned_points_bombus$species%in%tree_spp_traits_bombus$species[which(tree_spp_traits_bombus$sociality_binary=="social")])
richness_social_per_area_bombus <- organize.bubble.plot2(points = social_subset_bombus, twgd_data) # 369
head(richness_social_per_area_bombus)
colnames(richness_social_per_area_bombus)[1] <- "social_rich"
head(richness_social_per_area_bombus)
all_rich_bombus <- merge(richness_social_per_area_bombus, richness_per_area_bombus, by = "one_area")
head(all_rich_bombus)
all_rich_bombus <- all_rich_bombus[, -c(6, 7)] 
colnames(all_rich_bombus)[3:4] <- c("lon", "lat") 
all_rich_bombus <- all_rich_bombus[, c("one_area", "lat", "lon", "social_rich", "spp_rich")] 
head(all_rich_bombus)
nrow(all_rich_bombus) # 298

all_rich_bombus$prop_social <- all_rich_bombus$social_rich / all_rich_bombus$spp_rich
head(all_rich_bombus)
all_rich_bombus <- all_rich_bombus[ , c("one_area", "prop_social")]
head(all_rich_bombus) 

bee_twgd_sociality_bombus <- merge(twgd_data, all_rich_bombus, by.x = "LEVEL3_COD", by.y = "one_area")
bee_twgd_sociality_bombus <- subset(bee_twgd_sociality_bombus, bee_twgd_sociality_bombus$LEVEL1_COD%in%c(7,8)) # subset to Americas
bee_twgd_sociality_bombus <- st_as_sf(bee_twgd_sociality_bombus)
bee_twgd_sociality_bombus

write.csv(bee_twgd_sociality_bombus, file.path(data_wd, "bee_twgd_sociality_bombus.csv"), row.names = FALSE)

prop_social_heatmap_bombus <- ggplot(data = bee_twgd_sociality_bombus) +
  geom_sf(aes(fill = prop_social)) +
  labs(x = "Latitude", y = "Longitude", fill = "Proportion Social") +  # Relabel axes and legend
  scale_fill_viridis_c(option = "viridis") +
  xlim(-175, 0) +
  theme_classic()

prop_social_heatmap_bombus # Greenland is still showing up as 1.0 proportion social.

# Combined plot
combined_plot_bombus <- prop_aboveground_heatmap_bombus + prop_social_heatmap_bombus +
  plot_annotation(tag_levels = 'A') # This will add "A" and "B" labels automatically

combined_plot_bombus # Display combined plot
