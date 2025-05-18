# ==============================================================================
# Heatmaps
# ==============================================================================
# Plots heatmaps of the Americas showing proportion of social and above-ground species
# And makes scatterplots of absolute latitude vs. proportion social/above-ground nesting
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load libraries
# ------------------------------------------------------------------------------

rm(list=ls())

# Libraries
library(sp)
library(sf)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(data.table)
library(patchwork)
library(lwgeom)
library(nlme)

# Directories
wd <- "/Users/lenarh/Desktop/bee_sociality_dist/proportion_heatmaps"
setwd(wd)
data_wd <- paste0(wd, "/data")
plots_wd <- paste0(wd, "/plots")

# Load TWGD shapefile
twgd_path <- file.path(wd, "TWDG/wgsrpd-master/level3/level3.shp")
twgd_data <- st_read(twgd_path) %>% as("Spatial")
proj4string(twgd_data) <- CRS("")

# Source code
# functions_path <- file.path(wd, "example code", "00_utility_functions_synthesis.R") 
# source(functions_path)
# Note 11/6/24 - function that is in this document works, but one in 00_utility_functions_synthesis.R doesn't


# ==============================================================================
# Quick load (if skipping full re-run)
# This is everything you need to skip right to plotting
bee_twgd <- st_read(file.path(data_wd, "bee_twgd.gpkg")) # Overall species richness per polygon
st_crs(bee_twgd) <- 4326

bee_twgd_sociality <- st_read(file.path(data_wd, "bee_twgd_sociality.gpkg")) # Proportion social per polygon
st_crs(bee_twgd_sociality) <- 4326

bee_twgd_nest <- st_read(file.path(data_wd, "bee_twgd_nest.gpkg")) # Proportion above-ground nesting per polygon
st_crs(bee_twgd_nest) <- 4326

scatterplot_data_clean <- read.csv(file.path(data_wd, "scatterplot_data_clean.csv")) # Data prepped for scatterplots
# ==============================================================================


# ------------------------------------------------------------------------------
# Load trait and occurrence data
# ------------------------------------------------------------------------------

# Load trait data
tree_spp_traits <- read.csv(file.path(data_wd, "bees_traits.csv"))
colnames(tree_spp_traits)[3] <- "species"
tree_spp_traits$species <- gsub("_", " ", tree_spp_traits$species)
# tree_spp_traits contains trait information (sociality, nesting) for phylogeny tip species

# Load GBIF occurrence points
load("thinned_points_res1.Rsave") # Use thinned points instead of all from GBIF
colnames(thinned_points) <- c("species", "lat", "lon")
thinned_points$species <- gsub("_", " ", thinned_points$species)
head(thinned_points) # filtered GBIF data for phylogeny tip species
# (only one occurrence point per cell, exotic species removed)
# Includes species not in tree_spp_traits that we haven't scored yet
# Do we need to harmonize names (GBIF vs. phylogeny tip species name formatting?)


# ------------------------------------------------------------------------------
# Load or create species richness per region
# ------------------------------------------------------------------------------

# Load function
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


# Generate species richness per region using organize.bubble.plot() and thinned GBIF points
richness_per_area <- organize.bubble.plot2(points = thinned_points, twgd_data) # 369
colnames(richness_per_area)[1] <- "spp_rich"
write.csv(richness_per_area, file.path(data_wd, "richness_per_area.csv"), row.names = FALSE)
richness_per_area <- read.csv(file.path(data_wd, "richness_per_area.csv"))
head(richness_per_area)


# ------------------------------------------------------------------------------
# Compute proportion social and above-ground nesting per polygon
# ------------------------------------------------------------------------------

# Social
social_subset <- thinned_points[thinned_points$species %in% tree_spp_traits$species[tree_spp_traits$sociality_binary == "social"], ]
richness_social_per_area <- organize.bubble.plot2(social_subset, twgd_data)
colnames(richness_social_per_area)[1] <- "social_rich"

# Nesting
abvgrnd_subset <- thinned_points[thinned_points$species %in% tree_spp_traits$species[tree_spp_traits$nest_binary == "aboveground"], ]
richness_abvgrnd_per_area <- organize.bubble.plot2(abvgrnd_subset, twgd_data)
colnames(richness_abvgrnd_per_area)[1] <- "aboveground_rich"

# Merge with total richness
all_rich_social<- merge(richness_social_per_area, richness_per_area, by = "one_area") %>%
  mutate(prop_social = social_rich / spp_rich) %>%
  select(one_area, prop_social)
write.csv(all_rich_social, file.path(data_wd, "all_rich_social.csv"), row.names = FALSE)

all_rich_nest <- merge(richness_abvgrnd_per_area, richness_per_area, by = "one_area") %>%
  mutate(prop_aboveground = aboveground_rich / spp_rich) %>%
  select(one_area, prop_aboveground)
write.csv(all_rich_nest, file.path(data_wd, "all_rich_nest.csv"), row.names = FALSE)


# ------------------------------------------------------------------------------
# Join trait proportion data to spatial polygons
# ------------------------------------------------------------------------------

twgd_data_sf <- st_as_sf(twgd_data)
st_crs(twgd_data_sf) <- 4326

# Join polygons with species richness
bee_twgd <- merge(twgd_data_sf, richness_per_area, by.x = "LEVEL3_COD", by.y = "one_area") %>%
  filter(LEVEL1_COD %in% c(7,8))
st_write(bee_twgd, file.path(data_wd, "bee_twgd.gpkg"), delete_dsn = TRUE)

# Join polygons with sociality proportion
bee_twgd_sociality <- merge(twgd_data_sf, all_rich_social, by.x = "LEVEL3_COD", by.y = "one_area") %>%
  filter(LEVEL1_COD %in% c(7,8)) %>% st_as_sf()
st_write(bee_twgd_sociality, file.path(data_wd, "bee_twgd_sociality.gpkg"), delete_dsn = TRUE)

# Join polygons with nesting proportion
bee_twgd_nest <- merge(twgd_data_sf, all_rich_nest, by.x = "LEVEL3_COD", by.y = "one_area") %>%
  filter(LEVEL1_COD %in% c(7,8)) %>% st_as_sf()
st_write(bee_twgd_nest, file.path(data_wd, "bee_twgd_nest.gpkg"), delete_dsn = TRUE)

# See which regions are in bee_twgd but not in bee_twgd_sociality
anti_join(bee_twgd %>% st_drop_geometry(), bee_twgd_sociality %>% st_drop_geometry(), by = "LEVEL3_COD") %>%
  select(LEVEL3_COD, LEVEL3_NAM)


# ------------------------------------------------------------------------------
# Plot heatmaps
# ------------------------------------------------------------------------------

plot_heatmap <- function(data, var, title, viridis_option = "viridis") {
  ggplot(data) +
    # Plot polygons filled by the specified variable
    geom_sf(aes(fill = !!sym(var))) +
    # Axis and legend labels
    labs(x = "Longitude", y = "Latitude", fill = title) +
    # Fill scale
    scale_fill_viridis_c(
      option = viridis_option,
      limits = c(0, 1)
    ) +
    # Map limits and projection
    coord_sf(
      xlim = c(-175, 0),
      ylim = c(-60, 90),
      expand = FALSE
    ) +
    # Base theme
    theme_bw() +
    # Theme customization
    theme(
      legend.position       = c(0.95, 0.80),
      legend.justification  = "left",
      legend.background     = element_rect(fill = "white", color = NA),
      legend.key.height     = unit(0.6, "cm"),
      legend.key.width      = unit(0.6, "cm"),
      legend.title          = element_text(size = 12, margin = margin(b = 6)),
      legend.text           = element_text(size = 10),
      panel.grid            = element_blank(),
      panel.border          = element_blank(),
      plot.background       = element_blank(),
      axis.text             = element_text(color = "black", size = 10),
      axis.title            = element_text(color = "black", size = 12),
      axis.ticks            = element_line(color = "black"),
      axis.ticks.length     = unit(0.15, "cm"),
      axis.line             = element_line(color = "black")
    )
}

prop_social_heatmap <- plot_heatmap(bee_twgd_sociality, "prop_social", "Proportion Social", viridis_option = "cividis")
prop_aboveground_heatmap <- plot_heatmap(bee_twgd_nest, "prop_aboveground", "Proportion\nAbove-Ground\nNesting", viridis_option = "cividis")

combined_heatmaps <- prop_aboveground_heatmap + prop_social_heatmap +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 16, face = "bold"))

quartz()

print(combined_heatmaps)

ggsave(file.path(plots_wd, "combined_heatmaps.pdf"), combined_heatmaps, width = 12, height = 12, dpi = 1500)
# 5/6/2025: Why does Greenland now have 4 species one of which is ground-nesting?

ggsave(file.path(plots_wd, "combined_heatmaps.png"), combined_heatmaps, width = 12, height = 12, dpi = 1500)


# ------------------------------------------------------------------------------
# Scatterplots: abs(latitude) vs trait proportions with spatial autocorrelation correction
# ------------------------------------------------------------------------------

# Drop geometry from both heatmap sf objects to get raw trait values
social_df <- bee_twgd_sociality %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, prop_social)
social_df

nesting_df <- bee_twgd_nest %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, prop_aboveground)
nesting_df

# Repair invalid geometries first
bee_twgd_valid <- st_make_valid(bee_twgd)

# Then compute centroids
centroid_coords <- bee_twgd_valid %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame()

# Combine with region ID and species richness
richness_clean <- bee_twgd_valid %>%
  st_drop_geometry() %>%
  select(one_area = LEVEL3_COD, spp_rich) %>%
  bind_cols(centroid_coords) %>%
  rename(lon = X, lat = Y)

head(richness_clean)

# Merge into dataframe for plotting
scatterplot_data <- social_df %>%
  inner_join(nesting_df, by = "one_area") %>%
  inner_join(richness_clean, by = "one_area")
head(scatterplot_data)
colnames(scatterplot_data)
nrow(scatterplot_data) # 113

# Calculate absolute latitude
scatterplot_data$abs_lat <- abs(scatterplot_data$lat)

# Remove NAs and inspect
scatterplot_data_clean <- scatterplot_data %>%
  filter(!is.na(abs_lat), !is.na(prop_social), !is.na(prop_aboveground), !is.na(lon), !is.na(lat))

nrow(scatterplot_data_clean) # 113

scatterplot_data_clean %>% # Check which polygons have prop_social == 1
  filter(prop_social == 1) # ALU, ARU, BER, GNL

bee_twgd %>% # Get country names for polygons which prop_social == 1
  filter(LEVEL3_COD %in% c("ALU", "ARU", "BER", "GNL")) %>%
  select(LEVEL3_COD, LEVEL3_NAM)

scatterplot_data_clean %>% # Check which polygons have prop_aboveground == 1
  filter(prop_aboveground == 1) # ARU, BER, NLA
  
bee_twgd %>% # Get country names for polygons which prop_aboveground == 1
  filter(LEVEL3_COD %in% c("ARU", "BER", "NLA")) %>%
  select(LEVEL3_COD, LEVEL3_NAM)

# Fit GLS models using cleaned dataset
# gls_social <- gls(prop_social ~ abs_lat, data = scatterplot_data_clean,
#                   correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
#                   method = "REML")
# 
# gls_aboveground <- gls(prop_aboveground ~ abs_lat, data = scatterplot_data_clean,
#                        correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
#                        method = "REML")

# Fit GLS models using cleaned dataset: quadratic version since data has multiple modes
gls_social <- gls(prop_social ~ abs_lat + I(abs_lat^2),
                  data = scatterplot_data_clean,
                  correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
                  method = "REML")

gls_aboveground <- gls(prop_aboveground ~ abs_lat + I(abs_lat^2),
                       data = scatterplot_data_clean,
                       correlation = corSpher(form = ~ lon + lat, nugget = TRUE),
                       method = "REML")


summary(gls_social)
summary(gls_aboveground)

# Save/read CSV
write.csv(scatterplot_data_clean, file.path(data_wd, "scatterplot_data_clean.csv"), row.names = FALSE)
scatterplot_data_clean <- read.csv(file.path(data_wd, "scatterplot_data_clean.csv"))

# Extract predicted values from GLS models
scatterplot_data_clean$gls_fit_social <- predict(gls_social)
scatterplot_data_clean$gls_fit_aboveground <- predict(gls_aboveground)

# ------------------------------------------------------------------------------
# Plot scatterplots
# ------------------------------------------------------------------------------

plot_trait_scatter <- function(data, trait_var, color_var, y_label, viridis_option = "viridis", show_size_legend = TRUE, gls_fit_var = NULL) {
  p <- ggplot(data, aes(x = abs_lat, y = !!sym(trait_var))) +
    geom_point(aes(size = spp_rich, color = !!sym(color_var)), alpha = 0.8) +
    # Loess line
    geom_smooth(mapping = aes(), 
                se = FALSE, 
                formula = y ~ x, 
                method = "loess", 
                span = 1.5,
                color = scales::alpha("darkgray", 0.8), 
                linewidth = 1)
  # Optional: include line from GLS models
  if (!is.null(gls_fit_var)) {
    p <- p + geom_line(aes(y = !!sym(gls_fit_var)), 
                       color = scales::alpha("gray", 0.8), 
                       linetype = "dashed", 
                       linewidth = 1)
  }
  
  p +
    # Color scale (no legend shown by default)
    scale_color_viridis_c(
      option = viridis_option,
      limits = c(0, 1),
      guide = "none"
    ) +
    # Size scale
    scale_size_continuous(
      name = if (show_size_legend) "Species\nRichness" else NULL,
      breaks = if (show_size_legend) c(5, 50, 150, 300, 500, 700) else waiver(),
      limits = range(data$spp_rich, na.rm = TRUE),
      guide = if (show_size_legend) "legend" else "none"
    ) +
    # Labels
    labs(x = "Absolute Latitude", y = y_label) +
    # Theme
    theme_classic() +
    theme(
      aspect.ratio      = 1,
      axis.title        = element_text(size = 12),
      axis.text         = element_text(size = 10),
      legend.title      = element_text(size = 12),
      legend.text       = element_text(size = 10),
      legend.key.height = unit(0.3, "in"),
      legend.key.width  = unit(0.3, "in")
    )
}

social_scatter <- plot_trait_scatter(
  data = scatterplot_data_clean,
  trait_var = "prop_social",
  color_var = "prop_social",
  y_label = "Proportion Social",
  viridis_option = "cividis",
  show_size_legend = TRUE,
  gls_fit_var = "gls_fit_social"
) +
  theme(
    legend.position = c(1.3, 0.7),        # top right
    legend.justification = "right",
    legend.background = element_rect(fill = "white", color = NA)
  )

nesting_scatter <- plot_trait_scatter(
  data = scatterplot_data_clean,
  trait_var = "prop_aboveground",
  color_var = "prop_aboveground",
  y_label = "Proportion Above-Ground Nesting",
  viridis_option = "cividis",
  show_size_legend = FALSE,
  gls_fit_var = "gls_fit_aboveground"
)


combined_scatterplots <- social_scatter + nesting_scatter +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

quartz()
print(combined_scatterplots)

ggsave(
  filename = file.path(wd, "plots", "combined_scatterplots.pdf"),
  plot = combined_scatterplots,
  width = 10,
  height = 5,
  dpi = 1500
)

ggsave(
  filename = file.path(wd, "plots", "combined_scatterplots.png"),
  plot = combined_scatterplots,
  width = 10,
  height = 5,
  dpi = 1500
)


# ------------------------------------------------------------------------------
# Assemble four-panel figure of heatmaps + scatterplots
# ------------------------------------------------------------------------------
combined_mapscatter <- (
  (prop_social_heatmap + social_scatter + plot_layout(widths = c(1, 1))) /
    (prop_aboveground_heatmap + nesting_scatter + plot_layout(widths = c(1, 1))))

quartz()
print(combined_mapscatter)

ggsave(
  filename = file.path(wd, "plots", "combined_mapscatter.png"),
  plot = combined_mapscatter,
  width = 12,
  height = 8,
  dpi = 1500
)

ggsave(
  filename = file.path(wd, "plots", "combined_mapscatter.pdf"),
  plot = combined_mapscatter,
  width = 12,
  height = 8,
  dpi = 1500
)

# ------------------------------------------------------------------------------
# Archive: unused code
# ------------------------------------------------------------------------------

# This chunk was originally where we load the bee trait dataset,
# but we weren't using it since we already had the thinned points Thais made
# gbif_occurrence_data <- fread(file.path(data_wd, "05_cleaned_database.csv"))
# # gbif_occurrence_data contains distribution data for all bees with data in GBIF
# 
# # Subset gbif occurrence data to only include species scored for nesting/sociality
# gbif_subset <- gbif_occurrence_data[gbif_occurrence_data$species %in% tree_spp_traits$species, 
#                                     c("species", "decimalLatitude", "decimalLongitude")]
# colnames(gbif_subset) <- c("species", "lat", "lon")
# write.csv(gbif_subset, file.path(data_wd, "gbif_subset.csv"), row.names = FALSE) 
# # gbif_subset contains occurrence (lat, lon) information for phylogeny tip species
# 
# # Checking data
# length(unique(gbif_subset$species)) # 3748 species (compared to 4293 in tree_spp_traits),
# # so some phylogeny tip species aren't in gbif_occurrence_data/lack occurrence points
# 
# length(unique(gbif_occurrence_data$species)) # was originally 11,607
# 
# # Once we finish scoring all bees from the Dorey et al. 2023 dataset, we can use gbif_occurrence_data instead
# # Now to use this subset to make plots.

