# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 PHYLOGENETIC GENERALISED LINEAR MODEL

# Setup
rm(list = ls())
setwd("~/Desktop/bee_sociality_dist")
library(phytools)
library(geiger)
library(dplyr)

# Loading traits, tree, and climatic data:
traits <- read.csv("curated_data/bees_traits.csv") # contains trait data for bees in tree
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
climate_files <- list.files("curated_data", "summstats.csv", full.names = TRUE) 

# Create separate traits dataset for use in loop to merge with climate data
traits_tmp <- traits

for (file_index in 1:length(climate_files)) {
  # Read the climate file
  temp_data <- read.csv(climate_files[file_index])
  
  # Extract relevant columns
  temp_data <- temp_data[, c(1, grep("mean", colnames(temp_data)))]
  
  # Remove NAs from the climate data
  temp_data <- temp_data[complete.cases(temp_data), ]
  
  # Merge with the traits data
  traits_tmp <- merge(traits_tmp, temp_data, by.x = "tips", by.y = "species", all = FALSE)
}

# Inspect the merged data which now contains trait data and mean climate values
head(traits_tmp)
colnames(traits_tmp)

# Fit phyloGLM for each climate variable
results <- list()  # Store results for each model
climate_vars <- grep("mean", colnames(filtered_traits), value = TRUE) # Identify climate variables
print(climate_vars)

for (climate_var in climate_vars) {
  # Create formula
  formula <- as.formula(paste(climate_var, "~ sociality_binary + nest_binary"))
  
  # Fit phyloGLM
  model <- phyloglm(formula, phy = tree, data = filtered_traits, method = "logistic_MPLE")
  
  # Store model summary in the results list
  results[[climate_var]] <- summary(model)
}

# Inspect results for the first climate variable
print(results[[climate_vars[1]]])

# Diagnostics (example for one model)
plot(model)





