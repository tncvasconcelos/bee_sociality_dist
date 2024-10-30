# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich & Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 PHYLOGENETIC GENERALISED LINEAR MODEL

# Setup
rm(list=ls())
setwd("~/Desktop/bee_sociality_dist")

library(phytools)
library(geiger)
library(dplyr)

# Loading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
climate <- list.files("curated_data", "summstats.csv", full.names = TRUE) # Selects summstats.csv files from curated_data

# Merge climate data to prepare for phyloGLM
traits_tmp <- traits 
# Loop through each file, read in the data, and extract the 'mean_bio_x' column
for (file_index in 1:length(climate)) {
  temp_data <- read.csv(climate[file_index])
  temp_data <- temp_data[,c(1, grep("mean",colnames(temp_data)))]
  traits_tmp <- merge(traits_tmp, temp_data, by.x="tips", by.y="species", all=F)
}

climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means

# Find species that are present in the trait data, climate data, and tree
sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label) 

# Create subsets of traits, climate data, and the tree that include only the sampled species
subset_traits <- subset(traits, traits$tips %in% sampled_species)
subset_climate <- subset(climate, climate$species %in% sampled_species)
subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])

# Combine traits and climate data into one data frame, now that they contain the same species
data <- cbind(subset_traits, subset_climate)

# Fit phyloGLM for each climate variable
for (climate_var in colnames(climate_data)) {
  model <- phyloGLM(as.formula(paste(climate_var, "~ sociality + nesting_biology")),
                    phy = tree, data = data)
}

# Summary of the model
summary(model)

# Check diagnostics 
plot(model)





