# rm(list=ls())
library(corHMM)
library(OUwie)
library(parallel)
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")

source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)

# Let's take bio1 and bio12 (temperature and precipitation) and seasonalities (bio4 and bio15)
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_1_","bio_4_","bio_12_","bio_15_"),collapse="|"), all_climatic_vars)]
climatic_list <- lapply(all_climatic_vars, read.csv)

# Now merge everything in one table
merged_climatic_vars <- climatic_list[[1]] 
for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}
# Select only mean columns
merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# And finally merge to the trait data
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

# log all continuous variables:
merged_traits$mean_bio_1 <- log((merged_traits$mean_bio_1)+273) # transform celcius to kelvin for temperature
merged_traits$mean_bio_12 <- log(merged_traits$mean_bio_12)
merged_traits$mean_bio_15 <- log(merged_traits$mean_bio_15)
merged_traits$mean_bio_4 <- log(merged_traits$mean_bio_4)

#--------------------------------------
# LIFE HISTORY TRAITS VS. MEAN ANNUAL TEMPERATURE (BIO1)
one.full.houwie.run(dat=merged_traits[,c("tips","sociality","mean_bio_1")], phy=phy, model_names="sociality_bio1_run1")
one.full.houwie.run(dat=merged_traits[,c("tips","nest","mean_bio_1")], phy=phy, model_names="nest_bio1_run1")
