# ==============================================================================
# 03. phylANOVA
# ==============================================================================
# Run phylogenetic ANOVA (phylANOVA) to test associations between 
# bee traits (sociality, nesting) and climatic variables.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load required library
# ------------------------------------------------------------------------------
rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)


# ------------------------------------------------------------------------------
# Load trait data, tree, and list of climate summary statistic files
# ------------------------------------------------------------------------------
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


# ------------------------------------------------------------------------------
# Run phylANOVA between sociality_binary and each climate variable
# ------------------------------------------------------------------------------
sink("results/phylANOVA_sociality_results.txt") # Output file

# Loop that iterates over each environmental variable in the all_climatic_vars vector
for(climate_index in 1:length(all_climatic_vars)) {
  
  # Load and clean climate data
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index])) # reading climate data corresponding to current index
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means
  
  # Identify species with complete data in all datasets
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Subset all datasets to sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Merge traits and climate data by species
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
    
  # Prepare named vectors for phylANOVA
  sociality <- merged_table$sociality_binary 
  names(sociality) <- merged_table$tips
    
  one_clim_var <- merged_table[,9] # column 9 is the mean climate value for that variable
  names(one_clim_var) <- merged_table$tips 
    
  # Print description of analysis being performed in each loop iteration
  label <- gsub("_climate_summstats.csv","", all_climatic_vars[climate_index])  # removing character string from names of climate variables
  print(paste0("sociality ~ ",label)) 
  
  # Run phylANOVA and save to results object
  results <- phylANOVA(subset_tree, sociality, one_clim_var, p.adj="bonferroni")
  print(results)
}

sink()


# ------------------------------------------------------------------------------
# Run phylANOVA between nesting_binary and each climate variable
# ------------------------------------------------------------------------------
sink("results/phylANOVA_nesting_results.txt")

for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  colnames(merged_table)
  
  boxplot(merged_table[,9]~merged_table$nest_binary, xlab="nest", ylab="env var")
  title(gsub("_climate_summstats.csv","", all_climatic_vars[climate_index]))
  
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[,9]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv","", all_climatic_vars[climate_index])
  print(paste0("nesting type ~ ",label))
  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj="bonferroni") 
  print(results)
}

sink()
dev.off()


# ------------------------------------------------------------------------------
# Run phylANOVA using 4-level combination of sociality and nesting traits
# ------------------------------------------------------------------------------
sink("results/phyloANOVA_combined_results.txt")

traits$comb_nest_soc <- paste(traits$sociality_binary, traits$nest_binary, sep="_")
table(traits$comb_nest_soc)

for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  colnames(merged_table)
  
  nests <- merged_table$comb_nest_soc
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[,9]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv","", all_climatic_vars[climate_index])
  print(paste0("nesting type ~ ",label))
  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj="bonferroni") 
  print(results)
}

sink()
