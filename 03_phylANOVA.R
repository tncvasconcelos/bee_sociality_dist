# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 PHYLANOVA

################################################################################

# Setup:

#rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)

# Loading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


################################################################################

# Run phylANOVA between sociality and climate variables using sociality_binary

# Open a sink for results
sink("results/phylANOVA_sociality_results.txt")

# Loop that iterates over each environmental variable in the all_climatic_vars vector
for(climate_index in 1:length(all_climatic_vars)) {
  
  # Read the climate data for the current variable
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index])) # reading climate data corresponding to current index
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label) # find species that are present in the trait data, climate data, and tree
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
    
  # Prepare data for phylANOVA
  sociality <- merged_table$sociality_binary 
  names(sociality) <- merged_table$tips # making named vector for phylANOVA
    
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

################################################################################

# phylANOVA between nesting and climate variables using nesting_binary:

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

################################################################################

# phylANOVA for the combination of nests and sociality:

sink("results/phyloANOVA_nesting*sociality_results.txt")

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
