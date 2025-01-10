# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 KRUSKAL-WALLIS TEST

################################################################################

# Setup:

rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")

# Loading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


################################################################################

# Kruskal-Wallis test (ignoring phylogeny) for nesting:

# Testing nesting type ~ all is significant
# Chi-squared = 3.8552, p-value = 0.04959 

pdf("plots/nesting_boxplots_kruskal.pdf")
sink("results/nesting_result_kruskal.txt")

# Loop through the climatic variables
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  colnames(merged_table)
  
  # Boxplot
  boxplot(merged_table[,9] ~ merged_table$nest_binary, xlab="Nest Type", ylab="Environmental Variable")
  title(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]))
  
  # Prepare data for Kruskal-Wallis test
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[,9]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label)) # print description of analysis being performed
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ nests)
  print(kruskal_result)
  
}

sink()
dev.off()


################################################################################

# Kruskal-Wallis test (ignoring phylogeny) for sociality:

# Testing sociality ~ all is significant
# Chi-squared = 14.431, p-value = 0.0001454

pdf("plots/sociality_boxplots_kruskal.pdf")
sink("results/sociality_result_kruskal.txt")

# Loop through the climatic variables
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  
  # Boxplot
  boxplot(merged_table[,9] ~ merged_table$sociality_binary, xlab="Sociality", ylab="Env. Var")
  title(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]))
  
  # Prepare data for Kruskal-Wallis test
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips
  
  one_clim_var <- merged_table[,9]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("sociality ~ ", label)) 
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ sociality)
  print(kruskal_result) 
}

sink()
dev.off()
