# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(phytools)

# Some exploratory analyses
# Reloading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv")

#traits <- subset(traits, traits$sociality!="parasite")
#traits <- subset(traits, grepl(paste(c("Lasioglossum"), collapse="|"), traits$tips))
#traits <- subset(traits, traits$tribe=="Xylocopini")

# Let's run a phylANOVA between sociality and env variables
# and plot boxplots for visualization:

climate_index <- 1
pdf("plots/sociality_boxplots_all_vars.pdf")
sink("sociality_result_phyloANOVA.txt")
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/",all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means (should remove these points before this step)
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # making sure every dataset include the same species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  #if(length(unique(merged_table$nest)) > 1) {
    boxplot(merged_table[,7]~merged_table$sociality, xlab="nest", ylab="env var")
    title(gsub("_climate_summstats.csv","",all_climatic_vars[climate_index]))
    
    sociality <- merged_table$sociality
    names(sociality) <- merged_table$tips
    
    one_clim_var <- merged_table[,7]
    names(one_clim_var) <- merged_table$tips
    
    label <- gsub("_climate_summstats.csv","",all_climatic_vars[climate_index])
    print(paste0("sociality ~ ",label))
    results <- phylANOVA(subset_tree, sociality, one_clim_var, p.adj="bonferroni") 
    print(results)
  #}
}
sink()
dev.off()

# nests
pdf("plots/nesting_boxplots_all_vars.pdf")
sink("nesting_result_phyloANOVA.txt")
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/",all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means (should remove these points before this step)
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # making sure every dataset include the same species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  #if(length(unique(merged_table$nest)) > 1) {
  boxplot(merged_table[,7]~merged_table$nest, xlab="nest", ylab="env var")
  title(gsub("_climate_summstats.csv","",all_climatic_vars[climate_index]))
  
  nests <- merged_table$nest
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[,7]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv","",all_climatic_vars[climate_index])
  print(paste0("nesting type ~ ",label))
  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj="bonferroni") 
  print(results)
  #}
}
sink()
dev.off()

# combination of nests and sociality
#pdf("plots/nesting_boxplots_all_vars.pdf")
#sink("nesting_result_phyloANOVA.txt")
traits$comb_nest_soc <- paste(traits$sociality, traits$nest, sep="_")
#table(traits$comb_nest_soc)
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/",all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means (should remove these points before this step)
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # making sure every dataset include the same species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  #if(length(unique(merged_table$nest)) > 1) {
  boxplot(merged_table[,7]~merged_table$comb_nest_soc, xlab="nest", ylab="env var")
  title(gsub("_climate_summstats.csv","",all_climatic_vars[climate_index]))
  
  nests <- merged_table$comb_nest_soc
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[,7]
  names(one_clim_var) <- merged_table$tips
  
  label <- gsub("_climate_summstats.csv","",all_climatic_vars[climate_index])
  print(paste0("nesting type ~ ",label))
  results <- phylANOVA(subset_tree, nests, one_clim_var, p.adj="bonferroni") 
  print(results)
  #}
}
sink()
dev.off()
