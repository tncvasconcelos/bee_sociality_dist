# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(phytools)

# Some exploratory analyses
# Reloading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_sociality.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv")

#traits <- subset(traits, traits$sociality!="parasite")
#traits <- subset(traits, grepl(paste(c("Lasioglossum"), collapse="|"), traits$tips))
#traits <- subset(traits, traits$tribe!="Centridini")

# Let's run a phylANOVA between sociality and env variables
# and plot boxplots for visualization:
pdf("plots/boxplots_all_vars.pdf")
for(climate_index in 1:length(all_climatic_vars)) {
  climate <- read.csv(paste0("curated_data/",all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3])) # removing species with NA means (should remove these points before this step)
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # making sure every dataset include the same species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips",by.y="species")
  if(length(unique(merged_table$sociality)) > 1) {
    boxplot(merged_table[,6]~merged_table$sociality, xlab="sociality", ylab="env var")
    title(gsub("_climate_summstats.csv","",all_climatic_vars[climate_index]))
    
    sociality <- merged_table$sociality
    names(sociality) <- merged_table$tips
    
    one_clim_var <- merged_table[,6]
    names(one_clim_var) <- merged_table$tips
    
    results <- phylANOVA(subset_tree, sociality, one_clim_var, p.adj="bonferroni")  
    print(results)
  }
}
dev.off()
