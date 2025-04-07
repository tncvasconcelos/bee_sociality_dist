
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(ape)
source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")

traits <- traits[,c("tips","sociality_binary","nest_binary")]
colnames(traits)[1] <- "species"

# table(traits$sociality_binary)
# table(traits$nest_binary)

library(ggtree)
library(ggplot2)
library(dplyr)

# Combine the two traits into one four-state category
traits <- traits %>%
  filter(species %in% phy$tip.label) %>%
  mutate(
    social_nest_combo = paste(sociality_binary, nest_binary, sep = "_")
  )

row.names(traits) <- traits$species

# Define nice colors for the four categories
combo_colors <- c(
  "social_aboveground" = "#F28E2B",  # orange
  "social_ground" = "#377EB8",       # blue
  "solitary_aboveground" = "#4DAF4A", # teal
  "solitary_ground" = "#9E1D42"      # deep pink
)

# Build and plot the tree
p <- ggtree(phy, layout = "fan", size = 0.2) %<+% traits
p +
  geom_tippoint(aes(color = social_nest_combo), 
                size = 1.6, 
                alpha = 0.9) +
  scale_color_manual(values = combo_colors) +
  theme(legend.position = "right") +
  ggtitle("Phylogeny with Combined Sociality and Nesting Traits")
