# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

################################################################################

# Setup
rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)
library(ape)

# Loading traits and tree
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
head(traits)
simple_traits <- traits[, c("tips", "sociality_binary", "nest_binary")]
head(simple_traits)

tree # 4293 tips

nrow(simple_traits) # 4293

# Check which species in the new_traits$tips are not in the tree
missing_in_tree <- setdiff(simple_traits$tips, tree$tip.label)
missing_in_tree # character(0)

# Check which species in the tree are not in the new_traits$tips
missing_in_traits <- setdiff(tree$tip.label, simple_traits$tips)
missing_in_traits # character(0)

# We want to:
# 1) Do ancestral state reconstruction of nesting and sociality using sociality_binary and nest_binary columns
#       Character states are social, solitary, and ground, aboveground
# 2) Create a circular phylogeny where the branches are colored by character state (social, solitary, ground, aboveground)
#       Which means 4 different colors/character state combinations 

# Convert sociality_binary into a factor
sociality_factor <- factor(simple_traits$sociality_binary, levels = c("solitary", "social"))

# Perform ancestral state reconstruction for sociality
sociality_ace <- ace(sociality_factor, tree, type = "discrete")

# Extract the ancestral states for the tips and internal nodes
ancestral_states <- sociality_ace$lik.anc

# Assign colors to the states
state_colors <- c("solitary" = "red", "social" = "blue")

# Color the tips based on the sociality states
tip_colors <- state_colors[as.character(simple_traits$sociality_binary)]

# Color the internal nodes based on the reconstructed ancestral states
node_colors <- state_colors[as.character(ancestral_states[, 1])]  # Using the first column of ancestral states for internal nodes

# Plot the tree in a circular (fan) style, using phytools' plotTree function
plotTree(tree, type = "fan", tip.color = tip_colors, node.color = node_colors, 
         main = "Circular Phylogeny with Ancestral State Reconstruction for Sociality", 
         show.tip.label = FALSE, show.node.label = FALSE, 
         cex = 0.8, # Increase the size of the tree
         edge.width = 2)  # Adjust edge width for better visibility

# Add a legend to the plot for clarity
legend("topright", legend = c("solitary", "social"), fill = c("red", "blue"))
