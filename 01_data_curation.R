# ==============================================================================
# 01. Data curation
# ==============================================================================
# - Loads raw trait data and simplifies into binary (and 3-state, separating parasites)
#   categories.
# - Removes species with unknown/polymorphic sociality/nesting behavior.
# - Saves curated trait table and a tip-matched, pruned tree.
# - Plots tip pies for each trait encoding (binary and three-state).
# ==============================================================================

#-------------------------------------------------------------------------------
# Setup: working directory and libraries
#-------------------------------------------------------------------------------
# rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")

library(ape)
library(phytools)

#-------------------------------------------------------------------------------
# Load raw trait data and keep only required columns
#-------------------------------------------------------------------------------
# Load raw trait data
bees <- read.csv("original_data/all_tips_bee_tree_11-25-24.csv")

# Extract columns we will need for the binary or 3-state scoring
bees <- bees[, c("family", "tribe", "tips",
                 "broad_sociality..solitary..social..kleptoparasite..social.parasite.",
                 "nesting_broad",
                 "parasite_nesting")]

colnames(bees) <- c("family", "tribe", "tips", "sociality", "nest", "parasite_nesting")

# Quick peek
head(bees)
table(bees$sociality) 

# Sample size before removing unknowns
nrow(bees) # 4651

#-------------------------------------------------------------------------------
# Sociality recode: Binary (solitary vs social)
#-------------------------------------------------------------------------------
# Initialize column
bees$sociality_binary <- NA

# Create objects that describe what counts as solitary, social, or drop
solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal","kleptoparasite")
social   <- c("social","communal/semisocial","communal","social parasite")
drop     <- c("unknown","polymorphic","")

# Assign "solitary" to all rows of sociality_binary where
# the sociality column matches any of the values in the solitary object
bees$sociality_binary[which(bees$sociality%in%solitary)] <- "solitary" 

# Do the same for social
bees$sociality_binary[which(bees$sociality%in%social)] <- "social"

# Mark unknowns/polymorphic species to drop
bees$sociality_binary[which(bees$sociality%in%drop)] <- "drop"

# Check counts
table(bees$sociality_binary)
# drop   social solitary 
# 113     1056     3482 

#-------------------------------------------------------------------------------
# Sociality recode: Three-state (solitary, social, parasite)
#-------------------------------------------------------------------------------
# Similar mapping but keep “parasite” as its own state

bees$sociality_three_states <- NA

solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal")
social <- c("social","communal/semisocial","communal")
drop <- c("unknown","polymorphic","")
parasite <- c("social parasite","kleptoparasite")

bees$sociality_three_states[which(bees$sociality%in%solitary)] <- "solitary"
bees$sociality_three_states[which(bees$sociality%in%social)] <- "social"
bees$sociality_three_states[which(bees$sociality%in%drop)] <- "drop"
bees$sociality_three_states[which(bees$sociality%in%parasite)] <- "parasite"

# Check counts
table(bees$sociality_three_states)
# drop  parasite   social solitary 
# 113      504     1025     3009 

#-------------------------------------------------------------------------------
# Remove “drop” rows (unknowns/polymorphic)
#-------------------------------------------------------------------------------
bees <- subset(bees, bees$sociality_three_states!="drop") 
bees <- subset(bees, bees$sociality_binary!="drop") 

#-------------------------------------------------------------------------------
# Nesting recode: Three-state (ground, aboveground, parasite)
#-------------------------------------------------------------------------------
# Lumping (simplifying) nesting categories into THREE groups

bees$nest_three_states <- NA

ground <- c("ground","ground ","ground surface")
aboveground <- c("above-ground","above ground")
drop <- c("unknown","possibly ground?","","above-ground or ground","variable")
parasite <- c("social parasite","kleptoparasite")

bees$nest_three_states[which(bees$nest%in%ground)] <- "ground"
bees$nest_three_states[which(bees$nest%in%aboveground)] <- "aboveground"
bees$nest_three_states[which(bees$nest%in%drop)] <- "drop"
bees$nest_three_states[which(bees$nest%in%parasite)] <- "parasite"

# Check counts
table(bees$nest_three_states)
# aboveground   drop      ground    parasite 
# 1331          140        2563         504 

#-------------------------------------------------------------------------------
# Nesting recode: Binary (ground vs aboveground)
#-------------------------------------------------------------------------------
# For parasites, replace their nesting state with that of their host
bees$nest[which(bees$nest=="kleptoparasite")] <- bees$parasite_nesting[bees$nest=="kleptoparasite"]
bees$nest[which(bees$nest=="social parasite")] <- bees$parasite_nesting[bees$nest=="social parasite"]

bees$nest_binary <- NA

ground <- c("ground","ground ","ground surface")
aboveground <- c("above-ground","communal/semisocial","communal")
drop <- c("unknown","possibly ground?","","above-ground or ground","variable")

bees$nest_binary[which(bees$nest%in%ground)] <- "ground"
bees$nest_binary[which(bees$nest%in%aboveground)] <- "aboveground"
bees$nest_binary[which(bees$nest%in%drop)] <- "drop"

# Check counts
table(bees$nest_binary)
#aboveground  drop      ground 
#1346         245        2947 

#-------------------------------------------------------------------------------
# Remove “drop” rows (unknown/ambiguous) for nesting
#-------------------------------------------------------------------------------
bees <- subset(bees, bees$nest_binary!="drop") 
bees <- subset(bees, bees$nest_three_states!="drop")

#-------------------------------------------------------------------------------
# Sample size after curation + drop removal
#-------------------------------------------------------------------------------
nrow(bees)  # expected: 4293

#-------------------------------------------------------------------------------
# Save curated traits (keep identifiers + recoded fields only)
#-------------------------------------------------------------------------------
colnames(bees)
# Remove raw columns that are no longer needed (sociality, nest, parasite_nesting)
bees <- bees[, -c(4:6)]
colnames(bees)

write.csv(bees, "curated_data/bees_traits.csv", row.names=F)

#-------------------------------------------------------------------------------
# Prune tree to match curated species list and save
#-------------------------------------------------------------------------------
bee_tree <- read.tree("original_data/ML_beetree.tre")
bee_tree <- keep.tip(bee_tree, bees$tips)  # drop non-overlapping tips
write.tree(bee_tree, "curated_data/ML_beetree_pruned.tre")


#-------------------------------------------------------------------------------
# Visualize trait distributions on the pruned tree (tip pies)
#-------------------------------------------------------------------------------
# Align trait table to pruned tree tip order
colnames(bees)
group_traits <- bees[bees[,"tips"] %in% bee_tree$tip.label,]  # filters the bees data frame to include only those rows with tips that match bee_tree$tip.label
group_traits <- group_traits[order(match(group_traits[,"tips"],bee_tree$tip.label)),] # orders this filtered data frame to align the tips with the order of bee_tree$tip.label

#-------------------------
# Sociality: three states
#-------------------------
mode <- group_traits[, "sociality_three_states"]
names(mode) <- group_traits[, "tips"]
colors_states <- c("midnightblue", "goldenrod", "darkred")  # solitary, social, parasite

pdf("plots/bee_sociality_three_states.pdf", width = 4, height = 45)
plot(bee_tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)
par(fg = "transparent")
tiplabels(pie = to.matrix(mode, sort(unique(mode))),
          piecol = colors_states, cex = 0.1, lwd = 0.2, frame = "n")
par(fg = "black")
legend("topleft", legend = sort(unique(mode)), pt.bg = colors_states, pch = 21, cex = 0.8)
axisPhylo()
dev.off()

#---------------------
# Sociality: binary
#---------------------
mode <- group_traits[, "sociality_binary"]
names(mode) <- group_traits[, "tips"]
colors_states <- c("midnightblue", "goldenrod")  # solitary, social

pdf("plots/bee_sociality_binary.pdf", width = 4, height = 45)
plot(bee_tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)
par(fg = "transparent")
tiplabels(pie = to.matrix(mode, sort(unique(mode))),
          piecol = colors_states, cex = 0.1, lwd = 0.2, frame = "n")
par(fg = "black")
legend("topleft", legend = sort(unique(mode)), pt.bg = colors_states, pch = 21, cex = 0.8)
axisPhylo()
dev.off()

#---------------------
# Nesting: binary
#---------------------
mode <- group_traits[, "nest_binary"]
names(mode) <- group_traits[, "tips"]
colors_states <- c("midnightblue", "goldenrod")  # ground, aboveground

pdf("plots/bee_nest_binary.pdf", width = 4, height = 45)
plot(bee_tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)
par(fg = "transparent")
tiplabels(pie = to.matrix(mode, sort(unique(mode))),
          piecol = colors_states, cex = 0.1, lwd = 0.2, frame = "n")
par(fg = "black")
legend("topleft", legend = sort(unique(mode)), pt.bg = colors_states, pch = 21, cex = 0.8)
axisPhylo()
dev.off()

#-------------------------
# Nesting: three states
#-------------------------
mode <- group_traits[, "nest_three_states"]
names(mode) <- group_traits[, "tips"]
colors_states <- c("midnightblue", "goldenrod", "darkred")  # ground, aboveground, parasite

pdf("plots/bee_nest_three_states.pdf", width = 4, height = 45)
plot(bee_tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)
par(fg = "transparent")
tiplabels(pie = to.matrix(mode, sort(unique(mode))),
          piecol = colors_states, cex = 0.1, lwd = 0.2, frame = "n")
par(fg = "black")
legend("topleft", legend = sort(unique(mode)), pt.bg = colors_states, pch = 21, cex = 0.8)
axisPhylo()
dev.off()

