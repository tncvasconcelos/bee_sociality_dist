# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(ape)
library(phytools)

#--------------------------------
# Loading trait data (sociality)
bees <- read.csv("original_data/all_tips_bee_tree_Dec4_2023.csv")
# Some data curation:
bees <- bees[,c(1,3,5,6)] 
# head(bees)
colnames(bees) <- c("family","tribe","tips","sociality")
# categories <- unique(bees$sociality)

# Lumping (simplifying) sociality categories into a few groups:
solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal")
social <- c("social","communal/semisocial","communal")
drop <- c("unknown","polymorphic","")
parasite <- c("social parasite","kleptoparasite")
bees$sociality[which(bees$sociality%in%solitary)] <- "solitary"
bees$sociality[which(bees$sociality%in%social)] <- "social"
bees$sociality[which(bees$sociality%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
bees$sociality[which(bees$sociality%in%parasite)] <- "parasite"

# Species in each category:
# table(bees$sociality)
# drop parasite   social solitary 
# 50      414     1029     3160

bees <- subset(bees, bees$sociality!="drop") # removing ones to drop
# Saving curated dataset:
write.csv(bees, "curated_data/bees_sociality.csv", row.names=F)

#--------------------------------
# Pruning tree so that it matches species with data:
bee_tree <- read.tree("original_data/ML_beetree.tre")
bee_tree <- keep.tip(bee_tree, bees$tips)
# Saving pruned tree:
write.tree(bee_tree, "curated_data/ML_beetree_pruned.tre")

#--------------------------------
# Plot to visualize sociality distribution in the tree:
group_traits <- bees[bees[,1] %in% bee_tree$tip.label,]
group_traits <- group_traits[order(match(group_traits[,1],bee_tree$tip.label)),]

mode <- group_traits[,2]
names(mode) <- group_traits[,1]
colors_states <- c("midnightblue", "goldenrod","darkred")

pdf("plots/bee_sociality.pdf", width= 4, height= 40)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.2,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()






