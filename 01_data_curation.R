# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(ape)
library(phytools)

#--------------------------------
# Loading trait data (sociality and nest strategy)
bees <- read.csv("original_data/all_tips_bee_tree_sep4.csv")
# Some data curation:
bees <- bees[,c("family","tribe","tips","broad_sociality..solitary..social..kleptoparasite..social.parasite.","nesting_broad", "parasite_nesting")]
# head(bees)
colnames(bees) <- c("family","tribe","tips","sociality","nest", "parasite_nesting")
# categories <- unique(bees$sociality)

# Lumping (simplifying) sociality categories into a few groups:
solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal","kleptoparasite")
social <- c("social","communal/semisocial","communal","social parasite")
drop <- c("unknown","polymorphic","")
#parasite <- c("social parasite","kleptoparasite")
bees$sociality[which(bees$sociality%in%solitary)] <- "solitary"
bees$sociality[which(bees$sociality%in%social)] <- "social"
bees$sociality[which(bees$sociality%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
#bees$sociality[which(bees$sociality%in%parasite)] <- "parasite"

# Species in each category:
table(bees$sociality)
# drop parasite   social solitary 
# 113      497     1029     3012

bees <- subset(bees, bees$sociality!="drop") # removing ones to drop

bees$nest[which(bees$nest=="kleptoparasite")] <- bees$parasite_nesting[bees$nest=="kleptoparasite"]
bees$nest[which(bees$nest=="social parasite")] <- bees$parasite_nesting[bees$nest=="social parasite"]

# Lumping (simplifying) nesting categories into a few groups:
ground <- c("ground","ground ","ground surface")
aboveground <- c("above-ground","communal/semisocial","communal")
drop <- c("unknown","possibly ground?","","above-ground or ground","variable")
#parasite <- c("social parasite","kleptoparasite")
bees$nest[which(bees$nest%in%ground)] <- "ground"
bees$nest[which(bees$nest%in%aboveground)] <- "aboveground"
bees$nest[which(bees$nest%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
#bees$nest[which(bees$nest%in%parasite)] <- "parasite"

# Species in each category:
 table(bees$nest)
#aboveground        drop      ground    parasite 
#     1521          45        2653         432 

bees <- subset(bees, bees$nest!="drop") # removing ones to drop

bees <- bees[,-6]
# Saving curated dataset:
write.csv(bees, "curated_data/bees_traits.csv", row.names=F)

#--------------------------------
# Pruning tree so that it matches species with data:
bee_tree <- read.tree("original_data/ML_beetree.tre")
bee_tree <- keep.tip(bee_tree, bees$tips)
# Saving pruned tree:
write.tree(bee_tree, "curated_data/ML_beetree_pruned.tre")

#--------------------------------
# Plot to visualize sociality and nesting distribution in the tree:
group_traits <- bees[bees[,3] %in% bee_tree$tip.label,]
group_traits <- group_traits[order(match(group_traits[,3],bee_tree$tip.label)),]

mode <- group_traits[,4]
names(mode) <- group_traits[,3]
colors_states <- c("midnightblue", "goldenrod")

pdf("plots/bee_sociality.pdf", width= 4, height= 40)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.2,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()

mode <- group_traits[,5]
names(mode) <- group_traits[,3]
colors_states <- c("midnightblue", "goldenrod")

pdf("plots/bee_nests.pdf", width= 4, height= 40)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.2,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()




