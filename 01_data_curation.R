# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(ape)
library(phytools)

#--------------------------------
# Loading trait data (sociality and nest strategy)
bees <- read.csv("original_data/all_tips_bee_tree_sep4.csv")
# Some data curation (extracting only columns we will need for the binary or three state scoring)
bees <- bees[,c("family","tribe","tips","broad_sociality..solitary..social..kleptoparasite..social.parasite.","nesting_broad", "parasite_nesting")]
colnames(bees) <- c("family","tribe","tips","sociality","nest", "parasite_nesting")

#--------------------------------
# sample size before removing unknowns
nrow(bees)
# [1] 4651

#--------------------------------
# Lumping (simplifying) sociality categories
#--------------------------------
# BINARY DATASET (TWO STATES)
solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal","kleptoparasite")
social <- c("social","communal/semisocial","communal","social parasite")
drop <- c("unknown","polymorphic","")
bees$sociality_binary <- NA
bees$sociality_binary[which(bees$sociality%in%solitary)] <- "solitary"
bees$sociality_binary[which(bees$sociality%in%social)] <- "social"
bees$sociality_binary[which(bees$sociality%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
# Species in each category:
table(bees$sociality_binary)
# drop   social solitary 
# 113     1056     3482 

# BROAD DATASET (THREE STATES)
# Lumping (simplifying) sociality categories into THREE groups:
solitary <- c("solitary","solitary?","solitary/semisocial","solitary/communal")
social <- c("social","communal/semisocial","communal")
drop <- c("unknown","polymorphic","")
parasite <- c("social parasite","kleptoparasite")
bees$sociality_three_states <- NA
bees$sociality_three_states[which(bees$sociality%in%solitary)] <- "solitary"
bees$sociality_three_states[which(bees$sociality%in%social)] <- "social"
bees$sociality_three_states[which(bees$sociality%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
bees$sociality_three_states[which(bees$sociality%in%parasite)] <- "parasite"
# Species in each category:
table(bees$sociality_three_states)
# drop parasite   social solitary 
# 113      504     1025     3009 

#--------------------------------
# Removing ones to drop
bees <- subset(bees, bees$sociality_three_states!="drop") 
bees <- subset(bees, bees$sociality_binary!="drop") 

#--------------------------------
# Lumping (simplifying) nesting strategy categories
#--------------------------------
# Lumping (simplifying) nesting categories into THREE groups:
ground <- c("ground","ground ","ground surface")
aboveground <- c("above-ground","communal/semisocial","communal")
drop <- c("unknown","possibly ground?","","above-ground or ground","variable")
parasite <- c("social parasite","kleptoparasite")
bees$nest_three_states <- NA
bees$nest_three_states[which(bees$nest%in%ground)] <- "ground"
bees$nest_three_states[which(bees$nest%in%aboveground)] <- "aboveground"
bees$nest_three_states[which(bees$nest%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
bees$nest_three_states[which(bees$nest%in%parasite)] <- "parasite"
# Species in each category:
table(bees$nest_three_states)
# aboveground        drop      ground    parasite 
# 1517                44        2593         497 

# Lumping (simplifying) nesting categories into TWO groups:
bees$nest_binary <- NA
bees$nest[which(bees$nest=="kleptoparasite")] <- bees$parasite_nesting[bees$nest=="kleptoparasite"]
bees$nest[which(bees$nest=="social parasite")] <- bees$parasite_nesting[bees$nest=="social parasite"]
ground <- c("ground","ground ","ground surface")
aboveground <- c("above-ground","communal/semisocial","communal")
drop <- c("unknown","possibly ground?","","above-ground or ground","variable")
bees$nest_binary[which(bees$nest%in%ground)] <- "ground"
bees$nest_binary[which(bees$nest%in%aboveground)] <- "aboveground"
bees$nest_binary[which(bees$nest%in%drop)] <- "drop" #dropping unknowns and polymorphic for now
# Species in each category:
table(bees$nest_binary)
# aboveground        drop      ground 
# 1575              104        2972 

#--------------------------------
# Removing ones to drop
bees <- subset(bees, bees$nest_binary!="drop") # removing ones to drop
bees <- subset(bees, bees$nest_three_states!="drop") # removing ones to drop

#--------------------------------
# sample size after removing unknowns
nrow(bees)
# [1] 4434
colnames(bees)
bees <- bees[,-c(4:6)]
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
colnames(bees)
group_traits <- bees[bees[,"tips"] %in% bee_tree$tip.label,]
group_traits <- group_traits[order(match(group_traits[,"tips"],bee_tree$tip.label)),]

# Sociality three states
mode <- group_traits[,"sociality_three_states"]
names(mode) <- group_traits[,"tips"]
colors_states <- c("midnightblue", "goldenrod", "darkred")

pdf("plots/bee_sociality_three_states.pdf", width= 4, height= 45)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.1,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()

# Sociality binary
mode <- group_traits[,"sociality_binary"]
names(mode) <- group_traits[,"tips"]
colors_states <- c("midnightblue", "goldenrod")

pdf("plots/bee_sociality_binary.pdf", width= 4, height= 45)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.1,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()

# Nesting binary
mode <- group_traits[,"nest_binary"]
names(mode) <- group_traits[,"tips"]
colors_states <- c("midnightblue", "goldenrod")

pdf("plots/bee_nest_binary.pdf", width= 4, height= 45)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.1,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()

# Nesting three states
mode <- group_traits[,"nest_three_states"]
names(mode) <- group_traits[,"tips"]
colors_states <- c("midnightblue", "goldenrod", "darkred")

pdf("plots/bee_nest_three_states.pdf", width= 4, height= 45)
plot(bee_tree, show.tip.label=T, edge.width=0.2, adj=1, cex=0.05)
par(fg="transparent")
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states,cex= 0.1,lwd=0.2, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
legend("topleft", legend=sort(unique(mode)), pt.bg = colors_states, pch=21, cex=0.8)
axisPhylo()
dev.off()

