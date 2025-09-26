
setwd("~/Desktop/bee_sociality_dist")
library(corHMM)
devtools::install_github("thej022214/corHMM")

#-------------------------------------------------------------------------------
# Organizing dataset
#-------------------------------------------------------------------------------
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")

# Drop tips in the tree that don’t match any in the trait data
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

# Find species shared by both datasets
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

# Reorder the data and tree to align species
dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

# Keep relevant trait columns
dat <- dat[,c("tips","sociality_binary","nest_binary")]
load("corHMMdredge_results/corhmm_dredge_binary.Rsave")


data <- dat

model <- dredge_sociality[[10]]$solution

##get simmap from corhmm solution
simmaps <- makeSimmap(tree=phy, data=dat, model=model, rate.cat=2, nSim=100, nCores=1)

# get the summary
simmap_summaries <- lapply(simmaps, summarize_single_simmap)
summary_df <- summarize_transition_stats(simmap_summaries)

print(summary_df)
write.csv(summary_df, file="summary_corhmm.csv", row.names = F)

plot_transition_summary(simmap_summaries)
