library(corHMM)

# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

# Some exploratory analyses
# Reloading traits, tree and climatic data:
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv")

traits <- traits[,c(3:5)]

# traits <- subset(traits, grepl("Xylocopa",traits$tips))
# tree <- keep.tip(tree, which(tree$tip.label%in%traits$tips))

# Correlation between nesting habit and sociality
corhmm_fits <- corHMM:::fitCorrelationTest(tree, traits) 
save(corhmm_fits, file = "results/corhmm_fits.Rsave")
load("results/corhmm_fits.Rsave")
corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)

# conducting a lrt
teststat <- -2 * (corhmm_tbl$lnLik[1] - corhmm_tbl$lnLik[2])
p.val50 <- pchisq(teststat, df = 1, lower.tail = FALSE)
print(p.val50)
