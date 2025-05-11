# rm(list=ls())
# setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(corHMM) 
source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits and tree
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
phy <- keep.tip(phy, which(phy$tip.label %in% traits$tips))

#------------------------------------
dat <- traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

sociality_scoring <- c("sociality_binary", "sociality_three_states")
nest_scoring <- c("nest_binary", "nest_three_states")

for(sociality_index in 1:length(sociality_scoring)) {
  for(nest_index in 1:length(nest_scoring)){
    dat1 <- dat[,c("tips",sociality_scoring[sociality_index],paste0(nest_scoring[nest_index]))]
    
    # Uncomment to run corHMM
    corhmm_fits <- corHMM:::fitCorrelationTest(phy, dat1)     
    corhmm_tbl <- corHMM:::getModelTable(corhmm_fits)
    
    write.csv(corhmm_tbl, file = paste0("correlation_test_results/corhmm_fits_cortest_table_",
                                        sociality_scoring[sociality_index],"_",
                                        nest_scoring[nest_index],".csv"))
    
    save(corhmm_fits, file = paste0("correlation_test_results/corhmm_fits_cortest_",
                                    sociality_scoring[sociality_index],"_",
                                    nest_scoring[nest_index],".Rsave"))
    
    

  }
}



#save(corhmm_fits, file = "corhmm_fits_cortest.Rsave")

#load("corhmm_fits_cortest.Rsave")


# LIKELIHOOD RATIO TEST
teststat <- -2 * (corhmm_tbl$lnLik[2] - corhmm_tbl$lnLik[4])
p.val <- pchisq(teststat, df = 8, lower.tail = FALSE)
