# rm(list=ls())
library(corHMM)
library(OUwie)
library(parallel)
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")

source("00_utility_functions.R")
#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)

# Let's take bio1 and bio12 (temperature and precipitation) and seasonalities (bio4 and bio15)
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_1_","bio_4_","ai"),collapse="|"), all_climatic_vars)]
climatic_list <- lapply(all_climatic_vars, read.csv)

# Now merge everything in one table
merged_climatic_vars <- climatic_list[[1]] 
for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}
# Select only mean columns
merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# And finally merge to the trait data
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

while(length(unique(merged_traits$sociality))<3) {
  merged_traits <- merged_traits[sample(1:nrow(merged_traits), 100),]
}

# log all continuous variables:
merged_traits$mean_bio_1 <- log((merged_traits$mean_bio_1)+273) # transform celcius to kelvin for temperature
merged_traits$mean_bio_4 <- log(merged_traits$mean_bio_4)
merged_traits$mean_awi_pm_sr_yr <- log(merged_traits$mean_awi_pm_sr_yr)

merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_1))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_4))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_awi_pm_sr_yr))

phy <- keep.tip(phy, which(phy$tip.label %in% merged_traits$tips))

#------------------------------------
# determining corHMM models with corhmm dredge
# dredge_sociality <- corHMM:::corHMMDredge(phy, merged_traits[,c("tips","sociality")],max.rate.cat=3)
# save(dredge_sociality, file="corhmm_dredge_sociality.Rsave")
# corhmm_tbl_sociality <- corHMM:::getModelTable(dredge_sociality)
# write.csv(corhmm_tbl_sociality, file="corhmm_tbl_sociality.csv")
# 
# dredge_nesting <- corHMM:::corHMMDredge(phy, merged_traits[,c("tips","nest")],max.rate.cat=3)
# save(dredge_nesting, file="corhmm_dredge_nesting.Rsave")
# corhmm_tbl_nesting <- corHMM:::getModelTable(dredge_nesting)
# write.csv(corhmm_tbl_nesting, file="corhmm_tbl_nesting.csv")
#--------------------------------------
# LIFE HISTORY TRAITS VS. MEAN ANNUAL TEMPERATURE (BIO1)

load("corhmm_dredge_sociality.Rsave")
corhmm_tbl_sociality <- read.csv("corhmm_tbl_sociality.csv")
disc_model_soc <- dredge_sociality[[which.min(corhmm_tbl_sociality$AIC)]]$index.mat

dat=merged_traits[,c("tips","sociality","mean_bio_1")]
phy=phy 
disc_model = disc_model_soc 
model_names="sociality_bio1_run3"

# organize data
shared_species <- intersect(dat$tips, phy$tip.label)
dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]
# model structure
#corhmm_set <- corhmm.model.setup(dat[,c(1,2)])
corhmm_set <- disc_model 

model_list <- houwie.model.setup(corhmm_set, model_names)
# setting up parallel

quickFunc <- function(model_list, model_name){
  res <- hOUwie(phy, dat, model_list[[1]], model_list[[2]], model_list[[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = FALSE, n_starts = 10, ncores = 10)
  file.name <- paste0("houwie_results/",model_name, ".Rsave")
  save(res, file=file.name)
}

mclapply(1:6, function(x) quickFunc(model_list[[x]], names(model_list)[x]), mc.cores = 10)



# load("corhmm_dredge_nesting.Rsave")
# corhmm_tbl_sociality <- read.csv("corhmm_tbl_nesting.csv")
# disc_model_nest <- dredge_nesting[[which.min(corhmm_tbl_nesting$AIC)]]$index.mat
# one.full.houwie.run(dat=merged_traits[,c("tips","nest","mean_bio_1")],
#                     phy=phy, disc_model = disc_model_nest, model_names="nest_bio1_run1")
