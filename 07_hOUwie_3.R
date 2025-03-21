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
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_1_","bio_12_","bio_4_",
                                                    "bio_15_","ai"),collapse="|"), all_climatic_vars)]
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

# merged_traits0 <- merged_traits[1,]
# while(length(unique(merged_traits0$sociality))<3) {
#   merged_traits0 <- merged_traits[sample(1:nrow(merged_traits), 100),]
# }
# merged_traits <- merged_traits0

# log all continuous variables:
merged_traits$mean_bio_1 <- log((merged_traits$mean_bio_1)+273) # transform celcius to kelvin for temperature
merged_traits$mean_bio_12 <- log(merged_traits$mean_bio_12)
merged_traits$mean_bio_15 <- log(merged_traits$mean_bio_15)
merged_traits$mean_bio_4 <- log(merged_traits$mean_bio_4)
merged_traits$mean_awi_pm_sr_yr <- log(merged_traits$mean_awi_pm_sr_yr)

merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_1))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_4))
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_awi_pm_sr_yr))

phy <- keep.tip(phy, which(phy$tip.label %in% merged_traits$tips))

#------------------------------------
dat <- merged_traits
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]
dat <- dat[,c("tips","sociality_binary","nest_binary","mean_awi_pm_sr_yr")]
dat <- dat[!dat$mean_awi_pm_sr_yr == -Inf,]

phy <- keep.tip(phy, dat$tips)
#--------------------------------------
# LIFE HISTORY TRAITS VS. ARIDITY INDEX

load("corhmm_dredge_binary.Rsave")
corhmm_tbl <- read.csv("corhmm_tbl_dredge.csv")
cid_disc_model <- dredge_sociality[[which.min(corhmm_tbl$AIC)]]$index.mat

# now setting up OU part
# Pure rate heterogeneity model for continuous and discrete
cid_oum_model <- nest_oum_model <- soc_oum_model<- full_oum_model<- getOUParamStructure("OUM", 4, 2, null.model = TRUE) # character independent model (i.e. all observed states have the same optima)
#               1                2                3                4 
#   [1]"social|aboveground"   "solitary|aboveground" "social|ground"        "solitary|ground"    

soc_oum_model[3,c(1,3,5,7)] <- 3
soc_oum_model[3,c(2,4,6,8)] <- 4
nest_oum_model[3,c(1,2,5,6)] <- 3
nest_oum_model[3,c(3,4,7,8)] <- 4
full_oum_model[3,c(1:8)] <- c(3:6)

# All models that were run
model_list <- list(list(2, cid_disc_model, "BM1"),
                   list(2, cid_disc_model, "OU1"),
                   list(2, cid_disc_model, soc_oum_model),
                   list(2, cid_disc_model, nest_oum_model),
                   list(2, cid_disc_model, full_oum_model),
                   list(2, cid_disc_model, cid_oum_model)) # the two is for two rate classes

names(model_list) <- c("bm1_8states_run2", "ou1_8states_run2", "oum_soc_8states_run2", "oum_nest_8states_run2", "oum_full_8states_run2", "oum_cid_8states_run2")

# model_list <- list(list(1, disc_model, oum_color),
#                   list(1, disc_model, oum_model))
# names(model_list) <- c("oum_col", "oum_full")


# quickFunc <- function(model_list, model_name){
#   res <- hOUwie(phy, dat, 2, model_list[[2]], model_list[[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = FALSE, n_starts = 10, ncores = 10)
#   file.name <- paste0("houwie_results/",model_name, ".Rsave")
#   save(res, file=file.name)
# }
# 
# mclapply(1:6, function(x) quickFunc(model_list[[x]], names(model_list)[x]), mc.cores = 80)

# 
# # all_model_res[[5]] <- hOUwie(phy, dat, model_set[[5]][[1]], model_set[[5]][[2]], model_set[[5]][[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 5, ncores = 5)
# # 
# oum_bm1_res <- hOUwie(phy, dat, 2, cid_disc_model, "BM1", FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
# save(oum_bm1_res, file="houwie_results/BM1.Rsave")
# oum_ou1_res <- hOUwie(phy, dat, 2, cid_disc_model, "OU1", FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
# save(oum_ou1_res, file="houwie_results/OU1.Rsave")
oum_soc_res <- hOUwie(phy, dat, 2, cid_disc_model, soc_oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
save(oum_soc_res, file="houwie_results/oum_soc.Rsave")
# oum_nest_res <- hOUwie(phy, dat, 2, cid_disc_model, nest_oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
# save(oum_nest_res, file="houwie_results/oum_nest.Rsave")
# oum_ful_res <- hOUwie(phy, dat, 2, cid_disc_model, full_oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
# save(oum_ful_res, file="houwie_results/oum_full.Rsave")
# oum_cid_res <- hOUwie(phy, dat, 2, cid_disc_model, cid_oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 10, ncores = 10)
# save(oum_cid_res, file="houwie_results/oum_cid.Rsave")
