# rm(list=ls())
library(corHMM)
library(OUwie)
library(parallel)
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")

#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)

# Let's take bio1 and bio12 (temperature and precipitation) and seasonalities (bio4 and bio15)
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_1_","bio_4_","bio_12_","bio_15_"),collapse="|"), all_climatic_vars)]
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

#--------------------------------------
# Life history traits vs. temperature (bio1)
trait_nest <- merged_traits[,c("tips","nest","mean_bio_1")]
trait_sociality <- merged_traits[,c("tips","sociality","mean_bio_1")]

# subset for speed
dat <- trait_sociality[sample(1:nrow(trait_sociality),50),] 

# organize data
shared_species <- intersect(dat$tips, phy$tip.label)

all(shared_species %in% dat$tips)
all(shared_species %in% phy$tip.label)

dat <- dat[match(shared_species, dat$tips),]
phy <- keep.tip(phy, shared_species)
dat <- dat[match(phy$tip.label, dat$tips),]

# log continuous data
dat[,3] <- log(dat[,3])

# getting model structure
# setting up both character independent and OUMV models 
cid_oumv_model <- full_oum_model <- full_ouv_model <- full_oumv_model <- getOUParamStructure("OUMV", 3, 2, null.model = TRUE) # character independent model (i.e. all observed states have the same optima)

# manually adjusting full OUMV model matrix (i.e. variable theta and sigma2 between observed characters)
full_oumv_model[2,c(1:6)] <- c(2:4)
full_oumv_model[3,c(1:6)] <- c(5:7)
# manually adjusting OUM model matrix
full_oum_model[3,c(1:6)] <- c(4:6)
# manually adjusting OUV model matrix
full_ouv_model[2,c(1:6)] <- c(2:4)

# All models that were run
model_list <- list(list(2, cid_oumv_model, "BM1"),
                   list(2, cid_oumv_model, "OU1"),
                   list(2, cid_oumv_model, full_ouv_model),
                   list(2, cid_oumv_model, full_oum_model),
                   list(2, cid_oumv_model, full_oumv_model)) # the two is for two rate classes

names(model_list) <- c("bm1_sociality_bio1_run1", "ou1_sociality_bio1_run1", "ouv_sociality_bio1_run1",
                       "oum_sociality_bio1_run1", "oumv_sociality_bio1_run1")

quickFunc <- function(model_list, model_name){
  res <- hOUwie(phy, dat, model_list[[1]], model_list[[2]], model_list[[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = FALSE, n_starts = 10, ncores = 4)
  file.name <- paste0("houwie_results/",model_name, ".Rsave")
  save(res, file=file.name)
}

mclapply(1:6, function(x) quickFunc(model_list[[x]], names(model_list)[x]), mc.cores = 4)



# all_model_res[[5]] <- hOUwie(phy, dat, model_set[[5]][[1]], model_set[[5]][[2]], model_set[[5]][[3]], nSim = 100, diagn_msg = TRUE, adaptive_sampling = TRUE, n_starts = 5, ncores = 5)
# 
# oum_bm1_res <- hOUwie(phy, dat, 1, disc_model, "BM1", FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)
# oum_ou1_res <- hOUwie(phy, dat, 1, disc_model, "OU1", FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)
# oum_col_res <- hOUwie(phy, dat, 1, disc_model, oum_color, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)
# oum_frt_res <- hOUwie(phy, dat, 1, disc_model, oum_fruit, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)
# oum_ful_res <- hOUwie(phy, dat, 1, disc_model, oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)
# oum_cid_res <- hOUwie(phy, dat, 2, cid_disc_model, cid_oum_model, FALSE, 100, diagn_msg = TRUE, adaptive_sampling = TRUE)


all_model_results <- list.files("houwie_results")
all_model_results <- subset(all_model_results, grepl("_8states_run2", all_model_results))
model_names <- gsub("_8states_run2.Rsave","",all_model_results)

all_results <- list()
for(i in 1:length(all_model_results)) {
  load(paste0("houwie_results/",all_model_results[i]))
  if(exists("res")) {
    all_results[[i]] <- res
    names(all_results)[i] <- model_names[i]
  }
  remove(res)
}

model_table <- getModelTable(all_results, type="AICc")
write.csv(model_table, file="model_table_results.csv")

average_pars <- getModelAvgParams(all_results, type="AICc")
write.csv(average_pars, file="average_pars.csv")

# what the tip values are expected to be when accounting for the models of evolution
# that is a way to integrate the other parameter of the model
boxplot(exp(average_pars$expected_mean)~average_pars$tip_state, xlab="tip state", ylab="expected mean (petal length)")





