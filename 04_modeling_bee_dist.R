# rm(list=ls())
library(raster)
library(maptools)
library(parallel)
library(data.table)
data("wrld_simpl")

source("00_utility_functions.R")

## Assemble bee sp rich map from rangers
load("curated_data/thinned_points_res1.Rsave")

#all_cleaned_points <- thinned_points
GetRanges(points=thinned_points, species="species", lat="decimalLatitude", lon="decimalLongitude", threshold=0.5, buffer=25, res=10)
# 
# species_list_per_clade <- matrix(nrow=0, ncol=2)
# to_delete <- c()
# for(i in 1:length(all_cleaned_points)) {
#   all_cleaned_points[[i]]$decimalLongitude <- as.numeric(all_cleaned_points[[i]]$decimalLongitude)
#   one_subset1 <- subset(all_cleaned_points[[i]], all_cleaned_points[[i]]$decimalLongitude > -140)
#   one_subset2 <- subset(one_subset1, one_subset1$decimalLongitude < -30)
#   if(nrow(one_subset2)>0) {
#     # unique(one_subset2$scientificName[grep("Ctenoplectra", one_subset2$scientificName)])
#     # one_record <- subset(one_subset2, one_subset2$scientificName %in% "Ctenoplectra terminalis Smith, 1879")
#     one_species_list <- as.data.frame(unique(one_subset2$scientificName))
#     one_species_list$clade <- names(all_cleaned_points)[i]
#     species_list_per_clade <- rbind(species_list_per_clade, one_species_list)    
#     all_cleaned_points[[i]] <- one_subset2
#   } else {
#     to_delete <- c(to_delete, i)
#   }
# }
# all_cleaned_points[to_delete] <- NULL
# # colnames(species_list_per_clade) <- c("species","clade")
# # write.csv(species_list_per_clade, file="species_list_per_clade_americas.csv", row.names=F)

# for(i in 1:length(all_cleaned_points)) {
#   all_cleaned_points[[i]]$decimalLatitude <- as.numeric(all_cleaned_points[[i]]$decimalLatitude)
#   all_cleaned_points[[i]]$decimalLongitude <- as.numeric(all_cleaned_points[[i]]$decimalLongitude)
# }

#mclapply(all_cleaned_points, function(x) GetRanges(x, species="scientificName", lat="decimalLatitude", 
#                                                   lon="decimalLongitude", threshold=0.5, buffer=25, res=10), mc.cores = 10)


