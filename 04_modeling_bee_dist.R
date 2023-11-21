# rm(list=ls())
library(raster)
library(maptools)
library(parallel)
library(data.table)
data("wrld_simpl")

source("R/00_utility_functions_rangers.R")
#setwd("~/Desktop/bee_tree")

## Assemble bee sp rich map from rangers

all_cleaned_points_files <- list.files("gbif_clean", full.names = T)
all_cleaned_points_files <- all_cleaned_points_files[grep("cleaned", all_cleaned_points_files)]
labels <- gsub(paste0(c("gbif_clean/","_cleaned.csv"), collapse="|"),"", all_cleaned_points_files)
all_cleaned_points <- lapply(all_cleaned_points_files, fread)
names(all_cleaned_points) <- labels

#all_all_points <- do.call(rbind, all_cleaned_points)

# length(unique(all_all_points$scientificName))
# [1] 9748

#all_all_points <- subset(all_all_points, !is.na(all_all_points$decimalLatitude))
#all_all_points <- subset(all_all_points, !is.na(all_all_points$decimalLongitude))

# species_fail <- read.csv("species_fail.csv")

#all_cleaned_points[which(unlist(lapply(all_cleaned_points, nrow))==0)]<-NULL
#all_all_points <- do.call(rbind, all_cleaned_points)
#GetRanges(points=all_all_points, species="scientificName", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10)

species_list_per_clade <- matrix(nrow=0, ncol=2)
to_delete <- c()
for(i in 1:length(all_cleaned_points)) {
  all_cleaned_points[[i]]$decimalLongitude <- as.numeric(all_cleaned_points[[i]]$decimalLongitude)
  one_subset1 <- subset(all_cleaned_points[[i]], all_cleaned_points[[i]]$decimalLongitude > -140)
  one_subset2 <- subset(one_subset1, one_subset1$decimalLongitude < -30)
  if(nrow(one_subset2)>0) {
    # unique(one_subset2$scientificName[grep("Ctenoplectra", one_subset2$scientificName)])
    # one_record <- subset(one_subset2, one_subset2$scientificName %in% "Ctenoplectra terminalis Smith, 1879")
    one_species_list <- as.data.frame(unique(one_subset2$scientificName))
    one_species_list$clade <- names(all_cleaned_points)[i]
    species_list_per_clade <- rbind(species_list_per_clade, one_species_list)    
    all_cleaned_points[[i]] <- one_subset2
  } else {
    to_delete <- c(to_delete, i)
  }
}
all_cleaned_points[to_delete] <- NULL
# colnames(species_list_per_clade) <- c("species","clade")
# write.csv(species_list_per_clade, file="species_list_per_clade_americas.csv", row.names=F)

for(i in 1:length(all_cleaned_points)) {
  all_cleaned_points[[i]]$decimalLatitude <- as.numeric(all_cleaned_points[[i]]$decimalLatitude)
  all_cleaned_points[[i]]$decimalLongitude <- as.numeric(all_cleaned_points[[i]]$decimalLongitude)
}

mclapply(all_cleaned_points, function(x) GetRanges(x, species="scientificName", lat="decimalLatitude", 
                                                   lon="decimalLongitude", threshold=0.5, buffer=25, res=10), mc.cores = 10)


