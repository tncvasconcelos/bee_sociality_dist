# rm(list=ls())
setwd("~/Desktop/bee_angio_mismatch/")

library(raster)
library(maptools)
library(parallel)
data("wrld_simpl")
source("00_utility_functions.R")

## Assemble bee sp rich map from rangers

all_shape_files <- list.files("shapefile_bees/")
labels <- gsub(".Rdata","",all_shape_files)
list_of_ranges <- lapply(paste0("shapefile_bees/",all_shape_files), readRDS)
names(list_of_ranges) <- labels

list_of_shapes <- list()
for(i in 1:length(list_of_ranges)){
  if(length(list_of_ranges[[i]])>1) {
    list_of_shapes[[i]] <- list_of_ranges[[i]]$range
    names(list_of_shapes)[i] <- names(list_of_ranges)[i]    
  }
}
species_fail <- names(list_of_ranges)[which(unlist(lapply(list_of_shapes, is.null)))]
write.csv(species_fail, file="species_fail.csv", row.names=F)
list_of_shapes[which(unlist(lapply(list_of_shapes, is.null)))] <- NULL

template.map <- readRDS("R/rangers/template.map.Rdata")
res(template.map)[] <- res(list_of_shapes[[1]])

#tmp.raster.list <- list()

if(length(list_of_shapes)>1){
  r0 <- list_of_shapes[[1]]
  r0 <- raster::resample(r0, template.map)
  cat(1, "\r")
  sum_raster <- r0
  for (j in 2:length(list_of_shapes)) {
    r1 <- list_of_shapes[[j]]
    r1 <- raster::resample(r1, template.map)
    sum_raster <- calc(stack(list(sum_raster, r1)), sum, na.rm=T)
    cat(j, "\r")
  }
  mean_raster <- sum_raster 
}

americas_raster <- crop(mean_raster, c(-140,-30,-60,60))
saveRDS(americas_raster, file="preliminary_sprich.Rdata")

pdf("ALL_BEES.pdf")
plot(americas_raster)
plot(wrld_simpl[wrld_simpl$REGION==19,], add=T)
dev.off()

####################################
## Per family
all_shape_files <- list.files("shapefile_bees/")
labels <- gsub(".Rdata","",all_shape_files)
list_of_ranges <- lapply(paste0("shapefile_bees/",all_shape_files), readRDS)
names(list_of_ranges) <- labels

list_of_shapes <- list()
for(i in 1:length(list_of_ranges)){
  if(length(list_of_ranges[[i]])>1) {
    list_of_shapes[[i]] <- list_of_ranges[[i]]$range
    names(list_of_shapes)[i] <- names(list_of_ranges)[i]    
  }
}
species_fail <- names(list_of_ranges)[which(unlist(lapply(list_of_shapes, is.null)))]
write.csv(species_fail, file="species_fail.csv", row.names=F)
list_of_shapes[which(unlist(lapply(list_of_shapes, is.null)))] <- NULL

template.map <- readRDS("R/rangers/template.map.Rdata")
res(template.map)[] <- res(list_of_shapes[[1]])


# all_cleaned_points_files <- list.files("gbif_clean", full.names = T)
# all_cleaned_points_files <- all_cleaned_points_files[grep("cleaned", all_cleaned_points_files)]
# labels <- gsub(paste0(c("gbif_clean/","_cleaned.csv"), collapse="|"),"", all_cleaned_points_files)
# all_cleaned_points <- lapply(all_cleaned_points_files, read.csv)
# names(all_cleaned_points) <- labels
# 
# species_list_per_clade <- matrix(nrow=0, ncol=2)
# for(i in 1:length(all_cleaned_points)) {
#   one_species_list <- as.data.frame(unique(all_cleaned_points[[i]]$scientificName))
#   one_species_list$clade <- names(all_cleaned_points)[i]
#   species_list_per_clade <- rbind(species_list_per_clade, one_species_list)
# }
# colnames(species_list_per_clade) <- c("species","clade")
# write.csv(species_list_per_clade, file="species_list_per_clade.csv", row.names=F)
#tmp.raster.list <- list()

species_list_per_clade <- read.csv("species_list_per_clade.csv")
for(label_index in 1:length(labels)) {
  species_in_one_clade <- subset(species_list_per_clade, species_list_per_clade$clade==labels[label_index])
  shapes_species_in_one_clade <- list_of_shapes[which(names(list_of_shapes)%in%species_in_one_clade$species)]
  
  if(length(shapes_species_in_one_clade)>1){
    r0 <- shapes_species_in_one_clade[[1]]
    r0 <- raster::resample(r0, template.map)
    cat(1, "\r")
    sum_raster <- r0
    for (j in 2:length(shapes_species_in_one_clade)) {
      r1 <- shapes_species_in_one_clade[[j]]
      r1 <- raster::resample(r1, template.map)
      sum_raster <- calc(stack(list(sum_raster, r1)), sum, na.rm=T)
      cat(j, "\r")
    }
    mean_raster <- sum_raster 
  }
  
  americas_raster <- crop(mean_raster, c(-140,-30,-60,60))
  saveRDS(americas_raster, file=paste0(labels[label_index],"_preliminary_sprich.Rdata"))
  
  pdf(paste0(labels[label_index],"_sprich.pdf"))
  plot(americas_raster)
  plot(wrld_simpl[wrld_simpl$REGION==19,], add=T)
  dev.off()
}



#------------------------------
# Shapefiles bees
all_shape_files_bees <- list.files("shapefile_bees/")
labels <- gsub(".Rdata","",all_shape_files_bees)
list_of_ranges_bees <- lapply(paste0("shapefile_bees/",all_shape_files_bees), readRDS)
names(list_of_ranges_bees) <- labels

#------------------------------
# Shapefiles angio
all_shape_files_angio <- list.files("shapefile_angio/")
labels <- gsub(".Rdata","",all_shape_files_angio)
list_of_ranges_angio <- lapply(paste0("shapefile_angio/",all_shape_files_angio), readRDS)
names(list_of_ranges_angio) <- labels

#------------------------------
ranges_bees <- FilterRanges(list_of_ranges_bees, label_fail="bee_species_range_fail")
ranges_angios <- FilterRanges(list_of_ranges_angio, label_fail="angio_species_range_fail")

#------------------------------
bee_richness <- GetSpRichness(ranges_bees)
saveRDS(bee_richness, "sp_rich_data/FullBee_preliminary_sprich.Rdata")
angio_richness <- GetSpRichness(ranges_angios)

