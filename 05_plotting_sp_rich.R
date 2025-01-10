# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(raster)
library(maptools)
library(parallel)
data("wrld_simpl")
source("00_utility_functions.R")

## Assemble bee sp rich map from rangers
#------------------------------
# Shapefiles bees
all_shape_files_bees <- list.files("shapefile_bees/")
labels <- gsub(".Rdata","",all_shape_files_bees)
list_of_ranges_bees <- lapply(paste0("shapefile_bees/",all_shape_files_bees), readRDS)
names(list_of_ranges_bees) <- labels
# Remove exotic and suspicious
exotic <- read.csv("original_data/exotic_species_to_remove.csv")[,1]
suspicious <- read.csv("original_data/species_with_suspicious_distribution.csv")[,1]
remove <- c(exotic, suspicious)
list_of_ranges_bees <- subset(list_of_ranges_bees, !names(list_of_ranges_bees)%in%remove)

#------------------------------
# Shapefiles angio
all_shape_files_angio <- list.files("shapefile_angiosperms/")
labels <- gsub(".Rdata","",all_shape_files_angio)
list_of_ranges_angio <- lapply(paste0("shapefile_angiosperms/",all_shape_files_angio), readRDS)
names(list_of_ranges_angio) <- labels

#------------------------------
# filter models that failed
ranges_bees <- FilterRanges(list_of_ranges_bees, label_fail="bee_species_range_fail")
ranges_angios <- FilterRanges(list_of_ranges_angio, label_fail="angio_species_range_fail")

#------------------------------
# Get species richness (it takes a while to run)
bee_richness <- GetSpRichness(ranges_bees)
save(bee_richness, "results/FullBee_preliminary_sprich.Rsave")
angio_richness <- GetSpRichness(ranges_angios)
save(angio_richness, "results/FullAngio_preliminary_sprich.Rsave")

#------------------------------
# plot species richness in the american continent
pdf("plots/ALL_BEES.pdf")
americas_raster <- crop(bee_richness, c(-200,-30,-60,80))
plot(americas_raster)
plot(wrld_simpl[wrld_simpl$REGION==19,], add=T)
dev.off()

#----------------
# Some tests:
# divided per sociality
bee_traits <- read.csv("curated_data/bees_sociality.csv")

socials <- which(names(list_of_ranges_bees)%in%bee_traits$tips[bee_traits$sociality=="social"])
#table(bee_traits$sociality)
social_shapes <- list_of_ranges_bees[socials]
#names(social_shapes)

solitaries <- which(names(list_of_ranges_bees)%in%bee_traits$tips[bee_traits$sociality=="solitary"])
#table(bee_traits$sociality)
solitaries_shapes <- list_of_ranges_bees[solitaries]
solitaries_raster <-  CalcSpRich(solitaries_shapes,template.map) 
save(solitaries_raster, file="solitaries_sprich.Rsave")
#load("solitaries_sprich.Rsave")

social_non_bombus_shapes <- social_shapes[which(!grepl("Bombus",names(social_shapes)))]
social_non_bombus_raster <- CalcSpRich(social_non_bombus_shapes,template.map) 
plot(social_non_bombus_raster)
save(social_non_bombus_raster, file="social_non_bombus_sprich.Rsave")

#-----------------------
# getting some proportion maps
pal <- hcl.colors(30, palette = "RdYlBu", alpha = 1)

plot(solitaries_raster, zlim=c(0,300))
plot(social_non_bombus_raster, zlim=c(0,300))

prop_solitaries <- solitaries_raster/mean_raster
prop_socials <- social_non_bombus_raster/mean_raster

prop_solitaries <- crop(prop_solitaries, c(-200,-30,-60,70))
prop_socials <- crop(prop_socials, c(-200,-30,-60,70))

plot(prop_solitaries, col=rev(pal))
plot(wrld_simpl[wrld_simpl$REGION==19,], add=T)

plot(prop_socials, col=rev(pal))
plot(wrld_simpl[wrld_simpl$REGION==19,], add=T)


####################################
## Per family
# all_shape_files <- list.files("shapefile_bees/")
# labels <- gsub(".Rdata","",all_shape_files)
# list_of_ranges <- lapply(paste0("shapefile_bees/",all_shape_files), readRDS)
# names(list_of_ranges) <- labels
# 
# list_of_shapes <- list()
# for(i in 1:length(list_of_ranges)){
#   if(length(list_of_ranges[[i]])>1) {
#     list_of_shapes[[i]] <- list_of_ranges[[i]]$range
#     names(list_of_shapes)[i] <- names(list_of_ranges)[i]    
#   }
# }
# species_fail <- names(list_of_ranges)[which(unlist(lapply(list_of_shapes, is.null)))]
# write.csv(species_fail, file="species_fail.csv", row.names=F)
# list_of_shapes[which(unlist(lapply(list_of_shapes, is.null)))] <- NULL
# 
# template.map <- readRDS("template.map.Rdata")
# res(template.map)[] <- res(list_of_shapes[[1]])
# 
# species_list_per_clade <- read.csv("original_data/all_tips_bee_tree.csv")
# taxon_rank="Tribe"
# 
# labels <- unique(species_list_per_clade[,taxon_rank])
# 
# for(label_index in 1:length(labels)) {
#   species_in_one_clade <- subset(species_list_per_clade, species_list_per_clade[,taxon_rank]==labels[label_index])
#   
#   shapes_species_in_one_clade <- list_of_shapes[which(names(list_of_shapes)%in%species_in_one_clade$tips)]
#   
#   if(length(shapes_species_in_one_clade)>1){
#     r0 <- shapes_species_in_one_clade[[1]]
#     r0 <- raster::resample(r0, template.map)
#     cat(1, "\r")
#     sum_raster <- r0
#     for (j in 2:length(shapes_species_in_one_clade)) {
#       r1 <- shapes_species_in_one_clade[[j]]
#       r1 <- raster::resample(r1, template.map)
#       sum_raster <- calc(stack(list(sum_raster, r1)), sum, na.rm=T)
#       cat(j, "\r")
#     }
#     mean_raster <- sum_raster 
#   }
#   
#   #americas_raster <- crop(mean_raster, c(-140,-30,-60,60))
#   save(mean_raster, file=paste0(labels[label_index],"_preliminary_sprich.Rsave"))
#   
#   pdf(paste0(labels[label_index],"_sprich.pdf"))
#   plot(mean_raster)
#   plot(wrld_simpl, add=T)
#   dev.off()
# }


