# rm(list=ls())
setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(raster)
library(maptools)
data("wrld_simpl")

#------------------------
load("curated_data/thinned_points.Rsave")

thinned_points
species <- unique(thinned_points$species)
species <- subset(species, species!="")

# Plotting to inspect distributions
pdf("plots/bee_maps.pdf")
for(spp_index in 1:length(species)){
  tmp_subset <- as.data.frame(thinned_points[thinned_points$species==species[spp_index],])
  coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
  coordinates(coord) <- ~ decimalLongitude + decimalLatitude
  plot(wrld_simpl)
  plot(coord, col="red", add=T)
  title(species[spp_index])
  cat(spp_index, "\r")
}
dev.off()
