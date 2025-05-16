# Climatic niche breadth
# rm(list=ls())
# setwd("/Users/tvasc/Desktop/bee_sociality_dist")

library(factoextra)
library(alphashape3d)
library(corrplot)
library(caret)

#--------------------------------------
# First organizing dataset:
# Reloading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)

# 1) Comparisons between life history and niche breadth
# Let's compare the extremes by taking the min max temperatures of the coldest and warmest months, 
# and the min max precipitations of the wettest and driest months
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_5_","bio_6_","bio_13_",
                                                    "bio_14_"),collapse="|"), all_climatic_vars)]
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
merged_traits$combined_trait <- paste(merged_traits$sociality_binary, merged_traits$nest_binary, sep="_")

# Removing NAs and NANs
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_5))
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_5))

#strategies <- unique(merged_traits$combined_trait)
head(merged_traits)
merged_traits$prec_breadth <- NA
merged_traits$temp_breadth <- NA
for(i in 1:nrow(merged_traits)) {
  
  merged_traits$prec_breadth[i] <- mean(merged_traits$mean_bio_13[i]) - mean(merged_traits$mean_bio_14[i])
  merged_traits$temp_breadth[i] <- mean(merged_traits$mean_bio_5[i]) - mean(merged_traits$mean_bio_6[i])
  
}  

#boxplot(merged_traits$prec_breadth~merged_traits$combined_trait)
#boxplot(merged_traits$temp_breadth~merged_traits$combined_trait)


#------------------------------------
#------------------------------------
# 2) Position of each species in the climatic multivariate space, 
# plus their volume and correlation between volume and life history. 
# 3) Then volume total for each grouping, and overlap between clouds of points (PERMANOVA)
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)
climatic_list <- lapply(all_climatic_vars, read.csv)

# Now merge everything in one table
merged_climatic_vars <- climatic_list[[1]] 
for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}

# Select only mean columns
#merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# And finally merge to the trait data
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")


# Assuming 'df' contains 19 numeric variables and a column named 'Group'
# Run PCA (scale. = TRUE standardizes the variables)

merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_1))
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_12))

cols <- grep("mean", colnames(merged_traits))
cols <- cols[2:(length(cols)-2)]

any(is.na(merged_traits[, cols]))
any(!is.finite(as.matrix(merged_traits[, cols])))

merged_traits <- merged_traits[complete.cases(merged_traits[, cols]), ]
merged_traits$combined_trait <- paste(merged_traits$sociality_binary, merged_traits$nest_binary, sep="_")

clim_vars <- merged_traits[cols]
clim_vars$tip <- merged_traits$tip
clim_vars$combined_trait <- merged_traits$combined_trait

#merged_traits <- merged_traits[,c("tips","combined_trait",cols)]
#merged_traits[,c("tips","combined_trait",cols)]
#merged_traits$combined_trait

# finite_rows <- apply(merged_traits, 1, function(row) all(is.finite(row)))
# clean_data <- merged_traits[finite_rows, ]
cols <- grep("mean", colnames(clim_vars))
cor_matrix <- cor(clim_vars[,cols])

corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)
cleaned_clim_vars <- clim_vars[,-highly_correlated]

pca_result <- prcomp(cleaned_clim_vars[,1:12], scale. = TRUE)

# Create a data frame with PCA scores
pca_df <- as.data.frame(pca_result$x)

# Add the grouping variable
pca_df$combined_trait <- cleaned_clim_vars$combined_trait

# # Plot the PCA (first two principal components)
# ggplot(pca_df, aes(x = PC1, y = PC2, color = combined_trait)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(title = "PCA of 19 Variables", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1")


# Assign a color to each group
group_colors <- as.factor(pca_df$combined_trait)
palette_colors <- rainbow(length(levels(group_colors)))

palette_colors=c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")
# Open 3D plot
plot3d(
  x = pca_df$PC1,
  y = pca_df$PC2,
  z = pca_df$PC3,
  col = palette_colors[group_colors],
  size = 0.5,
  type = "s",
  xlab = "PC1", ylab = "PC2", zlab = "PC3"
)

#Get contributions of the variables for the 3 axes
var <- get_pca_var(pca_result)
var$coord[,1:3]

#Variance explained by 3 axes
eig <- (pca_result$sdev)^2
variance <- eig*100/sum(eig)
sum(variance[1:3]) # 78.45% explained by first 3 axes

#Map the 3 first axes of the PCA
#PCs <- predict(predictors, pca_result, index = 1:3)
#Convert to Points
#PCs.points <- rasterToPoints(PCs)

strategies <- unique(pca_df$combined_trait)
volumes <- c()
for (i in 1:length(strategies)) {
  ashape3d.occ <- ashape3d(as.matrix(pca_df[pca_df$combined_trait==strategies[i],1:3]), alpha = 2)
  volumes <- c(volumes, volume_ashape3d(ashape3d.occ)) # niche size  
}
names(volumes) <- strategies

volumes
# solitary_ground solitary_aboveground   social_aboveground        social_ground 
# 242.35591            254.03667             68.07669            249.97488 
