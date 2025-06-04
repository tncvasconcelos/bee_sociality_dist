# ==============================================================================
# 09. Climatic niche breadth
# ==============================================================================
# Evaluates how combinations of sociality and nesting strategy influence
# bees’ climatic niche breadth, both univariately (temperature and precipitation) 
# and multivariately (using PCA and 3D volume estimation).
# ==============================================================================


#-------------------------------------------------------------------------------
# Setup: clear workspace, set wd, load packages
#-------------------------------------------------------------------------------
rm(list=ls())
# setwd("/Users/tvasc/Desktop/bee_sociality_dist")
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(factoextra)
library(alphashape3d)
library(corrplot)
library(caret)
library(vegan)
library(dplyr)
library(combinat)

#-------------------------------------------------------------------------------
# Load trait, tree, and climate summary data
#-------------------------------------------------------------------------------
traits <- read.csv("curated_data/bees_traits.csv")
phy <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)


# ==============================================================================
# 1) Univariate niche breadths (temp and precip)
# ==============================================================================
# Let's compare the extremes by taking the min max temperatures of the coldest and warmest months, 
# and the min max precipitations of the wettest and driest months
# ==============================================================================
# Filter for BIOCLIM variables corresponding to:
# - BIO5: Max temperature of warmest month
# - BIO6: Min temperature of coldest month
# - BIO13: Precipitation of wettest month
# - BIO14: Precipitation of driest month
all_climatic_vars <- all_climatic_vars[grep(paste(c("bio_5_","bio_6_","bio_13_",
                                                    "bio_14_"),collapse="|"), all_climatic_vars)]
climatic_list <- lapply(all_climatic_vars, read.csv)

# Merge climate means into one dataframe based on species
merged_climatic_vars <- climatic_list[[1]] 
for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}

# Keep only mean values (averaged across occurrence points)
merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# Merge climate data with trait data using species names
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

# Combine trait states into strategies: e.g., "social_aboveground"
merged_traits$combined_trait <- paste(merged_traits$sociality_binary, merged_traits$nest_binary, sep="_")
#strategies <- unique(merged_traits$combined_trait)

# Removing NAs and NANs
merged_traits <- subset(merged_traits, !is.nan(merged_traits$mean_bio_5))
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_5))

# Check data
head(merged_traits)

# Calculate univariate niche breadths for each species
# temp_breadth = BIO5 - BIO6 (temperature range)
# prec_breadth = BIO13 - BIO14 (precipitation range)
merged_traits$prec_breadth <- NA
merged_traits$temp_breadth <- NA

for(i in 1:nrow(merged_traits)) {
  merged_traits$prec_breadth[i] <- mean(merged_traits$mean_bio_13[i]) - mean(merged_traits$mean_bio_14[i])
  merged_traits$temp_breadth[i] <- mean(merged_traits$mean_bio_5[i]) - mean(merged_traits$mean_bio_6[i])
}  

# Visualize
#boxplot(merged_traits$prec_breadth~merged_traits$combined_trait)
#boxplot(merged_traits$temp_breadth~merged_traits$combined_trait)


# ==============================================================================
# 2) Multivariate climatic niche space and niche volume estimation
# ==============================================================================
# Position of each species in the climatic multivariate space, 
# plus their volume and correlation between volume and life history. 
# ==============================================================================
# Reload all climate summary stats for full PCA analysis
all_climatic_vars <- list.files("curated_data", "summstats.csv", full.names = T)
climatic_list <- lapply(all_climatic_vars, read.csv)

# Merge all climatic variables across species
merged_climatic_vars <- climatic_list[[1]] 

for(i in 2:length(climatic_list)) {
  one_climatic_var <- climatic_list[[i]]
  merged_climatic_vars <- merge(merged_climatic_vars, one_climatic_var, by="species") 
}

# Select only mean columns
#merged_climatic_vars <- merged_climatic_vars[,c(1, grep("mean", colnames(merged_climatic_vars)))]

# Merge with trait data again
merged_traits <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

# Remove rows with missing values
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_1))
merged_traits <- subset(merged_traits, !is.na(merged_traits$mean_bio_12))

# Extract only columns with climate mean values
cols <- grep("mean", colnames(merged_traits))
cols <- cols[2:(length(cols)-2)]

# Check for NA or non-finite values in climate columns
any(is.na(merged_traits[, cols]))
any(!is.finite(as.matrix(merged_traits[, cols])))

# Final filtering: keep only complete cases and re-define combined trait grouping
merged_traits <- merged_traits[complete.cases(merged_traits[, cols]), ]
merged_traits$combined_trait <- paste(merged_traits$sociality_binary, merged_traits$nest_binary, sep="_")

# Extract cleaned climate matrix + traits
clim_vars <- merged_traits[cols]
clim_vars$tip <- merged_traits$tip
clim_vars$combined_trait <- merged_traits$combined_trait

#merged_traits <- merged_traits[,c("tips","combined_trait",cols)]
#merged_traits[,c("tips","combined_trait",cols)]
#merged_traits$combined_trait

# finite_rows <- apply(merged_traits, 1, function(row) all(is.finite(row)))
# clean_data <- merged_traits[finite_rows, ]
cols <- grep("mean", colnames(clim_vars))

# Remove highly collinear variables (r > 0.9)
cor_matrix <- cor(clim_vars[,cols])
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)
cleaned_clim_vars <- clim_vars[,-highly_correlated]

# Principal component analysis (PCA)
# Assuming 'df' contains 19 numeric variables and a column named 'Group'
# Standardizes variables (mean=0, SD=1) using scale. = TRUE
pca_result <- prcomp(cleaned_clim_vars[,1:12], scale. = TRUE)

# Create a data frame with PCA scores for each species
pca_df <- as.data.frame(pca_result$x)

# Add the grouping variable
pca_df$combined_trait <- cleaned_clim_vars$combined_trait 

# # Plot 2D PCA (PC1 - PC2)
# ggplot(pca_df, aes(x = PC1, y = PC2, color = combined_trait)) +
#   geom_point(size = 3, alpha = 0.8) +
#   labs(title = "PCA of 19 Variables", x = "PC1", y = "PC2") +
#   theme_minimal() +
#   scale_color_brewer(palette = "Set1")

# Plot 3D PCA (PC1 - PC3)
# Assign a color to each group
group_colors <- as.factor(pca_df$combined_trait)
palette_colors <- rainbow(length(levels(group_colors)))
palette_colors = c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")

# 3D scatterplot of species in climate space
plot3d(
  x = pca_df$PC1,
  y = pca_df$PC2,
  z = pca_df$PC3,
  col = palette_colors[group_colors],
  size = 0.5,
  type = "s",
  xlab = "PC1", ylab = "PC2", zlab = "PC3"
)

# PCA variable loadings (i.e., contributions of bioclimatic variables to 3 PCs)
var <- get_pca_var(pca_result)
var$coord[,1:3]

# Cumulative variance explained by first 3 PCs (78.45% of total variation)
eig <- (pca_result$sdev)^2
variance <- eig*100/sum(eig)
sum(variance[1:3])

#Map the 3 first axes of the PCA
#PCs <- predict(predictors, pca_result, index = 1:3)
#Convert to Points
#PCs.points <- rasterToPoints(PCs)


# ==============================================================================
# 3) Niche volume calculation per life history syndrome (alpha shape volume)
# ==============================================================================
# Calculate the volume of the climatic space occupied by each life history syndrome.
# ==============================================================================
strategies <- unique(pca_df$combined_trait)

volumes <- c()

for (i in 1:length(strategies)) {
  ashape3d.occ <- ashape3d(as.matrix(pca_df[pca_df$combined_trait==strategies[i],1:3]), alpha = 2)
  volumes <- c(volumes, volume_ashape3d(ashape3d.occ)) # niche size  
}

names(volumes) <- strategies

# Final output: multivariate climatic niche volumes per trait group
volumes

# 6/4/25
# solitary_ground solitary_aboveground   social_aboveground        social_ground 
# 242.35591            254.03667             68.07669          249.97488 


# ==============================================================================
# 4) Statistical test: PERMANOVA on multivariate climatic niche space
# ==============================================================================
# Tests whether species grouped by trait combination occupy significantly 
# different regions of climate space (based on PCA scores).
# ==============================================================================

# Use the first 3 PCs for climate space (as above)
adonis_data <- pca_df[, c("PC1", "PC2", "PC3")]
adonis_group <- pca_df$combined_trait

# Run PERMANOVA
permanova_result <- adonis2(adonis_data ~ adonis_group, method = "euclidean", permutations = 999)

# Print results
print(permanova_result)

# 6/4/25
# R^2 = 0.11499
# This means 11.5% of the variation in species’ positions in climate space is explained by the trait group.
# Residual = 88.5% — what's left unexplained.

# p = 0.001 with 999 permutations
# With 999 permutations, the smallest possible p-value you can report is 0.001, 
# This means none of the 999 random groupings had as much between-group difference as your real data.

# ==============================================================================
# Now let's do a pairwise PERMANOVA to compare between trait groups
# ==============================================================================
# Get unique trait combinations
groups <- unique(pca_df$combined_trait)

# Store results
pairwise_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all unique pairs
for (i in 1:(length(groups)-1)) {
  for (j in (i+1):length(groups)) {
    g1 <- groups[i]
    g2 <- groups[j]
    
    # Subset data for just these two groups
    subset_df <- pca_df %>%
      filter(combined_trait %in% c(g1, g2))
    
    # Get PCA scores and group variable
    data <- subset_df[, c("PC1", "PC2", "PC3")]
    group <- subset_df$combined_trait
    
    # Run PERMANOVA
    result <- adonis2(data ~ group, method = "euclidean", permutations = 999)
    
    # Store in results table
    pairwise_results <- rbind(pairwise_results, data.frame(
      Group1 = g1,
      Group2 = g2,
      R2 = result$R2[1],
      F_value = result$F[1],
      p_value = result$`Pr(>F)`[1]
    ))
  }
}

# Print all pairwise results
print(pairwise_results)

# 6/4/25
# Group1               Group2         R2   F_value p_value
# 1      solitary_ground solitary_aboveground 0.03211794  92.78172   0.001
# 2      solitary_ground   social_aboveground 0.16466161 426.76404   0.001
# 3      solitary_ground        social_ground 0.03106720  80.25448   0.001
# 4 solitary_aboveground   social_aboveground 0.10204339 137.16295   0.001
# 5 solitary_aboveground        social_ground 0.03143050  50.13591   0.001
# 6   social_aboveground        social_ground 0.25537780 313.46810   0.001
