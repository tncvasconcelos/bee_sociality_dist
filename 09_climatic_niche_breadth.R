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
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(viridis)
library(colorspace)
library(gridExtra)
library(factoextra)
library(vegan)
library(caret)
library(corrplot)
library(phytools)
library(alphashape3d)
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

# Merge with trait data again
merged_traits_pca <- merge(traits, merged_climatic_vars, by.x="tips",by.y="species")

# Remove rows with missing values
merged_traits_pca <- subset(merged_traits_pca, !is.na(merged_traits_pca$mean_bio_1))
merged_traits_pca <- subset(merged_traits_pca, !is.na(merged_traits_pca$mean_bio_12))

# Extract only columns with climate mean values
cols <- grep("mean", colnames(merged_traits_pca))
cols <- cols[2:(length(cols)-2)]

# Check for NA or non-finite values in climate columns
any(is.na(merged_traits_pca[, cols]))
any(!is.finite(as.matrix(merged_traits_pca[, cols])))

# Keep only complete cases and re-define combined trait grouping
merged_traits_pca <- merged_traits_pca[complete.cases(merged_traits_pca[, cols]), ]
merged_traits_pca$combined_trait <- paste(merged_traits_pca$sociality_binary, merged_traits_pca$nest_binary, sep="_")

# Extract cleaned climate matrix + traits
clim_vars <- merged_traits_pca[cols]
clim_vars$tip <- merged_traits_pca$tip
clim_vars$combined_trait <- merged_traits_pca$combined_trait
cols <- grep("mean", colnames(clim_vars))

# Remove highly collinear variables (r > 0.9)
cor_matrix <- cor(clim_vars[,cols])
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.8)
highly_correlated <- findCorrelation(cor_matrix, cutoff = 0.9)
cleaned_clim_vars <- clim_vars[,-highly_correlated]

# Principal component analysis (PCA)
# Assuming 'df' contains 19 numeric variables and a column named 'Group'
# Standardizes variables (mean = 0, SD = 1) using scale. = TRUE
pca_result <- prcomp(cleaned_clim_vars[,1:12], scale. = TRUE)

# Create a data frame with PCA scores for each species
pca_df <- as.data.frame(pca_result$x)

# Add the grouping variable
pca_df$combined_trait <- cleaned_clim_vars$combined_trait 


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

# Use the first 3 PCs (as above)
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


# ==============================================================================
# 5) Plotting
# ==============================================================================
# Plot univariate niches as violin plots
#-------------------------------------------------------------------------------
# Set combined_trait factor level order
merged_traits$combined_trait <- factor(merged_traits$combined_trait,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))

# Color setup
violin_colors <- viridis::viridis(n = 4, option = "cividis")
names(violin_colors) <- levels(merged_traits$combined_trait)
point_colors <- sapply(violin_colors, function(x) darken(x, amount = 0.4))

# Labels for y-axis
trait_labels <- c("solitary_ground" = "Solitary/Ground", 
                  "solitary_aboveground" = "Solitary/Above-ground",
                  "social_ground" = "Social/Ground",
                  "social_aboveground" = "Social/Above-ground")

# Precipitation breadth violin plot
p1 <- ggplot(merged_traits, aes(x = prec_breadth, y = combined_trait, fill = combined_trait)) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(aes(color = combined_trait), height = 0.2, alpha = 0.7, size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  labs(x = "Precipitation niche breadth (BIO13 - BIO14)", y = "", title = "Precipitation Breadth") +
  scale_y_discrete(labels = trait_labels) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# Temperature breadth violin plot
p2 <- ggplot(merged_traits, aes(x = temp_breadth, y = combined_trait, fill = combined_trait)) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(aes(color = combined_trait), height = 0.2, alpha = 0.7, size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  labs(x = "Temperature niche breadth (BIO5 - BIO6)", y = "", title = "Temperature Breadth") +
  scale_y_discrete(labels = trait_labels) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

p1
p2

# Save plots
ggsave("plots/precipitation_breadth_violin.pdf", p1, width = 8, height = 6)
ggsave("plots/temperature_breadth_violin.pdf", p2, width = 8, height = 6)

# Combine into a two-panel figure (vertical layout)
combined_univariate_niches <- grid.arrange(p1, p2, ncol = 1)

# Save as PDF and PNG
ggsave("plots/niche_breadth_violin_combined.pdf", combined_univariate_niches, width = 10, height = 12)
ggsave("plots/niche_breadth_violin_combined.png", combined_univariate_niches, width = 10, height = 12)


#-------------------------------------------------------------------------------
# View PCA loadings and contributions before plotting
#-------------------------------------------------------------------------------
# Cumulative variance explained by first 2 PCs (62.44% of total variation)
eig <- (pca_result$sdev)^2
variance <- eig * 100 / sum(eig)
sum(variance[1:2])  # How much variation PC1 + PC2 account for - 62.44%

# Get PCA variable information
pca_vars <- get_pca_var(pca_result)

# Loadings: directions of variables in PC space
print(round(pca_vars$coord[, 1:2], 2))  # Loadings for PC1–PC3

#               Dim.1  Dim.2  Dim.3
# mean_alt.y  -0.347  0.166  0.728
# mean_bio_1   0.721  0.612 -0.228
# mean_bio_13  0.869 -0.084  0.187
# mean_bio_14  0.605 -0.667 -0.099
# mean_bio_15 -0.020  0.772  0.357
# mean_bio_18  0.738 -0.280 -0.023
# mean_bio_19  0.610 -0.344  0.248
# mean_bio_2  -0.656  0.546  0.004
# mean_bio_3   0.589  0.555  0.414
# mean_bio_4  -0.779 -0.320 -0.358
# mean_bio_5   0.079  0.695 -0.591
# mean_bio_8   0.663  0.379 -0.432

# Contributions: how much each variable contributes (%) to each PC
print(round(pca_vars$contrib[, 1:3], 2))  # Contributions for PC1–PC3

#             Dim.1 Dim.2 Dim.3
# mean_alt.y   2.67  0.93 32.08
# mean_bio_1  11.53 12.57  3.16
# mean_bio_13 16.75  0.24  2.12
# mean_bio_14  8.13 14.93  0.59
# mean_bio_15  0.01 19.96  7.73
# mean_bio_18 12.08  2.63  0.03
# mean_bio_19  8.25  3.97  3.74
# mean_bio_2   9.55  9.98  0.00
# mean_bio_3   7.70 10.34 10.36
# mean_bio_4  13.45  3.44  7.75
# mean_bio_5   0.14 16.19 21.17
# mean_bio_8   9.75  4.82 11.28

#-------------------------------------------------------------------------------
# Plot 2D PCA
#-------------------------------------------------------------------------------
# Set factor order and colors
pca_df$combined_trait <- factor(pca_df$combined_trait, levels = c("solitary_ground", "solitary_aboveground", 
                                                                  "social_ground", "social_aboveground"))
trait_colors <- viridis::viridis(n = 4, option = "cividis")
names(trait_colors) <- levels(pca_df$combined_trait)

# Plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = combined_trait, fill = combined_trait)) +
  stat_ellipse(aes(group = combined_trait), type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(group = combined_trait), type = "norm", size = 0.5, linetype = "solid") +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = trait_colors, labels = trait_labels) +
  scale_fill_manual(values = trait_colors, guide = "none") +
  theme_classic() +
  labs(
    title = "",
    color = "Trait Combination",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(size = 0),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = c(1.00, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  )

pca_plot

# Save PCA plot
ggsave("plots/pca_climatic_niche_space.pdf", pca_plot, width = 8, height = 8)
ggsave("plots/pca_climatic_niche_space.png", pca_plot, width = 8, height = 8)


#-------------------------------------------------------------------------------
# Plot PCA loading arrows (color by contribution, no arrow scaling)
#-------------------------------------------------------------------------------
# Extract PCA variable loadings and contributions
pca_vars <- get_pca_var(pca_result)
contributions <- pca_vars$contrib

loadings_df <- as.data.frame(pca_vars$coord)
loadings_df$variable <- rownames(loadings_df)

# Replace variable codes with descriptive labels
loadings_df$label <- recode(loadings_df$variable,
                            mean_bio_1 = "Annual Temp",
                            mean_bio_2 = "Diurnal Temp Range",
                            mean_bio_3 = "Isothermality",
                            mean_bio_4 = "Temp Seasonality",
                            mean_bio_5 = "Max Temp (Warmest Month)",
                            mean_bio_8 = "Mean Temp (Wettest Qtr)",
                            mean_bio_13 = "Precip (Wettest Month)",
                            mean_bio_14 = "Precip (Driest Month)",
                            mean_bio_15 = "Precip Seasonality",
                            mean_bio_18 = "Precip (Warmest Qtr)")

# Filter to top contributors (≥ 9% to PC1 or PC2)
top_vars <- contributions[, 1:2] >= 9
top_var_names <- rownames(contributions)[apply(top_vars, 1, any)]
loadings_df <- loadings_df[loadings_df$variable %in% top_var_names, ]

# Add contribution value for coloring
loadings_df$contrib <- apply(contributions[loadings_df$variable, 1:2], 1, max)

# Plot with unscaled arrow length and color-coded contribution
arrow_plot <- ggplot(loadings_df, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = Dim.1, yend = Dim.2, color = contrib),
               arrow = arrow(length = unit(0.25, "cm")),
               size = 1) +
  geom_text_repel(aes(x = Dim.1, y = Dim.2, label = label),
                  size = 4, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_viridis_c(name = "Max Contribution\n(PC1 or PC2)") +
  coord_equal() +
  theme_classic() +
  labs(
    title = "PCA Loadings",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = c(0.97, 0.03),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA)
  )

arrow_plot

# Arrow direction = the direction of the variable's influence on PC1 and PC2 (into which quadrant)

# Arrow length = the magnitude of the loading (i.e., how strongly the variable is associated with those axes)

# Arrow color = the contribution, i.e., how much that variable contributes to explaining variance along PC1 or PC2

# Save arrow plot
ggsave("plots/pca_loadings.pdf", arrow_plot, width = 8, height = 8)
ggsave("plots/pca_loadings.png", arrow_plot, width = 8, height = 8)


#-------------------------------------------------------------------------------
# Plot niche volume bar plots
#-------------------------------------------------------------------------------
# # Prepare volume data for plotting
# volume_df <- data.frame(
#   combined_trait = names(volumes),
#   volume = as.numeric(volumes)
# )
# 
# # Ensure consistent factor levels and colors
# volume_df$combined_trait <- factor(volume_df$combined_trait, levels = levels(pca_df$combined_trait))
# 
# volume_plot <- ggplot(volume_df, aes(x = combined_trait, y = volume, fill = combined_trait)) +
#   geom_bar(stat = "identity", width = 0.7, alpha = 0.8, show.legend = FALSE) +
#   scale_fill_manual(values = violin_colors) +
#   labs(title = "B. Climatic Niche Volume by Trait Group", x = "", y = "Niche Volume") +
#   scale_x_discrete(labels = c("solitary_ground" = "Solitary/Ground", 
#                               "solitary_aboveground" = "Solitary/Above-ground",
#                               "social_ground" = "Social/Ground",
#                               "social_aboveground" = "Social/Above-ground")) +
#   theme_classic() +
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title = element_text(size = 14),
#     axis.text.x = element_text(size = 12, angle = 20, hjust = 1),
#     axis.text.y = element_text(size = 12)
#   )

# #-------------------------------------------------------------------------------
# # Plot 3D PCA (PC1 - PC3)
# #-------------------------------------------------------------------------------
# # Assign a color to each group
# group_colors <- as.factor(pca_df$combined_trait)
# palette_colors <- rainbow(length(levels(group_colors)))
# palette_colors = c("#851170FF", "#040404FF", "#F36E35FF", "#FFFE9EFF")
# 
# # 3D scatterplot of species in climate space
# plot3d(
#   x = pca_df$PC1,
#   y = pca_df$PC2,
#   z = pca_df$PC3,
#   col = palette_colors[group_colors],
#   size = 0.5,
#   type = "s",
#   xlab = "PC1", ylab = "PC2", zlab = "PC3"
# )
# 
# # PCA variable loadings (i.e., contributions of bioclimatic variables to 3 PCs)
# var <- get_pca_var(pca_result)
# var$coord[,1:3]
# 
# # Cumulative variance explained by first 3 PCs (78.45% of total variation)
# eig <- (pca_result$sdev)^2
# variance <- eig*100/sum(eig)
# sum(variance[1:3])
# 
# #Map the 3 first axes of the PCA
# #PCs <- predict(predictors, pca_result, index = 1:3)
# #Convert to Points
# #PCs.points <- rasterToPoints(PCs)