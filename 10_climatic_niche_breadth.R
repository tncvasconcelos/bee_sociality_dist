# ==============================================================================
# 09. Climatic niche breadth
# ==============================================================================
# Evaluates how combinations of sociality and nesting strategy influence
# bees’ climatic niche breadth, both univariately (temperature and precipitation) 
# and multivariately (using PCA and 3D volume estimation).
# ==============================================================================

#-------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load data and packages
#-------------------------------------------------------------------------------
rm(list=ls())
wd <- "/Users/lenarh/Desktop/bee_sociality_dist"
setwd(wd)
results_wd <- file.path(wd, "results/PCA")

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

# Load trait, tree, and climate summary data
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

head(merged_traits$prec_breadth)
head(merged_traits$temp_breadth)
levels(merged_traits$combined_trait)
table(merged_traits$combined_trait)

# Save dataset with trait data, climatic data, and univariate niche breadths
saveRDS(merged_traits, file = file.path(results_wd, "merged_traits.rds"))


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

# Save PCA results and scores
saveRDS(pca_result, file = file.path(results_wd, "pca_result.rds"))
saveRDS(pca_df, file = file.path(results_wd, "pca_df.rds"))


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

# Convert to df and save
volumes_df <- data.frame(
  combined_trait = names(volumes),
  niche_volume = as.numeric(volumes)
)

head(volumes_df)

# Save niche volumes
write.csv(volumes_df, file = file.path(results_wd, "multivariate_niche_volumes.csv"), row.names = FALSE)

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

# Convert to df and save
permanova_df <- as.data.frame(permanova_result)
head(permanova_df)

# Save PERMANOVA
write.csv(permanova_df, file = file.path(results_wd, "permanova_main.csv"), row.names = TRUE)


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

# Save pairwise PERMANOVA
write.csv(pairwise_results, file = file.path(results_wd, "permanova_pairwise.csv"), row.names = FALSE)

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
# Setup required for ALL plotting
#-------------------------------------------------------------------------------
# Reload trait dataset, PCA dataframe and results
merged_traits <- readRDS("results/PCA/merged_traits.rds")
pca_result <- readRDS("results/PCA/pca_result.rds")
pca_df <- readRDS("results/PCA/pca_df.rds")

# Set combined_trait factor level order
merged_traits$combined_trait <- factor(merged_traits$combined_trait,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))

# Trait labels
trait_labels <- c("solitary_ground" = "Solitary/Ground", 
                  "solitary_aboveground" = "Solitary/Above-ground",
                  "social_ground" = "Social/Ground",
                  "social_aboveground" = "Social/Above-ground")

# Color setup for multipanel figure
trait_colors <- c(
  "solitary_ground"      = "#1E3D32",  # forest (dark)
  "solitary_aboveground" = "#7FAF9C",  # soft sage/forest light
  "social_ground"        = "#8B6C00",  # mustard brown (dark)
  "social_aboveground"   = "#F5D547"   # bright mustard yellow
)

sociality_colors <- c(
  "solitary" = "#2B6CA3",  # same as solitary_ground
  "social" = "#FDE725"     # same as social_ground
)

nesting_colors <- c(
  "ground" = "#505050",        # neutral dark gray for contrast
  "aboveground" = "#B0B0B0"    # reused from above-ground social
)


#-------------------------------------------------------------------------------
# Plot 2D PCA
#-------------------------------------------------------------------------------
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = combined_trait, fill = combined_trait)) +
  stat_ellipse(aes(group = combined_trait), type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(group = combined_trait), type = "norm", size = 0.5, linetype = "solid") +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = trait_colors , labels = trait_labels) +
  scale_fill_manual(values = trait_colors , guide = "none") +
  theme_classic() +
  labs(
    title = "",
    color = "Trait Combination",
    x = "PC1: Cool, dry, seasonal (T) ↔ Wet, warm, aseasonal (T)",
    y = ""
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    plot.title = element_text(size = 0),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA) 
  ) + guides(color = guide_legend(override.aes = list(size = 3)))# + theme(legend.position = "none")

pca_plot

# Save PCA plot of 4 trait combinations
ggsave(file.path(results_wd, "2D_pca.pdf"), pca_plot, width = 14, height = 8)
ggsave(file.path(results_wd, "2D_pca.png"), pca_plot, width = 14, height = 8)


# -------------------------------------------------------------------------------
# Plot separate PCAs for sociality and nesting strategy
# -------------------------------------------------------------------------------
# Make sure traits are aligned with PCA data
pca_df$sociality <- merged_traits$sociality_binary
pca_df$nesting <- merged_traits$nest_binary

# Sociality PCA
pca_sociality <- ggplot(pca_df, aes(x = PC1, y = PC2, color = sociality, fill = sociality)) +
  stat_ellipse(aes(group = sociality), type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(group = sociality), type = "norm", size = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = sociality_colors, labels = c("solitary" = "Solitary", "social" = "Social")) +
  scale_fill_manual(values = sociality_colors, guide = "none") +
  labs(
    title = "",
    color = "Sociality",
    x = "",
    y = "PC2: Mild summers (T), aseasonal (P) ↔c Extreme summers (T), seasonal (P)"
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  ) + guides(color = guide_legend(override.aes = list(size = 3)))

pca_sociality

# Nesting PCA
pca_nesting <- ggplot(pca_df, aes(x = PC1, y = PC2, color = nesting, fill = nesting)) +
  stat_ellipse(aes(group = nesting), type = "norm", geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(aes(group = nesting), type = "norm", size = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = nesting_colors, labels = c("ground" = "Ground", "aboveground" = "Above-ground")) +
  scale_fill_manual(values = nesting_colors, guide = "none") +
  labs(
    title = "",
    color = "Nesting Strategy",
    x = "",
    y = ""
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.position = "right",
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = NA)
  ) + guides(color = guide_legend(override.aes = list(size = 3)))

pca_nesting

# Combine and save 3 panel PCA plot
library(patchwork)
combined_pca_plot <- pca_nesting / pca_sociality / pca_plot
combined_pca_plot

ggsave(file.path(results_wd, "combined_pca_plot.pdf"), combined_pca_plot, width = 10, height = 15)
ggsave(file.path(results_wd, "combined_pca_plot.png"), combined_pca_plot, width = 10, height = 15)

#-------------------------------------------------------------------------------
# Plot univariate niches as violin plots
#-------------------------------------------------------------------------------
# Set combined_trait factor level order and rename
merged_traits$combined_trait <- factor(merged_traits$combined_trait,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))

trait_labels <- c("solitary_ground" = "Solitary/Ground", 
                  "solitary_aboveground" = "Solitary/Above-ground",
                  "social_ground" = "Social/Ground",
                  "social_aboveground" = "Social/Above-ground")

# Define colors
violin_colors <- viridis::viridis(n = 4, option = "cividis")
names(violin_colors) <- levels(merged_traits$combined_trait)
point_colors <- sapply(violin_colors, function(x) darken(x, amount = 0.4))

# Precipitation breadth violin plot
precip_breadth <- ggplot(merged_traits, aes(x = prec_breadth, y = combined_trait, fill = combined_trait)) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(aes(color = combined_trait), height = 0.2, alpha = 0.7, size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  labs(x = "BIO13 - BIO14", y = "", title = "Precipitation Niche Breadth") +
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
temp_breadth <- ggplot(merged_traits, aes(x = temp_breadth, y = combined_trait, fill = combined_trait)) +
  geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_jitter(aes(color = combined_trait), height = 0.2, alpha = 0.7, size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = violin_colors) +
  scale_color_manual(values = point_colors) +
  labs(x = "BIO5 - BIO6", y = "", title = "Temperature Niche Breadth") +
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

precip_breadth
temp_breadth

# Combine into a two-panel figure (vertical layout)
combined_univariate_niches <- grid.arrange(precip_breadth, temp_breadth, ncol = 1)

# Save niche breadth violin plots
ggsave(file.path(results_wd, "niche_breadth_violin_combined.pdf"), combined_univariate_niches, width = 10, height = 12)
ggsave(file.path(results_wd, "niche_breadth_violin_combined.png"), combined_univariate_niches, width = 10, height = 12)

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
    title = "",
    x = "PC1",
    y = "PC2"
  ) +
  theme(
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = c(1.2, 0.03),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = NA)
  )

arrow_plot

# Arrow direction = the direction of the variable's influence on PC1 and PC2 (into which quadrant)

# Arrow length = the magnitude of the loading (i.e., how strongly the variable is associated with those axes)

# Arrow color = the contribution, i.e., how much that variable contributes to explaining variance along PC1 or PC2

# Save arrow plot
ggsave(file.path(results_wd, "pca_loadings.pdf"), arrow_plot, width = 8, height = 8)
ggsave(file.path(results_wd, "pca_loadings.png"), arrow_plot, width = 8, height = 8)


#-------------------------------------------------------------------------------
# View PCA loadings and contributions
#-------------------------------------------------------------------------------
# Cumulative variance explained by first 2 PCs (62.44% of total variation)
eig <- (pca_result$sdev)^2
variance <- eig * 100 / sum(eig)
sum(variance[1:2])  # How much variation PC1 + PC2 account for - 62.44%

# Get PCA variable information
pca_vars <- get_pca_var(pca_result)

# Loadings: directions of variables in PC space
print(round(pca_vars$coord[, 1:2], 2))  # Loadings for PC1–PC2

#             Dim.1 Dim.2
# mean_alt.y  -0.35  0.17
# mean_bio_1   0.72  0.61
# mean_bio_13  0.87 -0.08
# mean_bio_14  0.61 -0.67
# mean_bio_15 -0.02  0.77
# mean_bio_18  0.74 -0.28
# mean_bio_19  0.61 -0.34
# mean_bio_2  -0.66  0.55
# mean_bio_3   0.59  0.56
# mean_bio_4  -0.78 -0.32
# mean_bio_5   0.08  0.69
# mean_bio_8   0.66  0.38

# Contributions: how much each variable contributes (%) to each PC
print(round(pca_vars$contrib[, 1:2], 2))  # Contributions for PC1–PC2

#             Dim.1 Dim.2
# mean_alt.y   2.67  0.93
# mean_bio_1  11.53 12.57
# mean_bio_13 16.75  0.24
# mean_bio_14  8.13 14.93
# mean_bio_15  0.01 19.96
# mean_bio_18 12.08  2.63
# mean_bio_19  8.25  3.97
# mean_bio_2   9.55  9.98
# mean_bio_3   7.70 10.34
# mean_bio_4  13.45  3.44
# mean_bio_5   0.14 16.19
# mean_bio_8   9.75  4.82

# Filter to top contributors (≥ 9% to PC1 or PC2)
top_vars <- contributions[, 1:2] >= 9
top_var_names <- rownames(contributions)[apply(top_vars, 1, any)]
top_var_names


# Save to CSV
# Extract loadings and contributions for top variables
top_loadings <- as.data.frame(pca_vars$coord[top_var_names, 1:2])
top_contribs <- as.data.frame(pca_vars$contrib[top_var_names, 1:2])

# Add variable names as a column
top_loadings$Variable <- rownames(top_loadings)
top_contribs$Variable  <- rownames(top_contribs)

# Merge both into one table
top_pca_summary <- merge(top_loadings, top_contribs, by = "Variable")

# Rename columns for clarity
colnames(top_pca_summary) <- c("Variable", "PC1_Loading", "PC2_Loading", "PC1_Contribution", "PC2_Contribution")

# Round values to 2 decimal places
top_pca_summary[, 2:5] <- round(top_pca_summary[, 2:5], 2)

head(top_pca_summary)

# 8/8/2025
#       Variable PC1_Loading PC2_Loading PC1_Contribution PC2_Contribution
# 1  mean_bio_1        0.72        0.61            11.53            12.57
# 2 mean_bio_13        0.87       -0.08            16.75             0.24
# 3 mean_bio_14        0.61       -0.67             8.13            14.93
# 4 mean_bio_15       -0.02        0.77             0.01            19.96
# 5 mean_bio_18        0.74       -0.28            12.08             2.63
# 6  mean_bio_2       -0.66        0.55             9.55             9.98

# Save PCA loadings/contributions
write.csv(top_pca_summary, file.path(results_wd, "pca_top_contributors_PC1_PC2.csv"), row.names = FALSE)





# ==============================================================================
# Code graveyard: optional things not currently used in manuscript
# ==============================================================================
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

#-------------------------------------------------------------------------------
# Plot 3D PCA (PC1 - PC3)
#-------------------------------------------------------------------------------
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

# # PCA variable loadings (i.e., contributions of bioclimatic variables to 3 PCs)
# var <- get_pca_var(pca_result)
# var$coord[,1:3]

# # Cumulative variance explained by first 3 PCs (78.45% of total variation)
# eig <- (pca_result$sdev)^2
# variance <- eig*100/sum(eig)
# sum(variance[1:3])

#Map the 3 first axes of the PCA
#PCs <- predict(predictors, pca_result, index = 1:3)
#Convert to Points
#PCs.points <- rasterToPoints(PCs)

