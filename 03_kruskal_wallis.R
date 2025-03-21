# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 KRUSKAL-WALLIS TEST

################################################################################

# Setup:

#rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")
library(phytools)
library(ggplot2)

# Loading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


################################################################################
# Kruskal-Wallis test (ignoring phylogeny) for nesting:

# Testing nesting type ~ all is significant
# Chi-squared = 3.8552, p-value = 0.04959 

# Open a PDF for saving the plots and sink for results
pdf("plots/kruskal_nesting_boxplots_all.pdf", width = 10, height = 25)
sink("results/kruskal_nesting_result_all.txt")

# Set up the plotting area for a multi-panel plot
n_vars <- length(all_climatic_vars)
ncols <- 3  # Set the number of columns
nrows <- ceiling(n_vars / ncols)  # Calculate the number of rows

# Adjust the size of the plotting area to fit all plots
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1))  # Increase the top margin for space

# Loop through the climatic variables for the nesting data
for(climate_index in 1:length(all_climatic_vars)) {
  # Read the climate data for the current variable
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))  # Remove rows with NA in the third column
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  # Prepare data for Kruskal-Wallis test
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ nests)
  p_value <- kruskal_result$p.value
  
  # Determine significance symbol
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""  # No significance symbol if p >= 0.05
  }
  
  # Change levels so that ground displays as first box in boxplot
  merged_table$nest_binary <- factor(merged_table$nest_binary, levels = c("ground", "aboveground"))

  # Boxplot for the current climate variable with significance level in title
  # boxplot(merged_table[, 9] ~ merged_table$nest_binary,
  #         xlab=" ", ylab="Env. Var",
  #         main=paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance),
  #         names = c("Ground","Above-ground"))

  colnames(merged_table)[9] <- "value"
  # Plot
  ggplot(merged_table, aes(x = nest_binary, y = value, fill = nest_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15)) +  # Half-width boxplot
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.1, color = "black") +  # Jittered points
    scale_fill_brewer(palette = "Set3") +  # Nice color palette
    theme_minimal() +
    labs(x = "Nest Binary", y = "", title = "") +
    theme(legend.position = "none")
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis Test result to file
  print(kruskal_result)
}

# Reset plotting layout to a single panel
par(mfrow = c(1, 1))

# Close the sink and PDF for all variables
sink()
dev.off()

################################################################################

# Now, create a separate PDF for bio_1, bio_4, bio_12, and bio_15
# bio_1 p < 2.2e-16
# bio_4 p < 2.2e-16
# bio_12 p < 2.2e-16
# bio_15 p = 0.0001512

# vs. phyloGLM:
# bio_1 p = 0.4253
# bio_4 p = 0.2193
# bio_12 p = 0.2016
# bio_15 p = 0.7411

# Open a PDF for saving the plots and sink for results for selected variables
pdf("plots/kruskal_nesting_boxplots_selected.pdf", width = 15, height = 15)
sink("results/kruskal_nesting_result_selected.txt")

# Set up the plotting area for a multi-panel plot
selected_vars <- c("bio_1_climate_summstats.csv", "bio_4_climate_summstats.csv", 
                   "bio_12_climate_summstats.csv", "bio_15_climate_summstats.csv")

# Mapping of selected variable names to descriptive labels
variable_labels <- c("bio_1_climate_summstats.csv" = "Mean annual temperature",
                     "bio_4_climate_summstats.csv" = "Temperature seasonality",
                     "bio_12_climate_summstats.csv" = "Annual precipitation",
                     "bio_15_climate_summstats.csv" = "Precipitation seasonality")

n_vars_selected <- length(selected_vars)
ncols <- 2  # Set the number of columns
nrows <- ceiling(n_vars_selected / ncols)  # Calculate the number of rows

# Adjust the size of the plotting area to fit all plots (same aspect ratio as before)
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1))  # Increase the top margin for space

# Loop through the selected climate variables for the nesting data
for(climate_index in 1:length(selected_vars)) {
  # Read the climate data for the current variable
  climate <- read.csv(paste0("curated_data/", selected_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))  # Remove rows with NA in the third column
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  # Prepare data for Kruskal-Wallis test
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ nests)
  p_value <- kruskal_result$p.value
  
  # Determine significance symbol
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""  # No significance symbol if p >= 0.05
  }
  
  # Set levels so box for ground nesters displays first
  merged_table$nest_binary <- factor(merged_table$nest_binary, levels = c("ground", "aboveground"))
  
  # Create a color vector based on the levels of 'nest_binary'
  box_colors <- ifelse(merged_table$nest_binary == "ground", "darkgreen", "lightgreen")
  
  # Boxplot for the current climate variable with significance level in title
  boxplot(merged_table[, 9] ~ merged_table$nest_binary, 
          xlab=" ", ylab=variable_labels[selected_vars[climate_index]], 
          main=paste(variable_labels[selected_vars[climate_index]], significance),
          names = c("Ground", "Above-ground"), col = box_colors)
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis Test result to file
  print(kruskal_result)
}

# Reset plotting layout to a single panel
par(mfrow = c(1, 1))

# Close the sink and PDF
sink()
dev.off()


################################################################################

# Kruskal-Wallis test (ignoring phylogeny) for sociality:

# Testing sociality ~ all is significant
# Chi-squared = 14.431, p-value = 0.0001454

# Open a PDF for saving the plots and sink for results
pdf("plots/kruskal_sociality_boxplots.pdf", width = 10, height = 25)
sink("results/kruskal_sociality_result.txt")

# Set up the plotting area for a multi-panel plot
n_vars <- length(all_climatic_vars)
ncols <- 3  # Set the number of columns
nrows <- ceiling(n_vars / ncols)  # Calculate the number of rows

# Adjust the size of the plotting area to fit all plots
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1))  # Increase the top margin for space

# Loop through the climatic variables
for(climate_index in 1:length(all_climatic_vars)) {
  # Read the climate data for the current variable
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))  # Remove rows with NA in the third column
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  # Prepare data for Kruskal-Wallis test
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ sociality)
  p_value <- kruskal_result$p.value
  
  # Determine significance symbol
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""  # No significance symbol if p >= 0.05
  }

  # Change levels so that solitary displays as first box in boxplot
  merged_table$sociality_binary <- factor(merged_table$sociality_binary, levels = c("solitary", "social"))
  
  # Boxplot for the current climate variable with significance level in title
  # Note: boxplot is generated within loop, which means contents of column 9 change with each new loop through to a new climate variable
  boxplot(merged_table[, 9] ~ merged_table$sociality_binary, 
          xlab=" ", ylab="Env. Var", 
          main=paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance),
          names = c("Solitary", "Social"))
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("sociality ~ ", label)) # Print description of analysis being performed
  
  # Print Kruskal-Wallis test result to file
  print(kruskal_result) 
}

# Reset plotting layout to a single panel
par(mfrow = c(1, 1))

# Close the sink and PDF
sink()
dev.off()

################################################################################

# Now, create a separate PDF for bio_1, bio_4, bio_12, and bio_15 for sociality
# bio_1 p = 0.05596
# bio_4 p = 5.498e-06
# bio_12 p < 2.2e-16
# bio_15 p = 4.943e-07

# vs. phyloGLM:
# bio_2 p = 0.9884
# bio_4 p = 0.2007
# bio_12 p = 0.5928
# bio_15 p = 0.1941

# Open a PDF for saving the plots and sink for results for selected variables
pdf("plots/kruskal_sociality_boxplots_selected.pdf", width = 15, height = 15)
sink("results/kruskal_sociality_result_selected.txt")

# Set up the plotting area for a multi-panel plot
selected_vars <- c("bio_1_climate_summstats.csv", "bio_4_climate_summstats.csv", 
                   "bio_12_climate_summstats.csv", "bio_15_climate_summstats.csv")

# Mapping of selected variable names to descriptive labels
variable_labels <- c("bio_1_climate_summstats.csv" = "Mean annual temperature",
                     "bio_4_climate_summstats.csv" = "Temperature seasonality",
                     "bio_12_climate_summstats.csv" = "Annual precipitation",
                     "bio_15_climate_summstats.csv" = "Precipitation seasonality")

n_vars_selected <- length(selected_vars)
ncols <- 2  # Set the number of columns
nrows <- ceiling(n_vars_selected / ncols)  # Calculate the number of rows

# Adjust the size of the plotting area to fit all plots
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1))  # Increase the top margin for space

# Loop through the selected climate variables for sociality data
for(climate_index in 1:length(selected_vars)) {
  # Read the climate data for the current variable
  climate <- read.csv(paste0("curated_data/", selected_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))  # Remove rows with NA in the third column
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Create subsets of traits, climate data, and the tree that include only the sampled species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create merged dataset
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  # Prepare data for Kruskal-Wallis test
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Kruskal-Wallis Test
  kruskal_result <- kruskal.test(one_clim_var ~ sociality)
  p_value <- kruskal_result$p.value
  
  # Determine significance symbol
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""  # No significance symbol if p >= 0.05
  }
  
  # Set levels so box for solitary displays first
  merged_table$sociality_binary <- factor(merged_table$sociality_binary, levels = c("solitary", "social"))
  
  # Create a color vector based on the levels of 'sociality_binary'
  box_colors <- ifelse(merged_table$sociality_binary == "solitary", "darkgreen", "lightgreen")
  
  # Boxplot for the current climate variable with significance level in title
  boxplot(merged_table[, 9] ~ merged_table$sociality_binary, 
          xlab=" ", ylab=variable_labels[selected_vars[climate_index]], 
          main=paste(variable_labels[selected_vars[climate_index]], significance),
          names = c("Solitary", "Social"), col = box_colors)
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("sociality ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis Test result to file
  print(kruskal_result)
}

# Reset plotting layout to a single panel
par(mfrow = c(1, 1))

# Close the sink and PDF for selected variables
sink()
dev.off()

