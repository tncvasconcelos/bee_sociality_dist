# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 KRUSKAL-WALLIS TEST

################################################################################

# Setup:

#rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)
library(ggplot2)

# Loading traits, tree and climatic data
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


################################################################################

# Kruskal-Wallis test for nesting:

# Open a sink for results
sink("results/kruskal_nesting_results_all.txt")

# Initialize an empty list to store the ggplot objects
plot_list <- list()

# Loop through the climatic variables
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
  
  # Prepare data for Kruskal-Wallis
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Kruskal-Wallis test
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
  
  # Create violin plot
  colnames(merged_table)[9] <- "value" # Rename the variable to match with ggplot syntax
  
  plot <- ggplot(merged_table, aes(x = nest_binary, y = value, fill = nest_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +  # Jittered points
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance)) +
    scale_x_discrete(labels = c("ground" = "Ground", "aboveground" = "Above-ground")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none", # Remove the legend
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),  # Get 6 nice breaks that are whole numbers
      limits = range(one_clim_var, na.rm = TRUE)  # Set y-axis limits to match the data range
    )
  
  # Print the plot
  print(plot)
  
  # Add the plot to the list
  plot_list[[climate_index]] <- plot
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis test result to file
  print(kruskal_result)
}

# Dynamically calculate number of rows and columns based on the number of plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)  # Calculate the required number of rows to fit all plots

# Arrange the plots into a grid
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save the plot grid to a file
ggsave("plots/kruskal_nesting_boxplots_all.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)  # Change dimensions to un-scrunch plots

# Close the sink
sink()

################################################################################

# Now, create a separate PDF for bio_1, bio_4, bio_12, and bio_15

# Summary of p-values:
# Kruskal-Wallis
# bio_1 p < 2.2e-16
# bio_4 p < 2.2e-16
# bio_12 p < 2.2e-16
# bio_15 p = 0.0001512

# vs. phyloGLM
# bio_1 p = 0.4253
# bio_4 p = 0.2193
# bio_12 p = 0.2016
# bio_15 p = 0.7411

# Open sink for results
sink("results/kruskal_nesting_results_selected.txt")

# Select variables
selected_vars <- c("bio_1_climate_summstats.csv", "bio_4_climate_summstats.csv", 
                   "bio_12_climate_summstats.csv", "bio_15_climate_summstats.csv")

# Rename variables
variable_labels <- c("bio_1_climate_summstats.csv" = "Mean annual temperature",
                     "bio_4_climate_summstats.csv" = "Temperature seasonality",
                     "bio_12_climate_summstats.csv" = "Annual precipitation",
                     "bio_15_climate_summstats.csv" = "Precipitation seasonality")

# Initialize an empty list to store the ggplot objects
plot_list <- list()

# Loop through the selected climate variables
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
  
  # Kruskal-Wallis test
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
  
  # Create violin plot
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = nest_binary, y = value, fill = nest_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +  # Jittered points
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(variable_labels[selected_vars[climate_index]], significance)) +
    scale_x_discrete(labels = c("ground" = "Ground", "aboveground" = "Above-ground")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",  # Remove the legend
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),  # Get 6 nice breaks that are whole numbers
      limits = range(one_clim_var, na.rm = TRUE)  # Set y-axis limits to match the data range
    )
  
  # Print the plot
  print(plot)
  
  # Add the plot to the list
  plot_list[[climate_index]] <- plot
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis result to file
  print(kruskal_result)
}

# Dynamically calculate the number of rows and columns based on the number of plots
n_vars_selected <- length(selected_vars)
ncols <- 2  # Set the number of columns
nrows <- ceiling(n_vars_selected / ncols)  # Calculate the required number of rows

# Arrange the plots into a grid
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save the plot grid to a file
ggsave("plots/kruskal_nesting_boxplots_selected.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 10, height = 10)  # Adjust dimensions to un-scrunch plots

# Close the sink
sink()


################################################################################

# Kruskal-Wallis test for sociality:

# Open a sink for results
sink("results/kruskal_sociality_results_all.txt")

# Initialize an empty list to store the ggplot objects
plot_list <- list()

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
  
  # Kruskal-Wallis test
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
  
  # Create a violin plot
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = sociality_binary, y = value, fill = sociality_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +  # Jittered points
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance)) +
    scale_x_discrete(labels = c("solitary" = "Solitary", "social" = "Social")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",  # Remove the legend
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),  # Get 6 nice breaks that are whole numbers
      limits = range(one_clim_var, na.rm = TRUE)  # Set y-axis limits to match the data range
    )
  
  # Print the plot
  print(plot)
  
  # Add the plot to the list
  plot_list[[climate_index]] <- plot
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("sociality ~ ", label)) # Print description of analysis being performed
  
  # Print Kruskal-Wallis test result to file
  print(kruskal_result) 
}

# Dynamically calculate the number of rows and columns based on the number of plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)  # Calculate the required number of rows to fit all plots

# Arrange the plots into a grid
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save the plot grid to a file
ggsave("plots/kruskal_sociality_boxplots_all.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)  # Adjust dimensions to un-scrunch plots

# Close the sink
sink()

################################################################################

# Now, create a separate PDF for bio_1, bio_4, bio_12, and bio_15 for sociality

# Summary of p-values:
# Kruskal-Wallis
# bio_1 p = 0.05596
# bio_4 p = 5.498e-06
# bio_12 p < 2.2e-16
# bio_15 p = 4.943e-07

# vs. phyloGLM
# bio_2 p = 0.9884
# bio_4 p = 0.2007
# bio_12 p = 0.5928
# bio_15 p = 0.1941

# Open a sink for results
sink("results/kruskal_sociality_results_selected.txt")

# Select variables
selected_vars <- c("bio_1_climate_summstats.csv", "bio_4_climate_summstats.csv", 
                   "bio_12_climate_summstats.csv", "bio_15_climate_summstats.csv")

# Adjust variable names
variable_labels <- c("bio_1_climate_summstats.csv" = "Mean annual temperature",
                     "bio_4_climate_summstats.csv" = "Temperature seasonality",
                     "bio_12_climate_summstats.csv" = "Annual precipitation",
                     "bio_15_climate_summstats.csv" = "Precipitation seasonality")

# Initialize an empty list to store the ggplot objects
plot_list <- list()

# Loop through the selected climate variables
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
  
  # Kruskal-Wallis test
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
  
  # Create violin plot
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = sociality_binary, y = value, fill = sociality_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +  # Jittered points
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(variable_labels[selected_vars[climate_index]], significance)) +
    scale_x_discrete(labels = c("solitary" = "Solitary", "social" = "Social")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none",  # Remove the legend
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),  # Get 6 nice breaks that are whole numbers
      limits = range(one_clim_var, na.rm = TRUE)  # Set y-axis limits to match the data range
    )
  
  # Print the plot
  print(plot)
  
  # Add the plot to the list
  plot_list[[climate_index]] <- plot
  
  # Print the label for the current climate variable
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("sociality ~ ", label))  # Print description of analysis being performed
  
  # Print Kruskal-Wallis test result to file
  print(kruskal_result)
}

# Dynamically calculate the number of rows and columns based on the number of plots
n_vars_selected <- length(selected_vars)
ncols <- 2  # Set the number of columns
nrows <- ceiling(n_vars_selected / ncols)  # Calculate the number of rows

# Arrange the plots into a grid
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save the plot grid to a file
ggsave("plots/kruskal_sociality_boxplots_selected.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 10, height = 10)  # Adjust dimensions to un-scrunch plots

# Reset plotting layout to a single panel
par(mfrow = c(1, 1))

# Close the sink
sink()


