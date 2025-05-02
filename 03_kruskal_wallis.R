# ==============================================================================
# 03. Kruskal-Wallis test
# ==============================================================================
# Runs Kruskal-Wallis tests between bee traits and climatic variables,
# and generates violin plots.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directory, load libraries
# ------------------------------------------------------------------------------

#rm(list=ls())
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
library(phytools)
library(ggplot2)
library(gridExtra)


# ------------------------------------------------------------------------------
# Load trait data, phylogeny, and list of climatic variable summary files
# ------------------------------------------------------------------------------

traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


# ------------------------------------------------------------------------------
# Kruskal-Wallis tests and violin plots for nesting_binary vs. ALL climate vars
# ------------------------------------------------------------------------------

sink("results/kruskal_nesting_results_all.txt") # Open a sink for results
plot_list <- list() # Initialize an empty list to store the ggplot objects

# Loop through the climatic variables
for(climate_index in 1:length(all_climatic_vars)) {
  
  # Load and clean climate data
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))  # Remove rows with NA in the third column
 
  # Identify species shared across datasets
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  # Subset datasets to overlapping species
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Merge datasets by species
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  # Prepare data for Kruskal-Wallis
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  # Run Kruskal-Wallis test
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
    significance <- "" # No significance symbol if p >= 0.05
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

sink() # Close the sink


# ------------------------------------------------------------------------------
# Kruskal-Wallis test and violin plots for nesting_binary vs. SELECTED climate variables
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------

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
  
  climate <- read.csv(paste0("curated_data/", selected_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  nests <- merged_table$nest_binary
  names(nests) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  kruskal_result <- kruskal.test(one_clim_var ~ nests)
  p_value <- kruskal_result$p.value
  
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  merged_table$nest_binary <- factor(merged_table$nest_binary, levels = c("ground", "aboveground"))
  
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = nest_binary, y = value, fill = nest_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value",
         title = paste(variable_labels[selected_vars[climate_index]], significance)) +
    scale_x_discrete(labels = c("ground" = "Ground", "aboveground" = "Above-ground")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),
      limits = range(one_clim_var, na.rm = TRUE)
    )
  
  print(plot)
  
  plot_list[[climate_index]] <- plot
  
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("Testing nesting type ~ ", label))
  
  print(kruskal_result)
}

n_vars_selected <- length(selected_vars)
ncols <- 2
nrows <- ceiling(n_vars_selected / ncols)

grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

ggsave("plots/kruskal_nesting_boxplots_selected.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 10, height = 10)

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis tests and violin plots for sociality_binary vs. ALL climate variables
# ------------------------------------------------------------------------------

sink("results/kruskal_sociality_results_all.txt")

plot_list <- list()

for(climate_index in 1:length(all_climatic_vars)) {
  
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[,3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  kruskal_result <- kruskal.test(one_clim_var ~ sociality)
  p_value <- kruskal_result$p.value
  
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }

  merged_table$sociality_binary <- factor(merged_table$sociality_binary, levels = c("solitary", "social"))
  
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = sociality_binary, y = value, fill = sociality_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +
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
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),
      limits = range(one_clim_var, na.rm = TRUE)
    )
  
  print(plot)
  
  plot_list[[climate_index]] <- plot
  
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("sociality ~ ", label))
  
  print(kruskal_result) 
}

n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)

grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

ggsave("plots/kruskal_sociality_boxplots_all.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

sink()

# ------------------------------------------------------------------------------
# Kruskal-Wallis test and violin plots for nesting_binary vs. SELECTED climate variables
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------

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
  
  climate <- read.csv(paste0("curated_data/", selected_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x="tips", by.y="species")
  
  sociality <- merged_table$sociality_binary
  names(sociality) <- merged_table$tips
  
  one_clim_var <- merged_table[, 9]
  names(one_clim_var) <- merged_table$tips
  
  kruskal_result <- kruskal.test(one_clim_var ~ sociality)
  p_value <- kruskal_result$p.value
  
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  merged_table$sociality_binary <- factor(merged_table$sociality_binary, levels = c("solitary", "social"))
  
  colnames(merged_table)[9] <- "value"
  
  plot <- ggplot(merged_table, aes(x = sociality_binary, y = value, fill = sociality_binary)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(variable_labels[selected_vars[climate_index]], significance)) +
    scale_x_discrete(labels = c("solitary" = "Solitary", "social" = "Social")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(
      breaks = pretty(one_clim_var, n = 6),
      limits = range(one_clim_var, na.rm = TRUE)
    )
  
  print(plot)
  
  plot_list[[climate_index]] <- plot
  
  label <- gsub("_climate_summstats.csv", "", selected_vars[climate_index])
  print(paste0("sociality ~ ", label))
  
  print(kruskal_result)
}

n_vars_selected <- length(selected_vars)
ncols <- 2
nrows <- ceiling(n_vars_selected / ncols)

grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

ggsave("plots/kruskal_sociality_boxplots_selected.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 10, height = 10)

par(mfrow = c(1, 1))

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis test for combined nesting and sociality (4 character states)
# ------------------------------------------------------------------------------
# trait combinations: solitary_ground, solitary_aboveground, 
# social_ground, social_aboveground.
# ------------------------------------------------------------------------------

sink("results/kruskal_combined_nest_sociality_results.txt")

plot_list <- list()

for (climate_index in 1:length(all_climatic_vars)) {
  
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  # Create combined trait factor: e.g., "social_aboveground"
  merged_table <- merge(subset_traits, subset_climate, by.x = "tips", by.y = "species")
  merged_table$comb_nest_soc <- interaction(merged_table$sociality_binary, merged_table$nest_binary, sep = "_")
  
  # Define test group and response variable
  group <- merged_table$comb_nest_soc
  one_clim_var <- merged_table[, 9]
  colnames(merged_table)[9] <- "value" # Rename for plotting
  
  kruskal_result <- kruskal.test(value ~ group, data = merged_table)
  p_value <- kruskal_result$p.value
  
  if (p_value < 0.0001) {
    significance <- "****"
  } else if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  # Order factor levels for consistent x-axis positioning
  merged_table$comb_nest_soc <- factor(merged_table$comb_nest_soc,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))
  
  # Add numeric x-axis position for jittering
  merged_table$x_pos <- as.numeric(merged_table$comb_nest_soc)
  
  # Create violin + jitter plot
  plot <- ggplot(merged_table, aes(x = x_pos, y = value, fill = comb_nest_soc)) +
    geom_violin(aes(x = x_pos), width = 0.4, show.legend = FALSE) +
    geom_jitter(aes(x = x_pos + 0.15), width = 0.15, alpha = 0.6, size = 0.3, color = "black") +
    scale_x_continuous(breaks = 1:4, labels = levels(merged_table$comb_nest_soc)) +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value",
         title = paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    scale_y_continuous(breaks = pretty(one_clim_var, n = 6),
                       limits = range(one_clim_var, na.rm = TRUE))
  
  # Print plot and test results
  print(plot)
  plot_list[[climate_index]] <- plot
  
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("Testing nest*sociality combination ~ ", label))
  print(kruskal_result)
}

# Arrange and save plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)

grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

ggsave("plots/kruskal_combined_nesting_sociality_boxplots_all.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

sink()