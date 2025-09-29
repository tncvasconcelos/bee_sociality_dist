# ==============================================================================
# 03. Kruskal-Wallis test
# ==============================================================================
# Runs Kruskal-Wallis tests between bee traits and climatic variables,
# and generates violin plots.
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup: clear environment, set working directories, load libraries
# ------------------------------------------------------------------------------
rm(list=ls())
wd <- "/Users/lenarh/Desktop/bee_sociality_dist"
setwd(wd)
results_wd <- file.path(wd, "results/Kruskal Wallis")

library(phytools)
library(ggplot2)
library(gridExtra)
library(FSA)
library(multcompView)


# ------------------------------------------------------------------------------
# Load trait data, phylogeny, and list of climatic variable summary files
# ------------------------------------------------------------------------------
traits <- read.csv("curated_data/bees_traits.csv")
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
all_climatic_vars <- list.files("curated_data", "summstats.csv") # Selects summstats.csv files from curated_data


# ------------------------------------------------------------------------------
# Kruskal-Wallis tests and violin plots for nesting_binary vs. ALL climate vars
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_nesting_results_all.txt")) # Open sink for results
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

# Arrange plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)  # Calculate the required number of rows to fit all plots

# Display plots
quartz()
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plot grid
ggsave(file.path(results_wd, "kruskal_nesting_boxplots_all.pdf"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis test and violin plots for nesting_binary vs. SELECTED climate variables
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_nesting_results_selected.txt")) 

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

# Arrange plots
n_vars_selected <- length(selected_vars)
ncols <- 2
nrows <- ceiling(n_vars_selected / ncols)

# Display plots
quartz()
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plot grid
ggsave(file.path(results_wd, "kruskal_nesting_boxplots_selected.pdf"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 10, height = 10)

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis tests and violin plots for sociality_binary vs. ALL climate variables
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_sociality_results_all.txt")) 

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

# Arrange plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)

# Display plots
quartz()
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plot grid
ggsave(file.path(results_wd, "kruskal_sociality_boxplots_all.pdf"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

sink()

# ------------------------------------------------------------------------------
# Kruskal-Wallis test and violin plots for sociality_binary vs. SELECTED climate variables
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_sociality_results_selected.txt")) 

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

# Arrange plots
n_vars_selected <- length(selected_vars)
ncols <- 2
nrows <- ceiling(n_vars_selected / ncols)

# Display plots
quartz()
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plot grid
ggsave(file.path(results_wd, "kruskal_sociality_boxplots_selected.pdf"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

par(mfrow = c(1, 1))

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis test for combined nesting and sociality (4 character states) vs. ALL climate variables
# ------------------------------------------------------------------------------
# trait combinations: solitary_ground, solitary_aboveground, 
# social_ground, social_aboveground.
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_combined_results_all.txt")) 

plot_list <- list()

for (climate_index in 1:length(all_climatic_vars)) {
  
  climate <- read.csv(paste0("curated_data/", all_climatic_vars[climate_index]))
  climate <- subset(climate, !is.na(climate[, 3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x = "tips", by.y = "species")
  
  # Create combined trait group
  merged_table$comb_nest_soc <- interaction(merged_table$sociality_binary, 
                                            merged_table$nest_binary, 
                                            sep = "_")
  
  one_clim_var <- merged_table[, 9]
  colnames(merged_table)[9] <- "value"
  merged_table <- merged_table[!is.na(merged_table$value), ]  # Remove NAs
  
  kruskal_result <- kruskal.test(value ~ comb_nest_soc, data = merged_table)
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
  
  # Set factor level order
  merged_table$comb_nest_soc <- factor(merged_table$comb_nest_soc,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))
  
  # Make plot
  plot <- ggplot(merged_table, aes(x = comb_nest_soc, y = value, fill = comb_nest_soc)) +
    geom_violin(width = 0.3, position = position_nudge(x = -0.15), show.legend = FALSE) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 0.3, color = "black") +
    scale_fill_brewer(palette = "Set3") +
    labs(x = "", y = "Climatic Variable Value", 
         title = paste(gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index]), significance)) +
    scale_x_discrete(labels = c("solitary_ground" = "Solitary/Ground", 
                                "solitary_aboveground" = "Solitary/Above-ground",
                                "social_ground" = "Social/Ground",
                                "social_aboveground" = "Social/Above-ground")) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    scale_y_continuous(breaks = pretty(merged_table$value, n = 6))
  
  print(plot)
  
  plot_list[[climate_index]] <- plot
  
  label <- gsub("_climate_summstats.csv", "", all_climatic_vars[climate_index])
  print(paste0("comb_nest_soc ~ ", label))
  print(kruskal_result)
}

# Arrange plots
n_plots <- length(plot_list)
ncols <- 3
nrows <- ceiling(n_plots / ncols)

# Display plots
quartz(width = 15, height = 35)
grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plot grid
ggsave(file.path(results_wd, "kruskal_combined_boxplots_all.pdf"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 15, height = 35)

par(mfrow = c(1, 1))

sink()


# ------------------------------------------------------------------------------
# Kruskal-Wallis test for combined nesting and sociality (4 character states) vs. SELECTED climate variables
# ------------------------------------------------------------------------------
# trait combinations: solitary_ground, solitary_aboveground, 
# social_ground, social_aboveground.
# ------------------------------------------------------------------------------
# bio_1, bio_4, bio_12, bio_15
# ------------------------------------------------------------------------------
sink(file.path(results_wd, "kruskal_combined_results_selected.txt")) 

library(ggplot2)
library(viridis)
library(dunn.test)
library(colorspace)
library(gridExtra)

# Select variables
selected_vars <- c("bio_1_climate_summstats.csv", "bio_4_climate_summstats.csv", 
                   "bio_12_climate_summstats.csv", "bio_15_climate_summstats.csv")

# Adjust variable names
variable_labels <- c("bio_1_climate_summstats.csv" = "Mean annual temperature",
                     "bio_4_climate_summstats.csv" = "Temperature seasonality",
                     "bio_12_climate_summstats.csv" = "Annual precipitation",
                     "bio_15_climate_summstats.csv" = "Precipitation seasonality")

# Define variable units
variable_units <- c(
  "bio_1_climate_summstats.csv" = "Mean annual temperature (°C × 0.1)",
  "bio_4_climate_summstats.csv" = "Temperature seasonality (SD of monthly means, °C × 0.1)",
  "bio_12_climate_summstats.csv" = "Annual precipitation (mm)",
  "bio_15_climate_summstats.csv" = "Precipitation seasonality (% CV)"
)

# Initialize an empty list to store the ggplot objects
plot_list <- list()

# Compute global min/max for each variable
climate_ranges <- list()

for (file in selected_vars) {
  climate <- read.csv(paste0("curated_data/", file))
  climate_values <- climate[, 3]
  climate_ranges[[file]] <- c(min(climate_values, na.rm = TRUE),
                              max(climate_values, na.rm = TRUE))
}

# Loop through the selected climate variables
for (climate_index in 1:length(selected_vars)) {
  
  climate_file <- selected_vars[climate_index]
  climate <- read.csv(paste0("curated_data/", climate_file))
  climate <- subset(climate, !is.na(climate[, 3]))
  sampled_species <- intersect(intersect(traits$tips, climate$species), tree$tip.label)
  
  subset_traits <- subset(traits, traits$tips %in% sampled_species)
  subset_climate <- subset(climate, climate$species %in% sampled_species)
  subset_tree <- keep.tip(tree, tree$tip.label[tree$tip.label %in% sampled_species])
  
  merged_table <- merge(subset_traits, subset_climate, by.x = "tips", by.y = "species")
  
  # Create combined trait group
  merged_table$comb_nest_soc <- interaction(merged_table$sociality_binary, 
                                            merged_table$nest_binary, 
                                            sep = "_")
  
  one_clim_var <- merged_table[, 9]
  colnames(merged_table)[9] <- "value"
  merged_table <- merged_table[!is.na(merged_table$value), ]  # Remove NAs
  
  kruskal_result <- kruskal.test(value ~ comb_nest_soc, data = merged_table)
  p_value <- kruskal_result$p.value
  
  # Set factor level order
  merged_table$comb_nest_soc <- factor(merged_table$comb_nest_soc,
                                       levels = c("solitary_ground", "solitary_aboveground", 
                                                  "social_ground", "social_aboveground"))
  
  # Run Dunn's post-hoc test if Kruskal is significant
  label <- variable_labels[climate_file]
  print(paste0("comb_nest_soc ~ ", label))
  print(kruskal_result)
  
  if (p_value < 0.05) {
    dunn_result <- dunnTest(value ~ comb_nest_soc, data = merged_table, method = "bonferroni")
    print("Dunn's post-hoc test (Bonferroni corrected):")
    print(dunn_result$res)
  }
  
  # Compute breaks that end on a tick
  x_limits <- climate_ranges[[climate_file]]
  x_min <- x_limits[1]
  x_max <- x_limits[2]
  tick_step <- pretty(x_limits, n = 6)[2] - pretty(x_limits, n = 6)[1]
  tick_min <- floor(x_min / tick_step) * tick_step
  tick_max <- ceiling(x_max / tick_step) * tick_step
  tick_breaks <- seq(tick_min, tick_max, by = tick_step)
  
  # Get colors for violins and points
  violin_colors <- viridis::viridis(n = 4, option = "cividis")
  names(violin_colors) <- c("solitary_ground", "solitary_aboveground", 
                            "social_ground", "social_aboveground")
  point_colors <- sapply(violin_colors, function(x) darken(x, amount = 0.4))
  
  # Make plot
  plot <- ggplot(merged_table, aes(x = value, y = comb_nest_soc, fill = comb_nest_soc)) +
    geom_violin(width = 0.9, adjust = 0.8, alpha = 0.6, show.legend = FALSE) +
    geom_jitter(aes(color = comb_nest_soc), height = 0.2, alpha = 0.7, size = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = violin_colors) +
    scale_color_manual(values = point_colors) +
    labs(
      y = "",
      x = variable_units[climate_file],  # show unit label on x-axis
      title = label
    ) +
    scale_y_discrete(labels = c("solitary_ground" = "Solitary/Ground", 
                                "solitary_aboveground" = "Solitary/Above-ground",
                                "social_ground" = "Social/Ground",
                                "social_aboveground" = "Social/Above-ground")) +
    scale_x_continuous(limits = c(tick_min, tick_max), breaks = tick_breaks) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.y = element_text(size = 24, color = "black"),
      axis.text.x = element_text(size = 20, color = "black"),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.ticks.x = element_line(color = "black", size = 0.5),
      plot.margin = margin(t = 10, r = 40, b = 10, l = 40)
    )
  
  print(plot)
  plot_list[[climate_index]] <- plot
}

# Arrange plots in grid
n_vars_selected <- length(selected_vars)
ncols <- 1
nrows <- ceiling(n_vars_selected / ncols)

grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows)

# Save plots
ggsave("results/Kruskal-Wallis/kruskal_combined_boxplots_selected.pdf", 
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 12, height = 24)

ggsave(file.path(results_wd, "kruskal_combined_boxplots_selected.png"),
       plot = grid.arrange(grobs = plot_list, ncol = ncols, nrow = nrows),
       width = 12, height = 24)

sink()
