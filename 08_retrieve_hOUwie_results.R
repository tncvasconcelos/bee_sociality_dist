# ==============================================================================
# 07.3 hOUwie Results Processing
# ==============================================================================
# Loop over the four climatic variables of interest (bio1, bio12, bio15, bio4)
# 
# For each variable, 
#   Load the relevant .Rsave files
#   Extract model comparison tables
#   Compute model-averaged parameters
#   Save CSVs of the above
#   Plot expected means and variances
# ==============================================================================

library(OUwie)
library(corHMM)
library(parallel)

setwd("/Users/lenarh/Desktop/bee_sociality_dist")
source("00_utility_functions.R")  # getModelTable(), getModelAvgParams()

# Define focal BIO variables of interest
bio_vars <- c("bio1", "bio4", "bio12", "bio15")  # Make sure names match actual files

# Create output folder if it doesn't exist
if (!dir.exists("results/hOUwie")) dir.create("results/hOUwie", recursive = TRUE)

# Loop through each climate variable
for (bio in bio_vars) {
  
  message("Processing: ", bio)
  
  # Define expected output file paths
  csv_table_path <- paste0("results/hOUwie/hOUwie_model_table_", bio, ".csv")
  csv_pars_path  <- paste0("results/hOUwie/hOUwie_average_pars_", bio, ".csv")
  plot_mean_path <- paste0("results/hOUwie/hOUwie_boxplot_", bio, "_mean.png")
  plot_var_path  <- paste0("results/hOUwie/hOUwie_boxplot_", bio, "_variance.png")
  
  # # Skip if outputs already exist (comment this out if you want to re-run all)
  # if (all(file.exists(csv_table_path, csv_pars_path, plot_mean_path, plot_var_path))) {
  #   message("Skipping ", bio, " (results already exist)")
  #   next
  # }
  
  # Subset relevant .Rsave files
  all_model_results <- list.files("houwie_results", full.names = TRUE)
  bio_model_files <- all_model_results[grepl(paste0("_", bio, "\\.Rsave$"), all_model_results)]
  model_names <- gsub(".Rsave", "", basename(bio_model_files))
  
  # Load each result into a list
  all_results <- list()
  for (i in seq_along(bio_model_files)) {
    tryCatch({
      load(bio_model_files[i])  # Loads 'res'
      if (exists("res") && !is.null(res)) {
        all_results[[i]] <- res
        names(all_results)[i] <- model_names[i]
      }
      rm(res)
    }, error = function(e) {
      warning("Failed to load or process: ", bio_model_files[i], "\n", e$message)
    })
  }
  
  # Skip if no models loaded
  if (length(all_results) == 0) {
    warning("No valid models found for ", bio)
    next
  }
  
  # Generate and save model comparison table
  model_table <- getModelTable(all_results, type = "AICc")
  print(model_table)
  write.csv(model_table, file = csv_table_path, row.names = TRUE)
  
  # Get model-averaged parameters
  average_pars <- getModelAvgParams(all_results, type = "AICc")
  
  # Convert expected means to original units
  if (bio == "bio1") {
    average_pars$expected_mean <- exp(average_pars$expected_mean) - 273
    ylab <- "Expected mean (°C)"
  } else if (bio == "bio4") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (temperature seasonality)"
  } else if (bio == "bio12") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation)"
  } else if (bio == "bio15") {
    average_pars$expected_mean <- exp(average_pars$expected_mean)
    ylab <- "Expected mean (precipitation seasonality)"
  } else {
    ylab <- "Expected mean (log-transformed units)"
  }
  
  # Round for readability
  average_pars$expected_mean <- round(average_pars$expected_mean, 2)
  print(range(average_pars$expected_mean))
  
  # Save average parameters (AFTER converting back to original units!)
  write.csv(average_pars, file = csv_pars_path, row.names = TRUE)
  
  # Save boxplots
  plot_prefix <- paste0("results/hOUwie/hOUwie_boxplot_", bio)
  
  png(filename = paste0(plot_prefix, "_mean.png"), width = 800, height = 600)
  boxplot(average_pars$expected_mean ~ average_pars$tip_state,
          xlab = "Tip state", ylab = ylab,
          main = paste("Expected Means -", bio))
  dev.off()
  
  png(filename = paste0(plot_prefix, "_variance.png"), width = 800, height = 600)
  boxplot(average_pars$expected_var ~ average_pars$tip_state,
          xlab = "Tip state", ylab = "Expected variance",
          main = paste("Expected Variances -", bio))
  dev.off()
}

#-------------------------------------------------------------------------------
# Explore data
#-------------------------------------------------------------------------------
library(dplyr)

bio1  <- read.csv("results/hOUwie/hOUwie_average_pars_bio1.csv")
bio4 <- read.csv("results/hOUwie/hOUwie_average_pars_bio4.csv")
bio12 <- read.csv("results/hOUwie/hOUwie_average_pars_bio12.csv")
bio15 <- read.csv("results/hOUwie/hOUwie_average_pars_bio15.csv")

head(bio1)

bio1_summary <- bio1 %>%
  group_by(tip_state) %>%
  summarize(
    mean_temp = round(mean(expected_mean), 2),
    .groups = "drop"
  )

print(bio1_summary)

# tip_state            mean_temp (C)
# social_aboveground        22.1
# social_ground             13.8
# solitary_aboveground      16.3
# solitary_ground           16.9

bio4_summary <- bio4 %>%
  group_by(tip_state) %>%
  summarize(
    mean_seasonality = round(mean(expected_mean), 2),
    .groups = "drop"
  )

print(bio4_summary)

# tip_state            mean_seasonality (SD x 100)
# social_aboveground               186
# social_ground                    618
# solitary_aboveground             351
# solitary_ground                  430

bio12_summary <- bio12 %>%
  group_by(tip_state) %>%
  summarize(
    mean_precip = round(mean(expected_mean), 1),
    .groups = "drop"
  )

print(bio12_summary)

# tip_state            mean_precip (mm)
# social_aboveground          665
# social_ground               727
# solitary_aboveground        729
# solitary_ground             500

bio15_summary <- bio15 %>%
  group_by(tip_state) %>%
  summarize(
    mean_precip_seasonality = round(mean(expected_mean), 2),
    .groups = "drop"
  )

print(bio15_summary)


#-------------------------------------------------------------------------------
# Plotting
#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(viridis)
library(colorspace)
library(patchwork)

# Load results
bio1  <- read.csv("results/hOUwie/hOUwie_average_pars_bio1.csv")
bio4  <- read.csv("results/hOUwie/hOUwie_average_pars_bio4.csv")
bio12 <- read.csv("results/hOUwie/hOUwie_average_pars_bio12.csv")
bio15 <- read.csv("results/hOUwie/hOUwie_average_pars_bio15.csv")

# Add labels
bio1$variable  <- "BIO1: Mean Annual Temperature (°C)"
bio4$variable  <- "BIO4: Temperature Seasonality (SD × 100)"
bio12$variable <- "BIO12: Annual Precipitation (mm)"
bio15$variable <- "BIO15: Precipitation Seasonality (CV)"

# Combine datasets
all_data <- bind_rows(bio1, bio4, bio12, bio15)

# Define factor levels and color mapping
trait_order <- c("solitary_ground", "solitary_aboveground", "social_ground", "social_aboveground")

trait_labels <- c(
  "solitary_ground"     = "Solitary/Ground",
  "solitary_aboveground" = "Solitary/Above-ground",
  "social_ground"       = "Social/Ground",
  "social_aboveground"  = "Social/Above-ground"
)

all_data$tip_state <- factor(all_data$tip_state, levels = trait_order)

box_colors <- viridis::viridis(length(trait_order), option = "cividis")
names(box_colors) <- trait_order

# Generate combined boxplots for expected means
variables <- unique(all_data$variable)
mean_plots <- list()

for (v in variables) {
  subset_data <- filter(all_data, variable == v)
  
  p <- ggplot(subset_data, aes(x = tip_state, y = expected_mean, fill = tip_state)) +
    geom_boxplot(width = 0.6, alpha = 0.6, outlier.size = 0.8, outlier.alpha = 0.4) +
    scale_fill_manual(values = box_colors) +
    scale_x_discrete(labels = trait_labels) +
    labs(x = "", y = "Expected Mean", title = v) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      axis.text.x = element_text(size = 14, color = "black", angle = 15, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  mean_plots[[v]] <- p
}

# Combine plots
combined_mean_plot <- wrap_plots(mean_plots, ncol = 2) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag = element_text(size = 18, face = "bold")
  )


ggsave("results/hOUwie/hOUwie_combined_means.png", combined_mean_plot, width = 12, height = 12)
ggsave("results/hOUwie/hOUwie_combined_means.pdf", combined_mean_plot, width = 12, height = 12)

































# ORIGINAL SCRIPT BELOW
# 
# #-------------------------------------------------------------------------------
# # Setup: Load libraries, set working directory, and source helper functions
# #-------------------------------------------------------------------------------
# #rm(list=ls())
# #install.packages("devtools")
# #devtools::install_github("thej022214/OUwie")
# library(OUwie)
# library(corHMM)
# library(parallel)
# setwd("/Users/lenarh/Desktop/bee_sociality_dist")
# #setwd("/Users/tvasc/Desktop/bee_sociality_dist")
# source("00_utility_functions.R") # provides getModelTable() and getModelAvgParams()?
# 

# #-------------------------------------------------------------------------------
# # Load hOUwie model results (use comments to select BIO4 or BIO1)
# #-------------------------------------------------------------------------------
# all_model_results <- list.files("houwie_results") # List result files
# 
# all_model_results <- subset(all_model_results, grepl("_run2", all_model_results)) # Keep BIO1 only
# #all_model_results <- subset(all_model_results, grepl("_bio4", all_model_results))  # Keep BIO4 only
# 
# model_names <- gsub(".Rsave", "", all_model_results) # Strip extensions for naming
# 
# #-------------------------------------------------------------------------------
# # Load each result into a list named by model
# #-------------------------------------------------------------------------------
# all_results <- list() # Initialize
# 
# for (i in seq_along(all_model_results)) {
#   load(paste0("houwie_results/", all_model_results[i])) # Loads object named 'res'
#   if (exists("res")) {
#     all_results[[i]] <- res
#     names(all_results)[i] <- model_names[i]
#   }
#   remove(res)
# }
# 
# #-------------------------------------------------------------------------------
# # Generate AIC model comparison table
# #-------------------------------------------------------------------------------
# model_table <- getModelTable(all_results, type = "AICc")
# print(model_table)
# 
# write.csv(model_table, file="model_table_results_bio1.csv") # Keep for BIO1
# #write.csv(model_table, file = "model_table_results_bio4.csv") # Keep for BIO4
# 
# 
# # 7/22/25
# # Model comparison table summarizing the performance of your six hOUwie models for temperature seasonality (BIO4)
# #                       np     lnLik   DiscLik   ContLik     AICc     dAICc        AICcwt
# # bm1_8states_run_bio4  12 -3665.862 -649.4574 -3011.799 7355.808 749.88796 1.458472e-163
# # ou1_8states_bio4      13 -3314.317 -657.8155 -2651.897 6654.733  48.81264  2.514601e-11
# # oum_cid_8states_bio4  14 -3325.218 -663.4269 -2657.396 6678.550  72.62935  1.693313e-16
# # oum_full_8states_bio4 16 -3286.887 -660.7572 -2621.426 6605.920   0.00000  1.000000e+00
# # oum_nest_8states_bio4 14 -3317.567 -662.8891 -2649.841 6663.247  57.32636  3.562379e-13
# # oum_soc_8states_bio4  14 -3307.305 -666.6537 -2634.220 6642.723  36.80296  1.019387e-08
# 
# # Best model: oum_full_8states_bio4
# #   - Each trait combination has its own optimum
# #   - AICc = 6605.920, ΔAICc = 0, AICc weight = 1.0
# # Supports hypothesis that trait combinations (sociality × nesting) are tied to distinct climatic niches
# # Simpler models (e.g., BM1, OU1, sociality-only, nesting-only) fit worse
# 
# #-------------------------------------------------------------------------------
# # Get model-averaged parameter estimates for BIO4 (temperature seasonality)
# # These reflect expected means and variances across species under the best-fitting models
# #-------------------------------------------------------------------------------
# average_pars <- getModelAvgParams(all_results, type = "AICc")
# head(average_pars)
# range(average_pars$expected_mean) 
#     # 5.166410 <--> 6.467966 - BIO4
# 
# write.csv(average_pars, file = "average_pars_bio1.csv") # Keep for BIO1
# #write.csv(average_pars, file="average_pars_bio4.csv") # Keep for BIO4
# 
# # These are the expected trait values at the tips (species), averaged across models
# # based on AICc weights. They account for the full hOUwie model:
# # - The history of discrete trait evolution (e.g., sociality × nesting)
# # - The process of continuous trait evolution (e.g., optima, variance, strength of selection)
# # The values are in log scale and need to be transformed back to their original units
# 
# #-------------------------------------------------------------------------------
# # Visualize expected trait values by tip state
# #-------------------------------------------------------------------------------
# # BIO4: Convert log-transformed values back to original scale (standard deviation × 100)
# # average_pars$expected_mean <- exp(average_pars$expected_mean)
# # average_pars$expected_mean <- round(average_pars$expected_mean, 3)
# # head(average_pars$expected_mean)
# # range(average_pars$expected_mean) 
# #     # 175.284 <--> 644.172 - BIO4
# 
# # BIO1: Convert Kelvins back to original scale (°C)
# average_pars$expected_mean <- exp(average_pars$expected_mean) - 273
# average_pars$expected_mean <- round(average_pars$expected_mean, 2)
# head(average_pars$expected_mean)
# range(average_pars$expected_mean) 
# 
# # Boxplots of expected values and variances across tip states (BIO1)
# boxplot(average_pars$expected_mean ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected mean (temperature)")
# 
# boxplot(average_pars$expected_var ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected variance (temperature)")
# 
# # Boxplots of expected values and variances across tip states (BIO4)
# # boxplot(average_pars$expected_mean ~ average_pars$tip_state,
# #         xlab = "Tip state", ylab = "Expected mean (temperature seasonality)")
# # 
# # boxplot(average_pars$expected_var ~ average_pars$tip_state,
# #         xlab = "Tip state", ylab = "Expected variance (temperature seasonality)")
