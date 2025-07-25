# ==============================================================================
# 07.3 hOUwie Results Processing
# ==============================================================================
# Loads and compares fitted hOUwie models using temperature seasonality (BIO4)
# or mean annual temperature (BIO1)

# You will need to manually comment and uncomment lines pertaining to which variable to select!

# Outputs a model comparison table and model-averaged parameter estimates
# ==============================================================================

#-------------------------------------------------------------------------------
# Setup: Load libraries, set working directory, and source helper functions
#-------------------------------------------------------------------------------
#rm(list=ls())
#install.packages("devtools")
#devtools::install_github("thej022214/OUwie")
library(OUwie)
library(corHMM)
library(parallel)
setwd("/Users/lenarh/Desktop/bee_sociality_dist")
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")
source("00_utility_functions.R") # provides getModelTable() and getModelAvgParams()?

#-------------------------------------------------------------------------------
# Load hOUwie model results (use comments to select BIO4 or BIO1
#-------------------------------------------------------------------------------
all_model_results <- list.files("houwie_results") # List result files

all_model_results <- subset(all_model_results, grepl("_run2", all_model_results)) # Keep BIO1 only
#all_model_results <- subset(all_model_results, grepl("_bio4", all_model_results))  # Keep BIO4 only

model_names <- gsub(".Rsave", "", all_model_results) # Strip extensions for naming

#-------------------------------------------------------------------------------
# Load each result into a list named by model
#-------------------------------------------------------------------------------
all_results <- list() # Initialize

for (i in seq_along(all_model_results)) {
  load(paste0("houwie_results/", all_model_results[i])) # Loads object named 'res'
  if (exists("res")) {
    all_results[[i]] <- res
    names(all_results)[i] <- model_names[i]
  }
  remove(res)
}

#-------------------------------------------------------------------------------
# Generate AIC model comparison table
#-------------------------------------------------------------------------------
model_table <- getModelTable(all_results, type = "AICc")
print(model_table)

write.csv(model_table, file="model_table_results_bio1.csv") # Keep for BIO1
#write.csv(model_table, file = "model_table_results_bio4.csv") # Keep for BIO4


# 7/22/25
# Model comparison table summarizing the performance of your six hOUwie models for temperature seasonality (BIO4)
#                       np     lnLik   DiscLik   ContLik     AICc     dAICc        AICcwt
# bm1_8states_run_bio4  12 -3665.862 -649.4574 -3011.799 7355.808 749.88796 1.458472e-163
# ou1_8states_bio4      13 -3314.317 -657.8155 -2651.897 6654.733  48.81264  2.514601e-11
# oum_cid_8states_bio4  14 -3325.218 -663.4269 -2657.396 6678.550  72.62935  1.693313e-16
# oum_full_8states_bio4 16 -3286.887 -660.7572 -2621.426 6605.920   0.00000  1.000000e+00
# oum_nest_8states_bio4 14 -3317.567 -662.8891 -2649.841 6663.247  57.32636  3.562379e-13
# oum_soc_8states_bio4  14 -3307.305 -666.6537 -2634.220 6642.723  36.80296  1.019387e-08

# Best model: oum_full_8states_bio4
#   - Each trait combination has its own optimum
#   - AICc = 6605.920, ΔAICc = 0, AICc weight = 1.0
# Supports hypothesis that trait combinations (sociality × nesting) are tied to distinct climatic niches
# Simpler models (e.g., BM1, OU1, sociality-only, nesting-only) fit worse

#-------------------------------------------------------------------------------
# Get model-averaged parameter estimates for BIO4 (temperature seasonality)
# These reflect expected means and variances across species under the best-fitting models
#-------------------------------------------------------------------------------
average_pars <- getModelAvgParams(all_results, type = "AICc")
head(average_pars)
range(average_pars$expected_mean) 
    # 5.166410 <--> 6.467966 - BIO4

write.csv(average_pars, file = "average_pars_bio1.csv") # Keep for BIO1
#write.csv(average_pars, file="average_pars_bio4.csv") # Keep for BIO4

# These are the expected trait values at the tips (species), averaged across models
# based on AICc weights. They account for the full hOUwie model:
# - The history of discrete trait evolution (e.g., sociality × nesting)
# - The process of continuous trait evolution (e.g., optima, variance, strength of selection)
# The values are in log scale and need to be transformed back to their original units

#-------------------------------------------------------------------------------
# Visualize expected trait values by tip state
#-------------------------------------------------------------------------------
# BIO4: Convert log-transformed values back to original scale (standard deviation × 100)
# average_pars$expected_mean <- exp(average_pars$expected_mean)
# average_pars$expected_mean <- round(average_pars$expected_mean, 3)
# head(average_pars$expected_mean)
# range(average_pars$expected_mean) 
#     # 175.284 <--> 644.172 - BIO4

# BIO1: Convert Kelvins back to original scale (°C)
average_pars$expected_mean <- exp(average_pars$expected_mean) - 273
average_pars$expected_mean <- round(average_pars$expected_mean, 2)
head(average_pars$expected_mean)
range(average_pars$expected_mean) 

# Boxplots of expected values and variances across tip states (BIO1)
boxplot(average_pars$expected_mean ~ average_pars$tip_state,
        xlab = "Tip state", ylab = "Expected mean (temperature)")

boxplot(average_pars$expected_var ~ average_pars$tip_state,
        xlab = "Tip state", ylab = "Expected variance (temperature)")

# Boxplots of expected values and variances across tip states (BIO4)
# boxplot(average_pars$expected_mean ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected mean (temperature seasonality)")
# 
# boxplot(average_pars$expected_var ~ average_pars$tip_state,
#         xlab = "Tip state", ylab = "Expected variance (temperature seasonality)")
