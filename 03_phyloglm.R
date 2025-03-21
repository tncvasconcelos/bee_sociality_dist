# Project: TROPICAL BEE SOCIALITY/NESTING
# Authors: Lena Heinrich, Aline Martins, Thais Vasconcelos
# University of Michigan, Ann Arbor

# 03 PHYLOGENETIC GENERALISED LINEAR MODEL

################################################################################

# Setup:

# Libraries
rm(list = ls())
setwd("~/Desktop/bee_sociality_dist")
library(phytools)
library(geiger)
library(dplyr)
library(phylolm)
library(tidyr)
library(ggplot2)

# Loading traits, tree, and climatic data
traits <- read.csv("curated_data/bees_traits.csv") # contains trait data for bees in tree
traits_tmp <- traits  # create separate object for manipulation
tree <- read.tree("curated_data/ML_beetree_pruned.tre")
tree_tmp <- tree # create separate object for manipulation
climate_files <- list.files("curated_data", "summstats.csv", full.names = TRUE) 

################################################################################

# 1/8/25 Exploring curated traits dataset:

total_species <- nrow(traits)

# Summarize sociality_three_states and nest_three_states with counts and percentages
sociality_summary <- traits %>%
  count(sociality_three_states, name = "Count") %>%
  mutate(Percentage = (Count / total_species) * 100)

nest_summary <- traits %>%
  count(nest_three_states, name = "Count") %>%
  mutate(Percentage = (Count / total_species) * 100)

# Print summaries
print(sociality_summary)
print(nest_summary)


# Repeat using sociality_binary and nest_binary
# Summarize sociality_three_states and nest_three_states with counts and percentages
sociality_binary_summary <- traits %>%
  count(sociality_binary, name = "Count") %>%
  mutate(Percentage = (Count / total_species) * 100)

nest_binary_summary <- traits %>%
  count(nest_binary, name = "Count") %>%
  mutate(Percentage = (Count / total_species) * 100)

# Print summaries
print(sociality_binary_summary)
print(nest_binary_summary)

# Summarize nest_binary and sociality_binary by family
family_summary <- traits %>%
  group_by(family, nest_binary, sociality_binary) %>%
  summarize(Count = n(), .groups = "drop")

# Print the summary table
print(family_summary)

# Confirm counts are correct
nrow(subset(traits, family == "Mellitidae")) # 82
nrow(subset(traits, family == "Stenotritidae")) # 6
nrow(subset(traits, family == "Andrenidae")) # 644

################################################################################

# Prepare data for phyloGLM:

# Extract mean climate values and join with trait data
for (file_index in 1:length(climate_files)) {
  
  # Read the climate file
  temp_data <- read.csv(climate_files[file_index])
  
  # Extract relevant columns
  temp_data <- temp_data[, c(1, grep("mean", colnames(temp_data)))]
  
  # Remove NAs from the climate data
  temp_data <- temp_data[complete.cases(temp_data), ]
  
  # Merge with the traits data
  traits_tmp <- merge(traits_tmp, temp_data, by.x = "tips", by.y = "species", all = FALSE)
}

# Ensure that tree and trait dataset contain the same species
traits_tmp <- traits_tmp[complete.cases(traits_tmp), ] # remove rows with NAs from traits_tmp
tree_tmp <- drop.tip(tree_tmp, setdiff(tree_tmp$tip.label, traits_tmp$tips)) # remove unmatched species
setdiff(tree_tmp$tip.label, traits_tmp$tips) # character(0)
# NOTE: this pruning step must happen after merging climate and trait data in loop above!

# Inspect the merged data which now contains trait data and mean climate values
head(traits_tmp)
colnames(traits_tmp)

# Make sociality and nesting binary before putting into phyloGLM
traits_tmp$sociality_binary[traits_tmp$sociality_binary=="solitary"] <- 0
traits_tmp$sociality_binary[traits_tmp$sociality_binary=="social"] <- 1

traits_tmp$nest_binary[traits_tmp$nest_binary=="ground"] <- 0
traits_tmp$nest_binary[traits_tmp$nest_binary=="aboveground"] <- 1

unique(traits_tmp$nest_binary)
unique(traits_tmp$sociality_binary)

# Make sociality and nesting numeric
traits_tmp$sociality_binary <- as.numeric(traits_tmp$sociality_binary)
traits_tmp$nest_binary <- as.numeric(traits_tmp$nest_binary)

# Create vector of climate variables to loop through
climate_vars <- grep("mean", colnames(traits_tmp), value = TRUE) # Identify climate variables
rownames(traits_tmp) <- traits_tmp$tips
print(climate_vars)

# Center and standardize predictor (climate) variables
for(climate_var in 1:length(climate_vars)){
  
  curr_var <- climate_vars[climate_var] # name of current climate variable being processed by loop
 
  curr_col <- which(colnames(traits_tmp) == curr_var) # index of column that matches curr_var
  
  traits_tmp[ , curr_col] <- (traits_tmp[ , curr_col] - mean(traits_tmp[ , curr_col])) / sd(traits_tmp[ , curr_col])

  } # result is columns where the values have a mean of 0 and sd of 1 
# standardized columns replace original ones in the dataframe

head(traits_tmp)

# Check for NAs in dataset
colSums(is.na(traits_tmp)) # No NAs in any columns

################################################################################

# Run phyloGLMs:

# phyloGLM for sociality
sociality_results <- list()
sink("results/phyloGLM_sociality_result.txt") # prints sociality results to text file

for (climate_var in climate_vars) {
  
  # Create formula
  formula <- formula(paste("sociality_binary ~", climate_var, collapse = " "))
  
  # Fit phyloGLM
  model <- phyloglm(formula, phy = tree_tmp, btol = 10, data = traits_tmp, method = "logistic_MPLE")
  
  # Store model summary in the sociality results list
  sociality_results[[climate_var]] <- summary(model)
  
  print(summary(model))
}

print(sociality_results)
# none significant

# bio_2 p = 0.9884
# bio_4 p = 0.2007
# bio_12 p = 0.5928
# bio_15 p = 0.1941

sink()

# phyloGLM for nesting
nesting_results <- list()
sink("results/phyloGLM_nesting_result.txt") # prints nesting results to text file

for (climate_var in climate_vars) {
  
  # Create formula
  formula <- formula(paste("nest_binary ~", climate_var, collapse = " "))
  # Fit phyloGLM
  model <- phyloglm(formula, phy = tree_tmp, btol = 10, data = traits_tmp, method = "logistic_MPLE")
  
  # Store model summary in the nesting results list
  nesting_results[[climate_var]] <- summary(model)
  
  print(summary(model))
}

print(nesting_results)
# two marginally significant: mean_bio_13 and mean_bio_16, p = 0.09

sink()

# bio_1 p = 0.4253
# bio_4 p = 0.2193
# bio_12 p = 0.2016
# bio_15 p = 0.7411


################################################################################

# Interpreting p-values from phyloGLM:

# p-values = tests of whether each climate variable has a 
#   statistically significant effect on the likelihood of the bees' sociality or nesting biology, 
#   while accounting for phylogenetic relationships.

################################################################################

# Visualizing results:

# Open a PDF for sociality boxplots
pdf("plots/phyloGLM_sociality_boxplots.pdf", width = 10, height = 25)

# Set up the plotting area for a multi-panel plot
n_vars <- length(climate_vars)
ncols <- 3  # Number of columns
nrows <- ceiling(n_vars / ncols)  # Number of rows

# Adjust margins
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1)) 

# Loop through the climate variables for sociality data
for (climate_var in climate_vars) {
  
  # Extract p-value from the phyloGLM results
  p_value <- sociality_results[[climate_var]]$coefficients[2, 4]  # p-value for the climate variable
  
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
  
  # Boxplot for sociality
  boxplot(traits_tmp[[climate_var]] ~ traits_tmp$sociality_binary, 
          xlab = " ", ylab = "Env. Var", 
          main = paste(climate_var, significance),  # Add significance symbol to title
          names = c("Solitary", "Social"))
}

# Close PDF
dev.off()


# Open a PDF for nesting boxplots
pdf("plots/phyloGLM_nesting_boxplots.pdf", width = 10, height = 25)

# Set up the plotting area for a multi-panel plot
par(mfrow = c(nrows, ncols), mar = c(5, 5, 3, 1))  # Adjusted for better spacing

# Loop through the climate variables for nesting data
for (climate_var in climate_vars) {
  
  # Extract p-value from the phyloGLM results
  p_value <- nesting_results[[climate_var]]$coefficients[2, 4]  # p-value for the climate variable
  
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
  
  # Boxplot for nesting
  boxplot(traits_tmp[[climate_var]] ~ traits_tmp$nest_binary, 
          xlab = " ", ylab = "Env. Var", 
          main = paste(climate_var, significance), # Add significance symbol to title
          names = c("Ground", "Above-ground"))
}

# Close PDF
dev.off()
