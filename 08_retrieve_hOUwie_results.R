# rm(list=ls())
library(corHMM)
library(OUwie)
#devtools::install_github("thej022214/corHMM")
library(parallel)
#setwd("/Users/tvasc/Desktop/bee_sociality_dist")

all_model_results <- list.files("houwie_results")
all_model_results <- subset(all_model_results, grepl("_sociality_bio1_run1", all_model_results))
model_names <- gsub("_sociality_bio1_run2.Rsave","",all_model_results)

all_results <- list()
for(i in 1:length(all_model_results)) {
  load(paste0("houwie_results/",all_model_results[i]))
  if(exists("res")) {
    all_results[[i]] <- res
    names(all_results)[i] <- model_names[i]
  }
  remove(res)
}

model_table <- getModelTable(all_results, type="AICc")
write.csv(model_table, file="model_table_results.csv")

average_pars <- getModelAvgParams(all_results, type="AICc")

write.csv(average_pars, file="average_pars.csv")

# what the tip values are expected to be when accounting for the models of evolution
# that is a way to integrate the other parameter of the model
average_pars$expected_mean <- exp(average_pars$expected_mean)
average_pars$expected_mean <- average_pars$expected_mean - 273
average_pars$expected_mean <- round(average_pars$expected_mean, 3)
boxplot(average_pars$expected_mean~average_pars$tip_state, xlab="tip state", ylab="expected mean (temperature)")
boxplot(average_pars$expected_var~average_pars$tip_state, xlab="tip state", ylab="expected mean (temperature)")
