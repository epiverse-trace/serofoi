rm(list = ls())

library(devtools)
library(dplyr)

# install.packages("github/serofoi")
source("R/modeling.R")
source("R/seroprevalence_data.R")
source("R/model_comparison.R")
source("R/visualisation.R")

#----- Read Data
data_test <- readRDS("data/data.RDS")

survey_test <- data_test$survey[1]
model_object_0 <- run_model(model_data = data_test,
                           survey = survey_test,
                           model_name = "constant_foi_Bi",
                           n_iters = 1000) #the default n_iters=500 yields to an error in a sampling size.

model_object_1 <- run_model(model_data = data_test,
                           survey = survey_test,
                           model_name = "continuous_foi_normal_Bi",
                           n_iters = 1000) #the default n_iters=500 yields to an error in a sampling size.

model_object_2 <- run_model(model_data = data_test,
                           survey = survey_test,
                           model_name = "continuous_foi_normal_log",
                           n_iters = 1000) #the default n_iters=500 yields to an error in a sampling size.

# plots_model_0 <- generate_combined_plots(model_object_0, data_test)


