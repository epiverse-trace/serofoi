remotes::install_github("TRACE-LAC/serofoi", ref = "dev-zulma", force = TRUE)

library(serofoi)

data_test <- prepare_data(mydata)
?run_model

# Testing model comparison
#I'm sourcing here because installing from my branch is not working

source("R/model_comparison.R")
source("R/modelling.R")

model_0 <- run_model(model_data = data_test,
                     model_name = "constant_foi_bi",
                     n_iters = 1000)

model_1 <- run_model(model_data = data_test,
                     model_name = "",
                     n_iters = 1000)

model_2 <- run_model(model_data = data_test,
                     model_name = "",
                     n_iters = 1000)

