rm(list = ls())

library(devtools)
library(dplyr)

source("R/modelling.R")
source("R/seroprevalence_data.R")
source("R/model_comparison.R")
source("R/visualization.R")
mydata <- readRDS("data/data.RDS")

# Modelling module functions
?prepare_data
data_test <- prepare_data(model_data = mydata, alpha = 0.05)

?get_exposure_years
# exposure_years <- get_exposure_years(data_test)

?get_exposure_matrix
# exposure_matrix <- get_exposure_matrix(model_data = data_test)

?save_or_load_model
stan_model <- save_or_load_model(model_name = "constant_foi_bi")

?fit_model
fit_model_test <- fit_model(model_data = data_test,
                            model_name = "constant_foi_bi",
                            n_iters = 1000,
                            n_thin = 2,
                            delta = 0.90,
                            m_treed = 10,
                            decades = 0)


?run_model
model_object <- run_model(model_data = data_test, model_name = "constant_foi_bi")

