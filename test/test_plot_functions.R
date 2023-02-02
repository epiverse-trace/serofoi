rm(list = ls())

library(devtools)
library(dplyr)

# remotes::install_github("TRACE-LAC/serofoi", ref = "dev")
# library(serofoi)
data_test <- prepare_data(mydata)

# source("R/modelling.R")
# source("R/seroprevalence_data.R")
# source("R/model_comparison.R")
# source("R/visualization.R")
# data_test <- readRDS("data/data.RDS") %>% prepare_data(alpha = 0.05)


model_0_object <- run_model(
  model_data = data_test,
  model_name = "constant_foi_bi",
  n_iters = 1000
)

model_0_plot <- plot_model(model_0_object, size_text = 6)

plot_seroprev_fitted(model_0_object, size_text = 15)
plot_seroprev(model_0_object, size_text = 15)
