remotes::install_github("TRACE-LAC/serofoi", ref = "dev-zulma")

library(serofoi)

data_test <- prepare_data(mydata)

# Testing model comparison


model_0 <- run_model(model_data = data_test,
                     model_name = "constant_foi_bi",
                     n_iters = 1000)

model_1 <- run_model(model_data = data_test,
                     model_name = "",
                     n_iters = 1000)

model_2 <- run_model(model_data = data_test,
                     model_name = "",
                     n_iters = 1000)

