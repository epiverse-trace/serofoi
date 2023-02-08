remotes::install_github("TRACE-LAC/serofoi", ref = "dev-docker-tests", force = TRUE)
set.seed(1234) # For reproducibility
library(serofoi)

data_test <- prepare_data(mydata)
# ?run_model

model_0 <- run_model(model_data = data_test,
                     model_name = "constant_foi_bi",
                     n_iters = 1000)

model_1 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_bi",
                     n_iters = 1000)

model_2 <- run_model(model_data = data_test,
                     model_name = "continuous_foi_normal_log",
                     n_iters = 1000)

comp_table <- get_comparison_table(
  model_objects_list = c(m0 = model_0,
                         m1 = model_1,
                         m2 = model_2))

print(comp_table)
