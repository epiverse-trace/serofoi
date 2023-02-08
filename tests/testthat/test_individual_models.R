# TODO Complete

test_that("individual models", {
  library(devtools)
  library(dplyr)

  #----- Read and prepare data
  test_data_path <- test_path(
    "extdata", "data.RDS"
  )
  test_data <- readRDS(test_data_path) %>% prepare_data(alpha = 0.05)

  #----- Test each model
  model_0_object <- run_model(
    model_data = test_data,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  model_1_object <- run_model(
    model_data = test_data,
    model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  model_2_object <- run_model(
    model_data = test_data,
    model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )

  #----- Generate all plots for each model

  # model_0_plot <- plot_model(model_0_object, size_text = 6)
  # model_1_plot <- plot_model(model_1_object, size_text = 6)
  # model_2_plot <- plot_model(model_2_object, size_text = 6)

  #----- Generate each individual plot

  plot_seroprev_fitted(model_0_object, size_text = 15)
  plot_foi(model_0_object, size_text = 15)
  plot_rhats(model_0_object, size_text = 15)

  # bayesplot::mcmc_trace(model_1_object$fit, pars="lambda0")
})
