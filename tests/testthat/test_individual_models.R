
test_that("individual models", {
  library(devtools)
  library(dplyr)
  library(vdiffr)

  #----- Read and prepare data
  data_test_path <- test_path(
    "extdata", "data.RDS"
  )
  data_test <- readRDS(data_test_path) %>% preprare_seroprev_data(alpha = 0.05)

  #----- Test each model
  model_0_object <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "constant_foi_bi",
    n_iters = 1000
  )

  model_1_object <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  model_2_object <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )

  #----- Generate all plots for each model

  model_0_plot <- plot_seroprev_model(model_0_object, size_text = 6)
  vdiffr::expect_doppelganger("model_0_plot", model_0_plot)
  model_1_plot <- plot_seroprev_model(model_1_object, size_text = 6)
  vdiffr::expect_doppelganger("model_1_plot", model_1_plot)
  model_2_plot <- plot_seroprev_model(model_2_object, size_text = 6)
  vdiffr::expect_doppelganger("model_2_plot", model_2_plot)

  #----- Generate each individual plot

  sp_fitted_individual_plot <- plot_seroprev_fitted(model_2_object, size_text = 15)
  vdiffr::expect_doppelganger("seroprev_fitted_individual_plot", sp_fitted_individual_plot)

  foi_individual_plot <- plot_foi(model_2_object, size_text = 15)
  vdiffr::expect_doppelganger("foi_individual_plot", foi_individual_plot)

  rhats_individual_plot <- plot_rhats(model_2_object, size_text = 15)
  vdiffr::expect_doppelganger("rhats_individual_plot", rhats_individual_plot)


  # bayesplot::mcmc_trace(model_1_object$fit, pars="lambda0")

})
