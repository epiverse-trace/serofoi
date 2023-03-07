test_that("plot_seroprev_fitted", {
  # So far we are skipping tests on these platforms until
  # we find an efficient way to update rstan testthat snapshots on all of them
  skip_on_os(c("windows", "mac"))
  skip_on_ci()
  print("*** Test info ****")
  print(R.Version())
  cat("Interactive: ", interactive())
  library(devtools)
  library(dplyr)
  library(vdiffr)
  set.seed(1234) # For reproducibility
  data_test <- readRDS(testthat::test_path("extdata", "data.RDS")) %>% prepare_seroprev_data()

  actual_plot_seroprev <- plot_seroprev(data_test, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev", actual_plot_seroprev)

  # Constant Model
  model_object_constant <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "constant_foi_bi",
    n_iters = 1000
  )

  plot_seroprev_model_constant <- plot_seroprev_model(model_object_constant, size_text = 6)

  vdiffr::expect_doppelganger("plot_seroprev_model_constant", plot_seroprev_model_constant)

  plot_seroprev_fitted_constant <- plot_seroprev_fitted(model_object_constant, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_constant", plot_seroprev_fitted_constant)

  # Normal Bi Model
  model_object_normalbi <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  plot_seroprev_model_normalbi <- plot_seroprev_model(model_object_normalbi, size_text = 6)

  vdiffr::expect_doppelganger("plot_seroprev_model_normalbi", plot_seroprev_model_normalbi)

  plot_seroprev_fitted_normalbi <- plot_seroprev_fitted(model_object_normalbi, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_normalbi", plot_seroprev_fitted_normalbi)

  # Normal Log Model
  model_object_normallog <- run_seroprev_model(
    seroprev_data = data_test,
    seroprev_model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )

  plot_seroprev_model_normallog <- plot_seroprev_model(model_object_normallog, size_text = 6)

  vdiffr::expect_doppelganger("plot_seroprev_model_normallog", plot_seroprev_model_normallog)

  plot_seroprev_fitted_normallog <- plot_seroprev_fitted(model_object_normallog, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_normallog", plot_seroprev_fitted_normallog)

  # Models Comparison Plot
  plot_arrange_models <- plot_seroprev_models_grid(plot_seroprev_model_constant, plot_seroprev_model_normalbi, plot_seroprev_model_normallog, n_col = 3)

  vdiffr::expect_doppelganger("plot_arrange_models", plot_arrange_models)
})
