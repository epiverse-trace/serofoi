test_that("plot_seroprev_fitted", {
  library(devtools)
  library(dplyr)
  library(vdiffr)
  set.seed(1234) # For reproducibility
  data_test <- readRDS(test_path("extdata", "data.RDS")) %>% prepare_data()

  actual_plot_seroprev <- plot_seroprev(data_test, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev", actual_plot_seroprev)

  # Constant Model
  model_object_constant <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  plot_model_constant <- plot_model(model_object_constant, size_text = 6)

  vdiffr::expect_doppelganger("plot_model_constant", plot_model_constant)

  plot_seroprev_fitted_constant <- plot_seroprev_fitted(model_object_constant, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_constant", plot_seroprev_fitted_constant)

  # Normal Bi Model
  model_object_normalbi <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  plot_model_normalbi <- plot_model(model_object_normalbi, size_text = 6)

  vdiffr::expect_doppelganger("plot_model_normalbi", plot_model_normalbi)

  plot_seroprev_fitted_normalbi <- plot_seroprev_fitted(model_object_normalbi, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_normalbi", plot_seroprev_fitted_normalbi)

  # Normal Log Model
  model_object_normallog <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )

  plot_model_normallog <- plot_model(model_object_normallog, size_text = 6)

  vdiffr::expect_doppelganger("plot_model_normallog", plot_model_normallog)

  plot_seroprev_fitted_normallog <- plot_seroprev_fitted(model_object_normallog, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted_normallog", plot_seroprev_fitted_normallog)

  # Models Comparison Plot
  plot_arrange_models <- plot_models_list(plot_model_constant, plot_model_normalbi, plot_model_normallog, n_col = 3)

  vdiffr::expect_doppelganger("plot_arrange_models", plot_arrange_models)
})
