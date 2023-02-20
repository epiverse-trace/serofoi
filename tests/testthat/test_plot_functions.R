test_that("plot_seroprev_fitted", {
  library(devtools)
  library(dplyr)
  library(vdiffr)
  set.seed(1234) # For reproducibility
  plot_model_data <- readRDS(test_path("extdata", "data.RDS"))


  plot_data_test <- prepare_data(plot_model_data)

  plot_model_object <- run_model(
    model_data = plot_data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  # model_0_plot <- plot_model(plot_model_object, size_text = 6)

  actual_plot_seroprev_fitted <- plot_seroprev_fitted(plot_model_object, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev_fitted", actual_plot_seroprev_fitted)


  actual_plot_seroprev <- plot_seroprev(plot_data_test, size_text = 15)

  vdiffr::expect_doppelganger("plot_seroprev", actual_plot_seroprev)
})
