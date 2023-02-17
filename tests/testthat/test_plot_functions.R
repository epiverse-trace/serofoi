test_that("plot_seroprev_fitted", {
  library(devtools)
  library(dplyr)
  .pardefault <- par()
  source("testing_utils.R")
  set.seed(1234) # For reproducibility
  plot_model_data <- readRDS(test_path("extdata", "data.RDS"))


  plot_data_test <- prepare_data(plot_model_data)

  plot_model_object <- run_model(
    model_data = plot_data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  # model_0_plot <- plot_model(plot_model_object, size_text = 6)

  actual_plot <- plot_seroprev_fitted(plot_model_object, size_text = 15)

  expect_true(compare_plots("plot_seroprev_fitted", actual_plot))
})


test_that("plot_seroprev", {
  library(devtools)
  library(dplyr)

  plot_model_data <- readRDS(test_path("extdata", "data.RDS"))

  plot_data_test <- prepare_data(plot_model_data)

  plot_model_object <- run_model(
    model_data = plot_data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  # model_0_plot <- plot_model(plot_model_object, size_text = 6)

  actual_plot <- plot_seroprev(plot_data_test, size_text = 15)

  expect_true(compare_plots("plot_seroprev", actual_plot))
})
