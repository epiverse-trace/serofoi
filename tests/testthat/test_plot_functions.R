

test_that("plot_seroprev_fitted", {
  library(devtools)
  library(dplyr)
  source("testing_utils.R")
  set.seed(1234) # For reproducibility


  data_test <- prepare_data(mydata)

  model_0_object <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  # model_0_plot <- plot_model(model_0_object, size_text = 6)

  actual_plot <- plot_seroprev_fitted(model_0_object, size_text = 15)

  expect_true(compare_plots("plot_seroprev_fitted", actual_plot))
})


test_that("plot_seroprev", {
  library(devtools)
  library(dplyr)

  data_test <- prepare_data(mydata)

  model_0_object <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  # model_0_plot <- plot_model(model_0_object, size_text = 6)

  actual_plot <- plot_seroprev(data_test, size_text = 15)

  expect_true(compare_plots("plot_seroprev", actual_plot))
})
