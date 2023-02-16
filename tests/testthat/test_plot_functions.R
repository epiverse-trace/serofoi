# TODO complete

test_that("minimal data", {
  library(devtools)
  library(dplyr)

  data_test <- prepare_data(mydata)

  model_0_object <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  model_0_plot <- plot_model(model_0_object, size_text = 6)

  plot_seroprev_fitted(model_0_object, size_text = 15)
  plot_seroprev(data_test, size_text = 15)

  # TODO Complete test ###
})
