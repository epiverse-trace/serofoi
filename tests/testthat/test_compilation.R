test_that("compilation", {
  source("testing_utils.R")

  set.seed(1234) # For reproducibility


  data_test <- prepare_data(mydata)

  model_0_object <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )
  model_0_plot <- plot_model(model_0_object, size_text = 6)

  # plot_seroprev_fitted(model_0_object, size_text = 15)
  # plot_foi(model_0_object, size_text = 15)
  # plot_rhats(model_0_object, size_text = 15)
  summary_model <- extract_summary_model(model_0_object)
  expected_summary_model <- read.csv(
    test_path("extdata", "expected_summary_model.csv")
  )


  column_comparation_functions <- list(
    model = equal_exact(),
    dataset = equal_exact(),
    country = equal_exact(),
    year = equal_exact(),
    test = equal_exact(),
    antibody = equal_exact(),
    n_sample = equal_exact(),
    n_agec = equal_exact(),
    n_iter = equal_exact(),
    elpd = equal_with_tolerance(),
    se = equal_with_tolerance(),
    converged = equal_exact()
  )

  expect_true(
    compare_dataframes(
      expected_summary_model, summary_model, column_comparation_functions
    )
  )
})
