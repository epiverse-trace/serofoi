test_that("comparison", {
  library(dplyr)
  source("testing_utils.R")

  set.seed(1234) # For reproducibility

  package <- "serofoi"

  expected_comp_table <- read.csv(
    test_path("extdata", "expected_comp_table.csv")
  )
  wrong_comp_table <- read.csv(
    test_path("extdata", "wrong_comp_table.csv")
  )

  data_test <- prepare_data(mydata)

  model_0 <- run_model(
    model_data = data_test,
    model_name = "constant_foi_bi",
    n_iters = 1000
  )

  model_1 <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  model_2 <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )

  comp_table <- get_comparison_table(
    model_objects_list = c(
      m0 = model_0,
      m1 = model_1,
      m2 = model_2
    )
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
    converged = equal_exact(),
    difference = equal_exact(),
    diff_se = equal_with_tolerance(),
    best_elpd = equal_exact(),
    pvalue = equal_with_tolerance()
  )

  expect_true(
    compare_dataframes(
      expected_comp_table, comp_table, column_comparation_functions
    )
  )

  expect_false(
    compare_dataframes(
      wrong_comp_table, comp_table, column_comparation_functions
    )
  )
})
