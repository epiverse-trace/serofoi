test_that("comparison", {
  library(dplyr)
  set.seed(1234) # For reproducibility
  print("test path")

  compare_with_tolerance <- function(tolerance = 1e-4) {
    function(a, b, tolerance = 1e-4) {
      x <- mapply(function(x, y) abs(x - y), a, b)
      base::all(x < tolerance)
    }
  }
  compare_exact <- function() {
    function(a, b) {
      x <- mapply(function(x, y) x == y, a, b)
      base::all(x == TRUE)
    }
  }
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
  print("model 0 ----------")
  model_1 <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_bi",
    n_iters = 1000
  )

  print("model 1 ----------")
  model_2 <- run_model(
    model_data = data_test,
    model_name = "continuous_foi_normal_log",
    n_iters = 1000
  )
  print("model 2 ----------")

  comp_table <- get_comparison_table(
    model_objects_list = c(
      m0 = model_0,
      m1 = model_1,
      m2 = model_2
    )
  )

  col_tests <- list(
    model = compare_exact(),
    dataset = compare_exact(),
    country = compare_exact(),
    year = compare_exact(),
    test = compare_exact(),
    antibody = compare_exact(),
    n_sample = compare_exact(),
    n_agec = compare_exact(),
    n_iter = compare_exact(),
    elpd = compare_with_tolerance(),
    se = compare_with_tolerance(),
    converged = compare_exact(),
    difference = compare_exact(),
    diff_se = compare_with_tolerance(),
    best_elpd = compare_exact(),
    pvalue = compare_with_tolerance()
  )

  print(comp_table)
  for (col in names(col_tests)) {
    test_fun <- col_tests[[col]]
    testthat::expect_true(
      test_fun(expected_comp_table[col], comp_table[col])
    )

  }
})
