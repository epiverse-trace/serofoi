test_that("issue 47", {
  skip_on_os(c("windows", "mac"))
  source("testing_utils.R")

  library(devtools)
  library(dplyr)

  # Load data
  ## This dataset is already prepared
  data_path <- testthat::test_path("extdata", "haiti_ssa_sample.RDS")
  data_issue <- readRDS(data_path)

  # Error reproduction
  model_test <- run_seromodel(data_issue, foi_model = "tv_normal", print_summary = FALSE)
  foi <- rstan::extract(model_test$seromodel_fit, "foi", inc_warmup = FALSE)[[1]]
  age_max <- max(data_issue$age_mean_f)
  prev_expanded <- get_prev_expanded(foi, serodata = data_issue)

  # Test
  expect_length(prev_expanded$age, n = age_max)
})
