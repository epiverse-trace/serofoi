test_that("issue 47", {
  skip_on_os(c("windows", "mac"))
  source("testing_utils.R")

  library(dplyr)

  # Load data
  ## This dataset is already prepared
  serodata_path <- testthat::test_path("extdata", "haiti_ssa_sample.RDS")
  serodata <- readRDS(serodata_path)

  # Error reproduction
  model_test <- fit_seromodel(
    serodata = serodata,
    foi_model = "tv_normal",
    print_summary = FALSE
  )
  foi <- rstan::extract(model_test, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = serodata)

  # Test
  age_max <- max(serodata$age_mean_f)
  expect_length(prev_expanded$age, n = age_max)
})
