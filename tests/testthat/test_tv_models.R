# load and prepare data
data(chagas2012)
serodata <- prepare_serodata(chagas2012, alpha = 0.05)

test_that("Check cohort ages", {
  cohort_ages <- get_cohort_ages(serodata = serodata)
  testthat::expect_equal(nrow(cohort_ages), max(unique(serodata$tsur)) - min(serodata$birth_year))
})

test_that("Test tv_normal model", {
  model_name <- "tv_normal"
  # read benchmark data
  data_compare_path <- testthat::test_path(
    "extdata",
    paste0("prev_expanded_", model_name, ".RDS")
  )
  prev_expanded_compare <- readRDS(data_compare_path) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (predicted_prev_upper - predicted_prev_lower)/2
    )

  # run tv_normal model
  model_object <- run_seromodel(
    serodata = serodata,
    foi_model = model_name,
    n_iters = 1000,
    print_summary = FALSE
  )

  foi <- rstan::extract(model_object, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = serodata)

  # compares expanded prevalence with benchmark
  expect_true(
    all(
      dplyr::near(
        prev_expanded$predicted_prev,
        prev_expanded_compare$predicted_prev,
        tol = prev_expanded_compare$tol
      )
    )
  )
})

test_that("Test tv_normal_log model", {
  model_name <- "tv_normal_log"
  # read benchmark data
  data_compare_path <- testthat::test_path(
    "extdata",
    paste0("prev_expanded_", model_name, ".RDS")
  )
  prev_expanded_compare <- readRDS(data_compare_path) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (predicted_prev_upper - predicted_prev_lower)/2
    )

  # run tv_normal_log model
  model_object <- run_seromodel(
    serodata = serodata,
    foi_model = model_name,
    n_iters = 1000,
    print_summary = FALSE
  )

  foi <- rstan::extract(model_object, "foi", inc_warmup = FALSE)[[1]]
  prev_expanded <- get_prev_expanded(foi, serodata = serodata)

  # compares expanded prevalence with benchmark
  # We use dplyr::near() rather than expect_equal() to allow passing a vector 
  # of tolerances.
  expect_true(
    all(
      dplyr::near(
        prev_expanded$predicted_prev,
        prev_expanded_compare$predicted_prev,
        tol = prev_expanded_compare$tol
      )
    )
  )
})
