test_that("Test tv_normal model", {
  model_name <- "tv_normal"

  tsur <- 2023
  n_years <- 75
  birth_years <- seq(tsur - n_years, tsur - 1)
  sample_size_by_age <- 100

  foi <- rep(
    c(0.01, 0.03, 0.06),
    c(25, 25, 25)
  )

  simdata <- generate_sim_data(
    sim_data = data.frame(
      age = seq(1:n_years),
      tsur = tsur
    ),
    foi = foi,
    sample_size_by_age = sample_size_by_age
  ) %>%
  prepare_serodata()

  model_object <- fit_seromodel(
    serodata = simdata,
    foi_model = model_name,
    chunk_size = 25,
    iter = 1500
  )

  cohort_ages <- get_cohort_ages(simdata)
  foi_central_estimates <- get_foi_central_estimates(
    seromodel_object = model_object,
    cohort_ages = cohort_ages
  ) %>%
    mutate(
      # calculates tolerance as half the confidence interval size
      tol = (upper - lower)/2
    )

  expect_true(
    all(
      dplyr::near(
        foi_central_estimates$medianv,
        foi,
        tol = foi_central_estimates$tol
      )
    )
  )
})


test_that("Test tv_normal_log model", {
  model_name <- "tv_normal_log"

  tsur <- 2023
  n_years <- 60
  birth_years <- seq(tsur - n_years, tsur - 1)
  sample_size_by_age <- 100


  foi <- rep(
    c(0.25, 0.1, 0.01, 0.0001),
    c(10, 10, 10, n_years-30)
  )

  chunks <- rep(
    c(1, 2, 3, 4),
    c(10, 10, 10, n_years-30)
  )

  simdata <- generate_sim_data(
    sim_data = data.frame(
      age = seq(1:n_years),
      tsur = tsur
    ),
    foi = foi,
    sample_size_by_age = sample_size_by_age
  ) %>%
  prepare_serodata()

  model_object <- fit_seromodel(
    serodata = simdata,
    foi_model = model_name,
    chunks = chunks,
    iter = 1500
  )

  cohort_ages <- get_cohort_ages(simdata)
  foi_central_estimates <- get_foi_central_estimates(
    seromodel_object = model_object,
    cohort_ages = cohort_ages
  ) %>%
  mutate(
    # calculates tolerance as half the confidence interval size
    tol = (upper - lower)/2
  )

  expect_true(
    all(
      dplyr::near(
        foi_central_estimates$medianv,
        foi,
        tol = foi_central_estimates$tol
      )
    )
  )
})
