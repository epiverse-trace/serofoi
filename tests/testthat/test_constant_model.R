test_that("Test constant model", {
  model_name <- "constant"

  tsur <- 2023
  n_years <- 80
  birth_years <- seq(tsur - n_years, tsur - 1)
  sample_size_by_age <- 100

  foi <- rep(0.03, n_years)

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
    iter = 1000
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
