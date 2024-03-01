library(dplyr)
library(purrr)
library(serofoi)

seed <- 1234
tsur <- 2050
birth_year_min <- 2000
n_years <- tsur - birth_year_min
sample_size_by_age <- rep(10^7, n_years)
tolerance = 10^-3

test_that("Test simulated data from constant force-of-infection", {
  foi_values <- c(0.001, 0.01, 0.1, 0.3, 0.4)
  for (foi_value in foi_values) {
    sim_data <- data.frame(
      age_mean_f = seq(1,n_years),
      tsur = tsur
    )
    foi_sim <- rep(foi_value, tsur - birth_year_min)
    #----- Test function generate_sim_data
    sim_data <- generate_sim_data(
      sim_data = sim_data,
      foi = foi_sim,
      sample_size_by_age = sample_size_by_age,
      survey_label = "foi_sim_constant",
      seed = seed
    ) %>%
    prepare_serodata()

    prev_exact <- 1 - exp(-foi_value * sim_data$age_mean_f)

    expect_s3_class(sim_data, "data.frame")
    expect_length(sim_data$birth_year, tsur - birth_year_min)
    expect_true(
      all(
        dplyr::near(
          sim_data$prev_obs,
          prev_exact,
          tol = tolerance
        )
      )
    )

    #----- Test function group_sim_data
    sim_data_grouped <- group_sim_data(sim_data = sim_data)
    expect_s3_class(sim_data_grouped, "data.frame")
    expect_s3_class(sim_data_grouped$age_group, "factor")
  }
})

test_that("Test simulated data from time-varying force-of-infection", {
  sim_data <- data.frame(
    age_mean_f = seq(1,n_years),
    tsur = tsur
  )

  no_transm <- 0.0000000001
  foi_sim <- c(rep(0.2, 25), rep(0.1, 10), rep(no_transm, 15))

  sim_data <- generate_sim_data(
    sim_data = sim_data,
    foi = foi_sim,
    sample_size_by_age = sample_size_by_age,
    survey_label = "sw_dec_foi",
    seed = seed
  ) %>%
  prepare_serodata()

  prev_exact <- 1 - exp(-cumsum(rev(foi_sim)))

  expect_s3_class(sim_data, "data.frame")
  expect_length(sim_data$birth_year, tsur - birth_year_min)
  expect_true(
    all(
      dplyr::near(
        sim_data$prev_obs,
        prev_exact,
        tol = tolerance
      )
    )
  )
})

test_that("Test simulated data from age-varying force-of-infection", {
  #----- Exact age-varying probability with seroreversion

  #----- Test without seroreversion
  mu <- 0.
  foi_sim <- rep(0.01, n_years)

  sim_data <- data.frame(
    age_mean_f = seq(1,n_years),
    tsur = tsur
  )

  sim_data <- generate_sim_data(
    sim_data = sim_data,
    foi = foi_sim,
    sample_size_by_age = sample_size_by_age,
    mu = mu,
    model_type = "age-varying",
    survey_label = "age-varying-foi",
    seed = seed
  ) %>%
  prepare_serodata()

  prev_exact <- (foi_sim/(foi_sim + mu)) * (1 - exp(-seq(1,n_years) * (foi_sim + mu)))

  expect_s3_class(sim_data, "data.frame")
  expect_length(sim_data$birth_year, tsur - birth_year_min)
  expect_true(
    all(
      dplyr::near(
        sim_data$prev_obs,
        prev_exact,
        tol = tolerance
      )
    )
  )

  #----- Test with seroreversion
  mu <- 0.1
  fois <- c(0.01, 0.05, 0.1, 0.03, 0.01)
  foi_sim <- unlist(lapply(fois, function(x) rep(x, n_years / length(fois))))

  sim_data <- data.frame(
    age_mean_f = seq(1,n_years),
    tsur = tsur
  )

  sim_data <- generate_sim_data(
    sim_data = sim_data,
    foi = foi_sim,
    sample_size_by_age = sample_size_by_age,
    mu = mu,
    model_type = "age-varying",
    survey_label = "age-varying-foi",
    seed = seed
  ) %>%
  prepare_serodata()

  prev_exact <- map_dbl(
    1:n_years,
    ~probability_exact_age_varying(., foi_sim, mu)
    )

  expect_s3_class(sim_data, "data.frame")
  expect_length(sim_data$birth_year, tsur - birth_year_min)
  expect_true(
    all(
      dplyr::near(
        sim_data$prev_obs,
        prev_exact,
        tol = tolerance
      )
    )
  )
})
