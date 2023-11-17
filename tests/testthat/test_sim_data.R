test_that("simulated data", {
  library(dplyr)
  library(serofoi)

  seed <- 1234
  sample_size_by_age <- 10^7
  tsur <- 2050
  birth_year_min <- 2000

  #----- Test for constant FoI

  foi_values <- c(0.001, 0.01, 0.1, 0.3, 0.4)
  for (foi_value in foi_values) {
    foi_sim <- rep(foi_value, tsur - birth_year_min)
    #----- Test function generate_sim_data
    sim_data <- generate_sim_data(
      foi = foi_sim,
      sample_size_by_age = sample_size_by_age,
      tsur = tsur,
      birth_year_min = birth_year_min,
      survey_label = "foi_sim_constant",
      seed = seed
    )
    prev_exact <- 1 - exp(-foi_value * sim_data$age_mean_f)

    expect_s3_class(sim_data, "data.frame")
    expect_length(sim_data$birth_year, tsur - birth_year_min)
    expect_equal(sim_data$prev_obs, prev_exact, tolerance = TRUE)

    #----- Test function group_sim_data
    sim_data <- sim_data %>% mutate(age_min = age_mean_f, age_max = age_mean_f)
    sim_data_grouped <- group_sim_data(sim_data = sim_data)
    expect_s3_class(sim_data_grouped, "data.frame")
    expect_s3_class(sim_data_grouped$age_group, "factor")
  }

  #----- Test for time-varying FoI
  no_transm <- 0.0000000001
  foi_sim <- c(rep(0.2, 25), rep(0.1, 10), rep(no_transm, 15))

  sim_data <- generate_sim_data(
    foi = foi_sim,
    sample_size_by_age = sample_size_by_age,
    tsur = tsur,
    birth_year_min = birth_year_min,
    survey_label = "sw_dec_foi",
    seed = seed
  )

  prev_exact <- 1 - exp(-cumsum(rev(foi_sim)))

  expect_s3_class(sim_data, "data.frame")
  expect_length(sim_data$birth_year, tsur - birth_year_min)
  expect_equal(sim_data$prev_obs, prev_exact, tolerance = TRUE)
})
