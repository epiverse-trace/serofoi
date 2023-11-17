test_that("test get_sim functionalities", {
  library(serofoi)
  library(dplyr)

  tsur <- 2050
  birth_year_min <- 2000
  survey_label <- "foi_sim"
  country <- "None"
  test <- "test"
  antibody <- "IgG"

  sim_data <- data.frame(birth_year = c(birth_year_min:(tsur - 1))) %>%
    mutate(
      tsur = tsur,
      country = country,
      test = test,
      antibody = antibody,
      survey = survey_label,
      age_mean_f = tsur - birth_year
    )

  #----- Test function generate_sim_probability
  n_years <- 50
  foi_sim <- rep(0.02, n_years)
  sim_probability <- get_sim_probability(
    sim_data = sim_data,
    foi = foi_sim
  )

  exposure_matrix <- matrix(1, n_years, n_years)
  exposure_matrix[lower.tri(exposure_matrix)] <- 0
  probabilities <- purrr::map_dbl(1:n_years, ~ 1 - exp(-pracma::dot(exposure_matrix[., ], foi_sim)))

  expect_s3_class(sim_probability, "data.frame")
  expect_type(sim_probability$age, "integer")
  expect_type(sim_probability$probability, "double")
  expect_true(all(probabilities == sim_probability$probability))

  #----- Test function generate_sim_n_seropositivity
  sample_size_by_age <- 5
  sim_n_seropositive <- get_sim_n_seropositive(
    sim_data = sim_data,
    foi = foi_sim,
    sample_size_by_age = sample_size_by_age
  )
  expect_s3_class(sim_n_seropositive, "data.frame")
  expect_type(sim_n_seropositive$age, "integer")
  expect_type(sim_n_seropositive$n_seropositive, "integer")
})
