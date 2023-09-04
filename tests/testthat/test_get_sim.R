test_that("test get_sim functionalities", {

  library(serofoi)

  tsur <- 2050
  birth_year_min <- 2000
  survey_label <- "foi_sim"
  country = 'None'
  test = "test"
  antibody = "IgG"

  sim_data <- data.frame(birth_year = c(birth_year_min:(tsur - 1))) %>%
    mutate(tsur = tsur,
           country = country,
           test = test,
           antibody = antibody,
           survey = survey_label,
           age_mean_f = tsur - birth_year)

  #----- Test function generate_sim_probability
  foi_sim <- rep(0.02, 50)
  sim_probability <- get_sim_probability(sim_data = sim_data,
                                         foi = foi_sim)
  expect_s3_class(sim_probability, "data.frame")
  expect_type(sim_probability$age, "integer")
  expect_type(sim_probability$probability, "double")

  #----- Test function generate_sim_n_seropositivity
  sample_size_by_age <- 5
  sim_n_seropositive <- get_sim_n_seropositive(sim_data = sim_data,
                                               foi = foi_sim,
                                               sample_size_by_age = sample_size_by_age)
  expect_s3_class(sim_n_seropositive, "data.frame")
  expect_type(sim_n_seropositive$age, "integer")
  expect_type(sim_n_seropositive$n_seropositive, "integer")
})
