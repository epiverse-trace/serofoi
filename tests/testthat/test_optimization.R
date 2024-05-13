test_that("Test constant model optimization", {
  
  # Characteristics of data for testing
  tsur <- 2024
  n_years <- 80
  birth_years <- seq(tsur - n_years, tsur - 1)
  true_foi <- rep(0.05, n_years)
  
  # Create simulated data
  sim_data <- data.frame(
    age = seq(1:n_years),
    tsur = tsur
  )
  serodata <- generate_sim_data(
    sim_data = sim_data,
    foi = true_foi,
    sample_size_by_age = 200) %>% prepare_serodata()
  
  cohort_ages <- get_cohort_ages(serodata)
  
  # Fit model with optimization
  model_object <- fit_seromodel_optimization(
    serodata = serodata,
    foi_model = "constant",
    iter = 1000
  )
  
  # Expect that the optimized fit is within 0.01 of the true FOI that generated
  # the data
  expect_true(
    all(
      dplyr::near(
        model_object$par$foi,
        true_foi,
        tol = 0.01
      )
    )
  )
  
})


test_that("Test tv_normal model optimization", {
  
  # Characteristics of data for testing
  tsur <- 2024
  n_years <- 80
  birth_years <- seq(tsur - n_years, tsur - 1)
  true_foi <- rep(
    c(0.01, 0.03, 0.06, 0.03),
    c(20, 20, 20, 20)
  )
  
  # Create simulated data
  sim_data <- data.frame(
    age = seq(1:n_years),
    tsur = tsur
  )
  serodata <- generate_sim_data(
    sim_data = sim_data,
    foi = true_foi,
    sample_size_by_age = 200) %>% prepare_serodata()
  
  cohort_ages <- get_cohort_ages(serodata)
  
  # Fit model with optimization
  model_object <- fit_seromodel_optimization(
    serodata = serodata,
    foi_model = "tv_normal",
    iter = 1000,
    chunk_size = 20
  )
  
  # Expect that the optimized fit is within 0.01 of the true FOI that generated
  # the data
  expect_true(
    all(
      dplyr::near(
        model_object$par$foi,
        true_foi,
        tol = 0.01
      )
    )
  )
  
})


test_that("Test tv_normal_log model optimization", {
  
  # Characteristics of data for testing
  tsur <- 2024
  n_years <- 80
  birth_years <- seq(tsur - n_years, tsur - 1)
  true_foi <- rep(
    c(0.01, 0.03, 0.06, 0.0001),
    c(20, 20, 20, 20)
  )
  
  # Create simulated data
  sim_data <- data.frame(
    age = seq(1:n_years),
    tsur = tsur
  )
  serodata <- generate_sim_data(
    sim_data = sim_data,
    foi = true_foi,
    sample_size_by_age = 200) %>% prepare_serodata()
  
  cohort_ages <- get_cohort_ages(serodata)
  
  # Fit model with optimization
  model_object <- fit_seromodel_optimization(
    serodata = serodata,
    foi_model = "tv_normal_log",
    iter = 1000,
    chunk_size = 20
  )
  
  # Expect that the optimized fit is within 0.01 of the true FOI that generated
  # the data
  expect_true(
    all(
      dplyr::near(
        model_object$par$foi,
        true_foi,
        tol = 0.01
      )
    )
  )
  
})

