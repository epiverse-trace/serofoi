
test_that("probability_exact_age_varying calculates probabilities correctly", {
  # Test with simple input
  ages <- c(1, 2, 3)
  foi <- 0.1
  fois <- rep(foi, length(ages))
  probabilities <- serofoi:::probability_exact_age_varying(ages, fois)

  exact_probability_constant <- function(age, foi) {
    1 - exp(-age * foi)
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant(., foi))
  expect_equal(probabilities, expected, tolerance = 1e-6) # TODO change to dplyr::near

  # Test if FOIs increase that this leads to increased seropositivity
  fois_delta <- runif(length(ages))
  fois_h <- fois + fois_delta
  probabilities_h <- serofoi:::probability_exact_age_varying(ages, fois_h)
  expect_true(all(probabilities_h > probabilities))

  # Test with seroreversion
  seroreversion_rate <- 0.05
  probabilities <- serofoi:::probability_exact_age_varying(ages, fois, seroreversion_rate)

  exact_probability_constant_seroreversion <- function(age, foi, seroreversion) {
    foi / (foi + seroreversion_rate) * (1 - exp(-(foi + seroreversion_rate) * age))
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant_seroreversion(., foi, seroreversion))

  expect_equal(probabilities, expected, tolerance = 1e-6)

  # Test if FOIs increase that this leads to increased seropositivity when seroreversion present
  probabilities_h <- serofoi:::probability_exact_age_varying(ages, fois_h, seroreversion_rate)
  expect_true(all(probabilities_h > probabilities))

  # Test with analytical solution for non-constant FOIs
  ages <- c(1, 2)
  fois <- c(0.1, 0.2)
  probabilities <- serofoi:::probability_exact_age_varying(ages, fois)
  expected <- c(1 - exp(-0.1), 1 - exp(-(0.1 + 0.2)))
  expect_true(
    all(
      dplyr::near(
        probabilities,
        expected,
        tol = 1e-6
      )
    )
  )
})

test_that("probability_exact_time_varying calculates probabilities correctly", {
  # Test with constant FOI
  years <- c(1, 2, 3)
  foi <- 0.1
  fois <- rep(foi, length(years))
  probabilities <- serofoi:::probability_exact_time_varying(years, fois)

  exact_probability_constant <- function(age, foi) {
    1 - exp(-age * foi)
  }
  ages <- seq_along(years)
  expected <- purrr::map_dbl(ages, ~exact_probability_constant(., foi))
  expect_true(
    all(
      dplyr::near(
        probabilities,
        expected,
        tol = 1e-6
      )
    )
  )

  # Test with analytical solution
  years <- c(1, 2)
  fois <- c(0.1, 0.2)
  probabilities <- serofoi:::probability_exact_time_varying(years, fois)
  expected <- c(1 - exp(-0.2), 1 - exp(-(0.1 + 0.2)))
  expect_true(
    all(
      dplyr::near(
        probabilities,
        expected,
        tol = 1e-6
      )
    )
  )

  # Test that time-varying model gives a different answer to age-varying
  ages <- seq_along(years)
  probabilities_age <- serofoi:::probability_exact_age_varying(ages, fois)
  expect_true(
    probabilities_age[1] != probabilities[1] # for youngest age group these differ
  )

})

test_that("probability_seropositive_time_model_by_age works", {

  foi <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.2, 0.01))

  seroreversion <- 0.0
  prob_df <- serofoi:::probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df), nrow(foi))
  ages <- seq(1, nrow(foi), 1)
  expect_equal(ages, prob_df$age)

  # checking monotonicity
  derivative_foi <- diff(prob_df$seropositivity)
  expect_true(all(derivative_foi > 0))

  seroreversion <- 0.1
  prob_df_1 <- serofoi:::probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df_1), nrow(foi))
  expect_equal(ages, prob_df_1$age)

  # check seropositivities always lower (due to seroreversion)
  expect_true(all(prob_df_1$seropositivity < prob_df$seropositivity))
})

test_that("probability_seropositive_age_model_by_age works", {

  foi <- data.frame(
    age=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.2, 0.01))

  seroreversion <- 0.0
  prob_df <- serofoi:::probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df), nrow(foi))
  ages <- seq(1, nrow(foi), 1)
  expect_equal(ages, prob_df$age)

  # checking monotonicity
  derivative_foi <- diff(prob_df$seropositivity)
  expect_true(all(derivative_foi > 0))

  seroreversion <- 0.1
  prob_df_1 <- serofoi:::probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df_1), nrow(foi))
  expect_equal(ages, prob_df_1$age)

  # check seropositivities always lower (due to seroreversion)
  expect_true(all(prob_df_1$seropositivity < prob_df$seropositivity))
})

test_that("probability_seropositive_age_and_time_model_by_age works as expected", {
  us <- c(0.1, 0.2, 0.3)
  vs <- c(1, 0.5, 0.2)
  foi <- tidyr::expand_grid(
    u=us,
    v=vs
  ) %>%
    mutate(foi=u * v) %>%
    pull(foi)

  foi_df <- tidyr::expand_grid(
    year=c(1990, 1991, 1992),
    age=c(1, 2, 3)
  ) %>%
    mutate(foi=foi) %>%
    arrange(year)

  prob_df <- serofoi:::probability_seropositive_age_and_time_model_by_age(
    foi = foi_df,
    seroreversion_rate = 0
  )

  foi_matrix <- foi_df %>%
    tidyr::pivot_wider(
      values_from = foi,
      names_from = c(year)) %>%
    tibble::column_to_rownames("age") %>%
    as.matrix()
  serop_age_1 <- 1 - exp(-foi_matrix[1, 3])
  serop_age_2 <- 1 - exp(-(foi_matrix[1, 2] + foi_matrix[2, 3]))
  serop_age_3 <- 1 - exp(-(foi_matrix[1, 1] + foi_matrix[2, 2] + foi_matrix[3, 3]))
  expected <- c(serop_age_1, serop_age_2, serop_age_3)

  expect_true(
    all(
      dplyr::near(
        prob_df$seropositivity,
        expected,
        tol = 1e-6
      )
    )
  )

  # now add seroreversion
  mu <- 0.1
  prob_df_sr <- serofoi:::probability_seropositive_age_and_time_model_by_age(
    foi = foi_df,
    seroreversion_rate = mu
  )
  expect_true(all(prob_df_sr$seropositivity < prob_df$seropositivity))
  lambda <- foi_matrix[1, 3]
  serop_age_1 <- lambda / (lambda + mu) * (1 - exp(-(lambda + mu)))
  expect_true(
      dplyr::near(
        prob_df_sr$seropositivity[1],
        serop_age_1,
        tol = 1e-6
      )
  )
})

test_that("add_age_bins function works as expected", {
  # Test case 1: Check if intervals are created correctly for a single row dataframe
  survey_features <- data.frame(age_min = 20, age_max = 30)
  expected_intervals <- "[20,30]"
  actual_survey_features <- serofoi:::add_age_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)

  # Test case 2: Check if intervals are created correctly for multiple rows dataframe
  survey_features <- data.frame(age_min = c(20, 31), age_max = c(30, 50))
  expected_intervals <- c("[20,30]", "[31,50]")
  actual_survey_features <- serofoi:::add_age_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)
})

test_that("survey_by_individual_age function works as expected", {
  # Test case 1: Check if overall sample size is calculated correctly for a single row dataframe
  age_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]")
  survey_features <- data.frame(group = "[20,30]", n_sample = 100)
  expected_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]", overall_sample_size = 100)
  actual_df <- serofoi:::survey_by_individual_age(survey_features, age_df)
  expect_equal(actual_df, expected_df)

  # Test case 2: Check if overall sample size is calculated correctly for multiple rows dataframe
  age_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "[31,50]"))
  survey_features <- data.frame(group = c("[20,30]", "[31,50]"), n_sample = c(100, 150))
  expected_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "[31,50]"), overall_sample_size = c(100, 150))
  actual_df <- serofoi:::survey_by_individual_age(survey_features, age_df)
  expect_equal(actual_df, expected_df)
})

test_that("multinomial_sampling_group function works as expected", {
  # Test case 1: Check if sample sizes are generated correctly for a sample size of 100 and 5 age groups
  n_sample <- 100
  n_ages <- 5
  expected_length <- n_ages
  actual_sample_sizes <- serofoi:::multinomial_sampling_group(n_sample, n_ages)
  expect_length(actual_sample_sizes, expected_length)
  expect_equal(sum(actual_sample_sizes), n_sample)

  # Test case 2: Check if sample sizes are generated correctly for a sample size of 200 and 10 age groups
  n_sample <- 200
  n_ages <- 10
  expected_length <- n_ages
  actual_sample_sizes <- serofoi:::multinomial_sampling_group(n_sample, n_ages)
  expect_length(actual_sample_sizes, expected_length)
  expect_equal(sum(actual_sample_sizes), n_sample)
})

test_that("generate_random_sample_sizes function works as expected", {
  # Test case 1: Check if random sample sizes are generated correctly for a single interval
  survey_df <- data.frame(
    age=seq(20, 30, 1),
    group = "[20,30]",
    overall_sample_size = 100)
  actual_df <- serofoi:::generate_random_sample_sizes(survey_df)
  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      n_sample = sum(n_sample)
    )
  expect_equal(
    group_df$overall_sample_size[1],
    group_df$n_sample[1]
  )

  # Test case 2: Check if random sample sizes are generated correctly for two intervals
  survey_df <- data.frame(
    age=seq(20, 50, 1),
    group = c(rep("[20,30]", 11), rep("[31, 50)", 20)),
    overall_sample_size = c(rep(100, 11), rep(27, 20))
  )
  actual_df <- serofoi:::generate_random_sample_sizes(survey_df)
  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      n_sample = sum(n_sample)
    )
  expect_equal(group_df$n_sample, group_df$overall_sample_size)
})

test_that("sample_size_by_individual_age_random returns correct dataframe structure", {

  # Test with sample survey_features data: contiguous age bins
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )
  actual_df <- serofoi:::sample_size_by_individual_age_random(survey_features)
  expect_equal(nrow(actual_df), max(survey_features$age_max))

  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      n_sample = sum(n_sample)
    )
  expect_equal(group_df$n_sample, group_df$overall_sample_size)

  # Test with sample survey_features data: non-contiguous age bins
  # TODO: doesn't work as age_bins construction too simple currently.
  # It may just be that cut won't work reliably here.
  survey_features <- data.frame(
    age_min = c(1, 7, 18),
    age_max = c(2, 16, 20),
    n_sample = c(1000, 2000, 1500)
  )
  actual_df <- serofoi:::sample_size_by_individual_age_random(survey_features)
  expect_equal(nrow(actual_df), 15)
})

test_that("simulate_serosurvey_time_model function works as expected", {
  # Test case 1: Check if the output dataframe has the correct structure
  n_samples <- c(1000, 2000, 1500)
  foi_df <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.1, 0.01))
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = n_samples)
  actual_df <- simulate_serosurvey_time_model(foi_df, survey_features)
  expect_true("age_min" %in% colnames(actual_df))
  expect_true("age_max" %in% colnames(actual_df))
  expect_true("n_sample" %in% colnames(actual_df))
  expect_true("n_seropositive" %in% colnames(actual_df))

  # Test case 2: Check if the output dataframe has the correct number of rows
  expected_rows <- nrow(survey_features)
  actual_rows <- nrow(actual_df)
  expect_equal(actual_rows, expected_rows)

  # Test case 3: try a much higher FOI which should result in a higher proportion seropositive
  foi_df_1 <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rep(10, 20))
  actual_df_1 <- simulate_serosurvey_time_model(foi_df_1, survey_features)
  expect_true(all(actual_df_1$n_seropositive >= actual_df$n_seropositive))

  # Test case 4: allow a high rate of seroreversion which should reduce the proportion seropositive
  actual_df_2 <- simulate_serosurvey_time_model(
    foi=foi_df,
    survey_features=survey_features,
    seroreversion_rate=10
    )
  expect_true(all(actual_df_2$n_seropositive <= actual_df$n_seropositive))
})

test_that("simulate_serosurvey_time_model input validation", {

  foi_df <- data.frame(
    year = seq(1990, 2009, 1),
    foi = rnorm(20, 0.1, 0.01)
  )

  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )

  # Test with valid inputs
  expect_silent(simulate_serosurvey_time_model(foi_df, survey_features))

  # Test with non-dataframe foi dataframe
  expect_error(simulate_serosurvey_time_model(list(), survey_features),
               "foi must be a dataframe with columns foi and year.")

  # Test with non-dataframe survey_features dataframe
  expect_error(simulate_serosurvey_time_model(foi_df, list()),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with misspelt columns in foi dataframe
  expect_error(simulate_serosurvey_time_model(data.frame(years = c(1990), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi and year.")

  # Test with too many columns in foi dataframe
  expect_error(simulate_serosurvey_time_model(data.frame(age = c(1), year = c(2), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi and year.")

  # Test with missing columns in survey_features dataframe
  expect_error(simulate_serosurvey_time_model(foi_df, data.frame(age_min = c(1))),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with non-numeric seroreversion_rate
  expect_error(simulate_serosurvey_time_model(foi_df, survey_features, "seroreversion"),
               "seroreversion_rate must be a non-negative numeric value.")

  # Test with negative seroreversion_rate
  expect_error(simulate_serosurvey_time_model(foi_df, survey_features, -1),
               "seroreversion_rate must be a non-negative numeric value.")
})

test_that("simulate_serosurvey_age_model function works as expected", {
  # Test case 1: Check if the output dataframe has the correct structure
  n_samples <- c(1000, 2000, 1500)
  foi_df <- data.frame(
    age=seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.1, 0.01))
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = n_samples)
  actual_df <- simulate_serosurvey_age_model(foi_df, survey_features)
  expect_true("age_min" %in% colnames(actual_df))
  expect_true("age_max" %in% colnames(actual_df))
  expect_true("n_sample" %in% colnames(actual_df))
  expect_true("n_seropositive" %in% colnames(actual_df))

  # Test case 2: Check if the output dataframe has the correct number of rows
  expected_rows <- nrow(survey_features)
  actual_rows <- nrow(actual_df)
  expect_equal(actual_rows, expected_rows)

  # Test case 3: try a much higher FOI which should result in a higher proportion seropositive
  foi_df_1 <- data.frame(
    age=seq(1, 20, 1)
  ) %>%
    mutate(foi=rep(10, 20))
  actual_df_1 <- simulate_serosurvey_age_model(foi_df_1, survey_features)
  expect_true(all(actual_df_1$n_seropositive >= actual_df$n_seropositive))

  # Test case 4: allow a high rate of seroreversion which should reduce the proportion seropositive
  actual_df_2 <- simulate_serosurvey_age_model(
    foi=foi_df,
    survey_features=survey_features,
    seroreversion_rate=10
  )
  expect_true(all(actual_df_2$n_seropositive <= actual_df$n_seropositive))
})

test_that("simulate_serosurvey_age_model input validation", {

  foi_df <- data.frame(
    age = seq(1, 20, 1),
    foi = rnorm(20, 0.1, 0.01)
  )

  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )

  # Test with valid inputs
  expect_silent(simulate_serosurvey_age_model(foi_df, survey_features))

  # Test with non-dataframe foi dataframe
  expect_error(simulate_serosurvey_age_model(list(), survey_features),
               "foi must be a dataframe with columns foi and age.")

  # Test with non-dataframe survey_features dataframe
  expect_error(simulate_serosurvey_age_model(foi_df, list()),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with misspelt columns in foi dataframe
  expect_error(simulate_serosurvey_age_model(data.frame(ages = c(1), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi and age.")

  # Test with too many columns in foi dataframe
  expect_error(simulate_serosurvey_age_model(data.frame(age = c(1), year = c(2), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi and age.")

  # Test with missing columns in survey_features dataframe
  expect_error(simulate_serosurvey_age_model(foi_df, data.frame(age_min = c(1))),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with non-numeric seroreversion_rate
  expect_error(simulate_serosurvey_age_model(foi_df, survey_features, "seroreversion"),
               "seroreversion_rate must be a non-negative numeric value.")

  # Test with negative seroreversion_rate
  expect_error(simulate_serosurvey_age_model(foi_df, survey_features, -1),
               "seroreversion_rate must be a non-negative numeric value.")
})

test_that("simulate_serosurvey_age_and_time_model function works as expected", {
  # Test case 1: Check if the output dataframe has the correct structure
  n_samples <- c(1000, 2000, 1500)
  foi_df <- tidyr::expand_grid(
    year = seq(1990, 2009, 1),
    age = seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20 * 20, 0.1, 0.001))
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = n_samples)
  actual_df <- simulate_serosurvey_age_and_time_model(foi_df, survey_features)
  expect_true("age_min" %in% colnames(actual_df))
  expect_true("age_max" %in% colnames(actual_df))
  expect_true("n_sample" %in% colnames(actual_df))
  expect_true("n_seropositive" %in% colnames(actual_df))

  # Test case 2: Check if the output dataframe has the correct number of rows
  expected_rows <- nrow(survey_features)
  actual_rows <- nrow(actual_df)
  expect_equal(actual_rows, expected_rows)

  # Test case 3: try a much higher FOI which should result in a higher proportion seropositive
  foi_df_1 <- tidyr::expand_grid(
    year = seq(1990, 2009, 1),
    age = seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20 * 20, 10.1, 0.001))
  actual_df_1 <- simulate_serosurvey_age_and_time_model(foi_df_1, survey_features)
  expect_true(all(actual_df_1$n_seropositive >= actual_df$n_seropositive))

  # Test case 4: allow a high rate of seroreversion which should reduce the proportion seropositive
  actual_df_2 <- simulate_serosurvey_age_and_time_model(
    foi=foi_df,
    survey_features=survey_features,
    seroreversion_rate=10
  )
  expect_true(all(actual_df_2$n_seropositive <= actual_df$n_seropositive))
})

test_that("simulate_serosurvey_age_and_time_model input validation", {

  foi_df <- tidyr::expand_grid(
    year = seq(1990, 2009, 1),
    age = seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20 * 20, 0.1, 0.001))

  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )

  # Test with valid inputs
  expect_silent(simulate_serosurvey_age_and_time_model(foi_df, survey_features))

  # Test with non-dataframe foi dataframe
  expect_error(simulate_serosurvey_age_and_time_model(list(), survey_features),
               "foi must be a dataframe with columns foi, age and year.")

  # Test with non-dataframe survey_features dataframe
  expect_error(simulate_serosurvey_age_and_time_model(foi_df, list()),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with misspelt columns in foi dataframe
  expect_error(simulate_serosurvey_age_and_time_model(data.frame(ages = c(1), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi, age and year.")

  # Test with missing columns in foi dataframe
  expect_error(simulate_serosurvey_age_and_time_model(data.frame(age = c(1), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi, age and year.")
  expect_error(simulate_serosurvey_age_and_time_model(data.frame(year = c(1), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi, age and year.")

  # Test with too many columns in foi dataframe
  expect_error(simulate_serosurvey_time_model(data.frame(age = c(1), year = c(2), sex = c(3), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns foi and year.")

  # Test with missing columns in survey_features dataframe
  expect_error(simulate_serosurvey_age_and_time_model(foi_df, data.frame(age_min = c(1))),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'.")

  # Test with non-numeric seroreversion_rate
  expect_error(simulate_serosurvey_age_and_time_model(foi_df, survey_features, "seroreversion"),
               "seroreversion_rate must be a non-negative numeric value.")

  # Test with negative seroreversion_rate
  expect_error(simulate_serosurvey_age_and_time_model(foi_df, survey_features, -1),
               "seroreversion_rate must be a non-negative numeric value.")
})

test_that("simulate_serosurvey returns serosurvey data based on specified model", {
  # Test with 'age' model
  foi_df <- data.frame(
    age = seq(1, 20, 1),
    foi = runif(20, 0.05, 0.15)
  )
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )
  serosurvey <- simulate_serosurvey("age", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "n_sample", "n_seropositive")))

  # Test with 'time' model
  foi_df <- data.frame(
    year = seq(1990, 2009, 1),
    foi = runif(20, 0.05, 0.15)
  )
  serosurvey <- simulate_serosurvey("time", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "n_sample", "n_seropositive")))

  # Test with 'age-time' model
  foi_df <- tidyr::expand_grid(
    year = seq(1990, 2009, 1),
    age = seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20 * 20, 0.1, 0.001))

  serosurvey <- simulate_serosurvey("age-time", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "n_sample", "n_seropositive")))
})

test_that("simulate_serosurvey handles invalid model inputs", {
  # Test with invalid model
  foi_df <- data.frame(
    age = seq(1, 20, 1),
    foi = runif(20, 0.05, 0.15)
  )
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    n_sample = c(1000, 2000, 1500)
  )
  expect_error(simulate_serosurvey("invalid_model", foi_df, survey_features),
               "model must be one of 'age', 'time', or 'age-time'.")
})

test_that("probability_seropositive_general_model_by_age reduces to time-varying model under
          appropriate limits", {

    # simple time-varying FOI model
    construct_A <- function(t, tau, lambda) {
      A <- matrix(0, ncol = 2, nrow = 2)
      A[1, 1] <- -lambda[t]
      A[2, 1] <- lambda[t]
      A
    }

    # determines the sum of seropositive compartments of those still alive
    calculate_seropositivity_function <- function(Y) {
      Y[2]
    }

    # initial conditions in 12D state vector
    initial_conditions <- rep(0, 2)
    initial_conditions[1] <- 1

    # random FOIs
    lambda <- runif(70, 0, 0.01)

    # solve linear system of ODEs
    seropositive_linear_system <- probability_seropositive_general_model_by_age(
      construct_A,
      calculate_seropositivity_function,
      initial_conditions,
      max_age=length(lambda),
      lambda
    )

    foi_df <- data.frame(
      year=seq_along(lambda),
      foi=lambda
    )

    seropositive_true <- probability_seropositive_by_age(
      model = "time",
      foi = foi_df,
      seroreversion_rate = 0
    )

    expect_equal(seropositive_true, seropositive_linear_system)

})

test_that("probability_seropositive_general_model_by_age reduces to age-varying model under
          appropriate limits", {

  # simple age-varying FOI model
  construct_A <- function(t, tau, lambda) {
    A <- matrix(0, ncol = 2, nrow = 2)
    A[1, 1] <- -lambda[t - tau]
    A[2, 1] <- lambda[t - tau]
    A
  }

  # determines the sum of seropositive compartments of those still alive
  calculate_seropositivity_function <- function(Y) {
    Y[2]
  }

  # initial conditions in 12D state vector
  initial_conditions <- rep(0, 2)
  initial_conditions[1] <- 1

  # random FOIs
  lambda <- runif(70, 0, 0.01)

  # solve linear system of ODEs
  seropositive_linear_system <- probability_seropositive_general_model_by_age(
    construct_A,
    calculate_seropositivity_function,
    initial_conditions,
    max_age=length(lambda),
    lambda
  )

  foi_df <- data.frame(
    age=seq_along(lambda),
    foi=lambda
  )

  seropositive_true <- probability_seropositive_by_age(
    model = "age",
    foi = foi_df,
    seroreversion_rate = 0
  )

  expect_equal(seropositive_true, seropositive_linear_system)

})

test_that("probability_seropositive_general_model_by_age reduces to age- and time-varying model under
          appropriate limits", {

  # age- and time-varying FOI model
  construct_A <- function(t, tau, u, v) {
    A <- matrix(0, ncol = 2, nrow = 2)
    u_bar <- u[t - tau]
    v_bar <- v[t]

    A[1, 1] <- -u_bar * v_bar
    A[2, 1] <- u_bar * v_bar
    A
  }

  # determines the sum of seropositive compartments of those still alive
  calculate_seropositivity_function <- function(Y) {
    Y[2]
  }

  # initial conditions in 12D state vector
  initial_conditions <- rep(0, 2)
  initial_conditions[1] <- 1

  # age and time-varying FOIs
  ages <- seq(1, 70, 1)
  foi_age <- 2 * dlnorm(
    ages, meanlog = 3.5, sdlog = 0.5)

  foi_df_age <- data.frame(
    age = ages,
    foi = foi_age
  )

  foi_time <- c(rep(0, 30), rep(1, 40))
  foi_df_time <- data.frame(
    year = seq(1956, 2025, 1),
    foi = foi_time
  )

  u <- foi_df_age$foi
  v <- foi_df_time$foi

  # solve linear system of ODEs
  seropositive_linear_system <- probability_seropositive_general_model_by_age(
    construct_A,
    calculate_seropositivity_function,
    initial_conditions,
    max_age=length(lambda),
    u,
    v
  )

  foi_df <- expand.grid(
    year=foi_df_time$year,
    age=foi_df_age$age
  ) %>%
    left_join(foi_df_age, by="age") %>%
    rename(foi_age=foi) %>%
    left_join(foi_df_time, by="year") %>%
    rename(foi_time=foi) %>%
    mutate(foi = foi_age * foi_time) %>%
    select(-c("foi_age", "foi_time")) %>%
    mutate(birth_year = year - age) %>%
    filter(birth_year >= 1955) %>%
    arrange(birth_year, age)

  seropositive_true <- probability_seropositive_by_age(
    model = "age-time",
    foi = foi_df %>% select(-birth_year),
    seroreversion_rate = 0
  )

  expect_equal(seropositive_true, seropositive_linear_system)

})
