
test_that("probability_exact_time_varying calculates probabilities correctly", {
  # Test with simple input
  ages <- c(1, 2, 3)
  foi <- 0.1
  fois <- rep(foi, length(ages))
  probabilities <- probability_exact_time_varying(ages, fois)

  exact_probability_constant <- function(age, foi) {
    1 - exp(-age * foi)
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant(., foi))
  expect_equal(probabilities, expected, tolerance = 1e-6)

  # Test if FOIs increase that this leads to increased seropositivity
  fois_delta <- runif(length(ages))
  fois_h <- fois + fois_delta
  probabilities_h <- probability_exact_time_varying(ages, fois_h)
  expect_true(all(probabilities_h > probabilities))

  # Test with seroreversion
  seroreversion_rate <- 0.05
  probabilities <- probability_exact_time_varying(ages, fois, seroreversion_rate)

  exact_probability_constant_seroreversion <- function(age, foi, seroreversion) {
    foi / (foi + seroreversion_rate) * (1 - exp(-(foi + seroreversion_rate) * age))
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant_seroreversion(., foi, seroreversion))

  expect_equal(probabilities, expected, tolerance = 1e-6)

  # Test if FOIs increase that this leads to increased seropositivity when seroreversion present
  probabilities_h <- probability_exact_time_varying(ages, fois_h, seroreversion_rate)
  expect_true(all(probabilities_h > probabilities))
})

test_that("probability_exact_age_varying calculates probabilities correctly", {
  # Test with simple input
  ages <- c(1, 2, 3)
  foi <- 0.1
  fois <- rep(foi, length(ages))
  probabilities <- probability_exact_age_varying(ages, fois)

  exact_probability_constant <- function(age, foi) {
    1 - exp(-age * foi)
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant(., foi))
  expect_equal(probabilities, expected, tolerance = 1e-6)

  # Test with seroreversion
  seroreversion_rate <- 0.05
  probabilities <- probability_exact_age_varying(ages, fois, seroreversion_rate)

  exact_probability_constant_seroreversion <- function(age, foi, seroreversion) {
    foi / (foi + seroreversion_rate) * (1 - exp(-(foi + seroreversion_rate) * age))
  }
  expected <- purrr::map_dbl(ages, ~exact_probability_constant_seroreversion(., foi, seroreversion))

  expect_equal(probabilities, expected, tolerance = 1e-6)
})

test_that("probability_seropositive_time_model_by_age works", {

  foi <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.2, 0.01))

  seroreversion <- 0.0
  prob_df <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df), nrow(foi))
  ages <- seq(1, nrow(foi), 1)
  expect_equal(ages, prob_df$age)

  # checking monotonicity
  derivative_foi <- diff(prob_df$seropositivity)
  expect_true(all(derivative_foi > 0))

  seroreversion <- 0.1
  prob_df_1 <- probability_seropositive_time_model_by_age(
    foi = foi,
    seroreversion = seroreversion
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
  prob_df <- probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df), nrow(foi))
  ages <- seq(1, nrow(foi), 1)
  expect_equal(ages, prob_df$age)

  # checking monotonicity
  derivative_foi <- diff(prob_df$seropositivity)
  expect_true(all(derivative_foi > 0))

  seroreversion <- 0.1
  prob_df_1 <- probability_seropositive_age_model_by_age(
    foi = foi,
    seroreversion_rate = seroreversion
  )

  # check output dimensions
  expect_equal(nrow(prob_df_1), nrow(foi))
  expect_equal(ages, prob_df_1$age)

  # check seropositivities always lower (due to seroreversion)
  expect_true(all(prob_df_1$seropositivity < prob_df$seropositivity))
})

test_that("create_group_interval function works as expected", {
  # Test case 1: Check if interval is created correctly for is_first_row = TRUE
  age_min <- 20
  age_max <- 30
  is_first_row <- TRUE
  expected_interval <- "[20,30]"
  actual_interval <- create_group_interval(age_min, age_max, is_first_row)
  expect_equal(actual_interval, expected_interval)

  # Test case 2: Check if interval is created correctly for is_first_row = FALSE
  age_min <- 20
  age_max <- 30
  is_first_row <- FALSE
  expected_interval <- "(19,30]"
  actual_interval <- create_group_interval(age_min, age_max, is_first_row)
  expect_equal(actual_interval, expected_interval)
})

test_that("add_age_bins function works as expected", {
  # Test case 1: Check if intervals are created correctly for a single row dataframe
  survey_features <- data.frame(age_min = 20, age_max = 30)
  expected_intervals <- "[20,30]"
  actual_survey_features <- add_age_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)

  # Test case 2: Check if intervals are created correctly for multiple rows dataframe
  survey_features <- data.frame(age_min = c(20, 31), age_max = c(30, 50))
  expected_intervals <- c("[20,30]", "(30,50]")
  actual_survey_features <- add_age_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)
})

test_that("survey_by_individual_age function works as expected", {
  # Test case 1: Check if overall sample size is calculated correctly for a single row dataframe
  age_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]")
  survey_features <- data.frame(group = "[20,30]", sample_size = 100)
  expected_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]", overall_sample_size = 100)
  actual_df <- survey_by_individual_age(survey_features, age_df)
  expect_equal(actual_df, expected_df)

  # Test case 2: Check if overall sample size is calculated correctly for multiple rows dataframe
  age_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "(30,50]"))
  survey_features <- data.frame(group = c("[20,30]", "(30,50]"), sample_size = c(100, 150))
  expected_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "(30,50]"), overall_sample_size = c(100, 150))
  actual_df <- survey_by_individual_age(survey_features, age_df)
  expect_equal(actual_df, expected_df)
})

test_that("multinomial_sampling_group function works as expected", {
  # Test case 1: Check if sample sizes are generated correctly for a sample size of 100 and 5 age groups
  sample_size <- 100
  n_ages <- 5
  expected_length <- n_ages
  actual_sample_sizes <- multinomial_sampling_group(sample_size, n_ages)
  expect_length(actual_sample_sizes, expected_length)
  expect_equal(sum(actual_sample_sizes), sample_size)

  # Test case 2: Check if sample sizes are generated correctly for a sample size of 200 and 10 age groups
  sample_size <- 200
  n_ages <- 10
  expected_length <- n_ages
  actual_sample_sizes <- multinomial_sampling_group(sample_size, n_ages)
  expect_length(actual_sample_sizes, expected_length)
  expect_equal(sum(actual_sample_sizes), sample_size)
})

test_that("generate_random_sample_sizes function works as expected", {
  # Test case 1: Check if random sample sizes are generated correctly for a single interval
  survey_df <- data.frame(
    age=seq(20, 30, 1),
    group = "[20,30]",
    overall_sample_size = 100)
  actual_df <- generate_random_sample_sizes(survey_df)
  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      sample_size = sum(sample_size)
    )
  expect_equal(
    group_df$overall_sample_size[1],
    group_df$sample_size[1]
  )

  # Test case 2: Check if random sample sizes are generated correctly for two intervals
  survey_df <- data.frame(
    age=seq(20, 50, 1),
    group = c(rep("[20,30]", 11), rep("(30, 50)", 20)),
    overall_sample_size = c(rep(100, 11), rep(27, 20))
  )
  actual_df <- generate_random_sample_sizes(survey_df)
  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      sample_size = sum(sample_size)
    )
  expect_equal(group_df$sample_size, group_df$overall_sample_size)
})

test_that("sample_size_by_individual_age_random returns correct dataframe structure", {

  # Test with sample survey_features data: contiguous age bins
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    sample_size = c(1000, 2000, 1500)
  )
  actual_df <- sample_size_by_individual_age_random(survey_features)
  expect_equal(nrow(actual_df), max(survey_features$age_max))

  group_df <- actual_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      overall_sample_size = overall_sample_size[1],
      sample_size = sum(sample_size)
    )
  expect_equal(group_df$sample_size, group_df$overall_sample_size)

  # Test with sample survey_features data: non-contiguous age bins
  # TODO: doesn't work as age_bins construction too simple currently.
  # It may just be that cut won't work reliably here.
  survey_features <- data.frame(
    age_min = c(1, 7, 18),
    age_max = c(2, 16, 20),
    sample_size = c(1000, 2000, 1500)
  )
  actual_df <- sample_size_by_individual_age_random(survey_features)
  expect_equal(nrow(actual_df), max(survey_features$age_max))
})


test_that("simulate_serosurvey_time_model function works as expected", {
  # Test case 1: Check if the output dataframe has the correct structure
  sample_sizes <- c(1000, 2000, 1500)
  foi_df <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.1, 0.01))
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    sample_size = sample_sizes)
  actual_df <- simulate_serosurvey_time_model(foi_df, survey_features)
  expect_true("age_min" %in% colnames(actual_df))
  expect_true("age_max" %in% colnames(actual_df))
  expect_true("sample_size" %in% colnames(actual_df))
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
    sample_size = c(1000, 2000, 1500)
  )

  # Test with valid inputs
  expect_silent(simulate_serosurvey_time_model(foi_df, survey_features))

  # Test with non-dataframe foi dataframe
  expect_error(simulate_serosurvey_time_model(list(), survey_features),
               "foi must be a dataframe with columns 'year' and 'foi'.")

  # Test with non-dataframe survey_features dataframe
  expect_error(simulate_serosurvey_time_model(foi_df, list()),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")

  # Test with misspelt columns in foi dataframe
  expect_error(simulate_serosurvey_time_model(data.frame(years = c(1990), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns 'year' and 'foi'.")

  # Test with missing columns in survey_features dataframe
  expect_error(simulate_serosurvey_time_model(foi_df, data.frame(age_min = c(1))),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")

  # Test with non-numeric seroreversion_rate
  expect_error(simulate_serosurvey_time_model(foi_df, survey_features, "seroreversion"),
               "seroreversion_rate must be a non-negative numeric value.")

  # Test with negative seroreversion_rate
  expect_error(simulate_serosurvey_time_model(foi_df, survey_features, -1),
               "seroreversion_rate must be a non-negative numeric value.")
})


test_that("simulate_serosurvey_age_model function works as expected", {
  # Test case 1: Check if the output dataframe has the correct structure
  sample_sizes <- c(1000, 2000, 1500)
  foi_df <- data.frame(
    age=seq(1, 20, 1)
  ) %>%
    mutate(foi=rnorm(20, 0.1, 0.01))
  survey_features <- data.frame(
    age_min = c(1, 3, 15),
    age_max = c(2, 14, 20),
    sample_size = sample_sizes)
  actual_df <- simulate_serosurvey_age_model(foi_df, survey_features)
  expect_true("age_min" %in% colnames(actual_df))
  expect_true("age_max" %in% colnames(actual_df))
  expect_true("sample_size" %in% colnames(actual_df))
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
    sample_size = c(1000, 2000, 1500)
  )

  # Test with valid inputs
  expect_silent(simulate_serosurvey_age_model(foi_df, survey_features))

  # Test with non-dataframe foi dataframe
  expect_error(simulate_serosurvey_age_model(list(), survey_features),
               "foi must be a dataframe with columns 'age' and 'foi'.")

  # Test with non-dataframe survey_features dataframe
  expect_error(simulate_serosurvey_age_model(foi_df, list()),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")

  # Test with misspelt columns in foi dataframe
  expect_error(simulate_serosurvey_age_model(data.frame(ages = c(1), foi = c(0.1)), survey_features),
               "foi must be a dataframe with columns 'age' and 'foi'.")

  # Test with missing columns in survey_features dataframe
  expect_error(simulate_serosurvey_age_model(foi_df, data.frame(age_min = c(1))),
               "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'sample_size'.")

  # Test with non-numeric seroreversion_rate
  expect_error(simulate_serosurvey_age_model(foi_df, survey_features, "seroreversion"),
               "seroreversion_rate must be a non-negative numeric value.")

  # Test with negative seroreversion_rate
  expect_error(simulate_serosurvey_age_model(foi_df, survey_features, -1),
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
    sample_size = c(1000, 2000, 1500)
  )
  serosurvey <- simulate_serosurvey("age", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "sample_size", "n_seropositive")))

  # Test with 'time' model
  foi_df <- data.frame(
    year = seq(1990, 2009, 1),
    foi = runif(20, 0.05, 0.15)
  )
  serosurvey <- simulate_serosurvey("time", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "sample_size", "n_seropositive")))

  # Test with 'age-time' model: TODO
  serosurvey <- simulate_serosurvey("age-time", foi_df, survey_features)
  expect_true(all(names(serosurvey) %in% c("age_min", "age_max", "sample_size", "n_seropositive")))
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
    sample_size = c(1000, 2000, 1500)
  )
  expect_error(simulate_serosurvey("invalid_model", foi_df, survey_features),
               "model must be one of 'age', 'time', or 'age-time'.")
})
