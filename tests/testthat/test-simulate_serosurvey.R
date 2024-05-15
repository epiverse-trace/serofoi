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

test_that("create_features_with_bins function works as expected", {
  # Test case 1: Check if intervals are created correctly for a single row dataframe
  survey_features <- data.frame(age_min = 20, age_max = 30)
  expected_intervals <- "[20,30]"
  actual_survey_features <- create_features_with_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)

  # Test case 2: Check if intervals are created correctly for multiple rows dataframe
  survey_features <- data.frame(age_min = c(20, 31), age_max = c(30, 50))
  expected_intervals <- c("[20,30]", "(30,50]")
  actual_survey_features <- create_features_with_bins(survey_features)
  actual_intervals <- actual_survey_features$group
  expect_equal(actual_intervals, expected_intervals)
})

test_that("overall_sample_size_by_group function works as expected", {
  # Test case 1: Check if overall sample size is calculated correctly for a single row dataframe
  age_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]")
  survey_features <- data.frame(group = "[20,30]", sample_size = 100)
  expected_df <- data.frame(age_min = 20, age_max = 30, group = "[20,30]", overall_sample_size = 100)
  actual_df <- overall_sample_size_by_group(survey_features, age_df)
  expect_equal(actual_df, expected_df)

  # Test case 2: Check if overall sample size is calculated correctly for multiple rows dataframe
  age_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "(30,50]"))
  survey_features <- data.frame(group = c("[20,30]", "(30,50]"), sample_size = c(100, 150))
  expected_df <- data.frame(age_min = c(20, 30), age_max = c(31, 50), group = c("[20,30]", "(30,50]"), overall_sample_size = c(100, 150))
  actual_df <- overall_sample_size_by_group(survey_features, age_df)
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
  foi_df <- data.frame(
    year=seq(1990, 2009, 1)
  ) %>%
    mutate(foi=rep(10, 20))
  actual_df_1 <- simulate_serosurvey_time_model(foi_df, survey_features)
  expect_true(all(actual_df_1$n_seropositive >= actual_df$n_seropositive))
})
