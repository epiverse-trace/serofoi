
# Test for validate_serosurvey ----
test_that("validate_serosurvey throws an error for invalid input", {
  # Define an invalid serosurvey data frame with a missing column
  missing_column_serosurvey <- data.frame(
    age_min = c(1, 10, 20),
    age_max = c(9, 19, 30),
    n_sample = c(100, 100, 100)
  )

  # Test that function errors with a missing column
  expect_error(
    validate_serosurvey(missing_column_serosurvey),
    "must include"
  )

  # Define an invalid serosurvey data frame with incorrect column types
  incorrect_type_serosurvey <-  dplyr::mutate(
    missing_column_serosurvey,
    n_seropositive = c("10", "30", "70")
    )

  # Test that function errors with incorrect column types
  expect_error(
    validate_serosurvey(incorrect_type_serosurvey),
    "`n_seropositive` must be of any of these types: `numeric`"
  )
})

# Test for validate_survey_features ----
test_that("validate_survey_features stops if survey_features is not a dataframe", {
  # Create non-data frame survey features
  non_df_input <- list(
    age_min = c(1, 5),
    age_max = c(6, 10),
    n_sample = c(100, 200)
    )

  # Validate exception
  expect_error(
    validate_survey_features(non_df_input),
    "survey_features must be a dataframe"
  )
})

test_that("validate_survey_features stops if required columns are missing", {
  # Case where required columns are missing
  missing_columns_df <- data.frame(
    age_min = c(1, 10, 20),
    age_max = c(9, 19, 30)
    )
  expect_error(
    validate_survey_features(missing_columns_df),
    "survey_features must be a dataframe with columns 'age_min', 'age_max', and 'n_sample'."
  )
})

test_that("validate_survey_features stops if age bins have overlapping bounds", {
  # Case where age bins overlap (age_max of one row equals age_min of another)
  overlapping_age_df <- data.frame(
    age_min = c(1, 10, 20),
    age_max = c(10, 19, 30),
    n_sample = c(100, 200, 300)
  )

  expect_error(
    validate_survey_features(overlapping_age_df),
    "Age bins in a survey are inclusive of both bounds, so the age_max of one bin cannot equal the age_min of another."
  )
})

# Test validate_foi_index ----
test_that("validate_foi_index throws an error for non-consecutive indexes", {
  # Sample survey features for testing
  survey_features <- data.frame(
    age_min = c(1, 6, 11, 16, 21),
    age_max = c(5, 10, 15, 20, 25),
    survey_year = 2025
  )

  # Test validation works for invalid sizes
  ## shorter
  foi_index <- data.frame(
    age = 1:20,
    foi_index = c(rep(1, 10), rep(2, 10))
  )
  expect_error(
    serofoi:::validate_foi_index(foi_index, survey_features, model_type = "age")
  )
  ## longer
  foi_index <- data.frame(
    age = 1:30,
    foi_index = c(rep(1, 10), rep(2, 10), rep(3, 10))
  )
  expect_error(
    serofoi:::validate_foi_index(foi_index, survey_features, model_type = "age")
  )

  # Test validation works for missing indexes
  foi_index <- data.frame(
    age = 1:25,
    foi_index = c(rep(1, 10), rep(3, 15))
  )
  expect_error(
    serofoi:::validate_foi_index(foi_index, survey_features, model_type = "age")
  )

  # Test that validation works decreasing indexes
  foi_index <- data.frame(
    age = 1:25,
    foi_index = c(rep(1, 10), rep(2, 10), rep(1, 5))
  )
  expect_error(
    serofoi:::validate_foi_index(foi_index, survey_features, model_type = "age")
  )
})

# Test for validate_seroreversion_rate ----
test_that("validate_seroreversion_rate stops if seroreversion_rate is negative", {
  # Case where seroreversion_rate is negative
  expect_error(
    validate_seroreversion_rate(-0.5),
    "seroreversion_rate must be a non-negative numeric value."
  )
})

# Test for validate_simulation_age ----
test_that("validate_simulation_age throws error if max age in foi_df exceeds max age in survey_features", {
  # Sample survey features for testing
  survey_features <- data.frame(
    age_min = c(1, 6, 11, 16, 21),
    age_max = c(5, 10, 15, 20, 25),
    n_sample = c(100, 150, 50, 100, 75)
  )

  # Case where foi_df has more rows than the max age in survey_features
  foi_df <- data.frame(foi = rep(0.1, 30))
  expect_error(
    validate_simulation_age(survey_features, foi_df),
    "maximum age implicit in foi_df should not exceed max age in survey_features"
  )
})

# Test for validate_simulation_age_time ----
test_that("validate_simulation_age_time throws error if max age in foi_df exceeds max age in survey_features", {
  # Sample survey features for testing
  survey_features <- data.frame(
    age_min = c(1, 6, 11, 16, 21),
    age_max = c(5, 10, 15, 20, 25),
    n_sample = c(100, 150, 50, 100, 75)
  )

  # Case where implicit age in foi_df exceeds max age in survey_features
  foi_df <- expand.grid(
    year = seq(1980, 2009, 1),
    age = seq(1, 30, 1)
  )
  foi_df$foi <- rnorm(30 * 30, 0.1, 0.01)

  # Case where foi_df has more rows than the max age in survey_features
  expect_error(
    validate_simulation_age_time(survey_features, foi_df),
    "maximum age implicit in foi_df should not exceed max age in survey_features"
  )
})

# Test for validate_plot_constant ----
test_that("validate_plot_constant works as expected", {
  # Valid cases
  expect_true(
    validate_plot_constant(
      plot_constant = TRUE,
      x_axis = "age",
      model_name = "constant_model",
      error_msg_x_axis = "x_axis must be either 'age' or 'time'."
    )
  )

  expect_true(
    validate_plot_constant(
      plot_constant = FALSE,
      x_axis = "time",
      model_name = "time_varying_model",
      error_msg_x_axis = "x_axis must be either 'age' or 'time'."
    )
  )

  # Invalid cases
  expect_error(
    validate_plot_constant(
      plot_constant = TRUE,
      x_axis = "invalid_axis",
      model_name = "constant_model",
      error_msg_x_axis = "x_axis must be either 'age' or 'time'."
    ),
    "x_axis must be either 'age' or 'time'."
  )

  expect_error(
    validate_plot_constant(
      plot_constant = TRUE,
      x_axis = "age",
      model_name = "time_varying_model",
      error_msg_x_axis = "x_axis must be either 'age' or 'time'."
    ),
    "plot_constant is only relevant when `seromodel@model_name == 'constant'`"
  )
})
