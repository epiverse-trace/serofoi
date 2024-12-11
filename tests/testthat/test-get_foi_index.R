# Setup for testing ----
# Sample survey features for testing
survey_features <- data.frame(
  age_min = c(1, 6, 11, 16, 21),
  age_max = c(5, 10, 15, 20, 25),
  survey_year = 2025
)

# Test get_foi_index ----
test_that("get_foi_index returns correct output for valid model types", {
  # Test for model_type "age"
  result_age <- get_foi_index(survey_features, group_size = 5, model_type = "age")

  # Check if the data frame has the correct structure for "age"
  checkmate::assert_names(
    names(result_age),
    must.include = c("age", "foi_index")
  )
  expect_equal(nrow(result_age), max(survey_features$age_max))

  # Test for model_type "time"
  result_time <- get_foi_index(survey_features, group_size = 5, model_type = "time")

  # Check if the data frame has the correct structure for "time"
  checkmate::assert_names(
    names(result_time),
    must.include = c("year", "foi_index")
  )
  expect_equal(nrow(result_time), max(survey_features$age_max))
})

test_that("get_foi_index returns an error for invalid model_type", {
  # Test for invalid model_type
  expect_error(
    get_foi_index(serosurvey, model_type = "constant"),
    regexp = "model_type must be either 'time' or 'age'"
  )
})

test_that("get_foi_index handles different group_size correctly", {
  # Test when group_size = 1 (edge case)
  result_1 <- get_foi_index(
    survey_features,
    group_size = 1,
    model_type = "age"
  )
  expected_1 <- seq(1, max(survey_features$age_max))

  expect_equal(nrow(result_1), max(survey_features$age_max))
  expect_equal(result_1$foi_index, expected_1)

  # Test when group_size equals the maximum age (edge case)
  result_max <- get_foi_index(
    survey_features,
    group_size = max(survey_features$age_max),
    model_type = "time"
  )
  expected_max <- rep(1, 25)
  expect_equal(nrow(result_max), max(survey_features$age_max))
  expect_equal(result_max$foi_index, expected_max)

  # Test when max age is not divisible by group_size
  result_no_div <- get_foi_index(
    survey_features,
    group_size = 7,
    model_type = "time"
  )
  # the remaining times are indexed in the last chunk
  expected_no_div <- c(rep(1, 7), rep(2, 7), rep(3, 7 + 4))
  expect_equal(nrow(result_no_div), max(survey_features$age_max))
  expect_equal(result_no_div$foi_index, expected_no_div)
})

test_that("get_foi_index throws an error for invalid group_size", {
  # Test for group_size > max age_max
  expect_error(
    get_foi_index(survey_features, group_size = 30, model_type = "age"),
    regexp = "Assertion on 'group_size' failed"
  )

  # Test for group_size less than 1
  expect_error(
    get_foi_index(survey_features, group_size = 0, model_type = "age"),
    regexp = "Assertion on 'group_size' failed"
  )
})
