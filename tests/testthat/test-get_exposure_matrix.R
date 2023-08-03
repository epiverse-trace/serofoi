library(testthat)

# Load the function to be tested
source("modelling.R")

# Define the test
test_that("The function get_exposure_matrix generates a correct exposure matrix", {
  # Define the test data frame
  serodata <- data.frame(
    age_mean_f = c(3, 2, 4),  # Example age groups (age classes)
    tsur = c(2023, 2022, 2019)  # Example reference years
  )

  # Call the function and store the result
  exposure_matrix <- get_exposure_matrix(serodata)

  # Check the dimensions of the exposure matrix
  expect_equal(nrow(exposure_matrix), length(serodata$age_mean_f))
  expect_equal(ncol(exposure_matrix), (serodata$tsur[1] - serodata$age_mean_f[1] + 1))

  # Check the correct values in the exposure matrix
  # For the first age group (3 years) and reference year 2023, the matrix should have 3 rows with 1s in the last 3 columns.
  expect_equal(exposure_matrix[1, ], c(0, 0, 0, 1, 1, 1))

  # For the second age group (2 years) and reference year 2022, the matrix should have 2 rows with 1s in the last 2 columns.
  expect_equal(exposure_matrix[2, ], c(0, 0, 1, 1, 1))

  # For the third age group (4 years) and reference year 2019, the matrix should have 4 rows with 1s in the last 4 columns.
  expect_equal(exposure_matrix[3, ], c(0, 1, 1, 1, 1))
})
