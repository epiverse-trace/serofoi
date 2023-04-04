library(testthat)

# Define a test context
context("Testing plot_rhats function")

# Create a mock seromodel_object
mock_seromodel_object <- list(fit = "model did not run")

# Define a test case for the else statement
test_that("plot_rhats function works for else statement", {

  # Call the function with the mock seromodel_object
  rhats_plot <- plot_rhats(mock_seromodel_object)

  # Expect the output to be a ggplot object
  expect_is(rhats_plot, "ggplot")

  # Expect the plot to have a single point
  expect_equal(length(rhats_plot$layers[[1]]$data), 1)

  # Expect the plot to have a single label
  expect_equal(length(rhats_plot$layers[[2]]$data), 1)

  # Expect the label to be "errors"
  expect_equal(rhats_plot$layers[[2]]$label, "errors")
})
